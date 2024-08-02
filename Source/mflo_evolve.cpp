#include <AMReX_ParallelDescriptor.H>
#include <AMReX_ParmParse.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_FillPatchUtil.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_VisMF.H>
#include <AMReX_PhysBCFunct.H>
#include <AMReX_TimeIntegrator.H>
#include <kernels_3d.H>

#ifdef AMREX_MEM_PROFILING
#include <AMReX_MemProfiler.H>
#endif

#include <mflo.H>

// advance solution to final time
void mflo::EvolveAMR(Real final_time, bool only_flow)
{
    Real cur_time = t_new[0];

    for (int step = istep[0]; step < max_step && cur_time < final_time; ++step) 
    {
        amrex::Print() << "\nCoarse STEP " << step + 1 << " starts ..."
            << std::endl;
        
        if (potential_solve == 1 && step % pot_solve_int == 0)
        {
            solve_potential(cur_time);
        }
        if (mag_potential_solve == 1 && step % mag_pot_solve_int == 0)
        {
            solve_magnetic_vecpot(cur_time,0);
            solve_magnetic_vecpot(cur_time,1);
            solve_magnetic_vecpot(cur_time,2);
            update_Bfields();
        }
        
        ComputeDt();

        int lev = 0;
        int iteration = 1;
        timeStep(lev, cur_time, iteration, only_flow);

        // update primitive variables
        for (int l = 0; l <= finest_level; l++) 
        {
            update_primitive_vars(l);
        }
        
        if(clip_species)
        {
            for (int l = 0; l <= finest_level; l++) 
            {
                clip_neg_speciesconc(l);
            }
        }

        //print residuals
        if(track_residual_norms == 1)
        {
            compute_residual_norms();

            for(int c=0;c<residual_norms.size();c++)
            {
                PrintToFile("resnorms")<<residual_norms[c]<<"\t";
            }

            PrintToFile("resnorms")<<"\n";
        }

        cur_time += dt[0];

        amrex::Print() << "Coarse STEP " << step + 1 << " ends."
        << " TIME = " << cur_time << " DT = " << dt[0]
        << std::endl;

        // sync up time
        for (lev = 0; lev <= finest_level; ++lev) 
        {
            t_new[lev] = cur_time;
        }

        if (plot_int > 0 && (step + 1) % plot_int == 0) 
        {
            WritePlotFile();
        }

        if (chk_int > 0 && (step + 1) % chk_int == 0) 
        {
            WriteCheckpointFile();
        }

        if (cur_time >= stop_time - 1.e-6 * dt[0]) break;
    }

    if (plot_int > 0) 
    {
        WritePlotFile();
    }
}

void mflo::compute_residual_norms()
{
    int lev=0;
    MultiFab& S_new = phi_new[lev]; 
    MultiFab& S_old = phi_old[lev]; 
    MultiFab difference(grids[lev], dmap[lev], S_new.nComp(), 0);

    MultiFab::LinComb(
        difference, 1.0/dt[lev], S_new, 0, 
        -1.0/dt[lev], S_old, 0, 0, S_new.nComp(), 0);

    //scale flow variables by vfrac
    for(int var=RHO_INDX;var<VFRAC_INDX;var++)
    {
        MultiFab::Multiply(difference,S_old,VFRAC_INDX,var,1,0);
    }

    for(int c=0;c<residual_norms.size();c++)
    {
        residual_norms[c]=difference.norm2(c);
    }
}

void mflo::update_primitive_vars(int lev)
{
    MultiFab& S_new = phi_new[lev];
    bool nsflag=(do_ns==1)?true:false;
    int chtflag=conj_ht;

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    {
        for (MFIter mfi(S_new, TilingIfNotGPU()); mfi.isValid(); ++mfi) 
        {
            const Box& bx = mfi.tilebox();
            Array4<Real> snew_arr = S_new.array(mfi);

            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {

                //cant be greater than 1
                if(snew_arr(i,j,k,VFRAC_INDX) == one)
                {
                    Real u[NCVARS+NUM_SPECIES], p[NCVARS];

                    for (int c = 0; c < NCVARS; c++) 
                    {
                        u[c + RHO_IND] = snew_arr(i, j, k, c + RHO_INDX);
                    }
                    for (int c = 0; c < NUM_SPECIES; c++) 
                    {
                        u[c + NCVARS] = snew_arr(i, j, k, c + FLO_NVARS);
                    }
                    cons_to_prim(u, p);
                    for (int c = 0; c < NCVARS; c++) 
                    {
                        snew_arr(i, j, k, c + DENS_INDX) = p[c + DENS_IND];
                    }
                    snew_arr(i, j, k, TEMP_INDX)=mflo_thermo::get_t_from_rpc(snew_arr(i, j, k, DENS_INDX),
                                                                             snew_arr(i, j, k, PRES_INDX),u+NCVARS);
                }
                else
                {
                    //there is conjugate heat transfer and I am in the solid..
                    if(chtflag && snew_arr(i,j,k,VFRAC_INDX)==zeroval)
                    {
                        snew_arr(i,j,k,TEMP_INDX)=mflo_thermo::get_solid_t_from_rhoe(snew_arr(i,j,k,RHOE_INDX));
                    }
                }
            });

            Real minpressure=S_new[mfi].min<RunOn::Device>(PRES_INDX);
            if(minpressure < 0 and nsflag)
            {
                Print() << "Minimum pressure in domain:" << minpressure << "\n";
                amrex::Abort("Pressure has gone negative"); 
            }
        }
    }

}

// advance a level by dt
// includes a recursive call for finer levels
void mflo::timeStep(int lev, Real time, int iteration, bool only_flow)
{
    if (regrid_int > 0) // We may need to regrid
    {

        // help keep track of whether a level was already regridded
        // from a coarser level call to regrid
        static Vector<int> last_regrid_step(max_level + 1, 0);

        // regrid changes level "lev+1" so we don't regrid on max_level
        // also make sure we don't regrid fine levels again if
        // it was taken care of during a coarser regrid
        if (lev < max_level && istep[lev] > last_regrid_step[lev]) 
        {
            if (istep[lev] % regrid_int == 0) 
            {
                // regrid could add newly refine levels (if finest_level <
                // max_level) so we save the previous finest level index
                int old_finest = finest_level;
                regrid(lev, time);

                // mark that we have regridded this level already
                for (int k = lev; k <= finest_level; ++k) 
                {
                    last_regrid_step[k] = istep[k];
                }

                // if there are newly created levels, set the time step
                for (int k = old_finest + 1; k <= finest_level; ++k) 
                {
                    dt[k] = dt[k - 1] / MaxRefRatio(k - 1);
                }
            }
        }
    }
    

    if (Verbose()) 
    {
        amrex::Print() << "[Level " << lev << " step " << istep[lev] + 1 << "] ";
        amrex::Print() << "ADVANCE with time = " << t_new[lev]
        << " dt = " << dt[lev] << std::endl;
    }


    // advance a single level for a single time step, updates flux registers
    Advance_coupled_strang(lev, time, dt[lev], iteration, nsubsteps[lev],only_flow);

    ++istep[lev];

    if (Verbose()) 
    {
        amrex::Print() << "[Level " << lev << " step " << istep[lev] << "] ";
        amrex::Print() << "Advanced " << CountCells(lev) << " cells"
        << std::endl;
    }

    if (lev < finest_level) 
    {
        // recursive call for next-finer level
        for (int i = 1; i <= nsubsteps[lev + 1]; ++i) 
        {
            timeStep(lev + 1, time + (i - 1) * dt[lev + 1], i);
        }

        if (do_reflux) 
        {
            // update lev based on coarse-fine flux mismatch
            flux_reg[lev + 1]->Reflux(
                phi_new[lev], 1.0, 0, 0, phi_new[lev].nComp(), geom[lev]);
        }

        AverageDownTo(lev); // average lev+1 down to lev
    }
}

void mflo::Advance_coupled_strang(int lev, Real time, Real dt_lev, int iteration, int ncycle,bool only_flow)
{
    constexpr int num_grow = 3;
    std::swap(phi_old[lev], phi_new[lev]); // old becomes new and new becomes
                                           // old
    t_old[lev] = t_new[lev];        // old time is now current time (time)
    t_new[lev] += dt_lev;           // new time is ahead
    MultiFab& S_new = phi_new[lev]; // this is the old value, beware!
    MultiFab& S_old = phi_old[lev]; // current value

    //RK3 TVD scheme
    // State with ghost cells
    MultiFab Sborder(grids[lev], dmap[lev], S_new.nComp(), num_grow);
    // source term
    MultiFab dsdt_flow(grids[lev], dmap[lev], S_new.nComp(), 0);

    if (do_reflux)
    {
        if (flux_reg[lev + 1])
        {
            flux_reg[lev+1]->setVal(Real(0.0));
        }
    }


    // stage 1
    // time is current time which is t_old
    FillPatch(lev, time, Sborder, 0, Sborder.nComp());
    // compute dsdt for 1/2 timestep
    update_cutcell_data(lev, num_grow, Sborder, time, dt_lev);
    compute_dsdt_flow(lev, num_grow, Sborder, dsdt_flow, time, dt_lev, sixth, true);
    // S_new=S_old+dt*dsdt //sold is the current value
    MultiFab::LinComb(S_new, one, Sborder, 0, dt_lev, dsdt_flow, 0, 0, S_new.nComp(), 0);


    // stage 2
    // time+dt_lev lets me pick S_new for sborder
    FillPatch(lev, time + dt_lev, Sborder, 0, Sborder.nComp());
    update_cutcell_data(lev,num_grow,Sborder,time,dt_lev);
    compute_dsdt_flow(lev, num_grow, Sborder, dsdt_flow, time + dt_lev, dt_lev, sixth,true);

    // S_new=3/4 S_old+1/4 S_new + 1/4 dt*dsdt
    // S_new = S_new + dt*dsdt
    // S_new = S_new + 3*S_old
    // S_new =S_new/4
    MultiFab::Saxpy(S_new, dt_lev, dsdt_flow, 0, 0, S_new.nComp(), 0);
    MultiFab::Saxpy(S_new, Real(3.0), S_old, 0, 0, S_new.nComp(), 0);
    S_new.mult(fourth);

    // stage 3
    // time+dt_lev lets me pick S_new for sborder
    FillPatch(lev, time + dt_lev, Sborder, 0, Sborder.nComp());
    update_cutcell_data(lev,num_grow,Sborder,time,dt_lev);
    compute_dsdt_flow(lev, num_grow, Sborder, dsdt_flow, time + dt_lev, dt_lev, two3rd,true);
    // S_new=1/3 S_old+2/3 S_new + 2/3 dt*dsdt
    // S_new = S_new + dt*dsdt
    // S_new = S_new + 2*S_old
    // S_new =S_new/3
    MultiFab::Saxpy(S_new, dt_lev, dsdt_flow, 0, 0, S_new.nComp(), 0);
    MultiFab::Xpay(S_new, two, S_old, 0, 0, S_new.nComp(), 0);
    S_new.mult(third);

    //Advance chemistry here
    if(!only_flow)
    {
        Advance_chemistry_implicit(lev, time, dt_lev);
    }
}

void mflo::Advance_coupled(int lev, Real time, Real dt_lev, int iteration, int ncycle,bool only_flow)
{
    constexpr int num_grow = 3;
    std::swap(phi_old[lev], phi_new[lev]); // old becomes new and new becomes
    // old
    t_old[lev] = t_new[lev];        // old time is now current time (time)
    t_new[lev] += dt_lev;           // new time is ahead
    MultiFab& S_new = phi_new[lev]; // this is the old value, beware!
    MultiFab& S_old = phi_old[lev]; // current value

    //RK3 TVD scheme
    // State with ghost cells
    MultiFab Sborder(grids[lev], dmap[lev], S_new.nComp(), num_grow);
    // source term
    MultiFab dsdt_flow(grids[lev], dmap[lev], S_new.nComp(), 0);
    MultiFab dsdt_chemistry(grids[lev], dmap[lev], S_new.nComp(), 0);

    if (do_reflux)
    {
        if (flux_reg[lev + 1]) 
        {
            flux_reg[lev+1]->setVal(Real(0.0)); 
        }
    } 


    // stage 1
    // time is current time which is t_old
    FillPatch(lev, time, Sborder, 0, Sborder.nComp());
    // compute dsdt for 1/2 timestep
    update_cutcell_data(lev, num_grow, Sborder, time, dt_lev);
    compute_dsdt_flow(lev, num_grow, Sborder, dsdt_flow, time, dt_lev, sixth, true);
    if(!only_flow)
    {
        compute_dsdt_chemistry(lev, num_grow, Sborder, dsdt_chemistry, time);
    }
    // S_new=S_old+dt*dsdt //sold is the current value
    MultiFab::LinComb(S_new, one, Sborder, 0, dt_lev, dsdt_flow, 0, 0, S_new.nComp(), 0);

    if(!only_flow)
    {
        MultiFab::Saxpy(S_new, dt_lev, dsdt_chemistry, 0, 0, S_new.nComp(), 0);
    }

    // stage 2
    // time+dt_lev lets me pick S_new for sborder
    FillPatch(lev, time + dt_lev, Sborder, 0, Sborder.nComp());
    update_cutcell_data(lev,num_grow,Sborder,time,dt_lev);
    compute_dsdt_flow(lev, num_grow, Sborder, dsdt_flow, time + dt_lev, dt_lev, sixth,true);
    if(!only_flow)
    {
        compute_dsdt_chemistry(lev, num_grow, Sborder, dsdt_chemistry, time);
    }

    // S_new=3/4 S_old+1/4 S_new + 1/4 dt*dsdt
    // S_new = S_new + dt*dsdt
    // S_new = S_new + 3*S_old
    // S_new =S_new/4
    MultiFab::Saxpy(S_new, dt_lev, dsdt_flow, 0, 0, S_new.nComp(), 0);
    if(!only_flow)
    {
        MultiFab::Saxpy(S_new, dt_lev, dsdt_chemistry, 0, 0, S_new.nComp(), 0);
    }
    MultiFab::Saxpy(S_new, Real(3.0), S_old, 0, 0, S_new.nComp(), 0);
    S_new.mult(fourth);

    // stage 3
    // time+dt_lev lets me pick S_new for sborder
    FillPatch(lev, time + dt_lev, Sborder, 0, Sborder.nComp());
    update_cutcell_data(lev,num_grow,Sborder,time,dt_lev);
    //NOTE: Using time + dt_lev might be incorrect here, need to verify using beta from the RK3 coefficients
    //According to one source, this might need to be t + 0.5*dt_lev
    compute_dsdt_flow(lev, num_grow, Sborder, dsdt_flow, time + dt_lev, dt_lev, two3rd,true);
    if(!only_flow)
    {
        compute_dsdt_chemistry(lev, num_grow, Sborder, dsdt_chemistry, time);
    }
    //note: use Xpay instead of saxpy to reduce complications
    // S_new=1/3 S_old+2/3 S_new + 2/3 dt*dsdt
    // S_new = S_new + dt*dsdt
    // S_new = S_new*2
    // S_new = S_new + S_old
    // S_new =S_new/3
    MultiFab::Saxpy(S_new, dt_lev, dsdt_flow, 0, 0, S_new.nComp(), 0);
    if(!only_flow)
    {
        MultiFab::Saxpy(S_new, dt_lev, dsdt_chemistry, 0, 0, S_new.nComp(), 0);
    }
    S_new.mult(two);
    MultiFab::Saxpy(S_new, Real(1.0), S_old, 0, 0, S_new.nComp(), 0);
    S_new.mult(third);

    /*old RK2 scheme
    // State with ghost cells
    MultiFab Sborder(grids[lev], dmap[lev], S_new.nComp(), num_grow);
    // source term
    MultiFab dsdt_flow(grids[lev], dmap[lev], S_new.nComp(), 0);
    MultiFab dsdt_chemistry(grids[lev], dmap[lev], S_new.nComp(), 0);


    // stage 1
    // time is current time which is t_old
    FillPatch(lev, time, Sborder, 0, Sborder.nComp());
    // compute dsdt for 1/2 timestep
    update_cutcell_data(lev,num_grow,Sborder,time+0.5*dt_lev,dt_lev);
    compute_dsdt_flow(lev, num_grow, Sborder, dsdt_flow, time, half * dt_lev, false);
    if(!only_flow)
    {
    compute_dsdt_chemistry(lev, num_grow, Sborder, dsdt_chemistry, time);
    }
    // S_new=S_old+0.5*dt*dsdt //sold is the current value
    MultiFab::LinComb(S_new, one, Sborder, 0, half * dt_lev, dsdt_flow, 0, 0, S_new.nComp(), 0);

    if(!only_flow)
    {
    MultiFab::Saxpy(S_new, half*dt_lev, dsdt_chemistry, 0, 0, S_new.nComp(), 0);
    }

    // stage 2
    // time+dt_lev lets me pick S_new for sborder
    FillPatch(lev, time + dt_lev, Sborder, 0, Sborder.nComp());
    // dsdt for full time-step
    update_cutcell_data(lev,num_grow,Sborder,time+0.5*dt_lev,dt_lev);
    compute_dsdt_flow(lev, num_grow, Sborder, dsdt_flow, time + half * dt_lev, dt_lev, true);
    if(!only_flow)
    {
    compute_dsdt_chemistry(lev, num_grow, Sborder, dsdt_chemistry, time + half * dt_lev);
    }
    // S_new=S_old+dt*dsdt
    MultiFab::LinComb(S_new, one, S_old, 0, dt_lev, dsdt_flow, 0, 0, S_new.nComp(), 0);
    if(!only_flow)
    {
    MultiFab::Saxpy(S_new, dt_lev, dsdt_chemistry, 0, 0, S_new.nComp(), 0);
    }*/
}
/*
Advances the chemisty state from time -> time + dt_lev
This is done using TimeIntegrator from AMReX, which can
be either implicit or explicit depending on the parameters
in the input file. Implicit currently requires an existing
SUNDIALS installation, with USE_SUNDIALS=TRUE at compile time,
along with setting the following in the input file

Explicit (no sundials):
   integration.type = RungeKutta
   integration.rk.type = 3 (for 3rd order)
Implicit (with sundials):
   integration.type = SUNDIALS
   integration.sundials.strategy = CVODE
 */
void mflo::Advance_chemistry_implicit(int lev, Real time, Real dt_lev)
{
    constexpr int num_grow = 3;
    std::swap(phi_old[lev], phi_new[lev]); // old becomes new and new becomes old
    MultiFab& S_new = phi_new[lev]; // old value
    MultiFab& S_old = phi_old[lev]; // current value

    auto rhs_function = [&] ( Vector<MultiFab> & dSdt_vec, const Vector<MultiFab>& S_vec, const Real time) {
        auto & dSdt = dSdt_vec[0];
        MultiFab S(S_vec[0], amrex::make_alias, 0, S_vec[0].nComp());
        compute_dsdt_chemistry(lev, num_grow, S, dSdt, time);
    };
    Vector<MultiFab> state_old, state_new;
    // This term has the current state
    state_old.push_back(MultiFab(S_old, amrex::make_alias, 0, S_new.nComp()));
    // This is where the integrator puts the new state, hence aliased to S_new
    state_new.push_back(MultiFab(S_new, amrex::make_alias, 0, S_new.nComp()));
    // Define the integrator
    TimeIntegrator<Vector<MultiFab>> integrator(state_old);
    integrator.set_rhs(rhs_function);
    // Advance from time to time + dt_lev
    //S_new/phi_new should have the new state
    integrator.advance(state_old, state_new, time, dt_lev); 
}

void mflo::update_cutcell_data(
    int lev,
    const int num_grow,
    MultiFab& Sborder,
    Real time,
    Real dtstep)
{
    int ncomp = Sborder.nComp();
    auto prob_lo = geom[lev].ProbLoArray();
    const auto dx = geom[lev].CellSizeArray();
    int using_inert_gas=using_bg_inertgas;
    int spec_in_solid=species_in_solid;
    int chtflag=conj_ht;

    for (MFIter mfi(Sborder, TilingIfNotGPU()); mfi.isValid(); ++mfi) 
    {
        const Box& bx = mfi.tilebox();
        Array4<Real> sborder_arr = Sborder.array(mfi);

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
            update_cutcells(i,j,k,sborder_arr,dx,prob_lo,time,
                            spec_in_solid,using_inert_gas,chtflag);
        });
    }
}

void mflo::clip_neg_speciesconc(int lev)
{
    MultiFab& Snew = phi_new[lev];

    for (MFIter mfi(Snew); mfi.isValid(); ++mfi) 
    {
        const Box& bx = mfi.tilebox();
        Array4<Real> snew_arr = Snew.array(mfi);

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
            for(int sp=0;sp<NUM_SPECIES;sp++)
            {
                if(snew_arr(i,j,k,FLO_NVARS+sp) < 0.0)
                {
                    snew_arr(i,j,k,FLO_NVARS+sp)=0.0;
                }
            }
        });
    }
}

// advance a single level for a single time step, updates flux registers
void mflo::compute_dsdt_flow(
    int lev,
    const int num_grow,
    MultiFab& Sborder,
    MultiFab& dsdt,
    Real time,
    Real tstep,
    Real fluxfactor,
    bool reflux_this_stage)
{

    const auto dx = geom[lev].CellSizeArray();
    const auto prob_lo = geom[lev].ProbLoArray();
    bool nsflag=(do_ns==1)?true:false;
    int conjhtflag=conj_ht;

    int ncomp = Sborder.nComp();
    int hyperbolics_order = order_hyp;
    Real hyperbolics_dissfactor = dissfactor;
    int spec_in_solid=species_in_solid;
    dsdt.setVal(zeroval);

    // Build temporary multiFabs to work on.
    Array<MultiFab, AMREX_SPACEDIM> flux;
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) 
    {
        BoxArray ba =
        amrex::convert(dsdt.boxArray(), IntVect::TheDimensionVector(idim));
        flux[idim].define(ba, dsdt.DistributionMap(), ncomp, 0);
        flux[idim].setVal(zeroval);
    }

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    {
        for (MFIter mfi(dsdt, TilingIfNotGPU()); mfi.isValid(); ++mfi) 
        {
            const Box& bx = mfi.tilebox();
            const Box& gbx = amrex::grow(bx, num_grow);

            FArrayBox fluid_transport_fab(gbx,2); //viscosity and thermal-conductivity
            FArrayBox specdiff_fab(gbx,NUM_SPECIES); //species diffusion
            FArrayBox source_fab(bx,TOTAL_NVARS); //external sources
            source_fab.setVal<RunOn::Device>(zeroval);

            Elixir fluid_transport_fab_eli = fluid_transport_fab.elixir();
            Elixir specdiff_fab_eli = specdiff_fab.elixir();
            Elixir source_fab_eli = source_fab.elixir();

            Array4<Real> sborder_arr = Sborder.array(mfi);
            Array4<Real> fluid_transport_arr = fluid_transport_fab.array();
            Array4<Real> specdiff_arr = specdiff_fab.array();
            Array4<Real> dsdt_arr = dsdt.array(mfi);
            Array4<Real> source_arr = source_fab.array();

            GpuArray<Array4<Real>, AMREX_SPACEDIM> flux_arr{AMREX_D_DECL(
                    flux[0].array(mfi), flux[1].array(mfi), flux[2].array(mfi))};

            amrex::ParallelFor(gbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) { 
                mflo_user_funcs::compute_fluid_transport(
                    i, j, k, sborder_arr, fluid_transport_arr, prob_lo, dx, time);
            });

            amrex::ParallelFor(gbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) { 
                mflo_chem_transport::compute_spec_dcoeff(
                    i, j, k, sborder_arr, specdiff_arr, prob_lo, dx, time);
            });

            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                mflo_user_funcs::compute_fluid_source(
                    i, j, k, sborder_arr, source_arr, prob_lo, dx, time);
            });

            amrex::ParallelFor(
                amrex::growHi(bx, 0, 1),
                [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                    compute_flux(
                        i, j, k, XDIR, sborder_arr, fluid_transport_arr, 
                        specdiff_arr, flux_arr[0], 
                        dx, hyperbolics_order,hyperbolics_dissfactor,nsflag,spec_in_solid);
                });

            amrex::ParallelFor(
                amrex::growHi(bx, 1, 1),
                [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                    compute_flux(
                        i, j, k, YDIR, sborder_arr, fluid_transport_arr, 
                        specdiff_arr, flux_arr[1], 
                        dx, hyperbolics_order,hyperbolics_dissfactor,nsflag,spec_in_solid);
                });

            amrex::ParallelFor(
                amrex::growHi(bx, 2, 1),
                [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                    compute_flux(
                        i, j, k, ZDIR, sborder_arr, fluid_transport_arr, 
                        specdiff_arr, flux_arr[2], 
                        dx, hyperbolics_order,hyperbolics_dissfactor,nsflag,spec_in_solid);
                });

            // update residual
            amrex::ParallelFor(
                bx, ncomp, [=] AMREX_GPU_DEVICE(int i, int j, int k, int n) {
                    update_residual(
                        i, j, k, n, dsdt_arr, source_arr, sborder_arr,
                        AMREX_D_DECL(flux_arr[0], flux_arr[1], flux_arr[2]),
                        dx,spec_in_solid,conjhtflag);
                });
        }
    }

    if (do_reflux and reflux_this_stage) 
    {
        if (flux_reg[lev + 1]) 
        {
            for (int i = 0; i < BL_SPACEDIM; ++i) 
            {
                // update the lev+1/lev flux register (index lev+1)
                const Real dA =
                (i == 0) ? dx[1] * dx[2]
                : ((i == 1) ? dx[0] * dx[2] : dx[0] * dx[1]);
                const Real scale = -tstep * dA * fluxfactor;
                flux_reg[lev + 1]->CrseInit(flux[i], i, 0, 0, ncomp, scale, FluxRegister::ADD);
            }
        }
        if (flux_reg[lev]) 
        {
            for (int i = 0; i < BL_SPACEDIM; ++i) 
            {
                // update the lev/lev-1 flux register (index lev)
                const Real dA =
                (i == 0) ? dx[1] * dx[2]
                : ((i == 1) ? dx[0] * dx[2] : dx[0] * dx[1]);
                const Real scale = tstep * dA * fluxfactor;
                flux_reg[lev]->FineAdd(flux[i], i, 0, 0, ncomp, scale);
            }
        }
    }
}

void mflo::compute_dsdt_chemistry(
    int lev,
    const int num_grow,
    MultiFab& S,
    MultiFab& dsdt,
    Real time)
{
    dsdt.setVal(zeroval);

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    {
        for (MFIter mfi(dsdt, TilingIfNotGPU()); mfi.isValid(); ++mfi) 
        {
            const Box& bx = mfi.tilebox();

            FArrayBox source_fab(bx,TOTAL_NVARS); //external sources
            source_fab.setVal<RunOn::Device>(zeroval);
            Elixir source_fab_eli = source_fab.elixir();

            Array4<Real> s_arr = S.array(mfi);
            Array4<Real> dsdt_arr = dsdt.array(mfi);
            Array4<Real> source_arr = source_fab.array();

            auto prob_lo = geom[lev].ProbLoArray();
            const auto dx = geom[lev].CellSizeArray();

            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                mflo_chem_reactions::compute_spec_source(
                    i, j, k, s_arr, source_arr, prob_lo, dx, time);
            });

            dsdt[mfi].plus<RunOn::Device>(source_fab);
        }
    }
}

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

//Note: This split evolve is used when chemistry is happening 
//at a very slow timescale compared to flow. 
//the strategy is to solve the flow to steady-state, advance the
//chemistry, and then solve the flow to steady state again.
//the technique is better described in
//Sitaraman et al., Chem Engg Res. and Design, 197, 2023
//Sitaraman et al., Chemical Engineering Science 206 (2019): 348-360
void mflo::Evolve_split(Real t_ss,Real t_react,Real dt_react,Real dt_react_rk,Real dt_ss)
{
    Real cur_time = t_new[0];
    bool only_flow=true;
    Real react_time=0.0;
    Real flow_time;
    Real tchem;

    amrex::Print()<<"running flow to steady-state\n";
    EvolveAMR(t_ss,only_flow);
    flow_time=t_ss;

    while(react_time < t_react)
    {
        amrex::Print()<<"advance chemistry\n";
        //temporarily store inert gas concentration in dens index
        if(using_bg_inertgas)
        {
            for (int l = 0; l <= finest_level; l++) 
            {
                store_inertgas_conc(l);
            }
        }
        for (int l = 0; l <= finest_level; l++) 
        {
            tchem=0.0;
            while(tchem < dt_react)
            {
                Advance_chemistry(l, react_time+tchem, dt_react_rk);
                tchem += dt_react_rk;
            }
        }
        AverageDown();

        for (int l = 0; l <= finest_level; l++) 
        {
            //do update of flow variables
            update_vars_after_chemsolve(l);
            //update_primitive_vars(l);
        }

        amrex::Print()<<"advance flow\n";
        EvolveAMR(flow_time+dt_ss,only_flow);

        react_time += dt_react;
        flow_time += dt_ss;
    }

}

void mflo::update_vars_after_chemsolve(int lev)
{
    MultiFab& S_new = phi_new[lev];
    int using_inert_gas=using_bg_inertgas;
    bool nsflag=(do_ns==1)?true:false;

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    {
        for (MFIter mfi(S_new, TilingIfNotGPU()); mfi.isValid(); ++mfi) 
        {
            const Box& bx = mfi.tilebox();
            Array4<Real> snew_arr = S_new.array(mfi);

            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {

                if(snew_arr(i,j,k,VFRAC_INDX) >= one)
                {
                    Real u[NCVARS+NUM_SPECIES], p[NCVARS],spec[NUM_SPECIES];
                    for (int c = 0; c < NUM_SPECIES; c++) 
                    {
                        spec[c] = snew_arr(i, j, k, c + FLO_NVARS);
                    }

                    //also pass background gas concentration
                    if(using_inert_gas)
                    {
                        snew_arr(i,j,k,RHO_INDX)  = mflo_thermo::get_r_from_c(spec,snew_arr(i,j,k,DENS_INDX));
                        snew_arr(i,j,k,DENS_INDX) = snew_arr(i,j,k,RHO_INDX);
                    }
                    else
                    {
                        snew_arr(i,j,k,RHO_INDX)  = mflo_thermo::get_r_from_c(spec);
                        snew_arr(i,j,k,DENS_INDX) = snew_arr(i,j,k,RHO_INDX);
                    }
                    snew_arr(i,j,k,RHOU_INDX) = snew_arr(i,j,k,RHO_INDX)*snew_arr(i,j,k,VELX_INDX);
                    snew_arr(i,j,k,RHOV_INDX) = snew_arr(i,j,k,RHO_INDX)*snew_arr(i,j,k,VELY_INDX);
                    snew_arr(i,j,k,RHOW_INDX) = snew_arr(i,j,k,RHO_INDX)*snew_arr(i,j,k,VELZ_INDX);

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

void mflo::store_inertgas_conc(int lev)
{
    MultiFab& S_new = phi_new[lev];

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    {
        for (MFIter mfi(S_new, TilingIfNotGPU()); mfi.isValid(); ++mfi) 
        {
            const Box& bx = mfi.tilebox();
            Array4<Real> snew_arr = S_new.array(mfi);

            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {

                if(snew_arr(i,j,k,VFRAC_INDX) >= one)
                {
                    Real spec[NUM_SPECIES];
                    for (int c = 0; c < NUM_SPECIES; c++) 
                    {
                        spec[c] = snew_arr(i, j, k, c + FLO_NVARS);
                    }
                    snew_arr(i, j, k, DENS_INDX)=
                    mflo_thermo::get_bgasconc_from_rc(snew_arr(i,j,k,RHO_INDX),spec);
                }
            });
        }
    }
}

void mflo::Advance_chemistry(int lev, Real time, Real dt_lev)
{
    constexpr int num_grow = 3;
    std::swap(phi_old[lev], phi_new[lev]); // old becomes new and new becomes old
    MultiFab& S_new = phi_new[lev]; // this is the old value, beware!
    MultiFab& S_old = phi_old[lev]; // current value

    // source term
    MultiFab dsdt_chemistry(grids[lev], dmap[lev], S_new.nComp(), 0);


    // stage 1
    // compute dsdt for 1/2 timestep
    compute_dsdt_chemistry(lev, num_grow, S_new, dsdt_chemistry, time);
    // S_new=S_old+0.5*dt*dsdt //sold is the current value
    MultiFab::LinComb(S_new, one, S_old, 0, half * dt_lev, dsdt_chemistry, 0, 0, S_new.nComp(), 0);

    // stage 2
    // dsdt for full time-step
    compute_dsdt_chemistry(lev, num_grow, S_new, dsdt_chemistry, time + half * dt_lev);
    // S_new=S_old+dt*dsdt
    MultiFab::LinComb(S_new, one, S_old, 0, dt_lev, dsdt_chemistry, 0, 0, S_new.nComp(), 0);

}

#include <AMReX_ParallelDescriptor.H>
#include <AMReX_ParmParse.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_FillPatchUtil.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_VisMF.H>
#include <AMReX_PhysBCFunct.H>
#include <kernels_3d.H>

#ifdef AMREX_MEM_PROFILING
#include <AMReX_MemProfiler.H>
#endif

#include <mflo.H>
#include <userfuncs.H>

// a wrapper for EstTimeStep
void mflo::ComputeDt()
{
    Vector<Real> dt_tmp(finest_level + 1);

    for (int lev = 0; lev <= finest_level; ++lev) 
    {
        dt_tmp[lev] = EstTimeStep(lev);
    }
    ParallelDescriptor::ReduceRealMin(&dt_tmp[0], dt_tmp.size());

    constexpr Real change_max = 1.1;
    Real dt_0 = dt_tmp[0];
    int n_factor = 1;
    for (int lev = 0; lev <= finest_level; ++lev) 
    {
        dt_tmp[lev] = std::min(dt_tmp[lev], change_max * dt[lev]);
        n_factor *= nsubsteps[lev];
        dt_0 = std::min(dt_0, n_factor * dt_tmp[lev]);
    }

    // Limit dt's by the value of stop_time.
    const Real eps = 1.e-3 * dt_0;
    if (t_new[0] + dt_0 > stop_time - eps) 
    {
        dt_0 = stop_time - t_new[0];
    }

    dt[0] = dt_0;
    for (int lev = 1; lev <= finest_level; ++lev) 
    {
        dt[lev] = dt[lev - 1] / nsubsteps[lev];
    }
}

// compute dt from CFL considerations
Real mflo::EstTimeStep(int lev)
{
    BL_PROFILE("mflo::EstTimeStep()");

    Real dt_est = std::numeric_limits<Real>::max();

    const auto dx = geom[lev].CellSizeArray();
    //    const Real* prob_lo = geom[lev].ProbLo();
    const Real cur_time = t_new[lev];
    MultiFab& S_new = phi_new[lev];

    MultiFab wavespeed(S_new.boxArray(), S_new.DistributionMap(), 1, 0);
    MultiFab specdiffmax(S_new.boxArray(), S_new.DistributionMap(), 1, 0);
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    {
        for (MFIter mfi(S_new, TilingIfNotGPU()); mfi.isValid(); ++mfi) 
        {
            const Box& bx = mfi.tilebox();
            Array4<Real> state_array = S_new.array(mfi);
            Array4<Real> wavespeed_array = wavespeed.array(mfi);
            Array4<Real> specdiffmax_array = specdiffmax.array(mfi);
            auto prob_lo = geom[lev].ProbLoArray();
            
            FArrayBox specdiff_fab(bx,NUM_SPECIES); 
            Elixir specdiff_fab_eli   = specdiff_fab.elixir();
            Array4<Real> specdiff_arr = specdiff_fab.array();
            Real do_ns_flag=do_ns;

            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) 
            {
                wavespeed_array(i,j,k,0) = zeroval;
#ifdef IMM_BOUNDARY
                if(state_array(i,j,k,VFRAC_INDX) > zeroval)
#endif
                {
                    Real fluid_speed = sqrt(
                    state_array(i, j, k, VELX_INDX) *
                        state_array(i, j, k, VELX_INDX) +
                    state_array(i, j, k, VELY_INDX) *
                        state_array(i, j, k, VELY_INDX) +
                    state_array(i, j, k, VELZ_INDX) *
                        state_array(i, j, k, VELZ_INDX));

                    Real spec[NUM_SPECIES];

                    for(int sp=0;sp<NUM_SPECIES;sp++)
                    {
                        spec[sp]=state_array(i,j,k,FLO_NVARS+sp);
                    }
                    Real gama = mflo_thermo::get_gama_from_rpc(state_array(i,j,k,DENS_INDX),
                       state_array(i,j,k,PRES_INDX), spec);
                    Real sound_speed = sqrt(
                    gama * state_array(i, j, k, PRES_INDX) /
                    state_array(i, j, k, DENS_INDX));
    
                    //don't use sound speed if navier-stokes is off
                    wavespeed_array(i, j, k, 0) = fluid_speed + do_ns_flag*sound_speed;
                 }
            });

            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) 
                    {
                        mflo_chem_transport::compute_spec_dcoeff(
                            i, j, k, state_array, specdiff_arr, prob_lo, dx, cur_time);
                        Real maxdcoeff=specdiff_arr(i,j,k,0);
                        for(int sp=0;sp<NUM_SPECIES;sp++)
                        {
                            if(specdiff_arr(i,j,k,sp) > maxdcoeff)
                            {
                                maxdcoeff=specdiff_arr(i,j,k,sp);
                            }       
                        }
                   
                        specdiffmax_array(i,j,k,0)=maxdcoeff;
                    
                    });
        }
    }

    Real max_wavespeed = wavespeed.norm0(0,0,true);
    Real max_specdiff  = specdiffmax.norm0(0,0,true);
    if(max_wavespeed > 0)
    {
        for (int i = 0; i < AMREX_SPACEDIM; i++) 
        {
            dt_est = std::min(dt_est, dx[i] / max_wavespeed);
        }
    }
    if(max_specdiff > 0)
    {
        for (int i = 0; i < AMREX_SPACEDIM; i++) 
        {
            dt_est = std::min(dt_est, half*dx[i]*dx[i]/max_specdiff/AMREX_SPACEDIM);
        }
    }

    ParallelDescriptor::ReduceRealMin(dt_est);

    dt_est *= cfl;
    return dt_est;
}

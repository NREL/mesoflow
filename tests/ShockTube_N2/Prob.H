#ifndef _PROB_H_
#define _PROB_H_

#include <AMReX_Box.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_Geometry.H>
#include <userfuncs.H>
#include <globalDefines.H>
#include <mflo_constants.H>

using namespace amrex;

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void initdata(Box const& bx, Array4<Real> const& phi, GeometryData const& geomdata)
{
    const auto lo = lbound(bx);
    const auto hi = ubound(bx);    

    const Real* AMREX_RESTRICT prob_lo = geomdata.ProbLo();
    const Real* AMREX_RESTRICT dx      = geomdata.CellSize();

#ifdef _OPENMP
#pragma omp parallel for collapse(2) if (GPU::notInLaunchRegion)
#endif
    for (int k = lo.z; k <= hi.z; ++k) 
    {
        for (int j = lo.y; j <= hi.y; ++j) 
        {
            Real z = prob_lo[ZDIR] + (half+k) * dx[ZDIR];
            Real y = prob_lo[YDIR] + (half+j) * dx[YDIR];

            AMREX_PRAGMA_SIMD
            for (int i = lo.x; i <= hi.x; ++i) 
            {
                Real x = prob_lo[XDIR] + (half+i) * dx[XDIR];

                phi(i,j,k,VFRAC_INDX) = one;
                phi(i,j,k,VELX_INDX) = zeroval;
                phi(i,j,k,VELY_INDX) = zeroval;
                phi(i,j,k,VELZ_INDX) = zeroval;
                phi(i,j,k,RHOU_INDX) = zeroval;
                phi(i,j,k,RHOV_INDX) = zeroval;
                phi(i,j,k,RHOW_INDX) = zeroval;

                if(fabs(x) < half and fabs(y) < half and fabs(z) < half)
                {
                    phi(i,j,k,PRES_INDX) = mflo_user_funcs::p_l;
                    phi(i,j,k,DENS_INDX) = mflo_user_funcs::rho_l;
                    phi(i,j,k,TEMP_INDX) = mflo_user_funcs::T_l;

                    for(int sp=0;sp<NUM_SPECIES;sp++)
                    {
                        phi(i,j,k,FLO_NVARS+sp)=mflo_user_funcs::spec_l[sp];
                    }

                    phi(i,j,k,RHO_INDX)  = phi(i,j,k,DENS_INDX);
                    phi(i,j,k,RHOE_INDX) = phi(i,j,k,DENS_INDX)*mflo_thermo::get_e_from_rpc(mflo_user_funcs::rho_l,
                                                          mflo_user_funcs::p_l,
                                                          mflo_user_funcs::spec_l);
                
                }
                else
                {
                    phi(i,j,k,PRES_INDX) = mflo_user_funcs::p_r;
                    phi(i,j,k,DENS_INDX) = mflo_user_funcs::rho_r;
                    phi(i,j,k,TEMP_INDX) = mflo_user_funcs::T_r;

                    for(int sp=0;sp<NUM_SPECIES;sp++)
                    {
                        phi(i,j,k,FLO_NVARS+sp)=mflo_user_funcs::spec_r[sp];
                    }

                    phi(i,j,k,RHO_INDX)  = phi(i,j,k,DENS_INDX);
                    phi(i,j,k,RHOE_INDX) = phi(i,j,k,DENS_INDX)*mflo_thermo::get_e_from_rpc(mflo_user_funcs::rho_r,
                                                          mflo_user_funcs::p_r,
                                                          mflo_user_funcs::spec_r);
                }
            }
        }
    }
}
#endif

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
    const Real* AMREX_RESTRICT prob_hi = geomdata.ProbHi();
    const Real* AMREX_RESTRICT dx      = geomdata.CellSize();
    Real spec[NUM_SPECIES];

#ifdef _OPENMP
#pragma omp parallel for collapse(2) if (GPU::notInLaunchRegion)
#endif
    for (int k = lo.z; k <= hi.z; ++k) 
    {
        Real z = prob_lo[ZDIR] + (half+k) * dx[ZDIR];
        for (int j = lo.y; j <= hi.y; ++j) 
        {
            Real y = prob_lo[YDIR] + (half+j) * dx[YDIR];

            AMREX_PRAGMA_SIMD
            for (int i = lo.x; i <= hi.x; ++i) 
            {
                Real x    = prob_lo[XDIR] + (half+i) * dx[XDIR];

                Real v0   = mflo_user_funcs::v0;
                Real p0   = one;
                Real rho0 = mflo_user_funcs::rho0;

                
                phi(i,j,k,VFRAC_INDX) = one;
                phi(i,j,k,VELX_INDX) =  v0*sin(PI*x)*cos(PI*y)*cos(PI*z);
                phi(i,j,k,VELY_INDX) = -v0*cos(PI*x)*sin(PI*y)*cos(PI*z);
                phi(i,j,k,VELZ_INDX) =  zeroval;

                phi(i,j,k,PRES_INDX) = p0 + sixteenth*rho0*v0*v0
                    *(cos(two*PI*x) +cos(two*PI*y)) *
                     (cos(two*PI*z) + two);

                phi(i,j,k,DENS_INDX) = phi(i,j,k,PRES_INDX); //R and T are 1
                spec[AIR_ID] = phi(i,j,k,DENS_INDX)/mflo_species::molwts[AIR_ID];

                phi(i,j,k,TEMP_INDX) = mflo_thermo::get_t_from_rpc(
                        phi(i,j,k,DENS_INDX),phi(i,j,k,PRES_INDX),spec);

                phi(i,j,k,RHO_INDX)  = phi(i,j,k,DENS_INDX);
                phi(i,j,k,RHOU_INDX) = phi(i,j,k,DENS_INDX)*phi(i,j,k,VELX_INDX);
                phi(i,j,k,RHOV_INDX) = phi(i,j,k,DENS_INDX)*phi(i,j,k,VELY_INDX);
                phi(i,j,k,RHOW_INDX) = phi(i,j,k,DENS_INDX)*phi(i,j,k,VELZ_INDX);

                phi(i,j,k,RHOE_INDX)  = phi(i,j,k,DENS_INDX)*
                    mflo_thermo::get_e_from_rpc(phi(i,j,k,DENS_INDX),phi(i,j,k,PRES_INDX),spec);

                phi(i,j,k,RHOE_INDX) += half*phi(i,j,k,DENS_INDX)*
                                        ( phi(i,j,k,VELX_INDX)*phi(i,j,k,VELX_INDX)
                                        + phi(i,j,k,VELY_INDX)*phi(i,j,k,VELY_INDX)
                                        + phi(i,j,k,VELZ_INDX)*phi(i,j,k,VELZ_INDX) );

                for(int sp=0;sp<NUM_SPECIES;sp++)
                {
                    phi(i,j,k,FLO_NVARS+sp) = spec[sp];
                }
            }
        }
    }
}
#endif

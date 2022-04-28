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

    Real rad=mflo_user_funcs::cylinder_rad;

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
                Real x  = prob_lo[XDIR] + (half+i) * dx[XDIR];
                
                phi(i,j,k,PRES_INDX) = mflo_user_funcs::pres_fs; 
                phi(i,j,k,DENS_INDX) = mflo_user_funcs::dens_fs;
                phi(i,j,k,TEMP_INDX) = mflo_user_funcs::temp_fs;
                
                phi(i,j,k,VFRAC_INDX) = zeroval;
                for(int kk=0;kk<2;kk++)
                {
                    for(int jj=0;jj<2;jj++)
                    {
                        for(int ii=0;ii<2;ii++)
                        {
                           Real x  = prob_lo[XDIR] + (i+ii) * dx[XDIR];
                           Real y  = prob_lo[YDIR] + (j+jj) * dx[YDIR];
                           Real z  = prob_lo[ZDIR] + (k+kk) * dx[ZDIR];
                           if(x*x + y*y + z*z > rad*rad)
                           {
                                phi(i,j,k,VFRAC_INDX)+=one;
                           }
                                
                        }
                    }
                }
                phi(i,j,k,VFRAC_INDX)*=eighth;

               if(phi(i,j,k,VFRAC_INDX)<one)
               {
                    phi(i,j,k,VELX_INDX) =  zeroval;
                    phi(i,j,k,VELY_INDX) =  zeroval;
                    phi(i,j,k,VELZ_INDX) =  zeroval;
               } 
               else
               {
                    phi(i,j,k,VELX_INDX) =  mflo_user_funcs::vel_fs;
                    phi(i,j,k,VELY_INDX) =  zeroval;
                    phi(i,j,k,VELZ_INDX) =  zeroval;
                }

                for(int sp=0;sp<NUM_SPECIES;sp++)
                {
                    phi(i,j,k,FLO_NVARS+sp) = mflo_user_funcs::spec_fs[sp];
                }


                phi(i,j,k,RHO_INDX)  = phi(i,j,k,DENS_INDX);
                phi(i,j,k,RHOU_INDX) = phi(i,j,k,DENS_INDX)*phi(i,j,k,VELX_INDX);
                phi(i,j,k,RHOV_INDX) = phi(i,j,k,DENS_INDX)*phi(i,j,k,VELY_INDX);
                phi(i,j,k,RHOW_INDX) = phi(i,j,k,DENS_INDX)*phi(i,j,k,VELZ_INDX);
                phi(i,j,k,RHOE_INDX)  = phi(i,j,k,DENS_INDX)*
                    mflo_thermo::get_e_from_rpc(mflo_user_funcs::dens_fs,mflo_user_funcs::pres_fs,mflo_user_funcs::spec_fs);
                phi(i,j,k,RHOE_INDX) += half*phi(i,j,k,DENS_INDX)*
                                        ( phi(i,j,k,VELX_INDX)*phi(i,j,k,VELX_INDX)
                                        + phi(i,j,k,VELY_INDX)*phi(i,j,k,VELY_INDX)
                                        + phi(i,j,k,VELZ_INDX)*phi(i,j,k,VELZ_INDX) );
            }
        }
    }
}
#endif

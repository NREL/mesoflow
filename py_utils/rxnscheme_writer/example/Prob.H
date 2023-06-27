#ifndef _PROB_H_
#define _PROB_H_

#include <AMReX_Box.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_Geometry.H>
#include <userfuncs.H>
#include <globalDefines.H>
#include <cflo_constants.H>

using namespace amrex;

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void initdata(Box const& bx, Array4<Real> const& phi, GeometryData const& geomdata)
{
    int i_mrc,j_mrc,k_mrc;
    const auto lo = lbound(bx);
    const auto hi = ubound(bx);   
    Real zerospec[NUM_SPECIES]={zero}; 
    Real rho_bg=cflo_thermo::get_r_from_tpc(cflo_user_funcs::fs_temp,
            cflo_user_funcs::fs_p,zerospec);

    const Real* AMREX_RESTRICT prob_lo = geomdata.ProbLo();
    const Real* AMREX_RESTRICT prob_hi = geomdata.ProbHi();
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
                
                //defaults
                phi(i,j,k,PRES_INDX) = cflo_user_funcs::fs_p;
                phi(i,j,k,DENS_INDX) = cflo_user_funcs::fs_rho;
                phi(i,j,k,VELX_INDX) = cflo_user_funcs::fs_vel;
                phi(i,j,k,VELY_INDX) = zero;
                phi(i,j,k,VELZ_INDX) = zero;
                phi(i,j,k,VFRAC_INDX) = one;
                phi(i,j,k,TEMP_INDX)  = cflo_user_funcs::fs_temp;

                for(int sp=0;sp<NUM_SPECIES;sp++)
                {
                    phi(i,j,k,CFLO_NVARS+sp) = zerospec[sp];
                }
                
                phi(i,j,k,RHOE_INDX) = phi(i,j,k,DENS_INDX) *
                    cflo_thermo::get_e_from_rpc(phi(i,j,k,DENS_INDX),
                            phi(i,j,k,PRES_INDX),
                            zerospec);
                            //cflo_user_funcs::fs_spec);
            
                i_mrc = floor((x-cflo_user_funcs::lo_mrc_x)/cflo_user_funcs::dx_mrc_x);
                j_mrc = floor((y-cflo_user_funcs::lo_mrc_y)/cflo_user_funcs::dx_mrc_y);
                k_mrc = floor((z-cflo_user_funcs::lo_mrc_z)/cflo_user_funcs::dx_mrc_z);
            
                if( (i_mrc >= 0) and (i_mrc < cflo_user_funcs::nx_mrc)
                and (j_mrc >= 0) and (j_mrc < cflo_user_funcs::ny_mrc) 
                and (k_mrc >= 0) and (k_mrc < cflo_user_funcs::nz_mrc)  )
                {
                    int index=k_mrc*(cflo_user_funcs::nx_mrc*cflo_user_funcs::ny_mrc)
                              + j_mrc*cflo_user_funcs::nx_mrc + i_mrc;
                   
                   //zero inside box and let it evolve
                   phi(i,j,k,VELX_INDX) = zero;
                   phi(i,j,k,TEMP_INDX)  = cflo_user_funcs::fs_temp;
                   for(int sp=0;sp<NUM_SPECIES;sp++)
                   {
                       phi(i,j,k,CFLO_NVARS+sp) = zerospec[sp];
                   }
                   phi(i,j,k,DENS_INDX) = rho_bg;
                
                   phi(i,j,k,RHOE_INDX) = phi(i,j,k,DENS_INDX) *
                    cflo_thermo::get_e_from_rpc(phi(i,j,k,DENS_INDX),
                            phi(i,j,k,PRES_INDX),
                            zerospec);

                   if(cflo_user_funcs::mrcdata[index] > cflo_user_funcs::mrc_threshold)
                   { 
                       phi(i,j,k,VFRAC_INDX) = zero;
                       phi(i,j,k,CFLO_NVARS+S1_ID) = cflo_user_funcs::catalyst_sites;
                   }
                   else
                   {
                       phi(i,j,k,VFRAC_INDX) = one;
                       phi(i,j,k,CFLO_NVARS+S1_ID) = zero;
                   }
                }
                else
                {
                    phi(i,j,k,VELX_INDX) = cflo_user_funcs::fs_vel;
                    phi(i,j,k,VFRAC_INDX) = one;
                    phi(i,j,k,CFLO_NVARS+S1_ID) = zero;
                }


                phi(i,j,k,RHO_INDX)  = phi(i,j,k,DENS_INDX);
                phi(i,j,k,RHOU_INDX) = phi(i,j,k,DENS_INDX)*phi(i,j,k,VELX_INDX);
                phi(i,j,k,RHOV_INDX) = phi(i,j,k,DENS_INDX)*phi(i,j,k,VELY_INDX);
                phi(i,j,k,RHOW_INDX) = phi(i,j,k,DENS_INDX)*phi(i,j,k,VELZ_INDX);

                phi(i,j,k,RHOE_INDX) +=  half*phi(i,j,k,DENS_INDX)  *
                    ( phi(i,j,k,VELX_INDX)*phi(i,j,k,VELX_INDX)
                      + phi(i,j,k,VELY_INDX)*phi(i,j,k,VELY_INDX)
                      + phi(i,j,k,VELZ_INDX)*phi(i,j,k,VELZ_INDX) );

            }
        }
    }
}
#endif

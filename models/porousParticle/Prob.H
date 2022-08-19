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

    bool inside;
    int i_mrc,j_mrc,k_mrc,index;

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

                    phi(i,j,k,PRES_INDX) = mflo_user_funcs::pres_left; 
                    phi(i,j,k,DENS_INDX) = mflo_user_funcs::dens_left;
                    phi(i,j,k,TEMP_INDX) = mflo_user_funcs::temp_left;
                    phi(i,j,k,VELX_INDX) =  zeroval;
                    phi(i,j,k,VELY_INDX) =  zeroval;
                    phi(i,j,k,VELZ_INDX) =  zeroval;


                    Real vfrac=zeroval;
                    for(int kk=0;kk<2;kk++)
                    {
                        for(int jj=0;jj<2;jj++)
                        {
                            for(int ii=0;ii<2;ii++)
                            {
                                Real xx[3]={prob_lo[0]+(i+ii)*dx[0],
                                            prob_lo[1]+(j+jj)*dx[1],           
                                            prob_lo[2]+(k+kk)*dx[2]};

		                i_mrc = floor((xx[0]-mflo_user_funcs::lo_mrc_x)/mflo_user_funcs::dx_mrc_x);
		                j_mrc = floor((xx[1]-mflo_user_funcs::lo_mrc_y)/mflo_user_funcs::dx_mrc_y);
		                k_mrc = floor((xx[2]-mflo_user_funcs::lo_mrc_z)/mflo_user_funcs::dx_mrc_z);

		                 if( (i_mrc >= 0) && (i_mrc < mflo_user_funcs::nx_mrc)
			         && (j_mrc >= 0)  && (j_mrc < mflo_user_funcs::ny_mrc)
			         && (k_mrc >= 0)  && (k_mrc < mflo_user_funcs::nz_mrc)  )
		                 {
			            int index=k_mrc*(mflo_user_funcs::nx_mrc*mflo_user_funcs::ny_mrc)
			            + j_mrc*mflo_user_funcs::nx_mrc + i_mrc;
			
                                    if(mflo_user_funcs::mrcdata[index] < mflo_user_funcs::mrc_threshold)
			            {
		                        vfrac+=one;	                
			            }
		                }
                             }
                         }
                     }

                    vfrac*=eighth;
                    phi(i,j,k,VFRAC_INDX) = vfrac;

                    for(int sp=0;sp<NUM_SPECIES;sp++)
                    {
                        phi(i,j,k,FLO_NVARS+sp) = mflo_user_funcs::spec_left[sp];
                    }


                    phi(i,j,k,RHO_INDX)  = phi(i,j,k,DENS_INDX);
                    phi(i,j,k,RHOU_INDX) = phi(i,j,k,DENS_INDX)*phi(i,j,k,VELX_INDX);
                    phi(i,j,k,RHOV_INDX) = phi(i,j,k,DENS_INDX)*phi(i,j,k,VELY_INDX);
                    phi(i,j,k,RHOW_INDX) = phi(i,j,k,DENS_INDX)*phi(i,j,k,VELZ_INDX);
                    phi(i,j,k,RHOE_INDX)  = phi(i,j,k,DENS_INDX)*
                        mflo_thermo::get_e_from_rpc(mflo_user_funcs::dens_left,mflo_user_funcs::pres_left,mflo_user_funcs::spec_left);
                    phi(i,j,k,RHOE_INDX) += half*phi(i,j,k,DENS_INDX)*
                        ( phi(i,j,k,VELX_INDX)*phi(i,j,k,VELX_INDX)
                        + phi(i,j,k,VELY_INDX)*phi(i,j,k,VELY_INDX)
                        + phi(i,j,k,VELZ_INDX)*phi(i,j,k,VELZ_INDX) );
                }
        }
    }
}
#endif

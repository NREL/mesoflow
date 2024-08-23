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
    Real spec[NUM_SPECIES]={zeroval};
    Real dummyspec[NUM_SPECIES]={zeroval};
        
    Real mwbar = mflo_user_funcs::H2molfrac*mflo_species::molwts[H2_ID] +   //-v effective molecular weight (H, H2, Ar) applied
                     mflo_user_funcs::Hmolfrac*mflo_species::molwts[H_ID] +                //everywhere in the domain
                     (1.0-mflo_user_funcs::H2molfrac-mflo_user_funcs::Hmolfrac)
                    *mflo_species::molwts[AR_ID];

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

                for(int var=0;var<TOTAL_NVARS;var++)
                {
                    phi(i,j,k,var)=0.0;
                }
                
                Real vfrac=0.0;
                for(int kk=0;kk<2;kk++)
                {
                   for(int jj=0;jj<2;jj++)
                   {
                        for(int ii=0;ii<2;ii++)
                        {
                            Real x=prob_lo[0]+(i+ii)*dx[0];             
                            Real y=prob_lo[1]+(j+jj)*dx[1];             
                            Real z=prob_lo[2]+(k+kk)*dx[2];

                            Real yfrac=(y-prob_lo[1])/(prob_hi[1]-prob_lo[1]);

                            Real interfaceloc=(prob_hi[XDIR]-mflo_user_funcs::catblocksize);
                                                
                
                            if(x < interfaceloc)
                            {
                                vfrac+=one;
                            }
                        }
                   }
                }

                vfrac*=eighth;
                phi(i,j,k,VFRAC_INDX) = vfrac;
                phi(i,j,k,VELX_INDX) =  zeroval;
                phi(i,j,k,VELY_INDX) =  zeroval;
                phi(i,j,k,VELZ_INDX) =  zeroval;

                phi(i,j,k,PRES_INDX) = mflo_user_funcs::p0;
                if(vfrac>zeroval) 
                {
				
					//phi(i,j,k,TEMP_INDX) = ( (y*y + z*z) < (mflo_user_funcs::elec_rad*mflo_user_funcs::elec_rad) && x>0.005 && x<0.0125 ) ? mflo_user_funcs::T_init:
																												  //mflo_user_funcs::T_surr;
					phi(i,j,k,TEMP_INDX) = ( (y*y + z*z) < (mflo_user_funcs::elec_rad*mflo_user_funcs::elec_rad)  ) ? mflo_user_funcs::T_ign:
																												  mflo_user_funcs::T_surr;																							  
                }
                else
                {
                    phi(i,j,k,TEMP_INDX)=mflo_user_funcs::T_solid;
                }

                //for solid, ideally should change
                phi(i,j,k,DENS_INDX) = phi(i,j,k,PRES_INDX)/RU/
                phi(i,j,k,TEMP_INDX)*mwbar;

                Real totalconc=phi(i,j,k,DENS_INDX)/mwbar;

                spec[H_ID]   = vfrac*totalconc*mflo_user_funcs::Hmolfrac;
                spec[H2_ID]  = vfrac*totalconc*mflo_user_funcs::H2molfrac;
                spec[H2O_ID] = zeroval;

                //pure Ar in solid, strictly mwbar should be mw_argon
                //but it wouldnt matter as we dont solve any fluid equations
                spec[AR_ID]  = vfrac*totalconc*(1.0 - mflo_user_funcs::H2molfrac                        
                                                - mflo_user_funcs::Hmolfrac);
                spec[FEO_ID] = zeroval;
                spec[FE2O3_ID] = zeroval; 
                spec[FE_ID]  = zeroval;

                phi(i,j,k,RHO_INDX)  = phi(i,j,k,DENS_INDX);
                phi(i,j,k,RHOU_INDX) = phi(i,j,k,DENS_INDX)*phi(i,j,k,VELX_INDX);
                phi(i,j,k,RHOV_INDX) = phi(i,j,k,DENS_INDX)*phi(i,j,k,VELY_INDX);
                phi(i,j,k,RHOW_INDX) = phi(i,j,k,DENS_INDX)*phi(i,j,k,VELZ_INDX);

                if(vfrac>zeroval)   
                {                                                      
                    phi(i,j,k,RHOE_INDX)  = phi(i,j,k,DENS_INDX)*
                    mflo_thermo::get_e_from_rpc(phi(i,j,k,DENS_INDX),phi(i,j,k,PRES_INDX),spec);

                    phi(i,j,k,RHOE_INDX) += half*phi(i,j,k,DENS_INDX)*
                    ( phi(i,j,k,VELX_INDX)*phi(i,j,k,VELX_INDX)
                     + phi(i,j,k,VELY_INDX)*phi(i,j,k,VELY_INDX)
                     + phi(i,j,k,VELZ_INDX)*phi(i,j,k,VELZ_INDX) );
                }
                else
                {
                    phi(i,j,k,RHOE_INDX)  = mflo_thermo::get_solid_rhoe_from_t(phi(i,j,k,TEMP_INDX));
                }

                for(int sp=0;sp<NUM_SPECIES;sp++)
                {
                    phi(i,j,k,FLO_NVARS+sp) = spec[sp];
                }
                phi(i,j,k,SIG_INDX) = 0.1;
            }
        }
    }
}
#endif

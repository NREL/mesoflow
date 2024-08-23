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
    
    Real spec[NUM_SPECIES]={zeroval};
    Real dummyspec[NUM_SPECIES]={zeroval};
    
    Real mwbar = mflo_user_funcs::H2molfrac*mflo_species::molwts[H2_ID] +
        mflo_user_funcs::H2Omolfrac*mflo_species::molwts[H2O_ID] +  mflo_user_funcs::COmolfrac*mflo_species::molwts[CO_ID] +  mflo_user_funcs::N2molfrac*mflo_species::molwts[N2_ID];
    

#ifdef _OPENMP
#pragma omp parallel for collapse(2) if (GPU::notInLaunchRegion)
#endif
    
    //Looping Over Grid Points
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
                    
                    //Initialization of Physical Variables
                    phi(i,j,k,PRES_INDX) = mflo_user_funcs::pres_left; 
                    phi(i,j,k,DENS_INDX) = mflo_user_funcs::dens_left;
                    phi(i,j,k,TEMP_INDX) = mflo_user_funcs::temp_left;
                    phi(i,j,k,VELX_INDX) =  zeroval;
                    phi(i,j,k,VELY_INDX) =  mflo_user_funcs::fs_vel;
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
                        //
		                i_mrc = floor((xx[0]-mflo_user_funcs::lo_mrc_x)/mflo_user_funcs::dx_mrc_x);
		                j_mrc = floor((xx[1]-mflo_user_funcs::lo_mrc_y)/mflo_user_funcs::dx_mrc_y);
		                k_mrc = floor((xx[2]-mflo_user_funcs::lo_mrc_z)/mflo_user_funcs::dx_mrc_z);

                     //
		                 if( (i_mrc >= 0) && (i_mrc < mflo_user_funcs::nx_mrc)
			         && (j_mrc >= 0)  && (j_mrc < mflo_user_funcs::ny_mrc)
			         && (k_mrc >= 0)  && (k_mrc < mflo_user_funcs::nz_mrc)  )
		                 {
			            int index=k_mrc*(mflo_user_funcs::nx_mrc*mflo_user_funcs::ny_mrc)
			            + j_mrc*mflo_user_funcs::nx_mrc + i_mrc;
                    //Threshold
                                    if(mflo_user_funcs::mrcdata[index] < mflo_user_funcs::mrc_threshold)
			            {
		                        vfrac+=one;	                
			            }
		                }
                             }
                         }
                     }
                    //Final Volume Fraction Calculation, turn counting into fraction
                    vfrac*=eighth;
                    phi(i,j,k,VFRAC_INDX) = vfrac;
                    phi(i,j,k,VELX_INDX) =  zeroval;
                    phi(i,j,k,VELY_INDX) =  mflo_user_funcs::fs_vel;
                    phi(i,j,k,VELZ_INDX) =  zeroval;

                  //Additional Initialization (Species and Momentum)
                // rho; rhou;rhov;rhow ;rho e = (internal e +1/2v^2)
             
                    
                    phi(i,j,k,PRES_INDX) = mflo_user_funcs::pres_left;
                                    if(vfrac>zeroval)
                                    {
                                        phi(i,j,k,TEMP_INDX) = mflo_user_funcs::temp_left;
                                    }
                                    else
                                    {
                                        phi(i,j,k,TEMP_INDX)=mflo_user_funcs::temp_left;
                                    }

                                    //for solid, ideally should change
                                    phi(i,j,k,DENS_INDX) = phi(i,j,k,PRES_INDX)/RU/phi(i,j,k,TEMP_INDX)*mwbar;

                                    Real totalconc=phi(i,j,k,DENS_INDX)/mwbar;

//                                    spec[CO_ID]   = vfrac*totalconc*mflo_user_funcs::Hmolfrac;
//                                    spec[H2_ID]  = vfrac*totalconc*mflo_user_funcs::H2molfrac;
//                                    spec[H2O_ID] = zeroval;
//                                    spec[H2_ID]  = vfrac*totalconc*mflo_user_funcs::H2molfrac;
//                                    //pure Ar in solid, strictly mwbar should be mw_argon
//                                    //but it wouldnt matter as we dont solve any fluid equations
//                                    spec[AR_ID]  = vfrac*totalconc*(1.0 - mflo_user_funcs::H2molfrac                         //-v AR species concentration
//                                                                    - mflo_user_funcs::Hmolfrac)+(1-vfrac)*totalconc;
                    
                    spec[H2_ID] = vfrac*totalconc*mflo_user_funcs::H2molfrac;
                    spec[H2O_ID] = vfrac*totalconc*mflo_user_funcs::H2Omolfrac;
                    spec[CO_ID] = vfrac*totalconc*mflo_user_funcs::COmolfrac;
                    spec[N2_ID] = vfrac*totalconc*(1.0 - mflo_user_funcs::H2molfrac
                                                   - mflo_user_funcs::H2Omolfrac- mflo_user_funcs::COmolfrac)+(1-vfrac)*totalconc;
   
                    
                                    spec[FeO_ID] = zeroval;
                                    spec[Fe2O3_ID] = (1.0-vfrac)*mflo_user_funcs::siteconc;
                                    spec[Fe_ID]  = zeroval;
                    
                    
                    


                    phi(i,j,k,RHO_INDX)  = phi(i,j,k,DENS_INDX);
                    phi(i,j,k,RHOU_INDX) = phi(i,j,k,DENS_INDX)*phi(i,j,k,VELX_INDX);
                    phi(i,j,k,RHOV_INDX) = phi(i,j,k,DENS_INDX)*phi(i,j,k,VELY_INDX);
                    phi(i,j,k,RHOW_INDX) = phi(i,j,k,DENS_INDX)*phi(i,j,k,VELZ_INDX);
//                    phi(i,j,k,RHOE_INDX)  = phi(i,j,k,DENS_INDX)*
//                        mflo_thermo::get_e_from_rpc(mflo_user_funcs::dens_left,mflo_user_funcs::pres_left,mflo_user_funcs::spec);
//                    phi(i,j,k,RHOE_INDX) += half*phi(i,j,k,DENS_INDX)*
//                        ( phi(i,j,k,VELX_INDX)*phi(i,j,k,VELX_INDX)
//                        + phi(i,j,k,VELY_INDX)*phi(i,j,k,VELY_INDX)
//                        + phi(i,j,k,VELZ_INDX)*phi(i,j,k,VELZ_INDX) );
                    
                    
                    
                    if(vfrac>zeroval)   //-v average conserved quantities for the cell are obtained since the conc of spices
                                {                                                      //multiplied by volfr
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
                    
                }
        }
    }
}
#endif

#ifndef _EXTERNBC_H_
#define _EXTERNBC_H_

#include <AMReX_Box.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_Geometry.H>
#include <userfuncs.H>
#include <globalDefines.H>
#include <mflo_constants.H>

using namespace amrex;

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void externalbc(const amrex::Real x[AMREX_SPACEDIM],
        const amrex::Real s_int[],
        amrex::Real s_ext[],
        const int idir,
        const int sgn,
        const amrex::Real time,
        amrex::GeometryData const& geomdata)
{
    
    for(int n=0;n<TOTAL_NVARS;n++)
    {
        s_ext[n]=s_int[n];
    }
    
    
 
    if(idir==1)// Check if the boundary is in the y-direction
    {
        if(sgn > 0.0) // Positive direction (e.g., top boundary)
        {
            Real spec[NUM_SPECIES]={zeroval};
            
       
            s_ext[PRES_INDX] = mflo_user_funcs::pres_left;
            s_ext[TEMP_INDX] = mflo_user_funcs::temp_left;
            s_ext[VELX_INDX] = zeroval;
            s_ext[VELY_INDX] = mflo_user_funcs::fs_vel;
            s_ext[VELZ_INDX] = zeroval;
            
            
            Real mwbar = mflo_user_funcs::H2molfrac*mflo_species::molwts[H2_ID] +
                mflo_user_funcs::H2Omolfrac*mflo_species::molwts[H2O_ID] +  mflo_user_funcs::COmolfrac*mflo_species::molwts[CO_ID] +  mflo_user_funcs::N2molfrac*mflo_species::molwts[N2_ID]+  mflo_user_funcs::CO2molfrac*mflo_species::molwts[CO2_ID];
            
            
            //C_total= # of moles of gas/ Volume = P / RT
            Real totalconc=s_ext[PRES_INDX]/RU/s_ext[TEMP_INDX];
            
            // define T, solve rho
            s_ext[DENS_INDX] =  totalconc*mwbar;
            
            spec[H2_ID] = totalconc*mflo_user_funcs::H2molfrac;
            spec[N2_ID] = totalconc*mflo_user_funcs::N2molfrac;
            spec[H2O_ID] = totalconc*mflo_user_funcs::H2Omolfrac;
            spec[CO_ID] = totalconc*mflo_user_funcs::COmolfrac;
            spec[CO2_ID] = totalconc*mflo_user_funcs::CO2molfrac;
            
            s_ext[RHO_INDX]  = s_ext[DENS_INDX];
            s_ext[RHOU_INDX] = s_ext[DENS_INDX]*s_ext[VELX_INDX];
            s_ext[RHOV_INDX] = s_ext[DENS_INDX]*s_ext[VELY_INDX];
            s_ext[RHOW_INDX] = s_ext[DENS_INDX]*s_ext[VELZ_INDX];
            s_ext[RHOE_INDX]  = s_ext[DENS_INDX]*
            mflo_thermo::get_e_from_rpc(s_ext[DENS_INDX],s_ext[PRES_INDX],spec);

            s_ext[RHOE_INDX] += half*s_ext[DENS_INDX]*
            (s_ext[VELX_INDX]*s_ext[VELX_INDX] +
             s_ext[VELY_INDX]*s_ext[VELY_INDX] +
             s_ext[VELZ_INDX]*s_ext[VELZ_INDX]);
            
            

            for(int sp=0;sp<NUM_SPECIES;sp++)
                        {
                            s_ext[FLO_NVARS+sp]=spec[sp];
                        }
                        s_ext[VFRAC_INDX] = s_int[VFRAC_INDX];
            
            
            
        }
        else // Negative direction (e.g., bottom boundary)
        {
            
//            Real spec[NUM_SPECIES]={zeroval};
//
//           // s_ext[DENS_INDX] = s_int[DENS_INDX]; // 0 grad
//            s_ext[PRES_INDX] = s_int[DENS_INDX];
//            s_ext[VELX_INDX] = s_int[VELX_INDX];
//            s_ext[VELY_INDX] = s_int[VELY_INDX];
//            s_ext[VELZ_INDX] = s_int[VELZ_INDX];
//
//
//
//            Real mwbar = mflo_user_funcs::H2molfrac*mflo_species::molwts[H2_ID] +
//                mflo_user_funcs::H2Omolfrac*mflo_species::molwts[H2O_ID] +  mflo_user_funcs::COmolfrac*mflo_species::molwts[CO_ID] +  mflo_user_funcs::N2molfrac*mflo_species::molwts[N2_ID];
//
//            Real totalconc=s_ext[PRES_INDX]/RU/s_ext[TEMP_INDX];
//
//            // define T, solve rho
//            s_ext[DENS_INDX] =  totalconc*mwbar;
//
//            spec[H2_ID] = totalconc*mflo_user_funcs::H2molfrac;
//            spec[N2_ID] = totalconc*mflo_user_funcs::N2molfrac;
//            spec[H2O_ID] = totalconc*mflo_user_funcs::H2Omolfrac;
//            spec[CO_ID] = totalconc*mflo_user_funcs::COmolfrac;
//
//            s_ext[RHO_INDX]  = s_ext[DENS_INDX];
//                        s_ext[RHOU_INDX] = s_ext[DENS_INDX]*s_ext[VELX_INDX];
//                        s_ext[RHOV_INDX] = s_ext[DENS_INDX]*s_ext[VELY_INDX];
//                        s_ext[RHOW_INDX] = s_ext[DENS_INDX]*s_ext[VELZ_INDX];
//                        s_ext[RHOE_INDX]  = s_ext[DENS_INDX]*
//                        mflo_thermo::get_e_from_rpc(s_ext[DENS_INDX],s_ext[PRES_INDX],spec);
//
//            s_ext[RHOE_INDX] += half*s_ext[DENS_INDX]*
//                        (s_ext[VELX_INDX]*s_ext[VELX_INDX] +
//                         s_ext[VELY_INDX]*s_ext[VELY_INDX] +
//                         s_ext[VELZ_INDX]*s_ext[VELZ_INDX]);
//

            
            
            
            s_ext[PRES_INDX] = s_int[PRES_INDX];
            s_ext[DENS_INDX] = s_int[DENS_INDX];
            s_ext[TEMP_INDX] = s_int[TEMP_INDX];
            s_ext[VELX_INDX] = s_int[VELX_INDX];
            s_ext[VELY_INDX] = s_int[VELY_INDX];
            s_ext[VELZ_INDX] = s_int[VELZ_INDX];

            s_ext[RHO_INDX]  = s_int[RHO_INDX];
            s_ext[RHOU_INDX] = s_int[RHOU_INDX];
            s_ext[RHOV_INDX] = s_int[RHOV_INDX];
            s_ext[RHOW_INDX] = s_int[RHOW_INDX];
            s_ext[RHOE_INDX] = s_int[RHOE_INDX];
            s_ext[VFRAC_INDX] = s_int[VFRAC_INDX];
          //  s_ext[TEMP_INDX]=mflo_user_funcs::Tsolid;
            
//            for(int sp=0;sp<NUM_SPECIES;sp++)
//                        {
//                            s_ext[FLO_NVARS+sp]=spec[sp];
//                        }
                       
            
            
            
            
        }


    }
}
#endif

#ifndef _CHEMISTRY_H_
#define _CHEMISTRY_H_

#include <AMReX_Geometry.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_REAL.H>
#include <species.H>
#include <globalDefines.H>

using namespace amrex;
namespace mflo_chem_reactions
{
    AMREX_GPU_DEVICE AMREX_INLINE
        void compute_spec_source(int i, int j, int k,
                Array4<Real> const& phi,
                Array4<Real> const& specsource,
                GpuArray<Real, AMREX_SPACEDIM> prob_lo,
                GpuArray<Real, AMREX_SPACEDIM> dx,
                const Real time)
        {
            for(int sp=0;sp<NUM_SPECIES;sp++)
            {
                specsource(i,j,k,FLO_NVARS+sp) = zeroval;
            }
            amrex::Real kreac,reactrate;
            amrex::Real rg=1e-5; //grain radius
            amrex::Real vfrac=phi(i,j,k,VFRAC_INDX);
            amrex::Real twiddle=1e-10;
            amrex::Real totalconc = (1.0-vfrac)*mflo_user_funcs::siteconc+twiddle;
            
            //FeO + H2 -> Fe + H2O 
            //fix this to include area/vol scaling
            kreac=2857.34*std::exp(-14100.0/phi(i,j,k,TEMP_INDX))*(3.0/rg);
            reactrate=kreac*phi(i,j,k,FLO_NVARS+FEO_ID)/totalconc*phi(i,j,k,FLO_NVARS+H2_ID);
            specsource(i,j,k,FLO_NVARS+H2_ID)-=reactrate;
            specsource(i,j,k,FLO_NVARS+FEO_ID)-=reactrate;
            specsource(i,j,k,FLO_NVARS+FE_ID)+=reactrate;
            specsource(i,j,k,FLO_NVARS+H2O_ID)+=reactrate;
            
            //Fe2O3 + H2 -> 2FeO + H2O 
            kreac=80.0*std::exp(-8000.0/phi(i,j,k,TEMP_INDX))*(3.0/rg);
            reactrate=kreac*phi(i,j,k,FLO_NVARS+FE2O3_ID)/totalconc*phi(i,j,k,FLO_NVARS+H2_ID);
            specsource(i,j,k,FLO_NVARS+H2_ID)-=reactrate;
            specsource(i,j,k,FLO_NVARS+FE2O3_ID)-=reactrate;
            specsource(i,j,k,FLO_NVARS+FEO_ID)+=2.0*reactrate;
            specsource(i,j,k,FLO_NVARS+H2O_ID)+=reactrate;

            //Real gasfrac=phi(i,j,k,VFRAC_INDX);
            //Real solfrac=one-gasfrac;
            //specsource(i,j,k,FLO_NVARS+H2O_ID) += 50.0*specsource(i,j,k,FLO_NVARS+H_ID)*solfrac;
        }
}
#endif

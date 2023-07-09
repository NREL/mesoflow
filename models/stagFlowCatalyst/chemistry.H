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
            //2CH4 + O2 + S -> C2H4 + 2H2O + S
            for(int sp=0;sp<NUM_SPECIES;sp++)
            {
                specsource(i,j,k,FLO_NVARS+sp) = zeroval;
            }
            amrex::Real kreac=5.0;
            amrex::Real reactrate=kreac*phi(i,j,k,FLO_NVARS+S_ID)*phi(i,j,k,FLO_NVARS+CH4_ID)*phi(i,j,k,FLO_NVARS+O2_ID);
            specsource(i,j,k,FLO_NVARS+CH4_ID)-=2.0*reactrate;
            specsource(i,j,k,FLO_NVARS+O2_ID)-=reactrate;
            specsource(i,j,k,FLO_NVARS+C2H4_ID)+=reactrate;
            specsource(i,j,k,FLO_NVARS+H2O_ID)+=2.0*reactrate;
        }
}
#endif

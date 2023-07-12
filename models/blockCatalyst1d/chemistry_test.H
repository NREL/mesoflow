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
            amrex::Real kreac=5.0;
            amrex::Real siteconc=phi(i,j,k,FLO_NVARS+SITE1_ID);

            specsource(i,j,k,FLO_NVARS+S1_ID) 
            = -kreac*siteconc*phi(i,j,k,FLO_NVARS+S1_ID);

            specsource(i,j,k,FLO_NVARS+S2_ID) = kreac*siteconc*phi(i,j,k,FLO_NVARS+S1_ID);
            specsource(i,j,k,FLO_NVARS+SITE1_ID) = zeroval;
        }
}
#endif

#ifndef _TRANSPORT_H_
#define _TRANSPORT_H_

#include <AMReX_Geometry.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_REAL.H>
#include <species.H>

using namespace amrex;
namespace mflo_chem_transport
{
    AMREX_GPU_DEVICE AMREX_INLINE
        void compute_spec_dcoeff(int i, int j, int k,
                Array4<Real> const& phi,
                Array4<Real> const& dcoeff,
                GpuArray<Real, AMREX_SPACEDIM> prob_lo,
                GpuArray<Real, AMREX_SPACEDIM> dx,
                const Real time)
        {
            Real solid_dif=1.e-10;
            Real fluid_dif=4.27e-5;
            Real vfrac=phi(i,j,k,VFRAC_INDX);

            for(int sp=0;sp<NUM_SPECIES;sp++)
            {
                dcoeff(i,j,k,sp) = (one-vfrac)*solid_dif+vfrac*fluid_dif;
            }
            dcoeff(i,j,k,S1_ID)=zeroval;
       }
}
#endif
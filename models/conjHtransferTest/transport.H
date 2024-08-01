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
            dcoeff(i,j,k,AIR_ID)=1e-6;
            dcoeff(i,j,k,SOLSPEC_ID)=zeroval;
        }
}
#endif

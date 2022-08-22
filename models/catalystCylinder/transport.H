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
            Real gasfrac=phi(i,j,k,VFRAC_INDX);
            Real solfrac=one-gasfrac;
           
            dcoeff(i,j,k,AIR_ID)=1e-5;
            dcoeff(i,j,k,A_ID)=1e-5;
            dcoeff(i,j,k,B_ID)=1e-5;
            dcoeff(i,j,k,S_ID)=0.0;
        }
}
#endif

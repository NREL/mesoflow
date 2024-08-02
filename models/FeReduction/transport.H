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
            Real lscale=10*dx[1];
            Real y=prob_lo[1]+(j+0.5)*dx[1];
            Real dscale=5.0;
            Real eps=0.5;
            Real tort=5.0;

            Real dcoeff_H_gas=3e-4;
            Real dcoeff_H2_gas=1e-4;
            Real dcoeff_H2O_gas=1e-4;
            Real dcoeff_AR_gas=1e-4;
            
            Real dcoeff_H_solid=dcoeff_H_gas*eps/tort;
            Real dcoeff_H2_solid=dcoeff_H2_gas*eps/tort;
            Real dcoeff_H2O_solid=dcoeff_H2O_gas*eps/tort;
            Real dcoeff_AR_solid=dcoeff_AR_gas*eps/tort;

            Real sponge_fac=0.5*(1+std::tanh((amrex::Math::abs(y)-mflo_user_funcs::sponge_zone_dist)/lscale));
            
            dcoeff(i,j,k,H_ID)=gasfrac*dcoeff_H_gas*(1.0+sponge_fac*dscale)  
            + solfrac*dcoeff_H_solid;
            dcoeff(i,j,k,H2_ID)=gasfrac*dcoeff_H2_gas*(1.0+sponge_fac*dscale) 
            + solfrac*dcoeff_H2_solid;
            dcoeff(i,j,k,H2O_ID)=gasfrac*dcoeff_H2O_gas*(1.0+sponge_fac*dscale)
            + solfrac*dcoeff_H2O_solid;
            dcoeff(i,j,k,AR_ID)=gasfrac*dcoeff_AR_gas*(1.0+sponge_fac*dscale)
            + solfrac*dcoeff_AR_solid;

            dcoeff(i,j,k,FEO_ID)=solfrac*1e-12;
            dcoeff(i,j,k,FE2O3_ID)=solfrac*1e-12;
            dcoeff(i,j,k,FE_ID)=solfrac*1e-12;
        }
}
#endif

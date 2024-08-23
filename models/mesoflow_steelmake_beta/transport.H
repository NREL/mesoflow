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

            // Define diffusion coefficients for the reaction species (in m^2/s)
            
            
            Real gasfrac=phi(i,j,k,VFRAC_INDX);
            Real solfrac=one-gasfrac;
            Real lscale=10*dx[1];
            Real y=prob_lo[1]+(j+0.5)*dx[1];
            Real dscale=5.0;
            Real eps=0.5;
            Real tort=5.0;
            
            
            Real D_H2_gas    = 7.0e-5;  //
            Real D_CO_gas    = 1.6e-5;  //
            Real D_H2O_gas   = 2.5e-5;  //
            Real D_CO2_gas   = 1.39e-5;
            Real D_N2_gas   = 1.88e-5;
            
            Real D_Fe2O3 = 1.0e-9;
            Real D_FeO   = 1.0e-9;
            Real D_Fe    = 1.0e-9;
            
            Real D_H2_solid=D_H2_gas*eps/tort;
            Real D_CO_solid=D_CO_gas*eps/tort;
            Real D_H2O_solid=D_H2O_gas*eps/tort;
            Real D_CO2_solid=D_CO2_gas*eps/tort;
            Real D_N2_solid=D_N2_gas*eps/tort;
            
            
            Real sponge_fac=0.5*(1+std::tanh((amrex::Math::abs(y)-mflo_user_funcs::sponge_zone_dist)/lscale));
                        
                        dcoeff(i,j,k,H2_ID)=gasfrac*D_H2_gas*(1.0+sponge_fac*dscale)
                        + solfrac*D_H2_solid;
                        dcoeff(i,j,k,H2O_ID)=gasfrac*D_H2O_gas*(1.0+sponge_fac*dscale)
                        + solfrac*D_H2O_solid;
                        dcoeff(i,j,k,CO_ID)=gasfrac*D_CO_gas*(1.0+sponge_fac*dscale)
                        + solfrac*D_CO_solid;
                        dcoeff(i,j,k,CO2_ID)=gasfrac*D_CO2_gas*(1.0+sponge_fac*dscale)
                        + solfrac*D_CO2_solid;
                        dcoeff(i,j,k,N2_ID)=gasfrac*D_N2_gas*(1.0+sponge_fac*dscale)
                        + solfrac*D_N2_solid;

            // Assign diffusion coefficients to the appropriate species

            dcoeff(i,j,k,Fe2O3_ID) = solfrac*D_Fe2O3;
            dcoeff(i,j,k,FeO_ID)   = solfrac*D_FeO;
            dcoeff(i,j,k,Fe_ID)    = solfrac*D_Fe;
            
            
        }
}
#endif

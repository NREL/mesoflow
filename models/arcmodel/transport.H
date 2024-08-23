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
            Real eps=mflo_user_funcs::eps;
            Real tort=mflo_user_funcs::tort;

            Real thcond = (phi(i,j,k,FLO_NVARS+AR_ID)*thcond_AR_inter.interpolate(phi(i,j,k,TEMP_INDX))
                    + phi(i,j,k,FLO_NVARS+H2_ID)*thcond_H2_inter.interpolate(phi(i,j,k,TEMP_INDX)))
                / (phi(i,j,k,FLO_NVARS+AR_ID) + phi(i,j,k,FLO_NVARS+H2_ID));

            //FIX THIS, need a cp function thermo
            //assuming mostly argon and using a hack g R/(g-1)
            Real Cp=1.6*RU/0.6/0.04;
            Real thdiff=thcond/phi(i,j,k,DENS_INDX)/Cp;
            Real LewisNum_H2=0.3;
            Real LewisNum_Ar=1.0;
            Real LewisNum_H2O=1.0;

            Real dcoeff_H_gas=thdiff/LewisNum_H2;
            Real dcoeff_H2_gas=thdiff/LewisNum_H2;
            Real dcoeff_H2O_gas=thdiff/LewisNum_H2O;
            Real dcoeff_AR_gas=thdiff/LewisNum_Ar;

            //amrex::Print()<<"thdiff,dcoeff_H:"<<thdiff<<"\t"<<dcoeff_H_gas<<"\n";

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
                                                                                 

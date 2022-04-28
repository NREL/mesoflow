#ifndef _USERFUNCS_H_
#define _USERFUNCS_H_

#include<AMReX_REAL.H>
#include<AMReX.H>
#include<string>
#include<AMReX_Vector.H>
#include <AMReX_Geometry.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_Box.H>
#include <globalDefines.H>
#include <mflo_constants.H>

using namespace amrex;
namespace mflo_user_funcs
{
    extern AMREX_GPU_DEVICE_MANAGED Real Re; 
    extern AMREX_GPU_DEVICE_MANAGED Real v0; 
    extern AMREX_GPU_DEVICE_MANAGED Real rho0; 
    extern AMREX_GPU_DEVICE_MANAGED Real L; 
    extern AMREX_GPU_DEVICE_MANAGED Real gasconst; 

    void initialize_problem();

    AMREX_GPU_DEVICE AMREX_INLINE
        void compute_fluid_transport(int i, int j, int k,
                Array4<Real> const& phi,
                Array4<Real> const& transpcoeffs,
                GpuArray<Real, AMREX_SPACEDIM> prob_lo,
                GpuArray<Real, AMREX_SPACEDIM> dx,
                const Real time)
        {
            Real pr=0.7;
            Real gamma_air=1.4;
            Real R_air=one;
            Real Cp=one*gamma_air/(gamma_air-one);

            Real mu = rho0*v0*L/Re;
            transpcoeffs(i,j,k,VISC_INDX  ) = mu;       //viscosity
            transpcoeffs(i,j,k,THCOND_INDX) = mu*Cp/pr; //thermal conductivity
        }

    AMREX_GPU_DEVICE AMREX_INLINE
        void compute_fluid_source(int i, int j, int k,
                Array4<Real> const& phi,
                Array4<Real> const& source,
                GpuArray<Real, AMREX_SPACEDIM> prob_lo,
                GpuArray<Real, AMREX_SPACEDIM> dx,
                const Real time)
        {
            for(int nc=0;nc<FLO_NVARS;nc++)
            {
                source(i,j,k,nc)=zeroval;   
            }
        }
}
#endif

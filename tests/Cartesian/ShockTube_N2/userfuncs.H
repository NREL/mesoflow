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
#include <AMReX_ParmParse.H>
#include <thermo.H>

using namespace amrex;
namespace mflo_user_funcs
{
    extern AMREX_GPU_DEVICE_MANAGED Real p_l;
    extern AMREX_GPU_DEVICE_MANAGED Real rho_l;
    extern AMREX_GPU_DEVICE_MANAGED Real T_l;
    extern AMREX_GPU_DEVICE_MANAGED Real p_r;
    extern AMREX_GPU_DEVICE_MANAGED Real rho_r;
    extern AMREX_GPU_DEVICE_MANAGED Real T_r;
    extern AMREX_GPU_DEVICE_MANAGED Real spec_l[NUM_SPECIES];
    extern AMREX_GPU_DEVICE_MANAGED Real spec_r[NUM_SPECIES];

    void initialize_problem();

    AMREX_GPU_DEVICE AMREX_INLINE
        void compute_fluid_transport(int i, int j, int k,
                Array4<Real> const& phi,
                Array4<Real> const& transpcoeffs,
                GpuArray<Real, AMREX_SPACEDIM> prob_lo,
                GpuArray<Real, AMREX_SPACEDIM> dx,
                const Real time)
        {
            transpcoeffs(i,j,k,VISC_INDX  ) = zeroval; //viscosity
            transpcoeffs(i,j,k,THCOND_INDX) = zeroval; //thermal conductivity
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

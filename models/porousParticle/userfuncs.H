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
    extern AMREX_GPU_DEVICE_MANAGED int nx_mrc,ny_mrc,nz_mrc;
    extern Gpu::ManagedVector<Real>* mrcdatavec;
    extern AMREX_GPU_DEVICE_MANAGED Real *mrcdata;
    extern AMREX_GPU_DEVICE_MANAGED Real lo_mrc_x;
    extern AMREX_GPU_DEVICE_MANAGED Real lo_mrc_y;
    extern AMREX_GPU_DEVICE_MANAGED Real lo_mrc_z;
    extern AMREX_GPU_DEVICE_MANAGED Real dx_mrc_x;
    extern AMREX_GPU_DEVICE_MANAGED Real dx_mrc_y;
    extern AMREX_GPU_DEVICE_MANAGED Real dx_mrc_z;
    extern AMREX_GPU_DEVICE_MANAGED Real mrc_threshold;

    extern AMREX_GPU_DEVICE_MANAGED Real pres_left;
    extern AMREX_GPU_DEVICE_MANAGED Real temp_left;
    extern AMREX_GPU_DEVICE_MANAGED Real pres_right;
    extern AMREX_GPU_DEVICE_MANAGED Real dens_left;
    extern AMREX_GPU_DEVICE_MANAGED Real spec_left[NUM_SPECIES];

    void initialize_problem();

    AMREX_GPU_DEVICE AMREX_INLINE
    void compute_fluid_transport(int i, int j, int k,
                                 Array4<Real> const& phi,
                                 Array4<Real> const& transpcoeffs,
                                 GpuArray<Real, AMREX_SPACEDIM> prob_lo,
                                 GpuArray<Real, AMREX_SPACEDIM> dx,
                                 const Real time)
    {
        transpcoeffs(i,j,k,VISC_INDX  ) = 1e-5;   
        transpcoeffs(i,j,k,THCOND_INDX) = 0.014;
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

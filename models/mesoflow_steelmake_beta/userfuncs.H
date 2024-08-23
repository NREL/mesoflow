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
    extern AMREX_GPU_DEVICE_MANAGED Real siteconc;
    extern AMREX_GPU_DEVICE_MANAGED Real pres_left;
    extern AMREX_GPU_DEVICE_MANAGED Real temp_left;
    extern AMREX_GPU_DEVICE_MANAGED Real pres_right;
    extern AMREX_GPU_DEVICE_MANAGED Real dens_left;
    extern AMREX_GPU_DEVICE_MANAGED Real fs_vel;
    extern AMREX_GPU_DEVICE_MANAGED Real H2molfrac;
    extern AMREX_GPU_DEVICE_MANAGED Real H2Omolfrac;
    extern AMREX_GPU_DEVICE_MANAGED Real COmolfrac;
    extern AMREX_GPU_DEVICE_MANAGED Real N2molfrac;
    extern AMREX_GPU_DEVICE_MANAGED Real sponge_zone_dist;

    //extern AMREX_GPU_DEVICE_MANAGED Real spec_left[NUM_SPECIES];

    void initialize_problem();

    AMREX_GPU_DEVICE AMREX_INLINE
    void compute_fluid_transport(int i, int j, int k,
                                 Array4<Real> const& phi,
                                 Array4<Real> const& transpcoeffs,
                                 GpuArray<Real, AMREX_SPACEDIM> prob_lo,
                                 GpuArray<Real, AMREX_SPACEDIM> dx,
                                 const Real time)
    {
       // transpcoeffs(i,j,k,VISC_INDX  ) = 1e-5;
       // transpcoeffs(i,j,k,THCOND_INDX) = 0.014;
        
        // Dynamic viscosity and thermal conductivity values are fluid-specific.

           // Example values
           Real mu_H2 = 8.76e-6;   // Dynamic viscosity for H2 (Pa.s)
           Real mu_CO = 18.6e-6;   // Dynamic viscosity for CO (Pa.s)
           Real mu_H2O = 9.0e-6;   // Dynamic viscosity for H2O (Pa.s)
           Real mu_N2 = 1.78e-5;
           Real mu_CO2 = 1.48e-5;
 
           Real k_H2 = 0.182;      // Thermal conductivity for H2 (W/m.K)
           Real k_CO = 0.025;      // Thermal conductivity for CO (W/m.K)
           Real k_H2O = 0.016;     // Thermal conductivity for H2O (W/m.K)
           Real k_N2 = 0.026;
           Real k_CO2 = 0.0166;
        
        
        
           // fluid properties are computed as a weighted average based on species concentration
           Real mu = mu_H2 * phi(i,j,k,H2_ID) + mu_CO * phi(i,j,k,CO_ID) + mu_H2O * phi(i,j,k,H2O_ID)+ mu_CO2 * phi(i,j,k,CO2_ID)+mu_N2 * phi(i,j,k,N2_ID);
           Real k_mix = k_H2 * phi(i,j,k,H2_ID) + k_CO * phi(i,j,k,CO_ID) + k_H2O * phi(i,j,k,H2O_ID)+ k_CO2 * phi(i,j,k,CO2_ID)+k_N2 * phi(i,j,k,N2_ID);

           // Assign the computed values to the transport coefficients
           transpcoeffs(i,j,k,VISC_INDX) = mu;
           transpcoeffs(i,j,k,THCOND_INDX) = k_mix;
        
        
    }

    AMREX_GPU_DEVICE AMREX_INLINE
    void compute_fluid_source(int i, int j, int k,
                              Array4<Real> const& phi,
                              Array4<Real> const& source,
                              GpuArray<Real, AMREX_SPACEDIM> prob_lo,
                              GpuArray<Real, AMREX_SPACEDIM> dx,
                              const Real time)
    {
        // Initialize the source terms to zero
        // number of fluid variables
        for(int nc=0;nc<FLO_NVARS;nc++)
        {
            source(i,j,k,nc)=zeroval;   
        }
        // If there are any source terms to be added (e.g., chemical reactions, external forces),
            // they can be added here.
    }
}
#endif

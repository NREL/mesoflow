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
#include <thermo.H>

#include <AMReX_PlotFileUtil.H>


using namespace amrex;


namespace mflo_user_funcs
{
    extern AMREX_GPU_DEVICE_MANAGED Real cylrad;
    extern AMREX_GPU_DEVICE_MANAGED Real fs_vel;
    extern AMREX_GPU_DEVICE_MANAGED Real Tinit;
    extern AMREX_GPU_DEVICE_MANAGED Real Tjet;
    extern AMREX_GPU_DEVICE_MANAGED Real catblocksize;
    extern AMREX_GPU_DEVICE_MANAGED Real p0;
    extern AMREX_GPU_DEVICE_MANAGED Real H2molfrac;
    extern AMREX_GPU_DEVICE_MANAGED Real Hmolfrac;
    extern AMREX_GPU_DEVICE_MANAGED Real siteconc;
    extern AMREX_GPU_DEVICE_MANAGED Real jetrad;
    extern AMREX_GPU_DEVICE_MANAGED Real interface_ampl;
    extern AMREX_GPU_DEVICE_MANAGED Real interface_freq;
    extern AMREX_GPU_DEVICE_MANAGED Real T_coflow;
    extern AMREX_GPU_DEVICE_MANAGED Real vel_coflow;
    extern AMREX_GPU_DEVICE_MANAGED Real sponge_zone_dist;
    extern AMREX_GPU_DEVICE_MANAGED Real pin_voltage;
    extern AMREX_GPU_DEVICE_MANAGED Real pin_size;
    extern AMREX_GPU_DEVICE_MANAGED Real pin_currentden;
    extern AMREX_GPU_DEVICE_MANAGED int pin_VI_flag;
    extern AMREX_GPU_DEVICE_MANAGED Real Tsolid;
    



    
    void initialize_problem();

    AMREX_GPU_DEVICE AMREX_INLINE
        void compute_fluid_transport(int i, int j, int k,
                Array4<Real> const& phi,
                Array4<Real> const& transpcoeffs,
                GpuArray<Real, AMREX_SPACEDIM> prob_lo,
                GpuArray<Real, AMREX_SPACEDIM> dx,
                const Real time)
        {
            Real y=prob_lo[1]+(j+0.5)*dx[1];
            Real lscale=10*dx[1];
            Real trans_scale=20.0;
            Real sponge_fac=0.5*(1.0+std::tanh((amrex::Math::abs(y)-mflo_user_funcs::sponge_zone_dist)/lscale));
            Real visc, thcond;
			if(phi(i,j,k,TEMP_INDX) < 300 || phi(i,j,k,TEMP_INDX) > 30000){ 
				
					//from Boulos text book, page 1832
					//approximate fit for Ar+10% Fe
					Real scaled_temp=phi(i,j,k,TEMP_INDX)/1000.0;
					visc=(3.1*(scaled_temp)*exp(-std::pow((scaled_temp-6),2.0)/50.0)+2.5)*1e-5;
					thcond=0.006*std::pow(scaled_temp,2.0)+0.02;
				}
			else{
				
				visc =visc_inter.interpolate(phi(i,j,k,TEMP_INDX));
				thcond = thcond_inter.interpolate(phi(i,j,k,TEMP_INDX));
			}

            transpcoeffs(i,j,k,VISC_INDX  ) = visc*(1.0+sponge_fac*trans_scale);       //viscosity
            transpcoeffs(i,j,k,THCOND_INDX) = thcond*(1.0+sponge_fac*trans_scale); //thermal conductivity
        }

    AMREX_GPU_DEVICE AMREX_INLINE
    void compute_fluid_source(int i, int j, int k,
                              Array4<Real> const& phi,
                              Array4<Real> const& source,
                              GpuArray<Real, AMREX_SPACEDIM> prob_lo,
                              GpuArray<Real, AMREX_SPACEDIM> dx,
                              const Real time)
    {
        for(int nc=0;nc<TOTAL_NVARS;nc++)
        {
            source(i,j,k,nc)=zeroval;   
        }
        amrex::Real sig=mflo_thermo::conductivity(phi(i,j,k,TEMP_INDX),phi(i,j,k,VFRAC_INDX));
        amrex::Real Esq=phi(i,j,k,EFLDX_INDX)*phi(i,j,k,EFLDX_INDX) +
        phi(i,j,k,EFLDY_INDX)*phi(i,j,k,EFLDY_INDX) +
        phi(i,j,k,EFLDZ_INDX)*phi(i,j,k,EFLDZ_INDX);

        amrex::RealVect J(0.0, 0.0, 0.0);
        amrex::RealVect B(0.0, 0.0, 0.0);

        for(int d=0; d<3; d++)
        {
            B[d]=phi(i,j,k,BX_INDX+d);
            J[d]=sig*phi(i,j,k,EFLDX_INDX+d);
        }
        amrex::RealVect J_cross_B = J.crossProduct(B);
        
        //Lorentz force
        source(i,j,k,RHOU_INDX) += phi(i,j,k,VFRAC_INDX)*J_cross_B[0];
        source(i,j,k,RHOV_INDX) += phi(i,j,k,VFRAC_INDX)*J_cross_B[1];
        source(i,j,k,RHOW_INDX) += phi(i,j,k,VFRAC_INDX)*J_cross_B[2];

        //Joule heating 
        source(i,j,k,RHOE_INDX) += phi(i,j,k,VFRAC_INDX)*sig*Esq;
    }

    AMREX_GPU_DEVICE AMREX_INLINE
    void potential_bc(int i, int j, int k,
                      int dir, int sgn,
                      Array4<Real> const& phi,
                      Array4<Real> const& robin_a,
                      Array4<Real> const& robin_b,
                      Array4<Real> const& robin_f,
                      GpuArray<Real, AMREX_SPACEDIM> prob_lo,
                      GpuArray<Real, AMREX_SPACEDIM> prob_hi,
                      GpuArray<Real, AMREX_SPACEDIM> dx,
                      const Real time)
    {
        IntVect cell_int{AMREX_D_DECL(i,j,k)};
        IntVect ghost_cell{AMREX_D_DECL(i,j,k)};
        amrex::Real outward_normal[AMREX_SPACEDIM]={0.0};
        outward_normal[dir]=sgn;

        int gcell_adjust,cell_adjust;

        //ghost_cell is one behind
        gcell_adjust= (sgn==-1)?-1:0;
        //cell is one behind
        cell_adjust = (sgn==1)?-1:0;

        ghost_cell[dir]+=gcell_adjust;
        cell_int[dir]+=cell_adjust;

        if(sgn == 1 && dir==0) 
        { 

            robin_a(ghost_cell) = 1.0;
            robin_b(ghost_cell) = 0.0;
            robin_f(ghost_cell) = 0.0;
        }
        else if(sgn == -1 && dir==0)
        {
            // lo sides
            amrex::Real electrode_rad=mflo_user_funcs::pin_size;

            amrex::RealVect x;
            x[0]=prob_lo[0]+(i+0.5)*dx[0];
            x[1]=prob_lo[1]+(j+0.5)*dx[1];
            x[2]=prob_lo[2]+(k+0.5)*dx[2];

            amrex::Real loc_rad=x[0]*x[0]+x[1]*x[1]+x[2]*x[2];
            loc_rad -= x[dir]*x[dir];
            loc_rad=std::sqrt(loc_rad);

            if(loc_rad < electrode_rad)
            {
                if(mflo_user_funcs::pin_VI_flag==0)
                {
                    robin_a(ghost_cell) = 1.0;
                    robin_b(ghost_cell) = 0.0;
                    robin_f(ghost_cell) = mflo_user_funcs::pin_voltage;
                }
                else
                {
                    robin_a(ghost_cell) = 0.0;
                    robin_b(ghost_cell) = 1.0;
                    Real temp=phi(cell_int,TEMP_INDX);
                    Real sig=mflo_thermo::conductivity(temp,phi(cell_int,VFRAC_INDX));
                    robin_f(ghost_cell) = mflo_user_funcs::pin_currentden/sig;
                }
            }
            else
            {
                robin_a(ghost_cell) = 0.0;
                robin_b(ghost_cell) = 1.0;
                robin_f(ghost_cell) = 0.0;
            }
        }
        else
        {
            robin_a(ghost_cell) = 0.0;
            robin_b(ghost_cell) = 1.0;
            robin_f(ghost_cell) = 0.0;
        } 
    }

    AMREX_GPU_DEVICE AMREX_INLINE
    void add_potential_sources(int i, int j, int k,
                               Array4<Real> const& phi,
                               Array4<Real> const& rhs,
                               GpuArray<Real, AMREX_SPACEDIM> prob_lo,
                               GpuArray<Real, AMREX_SPACEDIM> prob_hi,
                               GpuArray<Real, AMREX_SPACEDIM> dx,
                               const Real time)
    {
        amrex::Real src=0.0;
        rhs(i,j,k)+=src; 
    }

    AMREX_GPU_DEVICE AMREX_INLINE
    void potential_dcoeff(int i, int j, int k,
                          Array4<Real> const& phi,
                          Array4<Real> const& dcoeff,
                          GpuArray<Real, AMREX_SPACEDIM> prob_lo,
                          GpuArray<Real, AMREX_SPACEDIM> prob_hi,
                          GpuArray<Real, AMREX_SPACEDIM> dx,
                          const Real time)
    {
        dcoeff(i,j,k)=mflo_thermo::conductivity(phi(i,j,k,TEMP_INDX),phi(i,j,k,VFRAC_INDX));
    }

    AMREX_GPU_DEVICE AMREX_INLINE
    void vec_potential_bc(int i, int j, int k, int Adir,
                          int dir, int sgn,
                          Array4<Real> const& phi,
                          Array4<Real> const& robin_a,
                          Array4<Real> const& robin_b,
                          Array4<Real> const& robin_f,
                          GpuArray<Real, AMREX_SPACEDIM> prob_lo,
                          GpuArray<Real, AMREX_SPACEDIM> prob_hi,
                          GpuArray<Real, AMREX_SPACEDIM> dx,
                          const Real time)
    {
        IntVect cell_int{AMREX_D_DECL(i,j,k)};
        IntVect ghost_cell{AMREX_D_DECL(i,j,k)};
        amrex::Real outward_normal[AMREX_SPACEDIM]={0.0};
        outward_normal[dir]=sgn;

        int gcell_adjust,cell_adjust;

        //ghost_cell is one behind
        gcell_adjust= (sgn==-1)?-1:0;
        //cell is one behind
        cell_adjust = (sgn==1)?-1:0;

        ghost_cell[dir]+=gcell_adjust;
        cell_int[dir]+=cell_adjust;
        robin_a(ghost_cell) = 1.0;
        robin_b(ghost_cell) = 0.0;
        robin_f(ghost_cell) = 0.0;
    }
}
#endif

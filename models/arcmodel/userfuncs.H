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

extern    AMREX_GPU_DEVICE_MANAGED Real fs_vel;
extern    AMREX_GPU_DEVICE_MANAGED Real Tinit;
extern    AMREX_GPU_DEVICE_MANAGED Real T_ign;
extern    AMREX_GPU_DEVICE_MANAGED Real T_surr;
extern    AMREX_GPU_DEVICE_MANAGED Real T_elec;
extern    AMREX_GPU_DEVICE_MANAGED Real T_plate;
extern    AMREX_GPU_DEVICE_MANAGED Real catblocksize;
extern    AMREX_GPU_DEVICE_MANAGED Real p0;
extern    AMREX_GPU_DEVICE_MANAGED Real H2molfrac;
extern    AMREX_GPU_DEVICE_MANAGED Real Hmolfrac;
extern    AMREX_GPU_DEVICE_MANAGED Real siteconc;
extern    AMREX_GPU_DEVICE_MANAGED Real sponge_zone_dist;
extern    AMREX_GPU_DEVICE_MANAGED Real eps;
extern    AMREX_GPU_DEVICE_MANAGED Real tort;
extern    AMREX_GPU_DEVICE_MANAGED Real elec_rad;
extern    AMREX_GPU_DEVICE_MANAGED Real elec_voltage;
extern    AMREX_GPU_DEVICE_MANAGED Real elec_currentden;
extern    AMREX_GPU_DEVICE_MANAGED Real T_solid;
extern    AMREX_GPU_DEVICE_MANAGED int use_V_flag;
extern    AMREX_GPU_DEVICE_MANAGED Real resist;
extern    AMREX_GPU_DEVICE_MANAGED Real T_coflow;
extern    AMREX_GPU_DEVICE_MANAGED Real vel_coflow;

    
    void initialize_problem();

    AMREX_GPU_DEVICE AMREX_INLINE
        void compute_fluid_transport(int i, int j, int k,
                Array4<Real> const& phi,
                Array4<Real> const& transpcoeffs,
                GpuArray<Real, AMREX_SPACEDIM> prob_lo,
                GpuArray<Real, AMREX_SPACEDIM> dx,
                const Real time)
        {
      		
			Real visc = (phi(i,j,k,FLO_NVARS+AR_ID)*visc_AR_inter.interpolate(phi(i,j,k,TEMP_INDX)) + phi(i,j,k,FLO_NVARS+H2_ID)*visc_H2_inter.interpolate(phi(i,j,k,TEMP_INDX)))/ ( phi(i,j,k,FLO_NVARS+AR_ID) + phi(i,j,k,FLO_NVARS+H2_ID));
			Real thcond = (phi(i,j,k,FLO_NVARS+AR_ID)*thcond_AR_inter.interpolate(phi(i,j,k,TEMP_INDX)) + phi(i,j,k,FLO_NVARS+H2_ID)*thcond_H2_inter.interpolate(phi(i,j,k,TEMP_INDX)))/ ( phi(i,j,k,FLO_NVARS+AR_ID) + phi(i,j,k,FLO_NVARS+H2_ID));
		

            transpcoeffs(i,j,k,VISC_INDX  ) = visc;     //viscosity
            transpcoeffs(i,j,k,THCOND_INDX) = thcond; //thermal conductivity
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
        amrex::Real sig=mflo_thermo::conductivity(phi(i,j,k,TEMP_INDX),phi(i,j,k,VFRAC_INDX),phi(i,j,k,FLO_NVARS+AR_ID),phi(i,j,k,FLO_NVARS+H2_ID));
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
        
        //Radiation loses
        source(i,j,k,RHOE_INDX) -= phi(i,j,k,VFRAC_INDX)*rad_inter.interpolate(phi(i,j,k,TEMP_INDX));
    }

    AMREX_GPU_DEVICE AMREX_INLINE
    void potential_bc(int i, int j, int k,
                      int dir, int sgn,
                      Array4<Real> const& phi,
                      Array4<Real> const& robin_a,
                      Array4<Real> const& robin_b,
                      Array4<Real> const& robin_f,
                      Real current,
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

        if(sgn == 1 && dir==0) // x_high dirichlet 0
        { 

            robin_a(ghost_cell) = 1.0;
            robin_b(ghost_cell) = 0.0;
            robin_f(ghost_cell) = 0.0;
        }
        else if(sgn == -1 && dir==0)  //  x_low robin bc
        {
            // lo sides
            amrex::Real electrode_rad=mflo_user_funcs::elec_rad;

            amrex::RealVect x;
            x[0]=prob_lo[0]+(i+0.5)*dx[0];
            x[1]=prob_lo[1]+(j+0.5)*dx[1];
            x[2]=prob_lo[2]+(k+0.5)*dx[2];

            amrex::Real loc_rad=x[0]*x[0]+x[1]*x[1]+x[2]*x[2];
            loc_rad -= x[dir]*x[dir];
            loc_rad=std::sqrt(loc_rad);   

            if(loc_rad < electrode_rad)     
            {
                if(mflo_user_funcs::use_V_flag==1)   //-v dirichlet pin_voltage
                {
                    robin_a(ghost_cell) = 1.0;
                    robin_b(ghost_cell) = 0.0;
                    robin_f(ghost_cell) = std::max(mflo_user_funcs::elec_voltage - current*resist, 0.0);
                }
                else
                {                                   //-v fixed gradient
                    robin_a(ghost_cell) = 0.0;
                    robin_b(ghost_cell) = 1.0;
                    Real temp=phi(cell_int,TEMP_INDX);
                    Real sig=mflo_thermo::conductivity(temp,phi(cell_int,VFRAC_INDX),phi(cell_int,FLO_NVARS+AR_ID),phi(cell_int,FLO_NVARS+H2_ID));
                    robin_f(ghost_cell) = mflo_user_funcs::elec_currentden/sig;
                }
            }
            else
            {
                robin_a(ghost_cell) = 0.0;
                robin_b(ghost_cell) = 1.0;
                robin_f(ghost_cell) = 0.0;
            }
        }
        else  // all other boundaries
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
                          Array4<Real> const& state,
                          GpuArray<Real, AMREX_SPACEDIM> prob_lo,
                          GpuArray<Real, AMREX_SPACEDIM> prob_hi,
                          GpuArray<Real, AMREX_SPACEDIM> dx,
                          const Real time)
    {
        dcoeff(i,j,k)=mflo_thermo::conductivity(phi(i,j,k,TEMP_INDX),phi(i,j,k,VFRAC_INDX),phi(i,j,k,FLO_NVARS+AR_ID),phi(i,j,k,FLO_NVARS+H2_ID));
        state(i,j,k,SIG_INDX)=dcoeff(i,j,k);
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
    
    AMREX_GPU_DEVICE AMREX_INLINE
    void update_bndry_face_bcoeff(int i, int j, int k,
                      int dir, int sgn,
                      Array4<Real> const& phi,
                      Array4<Real> const& face_bcoeff,
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
	   
	   if(sgn == 1 && dir==0) // x_high 
        { 
			Real face_temp = 0.5*( phi( cell_int,TEMP_INDX) + phi( ghost_cell,TEMP_INDX) );
			face_bcoeff(i,j,k) = mflo_thermo::conductivity(face_temp,phi( cell_int,VFRAC_INDX),phi( cell_int,FLO_NVARS+AR_ID),phi( cell_int,FLO_NVARS+H2_ID));
	   
		}
	   
	   else if(sgn == -1 && dir==0)  //  x_low 
        {
			Real face_temp = 0.5*( phi( cell_int,TEMP_INDX) + phi( ghost_cell,TEMP_INDX) );
			face_bcoeff(i,j,k) = mflo_thermo::conductivity(face_temp,phi( cell_int,VFRAC_INDX),phi( cell_int,FLO_NVARS+AR_ID),phi( cell_int,FLO_NVARS+H2_ID));
		}
		else  // everywhere else neumann conditions 
		{
			face_bcoeff(i,j,k) = mflo_thermo::conductivity(phi( cell_int,TEMP_INDX),phi( cell_int,VFRAC_INDX),phi( cell_int,FLO_NVARS+AR_ID),phi( cell_int,FLO_NVARS+H2_ID));
		}
	      
	   
   }
}
#endif

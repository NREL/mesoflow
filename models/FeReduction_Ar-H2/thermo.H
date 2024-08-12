#ifndef _THERMO_H_
#define _THERMO_H_

#include<AMReX_REAL.H>
#include<AMReX.H>
#include<string>
#include<AMReX_Vector.H>
#include <AMReX_Geometry.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_Box.H>
#include <globalDefines.H>
#include <mflo_constants.H>
#include <species.H>
#include <cubicSpline.H>

#include <AMReX_PlotFileUtil.H>
using namespace amrex;

extern AMREX_GPU_DEVICE_MANAGED CubicSplineInterpolator visc_AR_inter;
extern AMREX_GPU_DEVICE_MANAGED CubicSplineInterpolator thcond_AR_inter;
extern AMREX_GPU_DEVICE_MANAGED CubicSplineInterpolator sig_AR_inter;
extern AMREX_GPU_DEVICE_MANAGED CubicSplineInterpolator e_AR_inter;
extern AMREX_GPU_DEVICE_MANAGED CubicSplineInterpolator visc_H2_inter;
extern AMREX_GPU_DEVICE_MANAGED CubicSplineInterpolator thcond_H2_inter;
extern AMREX_GPU_DEVICE_MANAGED CubicSplineInterpolator sig_H2_inter;
extern AMREX_GPU_DEVICE_MANAGED CubicSplineInterpolator e_H2_inter;
extern AMREX_GPU_DEVICE_MANAGED CubicSplineInterpolator rad_inter;

namespace mflo_thermo
{
	

	
    AMREX_GPU_HOST_DEVICE AMREX_INLINE Real 
    get_solid_t_from_rhoe(Real rhoe)
    {
        //https://webbook.nist.gov/cgi/cbook.cgi?ID=C1317608&Type=JANAFS&Table=on
        // at 500 C, Cp is 131.7 J/mol/K = 131.7/160=0.82 J/g/K
        //160 is molwt of Fe2O3
        amrex::Real Cv_solid=820.0;
        amrex::Real rho_solid=5240.0;
        amrex::Real temp=rhoe/(rho_solid*Cv_solid);
        //amrex::Print()<<"rhoe,temp:"<<rhoe<<"\t"<<temp<<"\n";
        return(temp);
    }
    
    AMREX_GPU_HOST_DEVICE AMREX_INLINE Real 
    get_solid_rhoe_from_t(Real temp)
    {
        amrex::Real Cv_solid=820.0;
        amrex::Real rho_solid=5240.0;
        return(rho_solid*Cv_solid*temp);
    }
    
    AMREX_GPU_HOST_DEVICE AMREX_INLINE
    Real get_re_from_ct(Real spec[NUM_SPECIES],Real temp)
    {
        Real rhoe=0.0;
        Real gama_h2o=1.2;
        Real gama_h=1.66;
        Real tempk=temp/1000.0;
                                                   
	
		rhoe += spec[AR_ID]*mflo_species::molwts[AR_ID]*e_AR_inter.interpolate(temp);  //-v internal energy/kg basis
			
		rhoe += spec[H2_ID]*mflo_species::molwts[H2_ID]*e_H2_inter.interpolate(temp);
		  
        
        rhoe += spec[H_ID]*RU/(gama_h-1.0)*temp;
        rhoe += spec[H2O_ID]*RU/(gama_h2o-1.0)*temp;
        return(rhoe);
    }

    AMREX_GPU_HOST_DEVICE AMREX_INLINE
    Real get_t_from_cre(Real spec[NUM_SPECIES],Real re,Real tguess)
    {
        int nsteps=30;
        Real tn=tguess;
        Real tol=1e-7*tguess;
        Real fd_dt=1e-6*tguess; //finite diff dt
        Real delt_n=0.0;
        //Newton method
        for(int st=0;st<nsteps;st++)
        {
            Real re_t=get_re_from_ct(spec,tn);
            Real re_tpdt=get_re_from_ct(spec,tn+fd_dt);
            Real re_deriv=(re_tpdt-re_t)/fd_dt;

            Real func_T=re_t-re;
            Real funcderiv_T=re_deriv;
            delt_n=-func_T/funcderiv_T;
            tn+=delt_n;
            //amrex::Print()<<"st,delt_n,re_t,ret_tpdt,re:"<<st<<"\t"<<delt_n<<"\t"<<re_t<<"\t"<<re_tpdt<<"\t"<<re<<"\n";
            if(amrex::Math::abs(delt_n)<tol)
            {
                //amrex::Print()<<"breaking early\n";
                break;
            }
        }
        if(amrex::Math::abs(delt_n) > 1000*tol)
        {
            amrex::Print()<<"energy solve not converging\n";
        }
        //amrex::Print()<<"T:"<<tn<<"\n";
        //amrex::Print()<<"\n============\n";
        return(tn);
    }


    AMREX_GPU_HOST_DEVICE AMREX_INLINE
    Real get_r_from_c(Real spec[NUM_SPECIES],Real bgasconc=zeroval)
    {
        Real rho = zeroval;     
        for(int sp=0;sp<NUM_GAS_SPECIES;sp++)
        {
            rho += spec[sp]*mflo_species::molwts[sp];   
        }
        return(rho);
    }

    AMREX_GPU_HOST_DEVICE AMREX_INLINE
    Real get_bgasconc_from_rc(Real rho,Real spec[NUM_SPECIES])
    {
        //this function is used only when there is background species
        return(0.0);
    }

    AMREX_GPU_HOST_DEVICE AMREX_INLINE
    Real get_t_from_rpc(Real rho,Real pres,Real spec[NUM_SPECIES])
    {
        Real sum_ci   = zeroval;
        for(int sp=0;sp<NUM_GAS_SPECIES;sp++)
        {
            sum_ci   += spec[sp];   
        }
        return(pres/(sum_ci*RU));
    }

    AMREX_GPU_HOST_DEVICE AMREX_INLINE
    Real get_p_from_rec(Real rhoe,Real rho,Real spec[NUM_SPECIES])
    {
        Real T;
        Real tguess=1000.0;
        T = get_t_from_cre(spec,rhoe,tguess);
        Real p=zeroval;
        for(int sp=0;sp<NUM_GAS_SPECIES;sp++)
        {
            p += spec[sp]*RU*T;   
        }
        return(p);
    }

    AMREX_GPU_HOST_DEVICE AMREX_INLINE
    Real get_gama_from_rpc(Real rho,Real pres,Real spec[NUM_SPECIES])
    {
        Real gamma_spec[NUM_GAS_SPECIES]={zeroval};
        Real x_i[NUM_GAS_SPECIES]={zeroval};  //mole frac

        gamma_spec[H2_ID]  = 1.4;
        gamma_spec[H_ID]  = 1.66;
        gamma_spec[H2O_ID]  = 1.3;
        gamma_spec[AR_ID]  = 1.66;

        //find total number of moles
        Real sum_ci=zeroval;
        for(int sp=0;sp<NUM_GAS_SPECIES;sp++)
        {
            sum_ci += spec[sp];
        }

        //find molefraction
        for(int sp=0;sp<NUM_GAS_SPECIES;sp++)
        {
            x_i[sp]=spec[sp]/sum_ci;
        }

        Real term2=zeroval;
        for(int sp=0;sp<NUM_GAS_SPECIES;sp++)
        {
            term2 += x_i[sp]/(gamma_spec[sp]-one);
        }

        Real gama_mix=one+pow(term2,-one);

        return(gama_mix);
    }

    AMREX_GPU_HOST_DEVICE AMREX_INLINE
    Real get_e_from_rpc(Real rho,Real pres,Real spec[NUM_SPECIES])
    {
        Real temp=get_t_from_rpc(rho,pres,spec);
        Real rhoe=get_re_from_ct(spec,temp);
        return( rhoe/rho );
    }

    AMREX_GPU_HOST_DEVICE AMREX_INLINE
    Real conductivity(Real temp,Real vfrac,Real AR_conc, Real H2_conc)
    {
        Real sig_solid=1.0; //FIX: for FEO
        Real sig_gas=(1+std::tanh(0.0005*(temp-7000.0)))*1000.0+1.0;
        //Real sig_gas= (AR_conc*std::max(sig_AR_inter.interpolate(temp), 1e-1) + H2_conc*std::max(sig_H2_inter.interpolate(temp), 1e-1) )/ (AR_conc + H2_conc);      
        Real sig=vfrac*sig_gas+(1.0-vfrac)*sig_solid;
        return(sig);
    }
}
#endif
        

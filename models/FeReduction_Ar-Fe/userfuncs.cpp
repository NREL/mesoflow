#include<userfuncs.H>
#include <AMReX_ParmParse.H>

AMREX_GPU_DEVICE_MANAGED CubicSplineInterpolator visc_AR_inter;
AMREX_GPU_DEVICE_MANAGED CubicSplineInterpolator thcond_AR_inter;
AMREX_GPU_DEVICE_MANAGED CubicSplineInterpolator sig_AR_inter;
AMREX_GPU_DEVICE_MANAGED CubicSplineInterpolator e_AR_inter;
AMREX_GPU_DEVICE_MANAGED CubicSplineInterpolator visc_H2_inter;
AMREX_GPU_DEVICE_MANAGED CubicSplineInterpolator thcond_H2_inter;
AMREX_GPU_DEVICE_MANAGED CubicSplineInterpolator sig_H2_inter;
AMREX_GPU_DEVICE_MANAGED CubicSplineInterpolator e_H2_inter;
AMREX_GPU_DEVICE_MANAGED CubicSplineInterpolator visc_AR_Fe_inter;
AMREX_GPU_DEVICE_MANAGED CubicSplineInterpolator thcond_AR_Fe_inter;
AMREX_GPU_DEVICE_MANAGED CubicSplineInterpolator rad_inter;

namespace mflo_user_funcs
{
    AMREX_GPU_DEVICE_MANAGED Real fs_vel=1.0;
    AMREX_GPU_DEVICE_MANAGED Real Tinit=300.0;
    AMREX_GPU_DEVICE_MANAGED Real Tjet=1000.0;
    AMREX_GPU_DEVICE_MANAGED Real catblocksize=2e-3;
    AMREX_GPU_DEVICE_MANAGED Real p0=1e5;
    AMREX_GPU_DEVICE_MANAGED Real H2molfrac=0.05;
    AMREX_GPU_DEVICE_MANAGED Real Hmolfrac=0.05;
    AMREX_GPU_DEVICE_MANAGED Real siteconc=70000.0;
    AMREX_GPU_DEVICE_MANAGED Real jetrad=1.0e-3;
    AMREX_GPU_DEVICE_MANAGED Real interface_ampl=2e-4;
    AMREX_GPU_DEVICE_MANAGED Real interface_freq=10;
    AMREX_GPU_DEVICE_MANAGED Real T_coflow=300.0;
    AMREX_GPU_DEVICE_MANAGED Real vel_coflow=0.01;
    AMREX_GPU_DEVICE_MANAGED Real sponge_zone_dist=20e-3;
    AMREX_GPU_DEVICE_MANAGED Real pin_voltage=1000.0;
    AMREX_GPU_DEVICE_MANAGED Real pin_currentden=1e6;
    AMREX_GPU_DEVICE_MANAGED Real pin_size=1.25e-4;
    AMREX_GPU_DEVICE_MANAGED int pin_VI_flag=0;
    AMREX_GPU_DEVICE_MANAGED Real Tsolid=300.0;
    




    void initialize_problem()
    {
        Print()<<"Initializing problem\n";
    
        ParmParse pp("user");
        pp.query("fs_vel",fs_vel);
        pp.query("Tinit",Tinit);
        pp.query("Tjet",Tjet);
        pp.query("catblocksize",catblocksize);
        pp.query("interface_ampl",interface_ampl);
        pp.query("interface_freq",interface_freq);
        pp.query("jetrad",jetrad);
        pp.query("p0",p0);
        pp.query("H2molfrac",H2molfrac);
        pp.query("Hmolfrac",Hmolfrac);
        pp.query("siteconc",siteconc);
        pp.query("sponge_zone_dist",sponge_zone_dist);
        pp.query("pin_VI_flag",pin_VI_flag);
        pp.query("pin_voltage",pin_voltage);
        pp.query("pin_currentden",pin_currentden);
        pp.query("pin_size",pin_size);
        pp.query("Tsolid",Tsolid);

        T_coflow=Tinit;
        vel_coflow=0.01*fs_vel;
        
        pp.query("T_coflow",T_coflow);
        pp.query("vel_coflow",vel_coflow);
        

        
        
        visc_AR_inter.define("data_visc_AR.txt",1);
        thcond_AR_inter.define("data_thcond_AR.txt",1);
        sig_AR_inter.define("data_sig_AR.txt",1);  
        e_AR_inter.define("data_e_AR.txt",1);
        visc_H2_inter.define("data_visc_H2.txt",1);
        thcond_H2_inter.define("data_thcond_H2.txt",1);
        sig_H2_inter.define("data_sig_H2.txt",1);  
        e_H2_inter.define("data_e_H2.txt",1);
        visc_AR_Fe_inter.define("data_visc_Ar+Fe.txt",1);
        thcond_AR_Fe_inter.define("data_thcond_AR+Fe.txt",1);
        rad_inter.define("radiation.txt",1);       
		
	
    }
}

#include<userfuncs.H>
#include <AMReX_ParmParse.H>

namespace mflo_user_funcs
{
    AMREX_GPU_DEVICE_MANAGED Real cylrad=0.5;
    AMREX_GPU_DEVICE_MANAGED Real fs_vel=0.2366;
    AMREX_GPU_DEVICE_MANAGED Real rho0=1.0;
    AMREX_GPU_DEVICE_MANAGED Real p0=1.0;
    AMREX_GPU_DEVICE_MANAGED Real Aconc=0.0;
    AMREX_GPU_DEVICE_MANAGED Real siteconc=1.0;

    void initialize_problem()
    {
        Print()<<"Initializing problem\n";
    
        ParmParse pp("user");
        pp.query("fs_vel",fs_vel);
        pp.query("cylrad",cylrad);
        pp.query("rho0",rho0);
        pp.query("p0",p0);
        pp.query("Aconc",Aconc);
        pp.query("siteconc",siteconc);
    }
}

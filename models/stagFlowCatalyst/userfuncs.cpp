#include<userfuncs.H>
#include <AMReX_ParmParse.H>

namespace mflo_user_funcs
{
    AMREX_GPU_DEVICE_MANAGED Real fs_vel=0.2366;
    AMREX_GPU_DEVICE_MANAGED Real rho0=1.0;
    AMREX_GPU_DEVICE_MANAGED Real catblocksize=2e-3;
    AMREX_GPU_DEVICE_MANAGED Real p0=1e5;
    AMREX_GPU_DEVICE_MANAGED Real CH4conc=0.0;
    AMREX_GPU_DEVICE_MANAGED Real siteconc=1.0;
    AMREX_GPU_DEVICE_MANAGED Real jetrad=1.0e-3;
    AMREX_GPU_DEVICE_MANAGED Real Tsolid=600.0;

    void initialize_problem()
    {
        Print()<<"Initializing problem\n";
    
        ParmParse pp("user");
        pp.query("fs_vel",fs_vel);
        pp.query("rho0",rho0);
        pp.query("catblocksize",catblocksize);
        pp.query("jetrad",jetrad);
        pp.query("p0",p0);
        pp.query("CH4conc",CH4conc);
        pp.query("siteconc",siteconc);
        pp.query("Tsolid",Tsolid);
    }
}

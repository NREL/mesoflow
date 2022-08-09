#include<userfuncs.H>
#include <AMReX_ParmParse.H>

namespace mflo_user_funcs
{
    AMREX_GPU_DEVICE_MANAGED Real sphrad=0.5;
    AMREX_GPU_DEVICE_MANAGED Real fs_vel=0.2366;
    AMREX_GPU_DEVICE_MANAGED Real Re=300.0;
    AMREX_GPU_DEVICE_MANAGED Real rho0=1.0;
    AMREX_GPU_DEVICE_MANAGED Real p0=1.0;

    void initialize_problem()
    {
        Print()<<"Initializing problem\n";
    
        ParmParse pp("user");
        pp.query("fs_vel",fs_vel);
        pp.query("sphrad",sphrad);
        pp.query("Re",Re);
        pp.query("rho0",rho0);
        pp.query("p0",p0);
    }
}

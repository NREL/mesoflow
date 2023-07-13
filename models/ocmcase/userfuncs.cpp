
#include<userfuncs.H>
#include <AMReX_ParmParse.H>
#include <thermo.H>

namespace mflo_user_funcs
{
    AMREX_GPU_DEVICE_MANAGED int chemistry_on=0;
    AMREX_GPU_DEVICE_MANAGED Real catalyst_loc=0.7; //added from blockCatalyst
    AMREX_GPU_DEVICE_MANAGED Real siteconc=1.0;
    AMREX_GPU_DEVICE_MANAGED Real fs_vel=1.0;

    void initialize_problem()
    {
        ParmParse pp("problem");
        pp.query("catalyst_loc",catalyst_loc); //added from blockCatalyst
        pp.query("free_stream_vel",fs_vel); 
        pp.query("chemistry_on",chemistry_on);
        pp.query("siteconc",siteconc);
    }
}

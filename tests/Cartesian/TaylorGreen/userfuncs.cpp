#include<userfuncs.H>

namespace mflo_user_funcs
{
    AMREX_GPU_DEVICE_MANAGED Real Re=1600.0;
    AMREX_GPU_DEVICE_MANAGED Real v0=tenth;
    AMREX_GPU_DEVICE_MANAGED Real rho0=one;
    AMREX_GPU_DEVICE_MANAGED Real L=one/PI;

    void initialize_problem()
    {
        Print()<<"Initializing problem\n";

        //This is a good place to  have some global
        //parameters defined for initializing solution 
        //vector
    }
}

#include<userfuncs.H>
#include <AMReX_ParmParse.H>

namespace mflo_user_funcs
{
    AMREX_GPU_DEVICE_MANAGED Real iboffset=0.7;
    void initialize_problem()
    {
        ParmParse pp("user");
        pp.query("iboffset",iboffset);
        
        Print()<<"Initializing problem\n";

        //This is a good place to  have some global
        //parameters defined for initializing solution 
        //vector
    }
}

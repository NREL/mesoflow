#include<userfuncs.H>
#include <AMReX_ParmParse.H>

namespace mflo_user_funcs
{
    AMREX_GPU_DEVICE_MANAGED Real iboffset=0.7;
    AMREX_GPU_DEVICE_MANAGED Real specS1=0.7; //concentration of reactant gas
    void initialize_problem()
    {
        ParmParse pp("user");
        pp.query("iboffset",iboffset);
        pp.query("spec_S1",specS1);
        Print()<<"Initializing problem\n";

        //This is a good place to  have some global
        //parameters defined for initializing solution 
        //vector
    }
}

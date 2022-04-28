#include<userfuncs.H>

namespace mflo_user_funcs
{
    AMREX_GPU_DEVICE_MANAGED Real p_l=one;
    AMREX_GPU_DEVICE_MANAGED Real rho_l=one;
    AMREX_GPU_DEVICE_MANAGED Real T_l=one;
    AMREX_GPU_DEVICE_MANAGED Real p_r=tenth;
    AMREX_GPU_DEVICE_MANAGED Real rho_r=eighth;
    AMREX_GPU_DEVICE_MANAGED Real T_r=one;
    AMREX_GPU_DEVICE_MANAGED Real spec_l[NUM_SPECIES]={zeroval};
    AMREX_GPU_DEVICE_MANAGED Real spec_r[NUM_SPECIES]={zeroval};

    void initialize_problem()
    {
        //This is a good place to  have some global
        //parameters defined for initializing solution 
        //vector
        Print()<<"Initializing problem\n";

        amrex::ParmParse pp("prob");
        pp.query("p_l", p_l);
        pp.query("rho_l", rho_l);
        pp.query("p_r", p_r);
        pp.query("rho_r", rho_r);

        spec_l[N2_ID]=rho_l/mflo_species::molwts[N2_ID];
        spec_l[HE_ID]=zeroval;

        spec_r[N2_ID]=zeroval;
        spec_r[HE_ID]=rho_r/mflo_species::molwts[HE_ID];

        T_l=mflo_thermo::get_t_from_rpc(rho_l,p_l,spec_l);
        T_r=mflo_thermo::get_t_from_rpc(rho_r,p_r,spec_r);
    }
}

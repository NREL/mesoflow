#include<userfuncs.H>


namespace mflo_user_funcs
{
    AMREX_GPU_DEVICE_MANAGED Real cylinder_rad      = tenth;
    AMREX_GPU_DEVICE_MANAGED Real pres_fs         = one;
    AMREX_GPU_DEVICE_MANAGED Real temp_fs         = one;
    AMREX_GPU_DEVICE_MANAGED Real dens_fs         = one;
    AMREX_GPU_DEVICE_MANAGED Real vel_fs          = one;
    AMREX_GPU_DEVICE_MANAGED Real spec_fs[NUM_SPECIES] = {one};

    void initialize_problem()
    {
        ParmParse pp("prob");

        pp.query("cylinder_rad",cylinder_rad);
        pp.query("freestream_vel",vel_fs);

        spec_fs[AIR_ID]=dens_fs/mflo_species::molwts[AIR_ID];
        temp_fs=mflo_thermo::get_t_from_rpc(dens_fs,pres_fs,spec_fs);

    }
}

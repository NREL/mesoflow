#include<userfuncs.H>


namespace mflo_user_funcs
{
    AMREX_GPU_DEVICE_MANAGED Real sphere_rad      = tenth;
    AMREX_GPU_DEVICE_MANAGED Real Mach_fs         = two;
    AMREX_GPU_DEVICE_MANAGED Real pres_fs         = P_NTP;
    AMREX_GPU_DEVICE_MANAGED Real temp_fs         = T_NTP;
    AMREX_GPU_DEVICE_MANAGED Real dens_fs         = one;
    AMREX_GPU_DEVICE_MANAGED Real vel_fs          = one;
    AMREX_GPU_DEVICE_MANAGED Real spec_fs[NUM_SPECIES] = {one};

    void initialize_problem()
    {
        ParmParse pp("prob");

        pp.query("sphere_rad",sphere_rad);
        pp.query("freestream_Mach",Mach_fs);
        pp.query("freestream_pressure",pres_fs);
        pp.query("freestream_temperature",temp_fs);

        Real R_air=RU/mflo_species::molwts[AIR_ID];
        dens_fs=pres_fs/R_air/temp_fs;

        Real gamma_air=1.4;
        vel_fs=Mach_fs*sqrt(gamma_air*R_air*temp_fs);

        spec_fs[AIR_ID]=dens_fs/mflo_species::molwts[AIR_ID];
    }
}

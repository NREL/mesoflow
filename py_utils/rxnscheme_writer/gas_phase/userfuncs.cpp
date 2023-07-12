
#include<userfuncs.H>
#include <AMReX_ParmParse.H>

namespace mflo_user_funcs
{
    AMREX_GPU_DEVICE_MANAGED int nx_mrc=0;
    AMREX_GPU_DEVICE_MANAGED int ny_mrc=0;
    AMREX_GPU_DEVICE_MANAGED int nz_mrc=0;
    Gpu::ManagedVector<Real>* mrcdatavec=NULL;
    AMREX_GPU_DEVICE_MANAGED Real *mrcdata=NULL;
    AMREX_GPU_DEVICE_MANAGED Real lo_mrc_x=-6e-5;
    AMREX_GPU_DEVICE_MANAGED Real lo_mrc_y=-6e-5;
    AMREX_GPU_DEVICE_MANAGED Real lo_mrc_z=-6e-5;
    AMREX_GPU_DEVICE_MANAGED Real dx_mrc_x=1e-6;
    AMREX_GPU_DEVICE_MANAGED Real dx_mrc_y=1e-6;
    AMREX_GPU_DEVICE_MANAGED Real dx_mrc_z=1e-6;
    AMREX_GPU_DEVICE_MANAGED Real mrc_threshold=1800.0;

    AMREX_GPU_DEVICE_MANAGED Real fs_vel=fourth;
    AMREX_GPU_DEVICE_MANAGED Real fs_p=2*P_NTP;
    AMREX_GPU_DEVICE_MANAGED Real fs_rho=2*one;
    AMREX_GPU_DEVICE_MANAGED Real fs_temp=800.0;
    AMREX_GPU_DEVICE_MANAGED Real fs_Ar=1.396;
    AMREX_GPU_DEVICE_MANAGED Real fs_spec[NUM_SPECIES]={zeroval};
    AMREX_GPU_DEVICE_MANAGED Real catalyst_sites=10;
    AMREX_GPU_DEVICE_MANAGED int chemistry_on=0;

    void initialize_problem()
    {
        ParmParse pp("problem");
        pp.query("free_stream_vel",fs_vel); 
        pp.query("free_stream_pres",fs_p);
        pp.query("free_stream_primaryvapor",fs_Ar);
        pp.query("free_stream_temp",fs_temp);
        pp.query("mrc_threshold",mrc_threshold);
        pp.query("mrc_dx",dx_mrc_x);
        dx_mrc_y=dx_mrc_x;
        dx_mrc_z=dx_mrc_x;
        pp.query("mrc_lo_x",lo_mrc_x);
        pp.query("mrc_lo_y",lo_mrc_y);
        pp.query("mrc_lo_z",lo_mrc_z);
        pp.query("chemistry_on",chemistry_on);
        
        for(int sp=0;sp<NUM_SPECIES;sp++)
        {
            fs_spec[sp]=zeroval;
        }
        fs_spec[Ar_ID]=fs_Ar;

        ParmParse pp1("amr");
        std::string restfile="";
        pp1.query("restart",restfile);
        if(restfile=="")
        {

            std::ifstream infile("spherical_particle_ep50.txt");

            infile>>nx_mrc>>ny_mrc>>nz_mrc;

            mrcdatavec = new Gpu::ManagedVector<Real>;
            mrcdatavec->resize(nx_mrc*ny_mrc*nz_mrc);

            for(int i=0;i<nx_mrc*ny_mrc*nz_mrc;i++)
            {
                infile>>(*mrcdatavec)[i];
            }
            infile.close();
            mrcdata=mrcdatavec->dataPtr();
        }
        fs_rho = mflo_thermo::get_r_from_tpc(fs_temp,fs_p,fs_spec);
    }
}

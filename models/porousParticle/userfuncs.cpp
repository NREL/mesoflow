#include<userfuncs.H>
#include <AMReX_ParmParse.H>
#include <thermo.H>
#include <species.H>

namespace mflo_user_funcs
{
    AMREX_GPU_DEVICE_MANAGED int nx_mrc=0;
    AMREX_GPU_DEVICE_MANAGED int ny_mrc=0;
    AMREX_GPU_DEVICE_MANAGED int nz_mrc=0;
    Gpu::ManagedVector<Real>* mrcdatavec=NULL;
    AMREX_GPU_DEVICE_MANAGED Real *mrcdata=NULL;
    AMREX_GPU_DEVICE_MANAGED Real lo_mrc_x=0.0;
    AMREX_GPU_DEVICE_MANAGED Real lo_mrc_y=0.0;
    AMREX_GPU_DEVICE_MANAGED Real lo_mrc_z=0.0;
    AMREX_GPU_DEVICE_MANAGED Real dx_mrc_x=3e-6;
    AMREX_GPU_DEVICE_MANAGED Real dx_mrc_y=3e-6;
    AMREX_GPU_DEVICE_MANAGED Real dx_mrc_z=3e-6;
    AMREX_GPU_DEVICE_MANAGED Real mrc_threshold=3300.0;
    AMREX_GPU_DEVICE_MANAGED Real pres_left       = 2e5;
    AMREX_GPU_DEVICE_MANAGED Real pres_right      = 1e5;
    AMREX_GPU_DEVICE_MANAGED Real dens_left       = 2.32;
    AMREX_GPU_DEVICE_MANAGED Real temp_left       = 300.0;
    AMREX_GPU_DEVICE_MANAGED Real spec_left[NUM_SPECIES] = {one};

    void initialize_problem()
    {
        Print()<<"Initializing problem\n";
        std::string mrcfile="pine_box.txt";

        ParmParse pp("user");

        pp.query("pres_left",pres_left);
        pp.query("pres_right",pres_right);
        pp.query("dens_left",dens_left);
        pp.query("mrcfile",mrcfile);
        pp.query("mrc_threshold",mrc_threshold);

        spec_left[AIR_ID]=dens_left/mflo_species::molwts[AIR_ID];
        temp_left=mflo_thermo::get_t_from_rpc(dens_left,pres_left,spec_left);

	std::ifstream infile(mrcfile.c_str());

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
}

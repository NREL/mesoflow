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
    AMREX_GPU_DEVICE_MANAGED Real dx_mrc_x=5e-5;
    AMREX_GPU_DEVICE_MANAGED Real dx_mrc_y=5e-5;
    AMREX_GPU_DEVICE_MANAGED Real dx_mrc_z=5e-5;
    AMREX_GPU_DEVICE_MANAGED Real siteconc=70000.0;
    AMREX_GPU_DEVICE_MANAGED Real mrc_threshold=0.5;
    AMREX_GPU_DEVICE_MANAGED Real sponge_zone_dist=20e-3;
    AMREX_GPU_DEVICE_MANAGED Real pres_left       = 2e5;
    AMREX_GPU_DEVICE_MANAGED Real pres_right      = 1e5;
   // AMREX_GPU_DEVICE_MANAGED Real dens_left       = 2.32;
    AMREX_GPU_DEVICE_MANAGED Real temp_left       = 0.0;
    AMREX_GPU_DEVICE_MANAGED Real temp_solid      = 0.0;
    AMREX_GPU_DEVICE_MANAGED Real fs_vel          = 0.0;
    AMREX_GPU_DEVICE_MANAGED Real H2molfrac       = 0.2;
    AMREX_GPU_DEVICE_MANAGED Real H2Omolfrac      = 0.0;
    AMREX_GPU_DEVICE_MANAGED Real COmolfrac       = 0.0;
    AMREX_GPU_DEVICE_MANAGED Real N2molfrac       = 0.0;
    AMREX_GPU_DEVICE_MANAGED Real CO2molfrac      = 0.0;


 //   AMREX_GPU_DEVICE_MANAGED Real spec_left[NUM_SPECIES] = {one};

//   // Species concentrations for the inflow
//   AMREX_GPU_DEVICE_MANAGED Real spec_left[NUM_GAS_SPECIES] = {0.5, 0.3, 0.0, 0.0, 0.2}; // Example values for H2, CO, H2O, CO2, N2


    void initialize_problem()
    {
        Print()<<"Initializing problem\n";
        std::string mrcfile="pine_box.txt";

        ParmParse pp("user");

        pp.query("pres_left",pres_left);
        pp.query("pres_right",pres_right);
        pp.query("temp_left",temp_left);
        pp.query("temp_solid",temp_solid);
      //  pp.query("dens_left",dens_left);
        pp.query("fs_vel",fs_vel);
        pp.query("H2molfrac",H2molfrac);
        pp.query("H2Omolfrac",H2Omolfrac);
        pp.query("COmolfrac",COmolfrac);
        pp.query("N2molfrac",N2molfrac);
        
        pp.query("mrcfile",mrcfile);
        pp.query("mrc_threshold",mrc_threshold);

        
        
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

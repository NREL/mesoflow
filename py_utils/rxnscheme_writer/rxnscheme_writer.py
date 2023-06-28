'''
Reads chemical species, reactions, and rates from csv and produces species.cpp  and transport.H files
Written by Meagan Crowley, NREL 2023
'''
import os
import pandas as pd

#########################
#  Function definitions
#########################

def write_speciescpp(species_df, path):
    '''
    Write species.cpp file using species and properties
    defined in species.csv
    '''
    # Number of rows in dataframe is number of species
    num_species = species_df.shape[0]
    # write header to species.cpp
    with open(f'{path}/species.cpp', 'w+') as f:
        f.write('''#include<species.H>

namespace cflo_species
{
    amrex::Vector<std::string> specnames(NUM_SPECIES);
    AMREX_GPU_DEVICE_MANAGED amrex::Real molwts[NUM_SPECIES]={one};
    AMREX_GPU_DEVICE_MANAGED amrex::Real advect_flags[NUM_SPECIES];
    AMREX_GPU_DEVICE_MANAGED amrex::Real rateconsts[NUM_REACTIONS];
    AMREX_GPU_DEVICE_MANAGED amrex::Real gamma_spec[NUM_GAS_SPECIES+1];
    AMREX_GPU_DEVICE_MANAGED amrex::Real rxnarray[NUM_REACTIONS][NUM_SPECIES];

    void init()
    {''')
        f.close()
    # loop through species and write names, formula to species.cpp
    with open(f'{path}/species.cpp', 'a+') as f:
        for i in range(0,num_species):
            formula = species_df['formula'].iloc[i]
            species_name = species_df['species_name'].iloc[i]
            f.write(f'\n        specnames[{formula}_ID] = "{species_name}";')
        f.close()
    # loop through species and write molecular weights to species.cpp
    with open(f'{path}/species.cpp', 'a+') as f:
        for i in range(0,num_species):
            if i==0:
                f.write('\n\n        //kg/mol')
            formula = species_df['formula'].iloc[i]
            molwt = species_df['molwt'].iloc[i]
            f.write(f"\n        molwts[{formula}"+'_ID] = ' + "%5.5f"%(molwt)+";")
        f.close()

    # loop through species and write gamma for gas species to species.cpp
    with open(f'{path}/species.cpp', 'a+') as f:
        for i in range(0,num_species):
            if i == 0:
                f.write('\n\n        //gamma for each species\n')
            formula = species_df['formula'].iloc[i]
            gamma = species_df['gamma_spec'].iloc[i]
            #print(pd.isnull(gamma))
            if not pd.isnull(gamma):
                f.write(f'        gamma_spec[{formula}_ID] = {gamma};\n')
        f.write('        //background gas\n        gamma_spec[NUM_GAS_SPECIES]  = GAMMA_BG_GAS;\n\n')
        f.close()

    with open(f'{path}/species.cpp', 'a+') as f:
        f.write('''\n\n        for(int sp=0;sp<NUM_SPECIES;sp++)
        {
            advect_flags[sp]=one;
        }''')
        f.close()
    # Loop through species and write advection flags if zero
    with open(f'{path}/species.cpp', 'a+') as f:
        for i in range(0,num_species):
            formula = species_df['formula'].iloc[i]
            advect_flag = species_df['advect_flag'].iloc[i]
            if advect_flag == 0:
                f.write(f'\n        advect_flags[{formula}_ID] = zero;')
        f.write('\n')
        f.close()


def write_rxn_arrays(reactions_df, path):
    # Write reactions to species.cpp as defined in reactions.csv
    # Number of rows in dataframe is number of reactions (forward & backward)
    num_rxns = reactions_df.shape[0]
    species_names =  list(reactions_df.columns)
    with open(f'{path}/species.cpp', 'a+') as f:
        #  write header
        f.write('''        for(int i=0;i<NUM_REACTIONS;i++)
        {
            for(int j=0;j<NUM_SPECIES;j++)
            {
                rxnarray[i][j]=zero;
            }
        }\n''')
        f.close()
    # loop through reaction matrix and write species, stoichiometric coeffs,
    # and rate constants to species.cpp
    with open(f'{path}/species.cpp', 'a+') as f:
        for i in range(0,num_rxns):
            num_species = reactions_df.shape[1] # number of columns
            f.write(f'        //Reaction{i}\n')
            for j in range(1,num_species-1):
                if reactions_df.iat[i,j] != 0:
                    spec_name = species_names[j]
                    stoich = reactions_df.iat[i,j]
                    f.write(f'        rxnarray[{i}][{spec_name}_ID] = {stoich};\n')
            rate = reactions_df['rate'].iloc[i]
            f.write(f'        rateconsts[{i}] = ' + "%5.5f"%(rate)+';\n\n')
        f.close()
    # write footer
    with open(f'{path}/species.cpp', 'a+') as f:
        f.write('''    }    
    void close()
    {
        specnames.clear();
    }
    int find_id(std::string specname)
    {
        int loc=-1;
        auto it=std::find(specnames.begin(),specnames.end(),specname);
        if(it != specnames.end())
        {
            loc=it-specnames.begin();
        }
        return(loc);
    }
}''')
    f.close()


def write_transport(species_df, soldiff, fludiff, path):
    '''
    Read species from species.csv
    Set solid species' transport coefficients to zero in transport.H
    '''
    num_species = species_df.shape[0]
    # Write header
    with open(f'{path}/transport.H', 'w+') as f:
        f.write('''#ifndef _TRANSPORT_H_
#define _TRANSPORT_H_

#include <AMReX_Geometry.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_REAL.H>
#include <species.H>

using namespace amrex;
namespace cflo_chem_transport
{
    AMREX_GPU_DEVICE AMREX_INLINE
        void compute_spec_dcoeff(int i, int j, int k,
                Array4<Real> const& phi,
                Array4<Real> const& dcoeff,
                GpuArray<Real, AMREX_SPACEDIM> prob_lo,
                GpuArray<Real, AMREX_SPACEDIM> dx,
                const Real time)
        {
            Real solid_dif=%s;
            Real fluid_dif=%s;
            Real vfrac=phi(i,j,k,VFRAC_INDX);

            for(int sp=0;sp<NUM_SPECIES;sp++)
            {
                dcoeff(i,j,k,sp) = (one-vfrac)*solid_dif+vfrac*fluid_dif;
            }'''%(soldiff, fludiff))
    f.close()
    with open(f'{path}/transport.H', 'a+') as f:
        for i in range(0,num_species):
            formula = species_df['formula'].iloc[i]
            advect_flag = species_df['advect_flag'].iloc[i]
            if advect_flag == 0:
                f.write(f'\n            dcoeff(i,j,k,{formula}_ID)=zero;')
        f.write('\n')
    f.close()

    # write footer
    with open(f'{path}/transport.H', 'a+') as f:
        f.write('''       }
}
#endif''')
    f.close()


def write_speciesH(species_df, reactions_df, path):
    '''
    Write species.H file using species and properties
    defined in species.csv
    '''
    # Number of rows in dataframe is number of species
    num_species = species_df.shape[0]
    num_rxns = reactions_df.shape[0]
    num_gas_spec = len(species_df[species_df['advect_flag'] != 0])
    first_spec = species_df['formula'].iloc[0]
    # write header to species.H
    with open(f'{path}/species.H', 'w+') as f:
        f.write('''#ifndef _SPECIES_H_
#define _SPECIES_H_

#include<AMReX_REAL.H>
#include<AMReX.H>
#include<string>
#include<AMReX_Vector.H>
#include <AMReX_Geometry.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_Box.H>
#include <cflo_constants.H>''')
        f.write(f"\n\n#define NUM_SPECIES {num_species}\n")
        f.close()

    with open(f'{path}/species.H', 'a+') as f:
        for i in range(0,num_species):
            formula = species_df['formula'].iloc[i]
            f.write(f"#define {formula}_ID {i}\n")
        f.write(f"\n#define NUM_GAS_SPECIES {num_gas_spec}")
        f.write(f"\n#define FIRST_SPEC {first_spec}_ID\n")
        f.write("\n#define BG_GAS_MWT 0.04\n#define GAMMA_BG_GAS 1.67\n\n")
        f.write(f"#define NUM_REACTIONS {num_rxns}\n")
        f.close()

    # write footer to species.H
    with open(f'{path}/species.H', 'a+') as f:
        f.write('''
namespace cflo_species
{
    extern amrex::Vector<std::string> specnames;
    extern AMREX_GPU_DEVICE_MANAGED  amrex::Real molwts[NUM_SPECIES];
    extern AMREX_GPU_DEVICE_MANAGED  amrex::Real gamma_spec[NUM_GAS_SPECIES+1];
    extern AMREX_GPU_DEVICE_MANAGED  amrex::Real rxnarray[NUM_REACTIONS][NUM_SPECIES];
    extern AMREX_GPU_DEVICE_MANAGED amrex::Real rateconsts[NUM_REACTIONS];
    extern AMREX_GPU_DEVICE_MANAGED  amrex::Real advect_flags[NUM_SPECIES];
    void init();
    void close();
    int find_id(std::string specname);
}

#endif''')
        f.close()


def write_thermo(species_df, path):
    first_spec = species_df['formula'].iloc[0]
    with open(f'{path}/thermo.H', 'w+') as f:
        f.write('''
        #ifndef _THERMO_H_
#define _THERMO_H_

#include<AMReX_REAL.H>
#include<AMReX.H>
#include<string>
#include<AMReX_Vector.H>
#include <AMReX_Geometry.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_Box.H>
#include <globalDefines.H>
#include <cflo_constants.H>
#include <species.H>

using namespace amrex;
namespace cflo_thermo
{
    AMREX_GPU_HOST_DEVICE AMREX_INLINE
        Real get_r_from_c(Real spec[NUM_SPECIES],Real bgasconc=zero)
        {

            Real rho = zero;     
            for(int sp=FIRST_SPEC;sp<NUM_GAS_SPECIES;sp++)
            {
                rho += spec[sp]*cflo_species::molwts[sp];   
            }
            rho += bgasconc*BG_GAS_MWT;

            return(rho);
        }
    AMREX_GPU_HOST_DEVICE AMREX_INLINE
        Real get_bgasconc_from_rc(Real rho,Real spec[NUM_SPECIES])
        {
            Real sum_ciMi=0.0;
            for(int sp=FIRST_SPEC;sp<NUM_GAS_SPECIES;sp++)
            {
                sum_ciMi += spec[sp]*cflo_species::molwts[sp];   
            }
            Real c_bg_gas = (rho-sum_ciMi)/BG_GAS_MWT;
            return(c_bg_gas);
        }
    AMREX_GPU_HOST_DEVICE AMREX_INLINE
        Real get_r_from_tpc(Real temp,Real pres,Real spec[NUM_SPECIES])
        {
            Real sum_ci   = zero;
            Real sum_ciMi = zero;
            for(int sp=FIRST_SPEC;sp<NUM_GAS_SPECIES;sp++)
            {
                sum_ci   += spec[sp];   
                sum_ciMi += spec[sp]*cflo_species::molwts[sp];   
            }

            Real c_total = pres/RU/temp;
            Real c_bg_gas = c_total-sum_ci;

            return(sum_ciMi+c_bg_gas*BG_GAS_MWT);
        }

    AMREX_GPU_HOST_DEVICE AMREX_INLINE
        Real get_t_from_rpc(Real rho,Real pres,Real spec[NUM_SPECIES])
        {
            Real sum_ci   = zero;
            Real sum_ciMi = zero;
            for(int sp=FIRST_SPEC;sp<NUM_GAS_SPECIES;sp++)
            {
                sum_ci   += spec[sp];   
                sum_ciMi += spec[sp]*cflo_species::molwts[sp];   
            }

            Real c_bg_gas = (rho-sum_ciMi)/BG_GAS_MWT;
            sum_ci += c_bg_gas;

            //return(pres*BG_GAS_MWT/rho/RU);
            return(pres/(sum_ci*RU));
        }

    AMREX_GPU_HOST_DEVICE AMREX_INLINE
        Real get_p_from_rec(Real rhoe,Real rho,Real spec[NUM_SPECIES])
        {
            //not using coke and sites
            //these are only present in solid

            //include background gas also
            Real cv_spec[NUM_GAS_SPECIES+1]={zero};
            Real T;

            //find background gas concentration
            Real sum_ciMi = zero;
            for(int sp=FIRST_SPEC;sp<NUM_GAS_SPECIES;sp++)
            {
                sum_ciMi += spec[sp]*cflo_species::molwts[sp];   
            }
            Real c_bg_gas = (rho-sum_ciMi)/BG_GAS_MWT;

            for(int sp=FIRST_SPEC;sp<(NUM_GAS_SPECIES+1);sp++)
            {
                cv_spec[sp]=RU/(cflo_species::gamma_spec[sp]-one); //J/K/mol
            }

            Real cv_times_spec=zero;
            for(int sp=FIRST_SPEC;sp<NUM_GAS_SPECIES;sp++)
            {
                cv_times_spec += cv_spec[sp]*spec[sp]; //J/K/m3
            }
            cv_times_spec += cv_spec[NUM_GAS_SPECIES]*c_bg_gas;


            T=rhoe/cv_times_spec;
            Real p=zero;
            for(int sp=FIRST_SPEC;sp<NUM_GAS_SPECIES;sp++)
            {
                p += spec[sp]*RU*T;   
            }
            p += c_bg_gas*RU*T;

            //return(rhoe*(GAMMA_BG_GAS-one));
            return(p);
        }

    AMREX_GPU_HOST_DEVICE AMREX_INLINE
        Real get_gama_from_rpc(Real rho,Real pres,Real spec[NUM_SPECIES])
        {
            Real x_i[NUM_GAS_SPECIES+1]={zero};  //mole frac

            //find background gas concentration
            Real sum_ciMi = zero;
            for(int sp=%s_ID;sp<NUM_GAS_SPECIES;sp++)
            {
                sum_ciMi += spec[sp]*cflo_species::molwts[sp];   
            }
            Real c_bg_gas = (rho-sum_ciMi)/BG_GAS_MWT;


            Real sum_ci=zero;
            for(int sp=FIRST_SPEC;sp<NUM_GAS_SPECIES;sp++)
            {
                sum_ci += spec[sp];
            }
            sum_ci += c_bg_gas;

            for(int sp=FIRST_SPEC;sp<NUM_GAS_SPECIES;sp++)
            {
                x_i[sp]=spec[sp]/sum_ci;
            }
            x_i[NUM_GAS_SPECIES] = c_bg_gas/sum_ci;

            Real term2=zero;
            for(int sp=FIRST_SPEC;sp<(NUM_GAS_SPECIES+1);sp++)
            {
                term2 += x_i[sp]/(cflo_species::gamma_spec[sp]-one);
            }

            Real gama_mix=one+pow(term2,-one);

            //return(GAMMA_BG_GAS);
            return(gama_mix);
        }

    AMREX_GPU_HOST_DEVICE AMREX_INLINE
        Real get_e_from_rpc(Real rho,Real pres,Real spec[NUM_SPECIES])
        {
            Real cv_spec[NUM_GAS_SPECIES+1]={zero};
            Real temp;

            //find background gas concentration
            Real sum_ciMi = zero;
            for(int sp=FIRST_SPEC;sp<NUM_GAS_SPECIES;sp++)
            {
                sum_ciMi += spec[sp]*cflo_species::molwts[sp];   
            }
            Real c_bg_gas = (rho-sum_ciMi)/BG_GAS_MWT;

            for(int sp=FIRST_SPEC;sp<NUM_GAS_SPECIES;sp++)
            {
                cv_spec[sp] = RU/(cflo_species::gamma_spec[sp]-one); //J/K/mol
            }
            cv_spec[NUM_GAS_SPECIES] = RU/(cflo_species::gamma_spec[NUM_GAS_SPECIES]-one); 

            temp=get_t_from_rpc(rho,pres,spec);

            Real sum_cvi_ci=zero;
            for(int sp=FIRST_SPEC;sp<NUM_GAS_SPECIES;sp++)
            {
                sum_cvi_ci += cv_spec[sp]*spec[sp];
            }
            sum_cvi_ci += cv_spec[NUM_GAS_SPECIES]*c_bg_gas;

            //return(pres/rho/(GAMMA_BG_GAS-one));
            return( sum_cvi_ci*temp/rho );
        }

}
#endif

        '''% (first_spec))



def write_userfuncs(species_df, path):
    # This file is highly customizable based on use case
    # new variables can be declared here and called in inputs file
    first_spec = species_df['formula'].iloc[0]
    with open(f'{path}/userfuncs.cpp', 'w+') as f:
        f.write('''
#include<userfuncs.H>
#include <AMReX_ParmParse.H>

namespace cflo_user_funcs
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
    AMREX_GPU_DEVICE_MANAGED Real fs_%s=1.396;
    AMREX_GPU_DEVICE_MANAGED Real fs_spec[NUM_SPECIES]={zero};
    AMREX_GPU_DEVICE_MANAGED Real catalyst_sites=10;
    AMREX_GPU_DEVICE_MANAGED int chemistry_on=0;

    void initialize_problem()
    {
        ParmParse pp("problem");
        pp.query("free_stream_vel",fs_vel); 
        pp.query("free_stream_pres",fs_p);
        pp.query("free_stream_primaryvapor",fs_%s);
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
            fs_spec[sp]=zero;
        }
        fs_spec[%s_ID]=fs_%s;

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
        fs_rho = cflo_thermo::get_r_from_tpc(fs_temp,fs_p,fs_spec);
    }
}
'''%(first_spec, first_spec, first_spec, first_spec))


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

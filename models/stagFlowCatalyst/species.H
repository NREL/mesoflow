#ifndef _SPECIES_H_
#define _SPECIES_H_

#include<AMReX_REAL.H>
#include<AMReX.H>
#include<string>
#include<AMReX_Vector.H>
#include <AMReX_Geometry.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_Box.H>
#include <mflo_constants.H>

#define NUM_SPECIES 6

#define CH4_ID 0
#define O2_ID 1
#define N2_ID 2
#define C2H4_ID 3
#define H2O_ID 4
#define NUM_GAS_SPECIES 5
#define S_ID 5
namespace mflo_species
{
    extern amrex::Vector<std::string> specnames;
    extern AMREX_GPU_DEVICE_MANAGED  amrex::Real molwts[NUM_SPECIES];
    extern AMREX_GPU_DEVICE_MANAGED  amrex::Real advect_flags[NUM_SPECIES];
    void init();
    void close();
    int find_id(std::string specname);
}
#endif

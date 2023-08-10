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

#define NUM_SPECIES 3
#define CH4_ID 0
#define CH3_ID 1
#define H_ID 2

#define NUM_GAS_SPECIES 3
#define FIRST_SPEC CH4_ID

#define BG_GAS_MWT 0.04
#define GAMMA_BG_GAS 1.67

#define NUM_REACTIONS 105

namespace mflo_species
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

#endif
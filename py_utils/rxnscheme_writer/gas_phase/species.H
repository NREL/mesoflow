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

#define NUM_SPECIES 25
#define Ar_ID 0
#define CH4_ID 1
#define C2H6_ID 2
#define C2H4_ID 3
#define C2H2_ID 4
#define C3H8_ID 5
#define C3H6_ID 6
#define O2_ID 7
#define H2O_ID 8
#define H2O2_ID 9
#define CH2O_ID 10
#define CO_ID 11
#define CO2_ID 12
#define H2_ID 13
#define HO2_ID 14
#define O_ID 15
#define OH_ID 16
#define CH3_ID 17
#define C2H5_ID 18
#define C2H3_ID 19
#define C3H7_ID 20
#define CH3O_ID 21
#define CHO_ID 22
#define H_ID 23
#define S1_ID 24

#define NUM_GAS_SPECIES 24
#define FIRST_SPEC Ar_ID

#define BG_GAS_MWT 0.04
#define GAMMA_BG_GAS 1.67

#define NUM_REACTIONS 78

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
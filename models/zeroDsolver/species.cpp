#include<species.H>

namespace mflo_species
{
    amrex::Vector<std::string> specnames(NUM_SPECIES);
    AMREX_GPU_DEVICE_MANAGED amrex::Real molwts[NUM_SPECIES]={one};
    AMREX_GPU_DEVICE_MANAGED amrex::Real advect_flags[NUM_SPECIES];
    AMREX_GPU_DEVICE_MANAGED amrex::Real rateconsts[NUM_REACTIONS];
    AMREX_GPU_DEVICE_MANAGED amrex::Real gamma_spec[NUM_GAS_SPECIES+1];
    AMREX_GPU_DEVICE_MANAGED amrex::Real rxnarray[NUM_REACTIONS][NUM_SPECIES];

    void init()
    {
        specnames[CH4_ID] = "CH4";
        specnames[CH3_ID] = "CH3";
        specnames[H_ID] = "H";

        //kg/mol
        molwts[CH4_ID] = 0.01604;
        molwts[CH3_ID] = 0.01503;
        molwts[H_ID] = 0.00101;

        //gamma for each species
        gamma_spec[CH4_ID] = 1.13;
        gamma_spec[CH3_ID] = 1.16;
        gamma_spec[H_ID] = 1.67;

        //background gas
        gamma_spec[NUM_GAS_SPECIES]  = GAMMA_BG_GAS;



        for(int sp=0;sp<NUM_SPECIES;sp++)
        {
            advect_flags[sp]=one;
        }
       
        for(int i=0;i<NUM_REACTIONS;i++)
        {
            for(int j=0;j<NUM_SPECIES;j++)
            {
                rxnarray[i][j]=zeroval;
            }
        }
        //Reaction0 CH4 --> CH3 + H
        rxnarray[0][CH4_ID] = -1;
        rxnarray[0][CH3_ID] = 1;
        rxnarray[0][H_ID] = 1;
        rateconsts[0] = 30;


    }    
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
}
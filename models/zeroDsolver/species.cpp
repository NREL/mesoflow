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
        specnames[O2_ID] = "O2";
        specnames[HO2_ID] = "HO2";

        //kg/mol
        molwts[CH4_ID] = 0.01604;
        molwts[CH3_ID] = 0.01503;
        molwts[H_ID] = 0.00101;
        molwts[O2_ID] = 0.03200;
        molwts[HO2_ID] = 0.03301;

        //gamma for each species
        gamma_spec[CH4_ID] = 1.13;
        gamma_spec[CH3_ID] = 1.16;
        gamma_spec[H_ID] = 1.67;
        gamma_spec[O2_ID] = 1.31;
        gamma_spec[HO2_ID] = 1.21;

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
        //Reaction1 Forward CH4 + O2 --> CH3 + HO2
        rxnarray[0][CH4_ID] = -1;
        rxnarray[0][O2_ID] = -1;
        rxnarray[0][CH3_ID] = 1;
        rxnarray[0][HO2_ID] = 1;
        rateconsts[0] = 0.00124;

        //Reaction1 Backward CH4 + O2 <-- CH3 + HO2
        rxnarray[1][CH4_ID] = 1;
        rxnarray[1][O2_ID] = 1;
        rxnarray[1][CH3_ID] = -1;
        rxnarray[1][HO2_ID] = -1;
        rateconsts[1] = 1.24E-03;
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
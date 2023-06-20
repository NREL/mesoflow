#include<species.H>

namespace cflo_species
{
    amrex::Vector<std::string> specnames(NUM_SPECIES);
    AMREX_GPU_DEVICE_MANAGED amrex::Real molwts[NUM_SPECIES]={one};
    AMREX_GPU_DEVICE_MANAGED amrex::Real advect_flags[NUM_SPECIES];
    AMREX_GPU_DEVICE_MANAGED amrex::Real rateconsts[NUM_REACTIONS];
    AMREX_GPU_DEVICE_MANAGED amrex::Real gamma_spec[NUM_GAS_SPECIES+1];
    AMREX_GPU_DEVICE_MANAGED amrex::Real rxnarray[NUM_REACTIONS][NUM_SPECIES];

    void init()
    {
        specnames[H2_ID] = "H2";
        specnames[H2ADS1_ID] = "H2ADS1";
        specnames[S1_ID] = "S1";

        //kg/mol
        molwts[H2_ID] = 0.00200;
        molwts[H2ADS1_ID] = 0.04400;
        molwts[S1_ID] = 0.00000;

        //gamma for each species
        gamma_spec[H2_ID] = 1.4;
        //background gas
        gamma_spec[NUM_GAS_SPECIES]  = GAMMA_BG_GAS;



        for(int sp=0;sp<NUM_SPECIES;sp++)
        {
            advect_flags[sp]=one;
        }
        advect_flags[H2ADS1_ID] = zero;
        advect_flags[S1_ID] = zero;
        for(int i=0;i<NUM_REACTIONS;i++)
        {
            for(int j=0;j<NUM_SPECIES;j++)
            {
                rxnarray[i][j]=zero;
            }
        }
        //Reaction0
        rxnarray[0][H2_ID] = -1;
        rxnarray[0][H2ADS1_ID] = 1;
        rxnarray[0][S1_ID] = -1;
        rateconsts[0] = 165675000000.00000;

        //Reaction1
        rxnarray[1][H2_ID] = 1;
        rxnarray[1][H2ADS1_ID] = -1;
        rxnarray[1][S1_ID] = 1;
        rateconsts[1] = 16567452390.00000;

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
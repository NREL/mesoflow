#include<species.H>

namespace mflo_species
{
    amrex::Vector<std::string> specnames(NUM_SPECIES);
    AMREX_GPU_DEVICE_MANAGED amrex::Real advect_flags[NUM_SPECIES]={one};
    AMREX_GPU_DEVICE_MANAGED amrex::Real molwts[NUM_SPECIES]={one};

    void init()
    {
        specnames[Fe2O3_ID] = "Fe2O3";
        specnames[H2_ID] = "H2";
        specnames[FeO_ID] = "FeO";
        specnames[Fe_ID] = "Fe";
        specnames[CO_ID] = "CO";
        specnames[H2O_ID] = "H2O";
        specnames[CO2_ID] = "CO2";
        specnames[N2_ID] = "N2";
        
        //kg/mol
        molwts[Fe2O3_ID] = 0.1597;
        molwts[H2_ID] = 0.002;
        molwts[FeO_ID] = 0.0718;
        molwts[Fe_ID] = 0.0558;
        molwts[CO_ID] = 0.028;
        molwts[H2O_ID] = 0.018;
        molwts[CO2_ID] = 0.044;
        molwts[N2_ID] = 0.028;


        for(int sp=0;sp<NUM_SPECIES;sp++)
        {
            advect_flags[sp]=one;
        }
        
        
        advect_flags[FeO_ID]=0;
        advect_flags[Fe_ID]=0;
        advect_flags[Fe2O3_ID]=0;

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

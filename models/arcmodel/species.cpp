#include<species.H>

namespace mflo_species
{
    amrex::Vector<std::string> specnames(NUM_SPECIES);
    AMREX_GPU_DEVICE_MANAGED amrex::Real advect_flags[NUM_SPECIES]={one};
    AMREX_GPU_DEVICE_MANAGED amrex::Real molwts[NUM_SPECIES]={one};

    void init()
    {
        specnames[H_ID]="H";
        specnames[H2_ID]="H2";
        specnames[H2O_ID]="H2O";
        specnames[AR_ID]="AR";
        specnames[FEO_ID]="FEO";
        specnames[FE2O3_ID]="FE2O3";
        specnames[FE_ID]="FE";
        
        //kg/mol
        molwts[H_ID]=0.001;
        molwts[H2_ID]=0.002;
        molwts[H2O_ID]=0.018;
        molwts[AR_ID]=0.04;
    
        molwts[FEO_ID]=0.072;
        molwts[FE2O3_ID]=0.15969;
        molwts[FE_ID]=0.056;

        for(int sp=0;sp<NUM_GAS_SPECIES;sp++)
        {
            advect_flags[sp]=0.0;
        }
        advect_flags[FEO_ID]=0;
        advect_flags[FE_ID]=0;
        advect_flags[FE2O3_ID]=0;
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

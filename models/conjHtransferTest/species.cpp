#include<species.H>
#include<mflo_constants.H>

namespace mflo_species
{
    amrex::Vector<std::string> specnames(NUM_SPECIES);
    AMREX_GPU_DEVICE_MANAGED amrex::Real molwts[NUM_SPECIES]={zeroval};
    AMREX_GPU_DEVICE_MANAGED amrex::Real advect_flags[NUM_SPECIES]={one};

    void init()
    {
        specnames[AIR_ID]="air";
        specnames[SOLSPEC_ID]="solidmat";

        for(int sp=0;sp<NUM_SPECIES;sp++)
        {
            advect_flags[sp]=one;
        }
        advect_flags[SOLSPEC_ID]=zeroval;
        advect_flags[AIR_ID]=zeroval;

        //kg/mol
        molwts[AIR_ID]=0.028; 
        molwts[SOLSPEC_ID]=0.07; 
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

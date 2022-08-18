#include<species.H>
#include<mflo_constants.H>

namespace mflo_species
{
    amrex::Vector<std::string> specnames(NUM_SPECIES);
    AMREX_GPU_DEVICE_MANAGED amrex::Real molwts[NUM_SPECIES]={zeroval};
    AMREX_GPU_DEVICE_MANAGED amrex::Real advect_flags[NUM_SPECIES]={one};

    void init()
    {
        specnames[S1_ID]="S1";
        specnames[S2_ID]="S2";
        specnames[SITE1_ID]="site";

        for(int sp=0;sp<NUM_SPECIES;sp++)
        {
            advect_flags[sp]=one;
        }
        advect_flags[SITE1_ID]=zeroval;

        //kg/mol
        molwts[S1_ID]=0.032; 
        molwts[S2_ID]=0.028; 
        molwts[SITE1_ID]=0.001;
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

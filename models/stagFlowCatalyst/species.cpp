#include<species.H>

namespace mflo_species
{
    amrex::Vector<std::string> specnames(NUM_SPECIES);
    AMREX_GPU_DEVICE_MANAGED amrex::Real advect_flags[NUM_SPECIES]={one};
    AMREX_GPU_DEVICE_MANAGED amrex::Real molwts[NUM_SPECIES]={one};

    void init()
    {
        specnames[CH4_ID]="CH4";
        specnames[O2_ID]="O2";
        specnames[N2_ID]="N2";
        specnames[C2H4_ID]="C2H4";
        specnames[H2O_ID]="H2O";
        specnames[S_ID]="sites";
        
        //kg/mol
        molwts[CH4_ID]  =0.016;
        molwts[O2_ID]=0.032;
        molwts[N2_ID]=0.028;
        molwts[C2H4_ID]=0.028;
        molwts[H2O_ID]=0.018;

        molwts[S_ID]=0.001;

        advect_flags[S_ID]=0;
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

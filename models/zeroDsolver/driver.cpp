#include<species.H>
#include <AMReX_ParmParse.H>

AMREX_GPU_DEVICE AMREX_INLINE







void compute_spec_source(amrex::Real spec[NUM_SPECIES],
                         amrex::Real specsource[NUM_SPECIES],
                         amrex::Real time)
{
    for(int sp=0;sp<NUM_SPECIES;sp++)
    {
        specsource[sp]=0.0;
    }
    amrex::Real concmult;
    amrex::Real stoichcoeff;
    amrex::Real reactrate;

    for(int r=0;r<NUM_REACTIONS;r++)
    {
        concmult=1.0;
        for(int sp=0;sp<NUM_SPECIES;sp++)
        {
            if(mflo_species::rxnarray[r][sp] < 0.0)
            {
                stoichcoeff=amrex::Math::abs(mflo_species::rxnarray[r][sp]);      
                concmult*=std::pow(spec[sp],stoichcoeff);
            }
        }
        reactrate=mflo_species::rateconsts[r]*concmult;

        //update production/removal
        for(int sp=0;sp<NUM_SPECIES;sp++)
        {
            if(amrex::Math::abs(mflo_species::rxnarray[r][sp]) > 0.0)
            {
                specsource[sp] += mflo_species::rxnarray[r][sp]*reactrate;
            }
        }
	// keep primary vapor concentration constant, simulating continuous feed
	// comment out if you don't want continuous feed of primary vapor
<<<<<<< HEAD
     //   specsource[Ar_ID] = 0.0;
     //   specsource[CH4_ID] = 0.0;
     //   specsource[O2_ID] = 0.0;
=======
        //specsource[Ar_ID] = 0.0;
        //specsource[CH4_ID] = 0.0;
        //specsource[O2_ID] = 0.0;
>>>>>>> f2b75d105241e14be772bb8470b5c5cc53e0e625
    }

    /*for(int sp=0;sp<NUM_SPECIES;sp++)
    {
        amrex::Print()<<"specsource["<<sp<<"]:"<<specsource[sp]<<"\n";
    }*/
}

int main(int argc,char *argv[])
{
    amrex::Initialize(argc, argv);
    {
        mflo_species::init();
        amrex::Real finaltime=1.0;
        amrex::Real dt=1e-8;
        amrex::Real output_dt=1e-4;
        amrex::Real curtime=0.0;
        amrex::Real output_time=0.0;

        amrex::ParmParse pp;
        pp.query("finaltime",finaltime);
        pp.query("dt",dt);
        pp.query("output_dt",output_dt);

        amrex::Real spec_old[NUM_SPECIES]={0.0};
        amrex::Real spec_new[NUM_SPECIES]={0.0};
        amrex::Real specsource[NUM_SPECIES]={0.0};

        //initial conditions
        //spec_new[Ar_ID]=35.7; // primary vapor concentration
        spec_new[CH4_ID]=7.52; // primary vapor concentration
        spec_new[O2_ID]=3.76; // primary vapor concentration
        //spec_new[S1_ID]= 1000.0;//0.0000114; // active site concentration
        
	// This appends to quants.dat if it exists
	//
	// remove or rename quants.dat file if starting a new run
        std::string quantsFile = "quants.dat";
        amrex::PrintToFile(quantsFile)<<curtime<<"\t";
        for(int sp=0;sp<NUM_SPECIES;sp++)
        {
            amrex::PrintToFile(quantsFile)<<spec_new[sp]<<"\t";
        }
        amrex::PrintToFile(quantsFile)<<"\n";

        while(curtime<finaltime)
        {
            for(int sp=0;sp<NUM_SPECIES;sp++)
            {
                spec_old[sp]=spec_new[sp];
            }

            compute_spec_source(spec_old,specsource,curtime);

            //advance by 0.5*dt
            for(int sp=0;sp<NUM_SPECIES;sp++)
            {
                spec_new[sp]=spec_old[sp]+0.5*dt*specsource[sp];
            }

            compute_spec_source(spec_new,specsource,curtime+0.5*dt);

            //advance by dt
            for(int sp=0;sp<NUM_SPECIES;sp++)
            {
                spec_new[sp]=spec_old[sp]+dt*specsource[sp];
            }

            output_time+=dt;
            curtime+=dt;

            if(output_time > output_dt)
            {
                amrex::Print()<<"time:"<<curtime<<"\n";
                amrex::PrintToFile(quantsFile)<<curtime<<"\t";
                for(int sp=0;sp<NUM_SPECIES;sp++)
                {
                    amrex::PrintToFile(quantsFile)<<spec_new[sp]<<"\t";
                }
                amrex::PrintToFile(quantsFile)<<"\n";

                output_time=0.0;
            }
        }
    }

    amrex::Finalize();
}

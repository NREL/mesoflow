#include <iostream>

#include <AMReX.H>
#include <AMReX_BLProfiler.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_ParmParse.H>

#include <mflo.H>
#include <userfuncs.H>

using namespace amrex;


int main(int argc, char* argv[])
{
    amrex::Initialize(argc, argv);

    // timer for profiling
    BL_PROFILE_VAR("main()", pmain);

    // wallclock time
    const Real strt_total = amrex::second();
    

    {
        mflo_species::init();
        // constructor - reads in parameters from inputs file
        //             - sizes multilevel arrays and data structures

        mflo_user_funcs::initialize_problem();

        mflo mflo_obj;

        // initialize AMR data
        mflo_obj.InitData();
    
        ParmParse pp("mflo");
        int split_chemistry=0;
        pp.query("split_chemistry",split_chemistry);

        const Real strt_evolve = amrex::second();
        // advance solution to final time
        if(split_chemistry==0)
        {
            mflo_obj.EvolveAMR(mflo_obj.stop_time);
        }
        else
        {
            Real t_ss,t_react,dt_react,dt_react_rk,dt_ss;
            pp.get("steady_flow_time",t_ss);
            pp.get("react_final_time",t_react);
            pp.get("react_increment_time",dt_react);
            pp.get("react_time_step",dt_react_rk);
            pp.get("flow_coupling_time",dt_ss);
    
           mflo_obj.Evolve_split(t_ss,t_react,dt_react,dt_react_rk,dt_ss);
        }
        Real end_evolve = amrex::second() - strt_evolve;

        // wallclock time
        Real end_total = amrex::second() - strt_total;

        // print wallclock time
        ParallelDescriptor::ReduceRealMax(
            end_total, ParallelDescriptor::IOProcessorNumber());
        // print wallclock time
        ParallelDescriptor::ReduceRealMax(
            end_evolve, ParallelDescriptor::IOProcessorNumber());
        if (mflo_obj.Verbose()) 
        {
            amrex::Print() << "\nTotal Time: " << end_total << '\n';
            amrex::Print() << "\nEvolve Time: " << end_evolve << '\n';
        }
    }

    // destroy timer for profiling
    BL_PROFILE_VAR_STOP(pmain);

    amrex::Finalize();
}

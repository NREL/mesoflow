#ifndef _CHEMISTRY_H_
#define _CHEMISTRY_H_

#include <AMReX_Geometry.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_REAL.H>
#include <species.H>
#include <globalDefines.H>

using namespace amrex;
namespace mflo_chem_reactions
{
    AMREX_GPU_DEVICE AMREX_INLINE
        void compute_spec_source(int i, int j, int k,
                Array4<Real> const& phi,
                Array4<Real> const& specsource,
                GpuArray<Real, AMREX_SPACEDIM> prob_lo,
                GpuArray<Real, AMREX_SPACEDIM> dx,
                const Real time)
        {
         
            
                        for(int sp=0;sp<NUM_SPECIES;sp++)
                        {
                            specsource(i,j,k,FLO_NVARS+sp) = zeroval;
                        }
            
            
            
            amrex::Real rg=4e-4; //grain radius
            amrex::Real twiddle=1e-10;
            amrex::Real vfrac=phi(i,j,k,VFRAC_INDX);
            amrex::Real totalconc = (1.0-vfrac)*mflo_user_funcs::siteconc+twiddle;
            
            

            // Reaction 1: Fe2O3 + H2 -> 2FeO + H2
            amrex::Real kreac1 = 80.0*std::exp(-66516.0/phi(i,j,k,TEMP_INDX))*(3.0/rg); // Example reaction rate constant
            amrex::Real reactrate1 = kreac1 * phi(i, j, k, FLO_NVARS + Fe2O3_ID)/totalconc * phi(i, j, k, FLO_NVARS + H2_ID);
            specsource(i, j, k, FLO_NVARS + Fe2O3_ID) -= reactrate1;
            specsource(i, j, k, FLO_NVARS + H2_ID) -= reactrate1;
            specsource(i, j, k, FLO_NVARS + FeO_ID) += 2.0 * reactrate1;
            specsource(i, j, k, FLO_NVARS + H2_ID) += reactrate1;

            // Reaction 2: FeO + H2 -> Fe + H2O
            amrex::Real kreac2 = 2857.34*std::exp(-11723.0/phi(i,j,k,TEMP_INDX))*(3.0/rg); // Example reaction rate constant
            amrex::Real reactrate2 = kreac2 *phi(i,j,k,FLO_NVARS+FeO_ID)/totalconc*phi(i,j,k,FLO_NVARS+H2_ID);
            specsource(i, j, k, FLO_NVARS + FeO_ID) -= reactrate2;
            specsource(i, j, k, FLO_NVARS + H2_ID) -= reactrate2;
            specsource(i, j, k, FLO_NVARS + Fe_ID) += reactrate2;
            specsource(i, j, k, FLO_NVARS + H2O_ID) += reactrate2;

            // Reaction 3: Fe2O3 + CO -> 2FeO + CO2
            amrex::Real kreac3 = 25.0*std::exp(-73674.0/phi(i,j,k,TEMP_INDX))*(3.0/rg); // Example reaction rate constant
            amrex::Real reactrate3 = kreac3 * phi(i, j, k, FLO_NVARS + Fe2O3_ID)/totalconc * phi(i, j, k, FLO_NVARS + CO_ID);
            specsource(i, j, k, FLO_NVARS + Fe2O3_ID) -= reactrate3;
            specsource(i, j, k, FLO_NVARS + CO_ID) -= reactrate3;
            specsource(i, j, k, FLO_NVARS + FeO_ID) += 2.0 * reactrate3;
            specsource(i, j, k, FLO_NVARS + CO2_ID) += reactrate3;

            // Reaction 4: FeO + CO -> Fe + CO2
            amrex::Real kreac4 = 17.0*std::exp(-69488.0/phi(i,j,k,TEMP_INDX))*(3.0/rg); // Example reaction rate constant
            amrex::Real reactrate4 = kreac4 * phi(i, j, k, FLO_NVARS + FeO_ID)/totalconc * phi(i, j, k, FLO_NVARS + CO_ID);
            specsource(i, j, k, FLO_NVARS + FeO_ID) -= reactrate4;
            specsource(i, j, k, FLO_NVARS + CO_ID) -= reactrate4;
            specsource(i, j, k, FLO_NVARS + Fe_ID) += reactrate4;
            specsource(i, j, k, FLO_NVARS + CO2_ID) += reactrate4;

            // Reaction 5: CO + H2O -> CO2 + H2
            amrex::Real kreac5 = 1400.0*std::exp(-44895.0/phi(i,j,k,TEMP_INDX))*(3.0/rg); // Example reaction rate constant
            amrex::Real reactrate5 = kreac5 * phi(i, j, k, FLO_NVARS + CO_ID)/totalconc * phi(i, j, k, FLO_NVARS + H2O_ID);
            specsource(i, j, k, FLO_NVARS + CO_ID) -= reactrate5;
            specsource(i, j, k, FLO_NVARS + H2O_ID) -= reactrate5;
            specsource(i, j, k, FLO_NVARS + CO2_ID) += reactrate5;
            specsource(i, j, k, FLO_NVARS + H2_ID) += reactrate5;

            
            //Reaction 6:CO + H2O <-> CO2 + H2
            // Forward and Backward Reaction Rates
            
            amrex::Real kf = 93.32e6 * std::exp(-128200.0 / phi(i, j, k, TEMP_INDX)) * (3.0 / rg);
            amrex::Real kb = 18.27 * std::exp(137300.0 / phi(i, j, k, TEMP_INDX)) * (3.0 / rg);

            amrex::Real reactrate_f = kf * phi(i, j, k, FLO_NVARS + CO_ID) * phi(i, j, k, FLO_NVARS + H2O_ID) / totalconc;
            amrex::Real reactrate_b = kb * phi(i, j, k, FLO_NVARS + CO2_ID) * phi(i, j, k, FLO_NVARS + H2_ID) / totalconc;

            amrex::Real reactrate_net = reactrate_f - reactrate_b;

            specsource(i, j, k, FLO_NVARS + CO_ID) -= reactrate_net;
            specsource(i, j, k, FLO_NVARS + H2O_ID) -= reactrate_net;
            specsource(i, j, k, FLO_NVARS + CO2_ID) += reactrate_net;
            specsource(i, j, k, FLO_NVARS + H2_ID) += reactrate_net;
            
            
        }
}
#endif

#ifndef _CHEMISTRY_H_
#define _CHEMISTRY_H_

#include <AMReX_Geometry.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_REAL.H>
#include <species.H>
#include <globalDefines.H>
#include <userfuncs.H>

using namespace amrex;
namespace cflo_chem_reactions
{
    AMREX_GPU_DEVICE AMREX_INLINE
        void compute_spec_source(int i, int j, int k,
                Array4<Real> const& phi,
                Array4<Real> const& specsource,
                GpuArray<Real, AMREX_SPACEDIM> prob_lo,
                GpuArray<Real, AMREX_SPACEDIM> dx,
                const Real time)
        {
            Real wdot[NUM_SPECIES]={zero};
            amrex::Real concmult;
            amrex::Real stoichcoeff;
            amrex::Real reactrate;

	    if(cflo_user_funcs::chemistry_on){
		    for(int r=0;r<NUM_REACTIONS;r++)
		    {
			    concmult=1.0;
			    for(int sp=0;sp<NUM_SPECIES;sp++)
			    {
				    if(cflo_species::rxnarray[r][sp] < 0.0)
				    {
					    stoichcoeff=amrex::Math::abs(cflo_species::rxnarray[r][sp]);       
					    concmult*=std::pow(phi(i,j,k,CFLO_NVARS+sp),stoichcoeff);
				    }
			    }
			    reactrate=cflo_species::rateconsts[r]*concmult;

			    //update production/removal
			    for(int sp=0;sp<NUM_SPECIES;sp++)
			    {
				    if(amrex::Math::abs(cflo_species::rxnarray[r][sp]) > 0.0)
				    {
					    wdot[sp] += cflo_species::rxnarray[r][sp]*reactrate;
				    }
			    }
		    }
	    }
	    for(int sp=0;sp<NUM_SPECIES;sp++)
	    {
		    specsource(i,j,k,CFLO_NVARS+sp) += wdot[sp];
	    }
	}
}

#endif

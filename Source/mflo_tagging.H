#ifndef TAGGING_H
#define TAGGING_H

#include <AMReX_Array4.H>

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void state_based_refinement(
    int i,
    int j,
    int k,
    amrex::Array4<char> const& tag,
    amrex::Array4<amrex::Real const> const& state,
    amrex::Real* refine_phi,
    int* refine_phi_comps,
    int ntagvars,
    char tagval)
{

    for (int c = 0; c < ntagvars; c++) 
    {
        int comp = refine_phi_comps[c];
        if (fabs(state(i, j, k, comp)) > refine_phi[c])
        {
            tag(i, j, k) = tagval;
        }
    }
}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void stategrad_based_refinement(
    int i,
    int j,
    int k,
    amrex::Array4<char> const& tag,
    amrex::Array4<amrex::Real const> const& state,
    amrex::Real* refine_phigrad,
    int* refine_phi_comps,
    int ntagvars,
    char tagval)
{

    for (int c = 0; c < ntagvars; c++) 
    {
        int comp = refine_phi_comps[c];

        amrex::Real gradxplus =
            std::abs(state(i + 1, j, k, comp) - state(i, j, k, comp));
        amrex::Real gradxminus =
            std::abs(state(i, j, k, comp) - state(i - 1, j, k, comp));
        amrex::Real gradx = amrex::max(gradxplus, gradxminus);

        amrex::Real gradyplus =
            std::abs(state(i, j + 1, k, comp) - state(i, j, k, comp));
        amrex::Real gradyminus =
            std::abs(state(i, j, k, comp) - state(i, j - 1, k, comp));
        amrex::Real grady = amrex::max(gradyplus, gradyminus);

        amrex::Real gradzplus =
            std::abs(state(i, j, k + 1, comp) - state(i, j, k, comp));
        amrex::Real gradzminus =
            std::abs(state(i, j, k, comp) - state(i, j, k - 1, comp));
        amrex::Real gradz = amrex::max(gradzplus, gradzminus);

        if (amrex::max(gradx, grady, gradz) > refine_phigrad[c]) 
        {
            tag(i, j, k) = tagval;
        }
    }
}

#endif

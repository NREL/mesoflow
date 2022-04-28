#ifndef _Adv_3d_H_
#define _Adv_3d_H_

#include <AMReX_Box.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_Geometry.H>

using namespace amrex;

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void update_residual(
    int i,
    int j,
    int k,
    int n,
    Array4<Real> const& dsdt,
    Array4<Real> const& source,
    Array4<Real> const& phi,
    AMREX_D_DECL(
        Array4<Real> const& flxx,
        Array4<Real> const& flxy,
        Array4<Real> const& flxz),
    const GpuArray<Real, AMREX_SPACEDIM>& dx)
{
    // remember, we are solve dudt + del.F = 0
    dsdt(i, j, k, n) = (flxx(i, j, k, n) - flxx(i + 1, j, k, n)) / dx[0] +
                       (flxy(i, j, k, n) - flxy(i, j + 1, k, n)) / dx[1] +
                       (flxz(i, j, k, n) - flxz(i, j, k + 1, n)) / dx[2] +
                       source(i, j, k, n);

#ifdef IMM_BOUNDARY
    //only mask fluid quantities
    if(n < FLO_NVARS)
    {
        dsdt(i, j, k, n) *= phi(i, j, k, VFRAC_INDX);
    }
#endif
}
#endif

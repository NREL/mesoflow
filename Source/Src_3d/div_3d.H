#ifndef _Adv_3d_H_
#define _Adv_3d_H_

#include <AMReX_Box.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_Geometry.H>
#include <species.H>

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
    const GpuArray<Real, AMREX_SPACEDIM>& dx,
    int spec_in_solid,
    int conj_ht)
{
    // remember, we are solve dudt + del.F = 0
    dsdt(i, j, k, n) = (flxx(i, j, k, n) - flxx(i + 1, j, k, n)) / dx[0] +
                       (flxy(i, j, k, n) - flxy(i, j + 1, k, n)) / dx[1] +
                       (flxz(i, j, k, n) - flxz(i, j, k + 1, n)) / dx[2] +
                       source(i, j, k, n);

    int nvars_to_nullify=(spec_in_solid==1)?FLO_NVARS:FLO_NVARS+NUM_SPECIES;

    //check for rhoe because we may have conjugate heat transfer
    if(n!=RHOE_INDX)
    {
        if(n < nvars_to_nullify)
        {
            if(phi(i,j,k,VFRAC_INDX) < one)
            {    
                dsdt(i, j, k, n)=zeroval;
            }
        }
    }
    //if n is rhoe
    else
    {
        //if no conjugate heat transfer
        //null the dsdt
        if(!conj_ht)
        {
            if(phi(i,j,k,VFRAC_INDX) < one)
            {    
                dsdt(i, j, k, n)=zeroval;
            }
        }
        else
        {
            //when conjugate heat transfer is there
            //null the dsdt in fractional cells
            //as these are recalculated from cut cells
            if(phi(i,j,k,VFRAC_INDX) > zeroval && phi(i,j,k,VFRAC_INDX) < one)
            {
                dsdt(i, j, k, n)=zeroval;
            }
        }
    }

    if(dsdt(i,j,k,n)!=dsdt(i,j,k,n))
    {
        amrex::Print()<<"n="<<n<<"\t"<<RHOE_INDX<<"\t"<<i<<"\t"<<j<<"\t"<<k<<"\t"<<phi(i,j,k,VFRAC_INDX)<<"\n";
        amrex::Abort("NaN detected in dsdt\n");
    }
}
#endif

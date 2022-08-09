#ifndef _PROB_H_
#define _PROB_H_

#include <AMReX_Box.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_Geometry.H>
#include <userfuncs.H>
#include <globalDefines.H>
#include <mflo_constants.H>
#include <thermo.H>

using namespace amrex;

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void initdata(Box const& bx, Array4<Real> const& phi, GeometryData const& geomdata)
{
    const auto lo = lbound(bx);
    const auto hi = ubound(bx);    

    const Real* AMREX_RESTRICT prob_lo = geomdata.ProbLo();
    const Real* AMREX_RESTRICT prob_hi = geomdata.ProbHi();
    const Real* AMREX_RESTRICT dx      = geomdata.CellSize();

    Real xlen,ylen,zlen,maxlen;
    Real spec[NUM_SPECIES];
    spec[AIR_ID]=one/mflo_species::molwts[AIR_ID];
    int dir=0;

    xlen=prob_hi[XDIR]-prob_lo[XDIR];
    ylen=prob_hi[YDIR]-prob_lo[YDIR];
    zlen=prob_hi[ZDIR]-prob_lo[ZDIR];

    maxlen=std::max(xlen,std::max(ylen,zlen));

    if(maxlen==xlen) dir=0;
    if(maxlen==ylen) dir=1;
    if(maxlen==zlen) dir=2;


#ifdef _OPENMP
#pragma omp parallel for collapse(2) if (GPU::notInLaunchRegion)
#endif
    for (int k = lo.z; k <= hi.z; ++k) 
    {
        for (int j = lo.y; j <= hi.y; ++j) 
        {
            Real z = prob_lo[ZDIR] + (half+k) * dx[ZDIR];
            Real y = prob_lo[YDIR] + (half+j) * dx[YDIR];
            AMREX_PRAGMA_SIMD
            for (int i = lo.x; i <= hi.x; ++i) 
            {
                Real x = prob_lo[XDIR] + (half+i) * dx[XDIR]; 

                phi(i,j,k,PRES_INDX) = one;
                phi(i,j,k,DENS_INDX) = one;
                phi(i,j,k,VELX_INDX) = zeroval;
                phi(i,j,k,VELY_INDX) = zeroval;
                phi(i,j,k,VELZ_INDX) = zeroval;
                phi(i,j,k,TEMP_INDX) = mflo_thermo::get_t_from_rpc(one,one,spec);

                phi(i,j,k,VELX_INDX+dir)=tenth;

                phi(i,j,k,RHO_INDX)   = phi(i,j,k,DENS_INDX);
                phi(i,j,k,RHOU_INDX)  = phi(i,j,k,DENS_INDX)*phi(i,j,k,VELX_INDX);
                phi(i,j,k,RHOV_INDX)  = phi(i,j,k,DENS_INDX)*phi(i,j,k,VELY_INDX);
                phi(i,j,k,RHOW_INDX)  = phi(i,j,k,DENS_INDX)*phi(i,j,k,VELZ_INDX);
                phi(i,j,k,RHOE_INDX)  = phi(i,j,k,DENS_INDX)*mflo_thermo::get_e_from_rpc(one,one,spec);
                phi(i,j,k,RHOE_INDX) += half*phi(i,j,k,DENS_INDX)  *
                                            ( phi(i,j,k,VELX_INDX)*phi(i,j,k,VELX_INDX)
                                            + phi(i,j,k,VELY_INDX)*phi(i,j,k,VELY_INDX)
                                            + phi(i,j,k,VELZ_INDX)*phi(i,j,k,VELZ_INDX) );
                
                phi(i,j,k,VFRAC_INDX) = one;

                for(int sp=0;sp<NUM_SPECIES;sp++)
                {
                    phi(i,j,k,FLO_NVARS+sp)=one;
                }
            }
        }
    }
}
#endif

#ifndef _PROB_H_
#define _PROB_H_

#include <AMReX_Box.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_Geometry.H>
#include <userfuncs.H>
#include <globalDefines.H>
#include <mflo_constants.H>

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
    Real spec[NUM_SPECIES]={zeroval};
    int dir=0;

    xlen=prob_hi[0]-prob_lo[0];
    ylen=prob_hi[1]-prob_lo[1];
    zlen=prob_hi[2]-prob_lo[2];

    maxlen=std::max(xlen,std::max(ylen,zlen));

    if(maxlen==xlen) dir=0;
    if(maxlen==ylen) dir=1;
    if(maxlen==zlen) dir=2;

    Real xyz[3]={zeroval};

#ifdef _OPENMP
#pragma omp parallel for collapse(2) if (GPU::notInLaunchRegion)
#endif
    for (int k = lo.z; k <= hi.z; ++k) 
    {
        for (int j = lo.y; j <= hi.y; ++j) 
        {
            xyz[2] = prob_lo[2] + (0.5+k) * dx[2];
            xyz[1] = prob_lo[1] + (0.5+j) * dx[1];

            AMREX_PRAGMA_SIMD
            for (int i = lo.x; i <= hi.x; ++i) 
            {
                xyz[0] = prob_lo[0] + (0.5+i) * dx[0]; 

                
                Real vfrac=0.0;
                for(int kk=0;kk<2;kk++)
                {
                   for(int jj=0;jj<2;jj++)
                   {
                        for(int ii=0;ii<2;ii++)
                        {
                            Real xx[3]={prob_lo[0]+(i+ii)*dx[0],
                                prob_lo[1]+(j+jj)*dx[1],           
                                prob_lo[2]+(k+kk)*dx[2]};
                
                            if(xx[dir]<mflo_user_funcs::iboffset)
                            {
                                vfrac+=one;
                            }
                        }
                   }
                }
                
                vfrac*=eighth;
                phi(i,j,k,VFRAC_INDX) = vfrac;

                spec[SITE1_ID]=1.0-vfrac;
                
                phi(i,j,k,PRES_INDX) = one;
                phi(i,j,k,DENS_INDX) = one;
                phi(i,j,k,VELX_INDX) = zeroval;
                phi(i,j,k,VELY_INDX) = zeroval;
                phi(i,j,k,VELZ_INDX) = zeroval;
                phi(i,j,k,TEMP_INDX) = mflo_thermo::get_t_from_rpc(one,one,spec);


                if(phi(i,j,k,VFRAC_INDX) < one)
                {
                        phi(i,j,k,VELX_INDX+dir) = zeroval;
                }
                else
                {
                        phi(i,j,k,VELX_INDX+dir) = one;
                }


                phi(i,j,k,RHO_INDX)   = phi(i,j,k,DENS_INDX);
                phi(i,j,k,RHOU_INDX)  = phi(i,j,k,DENS_INDX)*phi(i,j,k,VELX_INDX);
                phi(i,j,k,RHOV_INDX)  = phi(i,j,k,DENS_INDX)*phi(i,j,k,VELY_INDX);
                phi(i,j,k,RHOW_INDX)  = phi(i,j,k,DENS_INDX)*phi(i,j,k,VELZ_INDX);
                phi(i,j,k,RHOE_INDX)  = phi(i,j,k,DENS_INDX)*mflo_thermo::get_e_from_rpc(one,one,spec);
                phi(i,j,k,RHOE_INDX) += half*phi(i,j,k,DENS_INDX)  *
                                            ( phi(i,j,k,VELX_INDX)*phi(i,j,k,VELX_INDX)
                                            + phi(i,j,k,VELY_INDX)*phi(i,j,k,VELY_INDX)
                                            + phi(i,j,k,VELZ_INDX)*phi(i,j,k,VELZ_INDX) );


                for(int sp=0;sp<NUM_SPECIES;sp++)
                {
                    phi(i,j,k,FLO_NVARS+sp)=spec[sp];
                }
            }
        }
    }
}
#endif

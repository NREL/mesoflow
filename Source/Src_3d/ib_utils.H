#ifndef _IB_UTILS_H_
#define _IB_UTILS_H_

#include <AMReX_BLFort.H>
#include <AMReX_Box.H>
#include <AMReX_Geometry.H>
#include <globalDefines.H>
#include <mflo_constants.H>
#include <consprim_funcs.H>
#include <species.H>

    
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
         amrex::Real basis1d(int i,amrex::Real zeta)
{
     amrex::Real val=0.0;
     if(i==-1)
     {
      val=half*zeta*(zeta-one);
     }
     else if(i==0)
     {
      val=(one-zeta*zeta);
     }
     else if(i==1)
     {
      val=half*zeta*(zeta+one);
     }
     else //this case shouldnt exist
     {
        val=0.0;
     }

     return(val);
}

    AMREX_GPU_DEVICE AMREX_FORCE_INLINE
         amrex::Real basisder1d(int i,amrex::Real zeta)
{
     amrex::Real val=0.0;
     if(i==-1)
     {
      val=half*(two*zeta-one);
     }
     else if(i==0)
     {
      val=-two*zeta;
     }
     else if(i==1)
     {
      val=half*(two*zeta+one);
     }
     else //this case shouldnt exist
     {
        val=0.0;
     }

     return(val);
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
     amrex::Real basisval(int l,int m,int n,
     amrex::Real xp[AMREX_SPACEDIM]) 
{
    amrex::Real lval,mval,nval;
    amrex::Real alpha,beta,gama;
    
   alpha = xp[XDIR];
    beta = xp[YDIR];
    gama = xp[ZDIR];

    lval=basis1d(l,alpha);
    mval=basis1d(m,beta);
    nval=basis1d(n,gama);

    return(lval*mval*nval);
}

    AMREX_GPU_DEVICE AMREX_FORCE_INLINE
     amrex::Real basisderval(int l,int m,int n,int dir,
     amrex::Real xp[AMREX_SPACEDIM],
     const GpuArray<Real, AMREX_SPACEDIM>& dx)
{
    amrex::Real lval,mval,nval;
    amrex::Real alpha,beta,gama;
    
   alpha = xp[XDIR]/dx[XDIR];
    beta = xp[YDIR]/dx[YDIR];
    gama = xp[ZDIR]/dx[ZDIR];

    lval=(dir==0)?basisder1d(l,alpha)/dx[XDIR]:basis1d(l,alpha);
    mval=(dir==1)?basisder1d(m,beta)/dx[YDIR]:basis1d(m,beta);
    nval=(dir==2)?basisder1d(n,gama)/dx[ZDIR]:basis1d(n,gama);

    return(lval*mval*nval);
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE amrex::Real 
do_interpolation_dirichlet(int i,int j,int k,int ngrow,
                           Real dist,Real dist_cutoff,
                           Real dircval,Array4<Real> phi,int comp,
                           const GpuArray<Real, AMREX_SPACEDIM>& dx)
{
    Real dxmin=std::min(dx[XDIR],std::min(dx[YDIR],dx[ZDIR]));
    Real dx_nd[3]={dx[0]/dxmin, dx[1]/dxmin, dx[2]/dxmin};
    Real dist_nd=amrex::Math::abs(dist)/dxmin;
    Real cellval;

    if(amrex::Math::abs(dist)<dist_cutoff)
    {
        cellval=dircval;
    }
    else
    {

        //IDW weighting for now
        Real sum=0.0;
        Real sumwts=0.0;
        for(int n=-ngrow;n<=ngrow;n++)
        {
            for(int m=-ngrow;m<=ngrow;m++)
            {
                for(int l=-ngrow;l<=ngrow;l++)
                {
                    if(l==0 && m==0 && n==0)
                    {  
                        //miss the center cell
                        //this is our output point
                        continue;
                    }

                    if(phi(i+l,j+m,k+n,VFRAC_INDX) == 1.0)
                    {
                        Real dist2 = std::pow(l*dx_nd[0],2.0)
                        + std::pow(m*dx_nd[1],2.0)
                        + std::pow(n*dx_nd[2],2.0);

                        sum+=phi(i+l,j+m,k+n,comp)/dist2;
                        sumwts+=1.0/dist2;

                    }
                }
            }
        }
        sum += dircval/(dist_nd*dist_nd);
        sumwts += 1.0/(dist_nd*dist_nd);

        cellval=sum/sumwts;
    }

    return(cellval);
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE amrex::Real 
do_interpolation_hom_neumann(int i,int j,int k,int ngrow,
                             Real xp[3],Array4<Real> phi,int comp,
                             const GpuArray<Real, AMREX_SPACEDIM>& dx)
{
    Real dxmin=std::min(dx[XDIR],std::min(dx[YDIR],dx[ZDIR]));
    Real dx_nd[3]={dx[0]/dxmin, dx[1]/dxmin, dx[2]/dxmin};
    Real xp_nd[3]={xp[0]/dxmin, xp[1]/dxmin, xp[2]/dxmin};

    //IDW weighting for now
    Real sum=0.0;
    Real sumwts=0.0;
    for(int n=-ngrow;n<=ngrow;n++)
    {
        for(int m=-ngrow;m<=ngrow;m++)
        {
            for(int l=-ngrow;l<=ngrow;l++)
            {
                if(phi(i+l,j+m,k+n,VFRAC_INDX) == 1.0)
                {
                    Real dist2 = std::pow((l*dx_nd[0]-xp_nd[0]),2.0)
                    + std::pow((m*dx_nd[1]-xp_nd[1]),2.0)
                    + std::pow((n*dx_nd[2]-xp_nd[2]),2.0);

                    sum+=phi(i+l,j+m,k+n,comp)/dist2;
                    sumwts+=1.0/dist2;
                }
            }
        }
    }

    Real cellval;
    //this can happen if the 
    //mesh is under-resolved
    if(sumwts>0.0)
    {
        cellval=sum/sumwts;
    }
    else
    {
        cellval=phi(i,j,k,comp);
    }


    return(cellval);
}


AMREX_GPU_DEVICE AMREX_FORCE_INLINE void 
update_cutcells(int i, int j, int k,
                Array4<Real> const& phi, //state variable
                const GpuArray<Real, AMREX_SPACEDIM>& dx,
                const GpuArray<Real, AMREX_SPACEDIM>& plo,
                Real time)
{

    Real lsgrad[3]={0.0};
    Real center[3]={0.0};
    Real tol=1e-5;
    Real primvars[NCVARS];
    Real consvars[NCVARS+NUM_SPECIES];
    Real dxmin=std::min(dx[XDIR],std::min(dx[YDIR],dx[ZDIR]));
    Real closest_dist=dxmin*tol;

    if(phi(i,j,k,VFRAC_INDX)>zeroval && phi(i,j,k,VFRAC_INDX)<one)
    {
        for(int n=-1;n<=1;n++)
        {
            for(int m=-1;m<=1;m++)
            {
                for(int l=-1;l<=1;l++)
                {
                    lsgrad[XDIR] += phi(i+l,j+m,k+n,VFRAC_INDX)*basisderval(l,m,n,XDIR,center,dx);
                    lsgrad[YDIR] += phi(i+l,j+m,k+n,VFRAC_INDX)*basisderval(l,m,n,YDIR,center,dx);
                    lsgrad[ZDIR] += phi(i+l,j+m,k+n,VFRAC_INDX)*basisderval(l,m,n,ZDIR,center,dx);
                }
            }
        }

        Real lsgrad_mag=std::sqrt( lsgrad[XDIR]*lsgrad[XDIR]
                                +  lsgrad[YDIR]*lsgrad[YDIR]
                                +  lsgrad[ZDIR]*lsgrad[ZDIR]);

        //absolute distance from center
        Real dist=(phi(i,j,k,VFRAC_INDX)-half)/lsgrad_mag; //can be negative for cells inside solid

        Real lsetnormal[3]={lsgrad[XDIR]/lsgrad_mag,lsgrad[YDIR]/lsgrad_mag,lsgrad[ZDIR]/lsgrad_mag};

        //minus sign to move back along the normal
        Real xp[3]={-dist*lsetnormal[0],-dist*lsetnormal[1],-dist*lsetnormal[2]};

        phi(i,j,k,VELX_INDX)=do_interpolation_dirichlet(i,j,k,2,dist,closest_dist,0.0,phi,VELX_INDX,dx);
        phi(i,j,k,VELY_INDX)=do_interpolation_dirichlet(i,j,k,2,dist,closest_dist,0.0,phi,VELY_INDX,dx);
        phi(i,j,k,VELZ_INDX)=do_interpolation_dirichlet(i,j,k,2,dist,closest_dist,0.0,phi,VELZ_INDX,dx);

        phi(i,j,k,PRES_INDX)=do_interpolation_hom_neumann(i,j,k,2,xp,phi,PRES_INDX,dx);
        phi(i,j,k,DENS_INDX)=do_interpolation_hom_neumann(i,j,k,2,xp,phi,DENS_INDX,dx);

        //adiabatic wall for now
        phi(i,j,k,TEMP_INDX)=do_interpolation_hom_neumann(i,j,k,2,xp,phi,TEMP_INDX,dx);
        //amrex::Print()<<"vars:"<<phi(i,j,k,VELX_INDX)<<"\t"<<phi(i,j,k,VELY_INDX)<<"\t"<<phi(i,j,k,VELZ_INDX)<<
        //"\t"<<phi(i,j,k,PRES_INDX)<<"\t"<<phi(i,j,k,DENS_INDX)<<"\n";
        //for species (zero gradient)
        //for(int sp=0;sp<NUM_SPECIES;sp++)
        //{
        //  phi(i,j,k,FLO_NVARS+sp)=do_interpolation_hom_neumann(i,j,k,2,xp,phi,FLO_NVARS+sp,dx);
        //  consvars[NCVARS+sp]=phi(i,j,k,FLO_NVARS+sp);
        //}

        primvars[DENS_IND] = phi(i,j,k,DENS_INDX);
        primvars[VELX_IND] = phi(i,j,k,VELX_INDX);
        primvars[VELY_IND] = phi(i,j,k,VELY_INDX);
        primvars[VELZ_IND] = phi(i,j,k,VELZ_INDX);
        primvars[PRES_IND] = phi(i,j,k,PRES_INDX);

        prim_to_cons(consvars,primvars);

        for(int c=0;c<NCVARS;c++)
        {
            phi(i,j,k,RHO_IND+c)=consvars[c];
        }

    }


}

#endif

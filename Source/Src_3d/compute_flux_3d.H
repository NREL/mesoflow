#ifndef _compute_hyp_flux_3d_H_
#define _compute_hyp_flux_3d_H_

#include <AMReX_BLFort.H>
#include <AMReX_Box.H>
#include <AMReX_Geometry.H>
#include <globalDefines.H>
#include <mflo_constants.H>
#include <consprim_funcs.H>
#include <flux_schemes.H>
#include <species.H>

using namespace amrex;

AMREX_GPU_DEVICE
    AMREX_FORCE_INLINE
Real minmod_flux_limiter(Real r)
{
    return(std::max(zeroval,std::min(one,r)));
}

AMREX_GPU_DEVICE
    AMREX_FORCE_INLINE
void get_higherorder_states(Real uL[NCVARS+NUM_SPECIES],Real uR[NCVARS+NUM_SPECIES],Real uim2[NCVARS+NUM_SPECIES],Real uip1[NCVARS+NUM_SPECIES])
{
    //refer to page 114,chap 4 in Blazek's textbook

    //@inproceedings{tatsumi1995new,
    //  title={A new high resolution scheme for compressible viscous flows with shocks},
    //    author={Tatsumi, S and Martinelli, L and Jameson, A},
    //    booktitle={33rd Aerospace Sciences Meeting and Exhibit},
    //    pages={466},
    //    year={1995}
    //    }

    Real delt1,delt2;
    Real L[NCVARS+NUM_SPECIES],psi;
    Real sigma=two;

    for(int c=0;c<(NCVARS+NUM_SPECIES);c++)
    {
        delt1 = uip1[c]   -   uR[c];
        delt2 =   uL[c]   - uim2[c];
        psi   = one - pow(fabs((delt1-delt2)/(fabs(delt1)+fabs(delt2)+TINYVAL)),sigma);
        L[c]  = half*psi*(delt1+delt2); 
    }

    for(int c=0;c<(NCVARS+NUM_SPECIES);c++)
    {
        uL[c] = uL[c] + half*L[c];    
        uR[c] = uR[c] - half*L[c];
    }   

}

AMREX_GPU_DEVICE
    AMREX_FORCE_INLINE
void compute_flux_x(int i, int j, int k,
        Array4<Real> const& phi, //state variable
        Array4<Real> const& fluid_transpcoeffs,
        Array4<Real> const& specdiff,
        Array4<Real> const& fx,
        const GpuArray<Real, AMREX_SPACEDIM>& dx,bool nsflag)
{

    Real ul[NCVARS+NUM_SPECIES],ur[NCVARS+NUM_SPECIES];
    Real pl[NCVARS],pr[NCVARS];
    Real normal[AMREX_SPACEDIM];
    Real dudx,dudy,dudz;
    Real dvdx,dvdy,dvdz;
    Real dwdx,dwdy,dwdz;
    Real delu;
    Real dTdx;
    Real tau[AMREX_SPACEDIM][AMREX_SPACEDIM];
    Real visc,thcond;

    normal[XDIR] = one;
    normal[YDIR] = zeroval;
    normal[ZDIR] = zeroval;

    Real fhalf[NCVARS+NUM_SPECIES]={zeroval};

    visc   = half*( fluid_transpcoeffs(i,j,k,VISC_INDX)   + fluid_transpcoeffs(i-1,j,k,VISC_INDX)   );
    thcond = half*( fluid_transpcoeffs(i,j,k,THCOND_INDX) + fluid_transpcoeffs(i-1,j,k,THCOND_INDX) );

#ifdef IMM_BOUNDARY
    Real vfrac_L=phi(i-1,j,k,VFRAC_INDX);
    Real vfrac_R=phi(i,j,k,VFRAC_INDX);
#endif

    //face velocity
    Real uface = half*(phi(i-1,j,k,VELX_INDX) + phi(i,j,k,VELX_INDX));
    Real vface = half*(phi(i-1,j,k,VELY_INDX) + phi(i,j,k,VELY_INDX));
    Real wface = half*(phi(i-1,j,k,VELZ_INDX) + phi(i,j,k,VELZ_INDX));

#ifdef IMM_BOUNDARY
        //vfrac is zeroval on left or right, then face vel is zeroval
        if(vfrac_L == zeroval or vfrac_R == zeroval)
        {
                uface=zeroval;
                vface=zeroval;
                wface=zeroval;
        }
#endif

    if(nsflag)
    {

        for(int c=0;c<NCVARS;c++)
        {
            ul[c+RHO_IND] = phi(i-1,j,k,c+RHO_INDX);
            ur[c+RHO_IND] = phi(i  ,j,k,c+RHO_INDX);
        }
        for(int c=0;c<NUM_SPECIES;c++)
        {
            ul[c+NCVARS] = phi(i-1,j,k,c+FLO_NVARS);
            ur[c+NCVARS] = phi(i  ,j,k,c+FLO_NVARS);
        }

#ifdef IMM_BOUNDARY
        if(vfrac_L > zeroval and vfrac_R > zeroval)
#endif
        {

#ifdef SEC_ORDER_FLUX

#ifdef IMM_BOUNDARY
            if(phi(i-2,j,k,VFRAC_INDX) > zeroval and phi(i+1,j,k,VFRAC_INDX) > zeroval)
#endif
            {

                Real ulm1[NCVARS+NUM_SPECIES],urp1[NCVARS+NUM_SPECIES];
                for(int c=0;c<NCVARS;c++)
                {
                    ulm1[c+RHO_IND] = phi(i-2,j,k,c+RHO_INDX);
                    urp1[c+RHO_IND] = phi(i+1,j,k,c+RHO_INDX);
                }
                for(int c=0;c<NUM_SPECIES;c++)
                {
                    ulm1[c+NCVARS] = phi(i-2,j,k,FLO_NVARS+c);
                    urp1[c+NCVARS] = phi(i+1,j,k,FLO_NVARS+c);
                }
                get_higherorder_states(ul,ur,ulm1,urp1);
            }
#endif
        }
#ifdef IMM_BOUNDARY
        else if(vfrac_R == zeroval)
        {
            cons_to_prim(ul,pl);

            pr[DENS_IND] =  pl[DENS_IND];
            pr[VELX_IND] = -pl[VELX_IND];
            pr[VELY_IND] = -pl[VELY_IND];
            pr[VELZ_IND] = -pl[VELZ_IND];
            pr[PRES_IND] =  pl[PRES_IND];

            //zeroval gradient for species
            //note: this is only to make sure
            //rho,rhou,rhov,rhow,rhoe fluxes are correct
            for(int c=0;c<NUM_SPECIES;c++)
            {
                ur[c+NCVARS] = ul[c+NCVARS];
            }

            prim_to_cons(ur,pr);
        }
        else
        {
            cons_to_prim(ur,pr);

            pl[DENS_IND] =  pr[DENS_IND];
            pl[VELX_IND] = -pr[VELX_IND];
            pl[VELY_IND] = -pr[VELY_IND];
            pl[VELZ_IND] = -pr[VELZ_IND];
            pl[PRES_IND] =  pr[PRES_IND];

            //zeroval gradient for species
            //note: this is only to make sure
            //rho,rhou,rhov,rhow,rhoe fluxes are correct
            for(int c=0;c<NUM_SPECIES;c++)
            {
                ul[c+NCVARS] = ur[c+NCVARS];
            }

            prim_to_cons(ul,pl);
        }
#endif

#ifdef IMM_BOUNDARY
        //if both vfracs are 0 then fhalf should be 0
        if(vfrac_L > zeroval or vfrac_R > zeroval)
#endif
        {
           ausmp_up_flux(ul,ur,fhalf,normal);
           //laxf_flux(ul,ur,fhalf,normal);
        }

        for(int c=0;c<NCVARS;c++)
        {
            fx(i,j,k,c+RHO_INDX) = fhalf[c+RHO_IND];
        }

#ifdef IMM_BOUNDARY
        //if both vfracs are 0 dont do anything
        if(vfrac_L > zeroval or vfrac_R > zeroval)
#endif
        {
            viscous_stencil_x(i,j,k,phi,VELX_INDX,dx,dudx,dudy,dudz);
            viscous_stencil_x(i,j,k,phi,VELY_INDX,dx,dvdx,dvdy,dvdz);
            viscous_stencil_x(i,j,k,phi,VELZ_INDX,dx,dwdx,dwdy,dwdz);

            dTdx = (phi(i,j,k,TEMP_INDX) - phi(i-1,j,k,TEMP_INDX))/dx[XDIR];

            delu = dudx + dvdy + dwdz;

            tau[XDIR][XDIR]= visc*(two*dudx - two3rd*delu);
            tau[XDIR][YDIR]= visc*(dudy + dvdx);
            tau[XDIR][ZDIR]= visc*(dudz + dwdx);

            tau[YDIR][XDIR]= visc*(dvdx + dudy);
            tau[YDIR][YDIR]= visc*(two*dvdy - two3rd*delu);
            tau[YDIR][ZDIR]= visc*(dvdz + dwdy);

            tau[ZDIR][XDIR]= visc*(dwdx + dudz);
            tau[ZDIR][YDIR]= visc*(dwdy + dvdz);
            tau[ZDIR][ZDIR]= visc*(two*dwdz - two3rd*delu);

            fx(i,j,k,RHOU_INDX) -= ( tau[XDIR][XDIR]*normal[XDIR] 
                    + tau[YDIR][XDIR]*normal[YDIR] 
                    + tau[ZDIR][XDIR]*normal[ZDIR]);

            fx(i,j,k,RHOV_INDX) -= ( tau[XDIR][YDIR]*normal[XDIR] 
                    + tau[YDIR][YDIR]*normal[YDIR] 
                    + tau[ZDIR][YDIR]*normal[ZDIR]);

            fx(i,j,k,RHOW_INDX) -= ( tau[XDIR][ZDIR]*normal[XDIR] 
                    + tau[YDIR][ZDIR]*normal[YDIR] 
                    + tau[ZDIR][ZDIR]*normal[ZDIR]);

            fx(i,j,k,RHOE_INDX) -= thcond*dTdx;

            //viscous dissipation term
            //this term is the work done on the control volume by viscous forces - tau_{ji} n_j v_i
            fx(i,j,k,RHOE_INDX) -= (tau[XDIR][XDIR]*normal[XDIR] + 
                    tau[YDIR][XDIR]*normal[YDIR] +
                    tau[ZDIR][XDIR]*normal[ZDIR] )*uface;

            fx(i,j,k,RHOE_INDX) -= (tau[XDIR][YDIR]*normal[XDIR] + 
                    tau[YDIR][YDIR]*normal[YDIR] +
                    tau[ZDIR][YDIR]*normal[ZDIR] )*vface;

            fx(i,j,k,RHOE_INDX) -= (tau[XDIR][ZDIR]*normal[XDIR] + 
                    tau[YDIR][ZDIR]*normal[YDIR] +
                    tau[ZDIR][ZDIR]*normal[ZDIR] )*wface;
        }

    }    
    //Species convection and diffusion
    //convection
    for(int sp=0;sp<NUM_SPECIES;sp++)
    {
        Real phi_im2   = phi(i-2,j,k,FLO_NVARS+sp);
        Real phi_im1   = phi(i-1,j,k,FLO_NVARS+sp);
        Real phi_i     = phi(i,j,k,FLO_NVARS+sp);
        Real phi_ip1   = phi(i+1,j,k,FLO_NVARS+sp);
        Real rL = (phi_im1-phi_im2)/(phi_i-phi_im1+SMALLVAL);
        Real rR = (phi_i-phi_im1)/(phi_ip1-phi_i+SMALLVAL);
        Real phiL = phi_im1 + minmod_flux_limiter(rL)*half*(phi_im1-phi_im2);
        Real phiR = phi_i   - minmod_flux_limiter(rR)*half*(phi_ip1-phi_i);

        fx(i,j,k,FLO_NVARS+sp) = mflo_species::advect_flags[sp]*
            (half*(uface+fabs(uface))*phiL +
             half*(uface-fabs(uface))*phiR);
    }

    //diffusion
    for(int sp=0;sp<NUM_SPECIES;sp++)
    {
        fx(i,j,k,FLO_NVARS+sp) -= specdiff(i,j,k,sp)*
            (phi(i,j,k,FLO_NVARS+sp)-phi(i-1,j,k,FLO_NVARS+sp))/dx[XDIR];
    }
}

AMREX_GPU_DEVICE
    AMREX_FORCE_INLINE
void compute_flux_y(int i, int j, int k,
        Array4<Real> const& phi, //state variable
        Array4<Real> const& fluid_transpcoeffs,
        Array4<Real> const& specdiff,
        Array4<Real> const& fy,
        const GpuArray<Real, AMREX_SPACEDIM>& dx,bool nsflag)
{
    Real ul[NCVARS+NUM_SPECIES],ur[NCVARS+NUM_SPECIES];
    Real pl[NCVARS],pr[NCVARS];
    Real normal[AMREX_SPACEDIM];
    Real dudx,dudy,dudz;
    Real dvdx,dvdy,dvdz;
    Real dwdx,dwdy,dwdz;
    Real delu;
    Real dTdy;
    Real tau[AMREX_SPACEDIM][AMREX_SPACEDIM];
    Real visc,thcond;

    normal[XDIR] = zeroval;
    normal[YDIR] = one;
    normal[ZDIR] = zeroval;

    Real fhalf[NCVARS+NUM_SPECIES]={zeroval};

    visc   = half*(fluid_transpcoeffs(i,j-1,k,VISC_INDX)   +   fluid_transpcoeffs(i,j,k,VISC_INDX));
    thcond = half*(fluid_transpcoeffs(i,j-1,k,THCOND_INDX) +   fluid_transpcoeffs(i,j,k,THCOND_INDX));

    //face velocity
    Real uface = half*(phi(i,j-1,k,VELX_INDX) + phi(i,j,k,VELX_INDX));
    Real vface = half*(phi(i,j-1,k,VELY_INDX) + phi(i,j,k,VELY_INDX));
    Real wface = half*(phi(i,j-1,k,VELZ_INDX) + phi(i,j,k,VELZ_INDX));

#ifdef IMM_BOUNDARY
    Real vfrac_L=phi(i,j-1,k,VFRAC_INDX);
    Real vfrac_R=phi(i,j,k,VFRAC_INDX);
#endif

#ifdef IMM_BOUNDARY
    //vfrac is zeroval on left or right, then face vel is zeroval
    if(vfrac_L == zeroval or vfrac_R == zeroval)
    {
        uface=zeroval;
        vface=zeroval;
        wface=zeroval;
    }
#endif

    if(nsflag)
    {
        for(int c=0;c<NCVARS;c++)
        {
            ul[c+RHO_IND]=phi(i,j-1,k,c+RHO_INDX);
            ur[c+RHO_IND]=phi(i,  j,k,c+RHO_INDX);
        }
        for(int c=0;c<NUM_SPECIES;c++)
        {
            ul[c+NCVARS] = phi(i,j-1,k,c+FLO_NVARS);
            ur[c+NCVARS] = phi(i,j  ,k,c+FLO_NVARS);
        }

#ifdef IMM_BOUNDARY
        if(vfrac_L > zeroval and vfrac_R > zeroval)
#endif
        {
#ifdef SEC_ORDER_FLUX

#ifdef IMM_BOUNDARY
            if(phi(i,j-2,k,VFRAC_INDX) > zeroval and phi(i,j+1,k,VFRAC_INDX) > zeroval)
#endif
            {
                Real ulm1[NCVARS+NUM_SPECIES],urp1[NCVARS+NUM_SPECIES];
                for(int c=0;c<NCVARS;c++)
                {
                    ulm1[c+RHO_IND] = phi(i,j-2,k,c+RHO_INDX);
                    urp1[c+RHO_IND] = phi(i,j+1,k,c+RHO_INDX);
                }
                for(int c=0;c<NUM_SPECIES;c++)
                {
                    ulm1[c+NCVARS] = phi(i,j-2,k,c+FLO_NVARS);
                    urp1[c+NCVARS] = phi(i,j+1,k,c+FLO_NVARS);
                }
                get_higherorder_states(ul,ur,ulm1,urp1);
            }
#endif
        }
#ifdef IMM_BOUNDARY
        else if(vfrac_R == zeroval)
        {
            cons_to_prim(ul,pl);

            pr[DENS_IND] =  pl[DENS_IND];
            pr[VELX_IND] = -pl[VELX_IND];
            pr[VELY_IND] = -pl[VELY_IND];
            pr[VELZ_IND] = -pl[VELZ_IND];
            pr[PRES_IND] =  pl[PRES_IND];

            //zeroval gradient for species
            //note: this is only to make sure
            //rho,rhou,rhov,rhow,rhoe fluxes are correct
            for(int c=0;c<NUM_SPECIES;c++)
            {
                ur[c+NCVARS] = ul[c+NCVARS];
            }

            prim_to_cons(ur,pr);
        }
        else
        {
            cons_to_prim(ur,pr);

            pl[DENS_IND] =  pr[DENS_IND];
            pl[VELX_IND] = -pr[VELX_IND];
            pl[VELY_IND] = -pr[VELY_IND];
            pl[VELZ_IND] = -pr[VELZ_IND];
            pl[PRES_IND] =  pr[PRES_IND];

            //zeroval gradient for species
            //note: this is only to make sure
            //rho,rhou,rhov,rhow,rhoe fluxes are correct
            for(int c=0;c<NUM_SPECIES;c++)
            {
                ul[c+NCVARS] = ur[c+NCVARS];
            }

            prim_to_cons(ul,pl);
        }
#endif

#ifdef IMM_BOUNDARY
        //if both vfracs are 0 then fhalf should be 0
        if(vfrac_L > zeroval or vfrac_R > zeroval)
#endif
        {
            ausmp_up_flux(ul,ur,fhalf,normal);
            //laxf_flux(ul,ur,fhalf,normal);
        }

        for(int c=0;c<NCVARS;c++)
        {
            fy(i,j,k,c+RHO_INDX) = fhalf[c+RHO_IND];
        } 

#ifdef IMM_BOUNDARY
        if(vfrac_L > zeroval or vfrac_R > zeroval)
#endif
        {
            viscous_stencil_y(i,j,k,phi,VELX_INDX,dx,dudx,dudy,dudz);
            viscous_stencil_y(i,j,k,phi,VELY_INDX,dx,dvdx,dvdy,dvdz);
            viscous_stencil_y(i,j,k,phi,VELZ_INDX,dx,dwdx,dwdy,dwdz);

            dTdy = (phi(i,j,k,TEMP_INDX) - phi(i,j-1,k,TEMP_INDX))/dx[YDIR];

            delu = dudx + dvdy + dwdz;

            tau[XDIR][XDIR]= visc*(two*dudx - two3rd*delu);
            tau[XDIR][YDIR]= visc*(dudy + dvdx);
            tau[XDIR][ZDIR]= visc*(dudz + dwdx);

            tau[YDIR][XDIR]= visc*(dvdx + dudy);
            tau[YDIR][YDIR]= visc*(two*dvdy - two3rd*delu);
            tau[YDIR][ZDIR]= visc*(dvdz + dwdy);

            tau[ZDIR][XDIR]= visc*(dwdx + dudz);
            tau[ZDIR][YDIR]= visc*(dwdy + dvdz);
            tau[ZDIR][ZDIR]= visc*(two*dwdz - two3rd*delu);

            fy(i,j,k,RHOU_INDX) -= ( tau[XDIR][XDIR]*normal[XDIR] 
                    + tau[YDIR][XDIR]*normal[YDIR] 
                    + tau[ZDIR][XDIR]*normal[ZDIR]);

            fy(i,j,k,RHOV_INDX) -= ( tau[XDIR][YDIR]*normal[XDIR] 
                    + tau[YDIR][YDIR]*normal[YDIR] 
                    + tau[ZDIR][YDIR]*normal[ZDIR]);

            fy(i,j,k,RHOW_INDX) -= ( tau[XDIR][ZDIR]*normal[XDIR] 
                    + tau[YDIR][ZDIR]*normal[YDIR] 
                    + tau[ZDIR][ZDIR]*normal[ZDIR]);

            fy(i,j,k,RHOE_INDX) -= thcond*dTdy;

            //viscous dissipation term

            //this term is the work done on the control volume by viscous forces - tau_{ji} n_j v_i
            fy(i,j,k,RHOE_INDX) -= (tau[XDIR][XDIR]*normal[XDIR] + 
                    tau[YDIR][XDIR]*normal[YDIR] +
                    tau[ZDIR][XDIR]*normal[ZDIR] )*uface;

            fy(i,j,k,RHOE_INDX) -= (tau[XDIR][YDIR]*normal[XDIR] + 
                    tau[YDIR][YDIR]*normal[YDIR] +
                    tau[ZDIR][YDIR]*normal[ZDIR] )*vface;

            fy(i,j,k,RHOE_INDX) -= (tau[XDIR][ZDIR]*normal[XDIR] + 
                    tau[YDIR][ZDIR]*normal[YDIR] +
                    tau[ZDIR][ZDIR]*normal[ZDIR] )*wface;
        }
    } 
    //Species convection and diffusion
    //convection
    for(int sp=0;sp<NUM_SPECIES;sp++)
    {
        Real phi_im2   = phi(i,j-2,k,FLO_NVARS+sp);
        Real phi_im1   = phi(i,j-1,k,FLO_NVARS+sp);
        Real phi_i     = phi(i,j,k,FLO_NVARS+sp);
        Real phi_ip1   = phi(i,j+1,k,FLO_NVARS+sp);
        Real rL = (phi_im1-phi_im2)/(phi_i-phi_im1+SMALLVAL);
        Real rR = (phi_i-phi_im1)/(phi_ip1-phi_i+SMALLVAL);
        Real phiL = phi_im1 + minmod_flux_limiter(rL)*half*(phi_im1-phi_im2);
        Real phiR = phi_i   - minmod_flux_limiter(rR)*half*(phi_ip1-phi_i);

        fy(i,j,k,FLO_NVARS+sp) = mflo_species::advect_flags[sp] * 
            (half*(vface+fabs(vface))*phiL +
             half*(vface-fabs(vface))*phiR);
    }

    //diffusion
    for(int sp=0;sp<NUM_SPECIES;sp++)
    {
        fy(i,j,k,FLO_NVARS+sp) -= specdiff(i,j,k,sp)*
            (phi(i,j,k,FLO_NVARS+sp)-phi(i,j-1,k,FLO_NVARS+sp))/dx[YDIR];
    }
}

AMREX_GPU_DEVICE
    AMREX_FORCE_INLINE
void compute_flux_z(int i, int j, int k,
        Array4<Real> const& phi, //state variable
        Array4<Real> const& fluid_transpcoeffs,
        Array4<Real> const& specdiff,
        Array4<Real> const& fz,
        const GpuArray<Real, AMREX_SPACEDIM>& dx,bool nsflag)
{
    Real ul[NCVARS+NUM_SPECIES],ur[NCVARS+NUM_SPECIES];
    Real pl[NCVARS],pr[NCVARS];
    Real normal[AMREX_SPACEDIM];
    Real dudx,dudy,dudz;
    Real dvdx,dvdy,dvdz;
    Real dwdx,dwdy,dwdz;
    Real delu;
    Real dTdz;
    Real tau[AMREX_SPACEDIM][AMREX_SPACEDIM];
    Real visc,thcond;

    normal[XDIR] = zeroval;
    normal[YDIR] = zeroval;
    normal[ZDIR] = one;

    Real fhalf[NCVARS+NUM_SPECIES]={zeroval};

    visc   = half*(fluid_transpcoeffs(i,j,k-1,VISC_INDX)   + fluid_transpcoeffs(i,j,k,VISC_INDX));
    thcond = half*(fluid_transpcoeffs(i,j,k-1,THCOND_INDX) + fluid_transpcoeffs(i,j,k,THCOND_INDX));

    //face velocity
    Real uface = half*(phi(i,j,k-1,VELX_INDX) + phi(i,j,k,VELX_INDX));
    Real vface = half*(phi(i,j,k-1,VELY_INDX) + phi(i,j,k,VELY_INDX));
    Real wface = half*(phi(i,j,k-1,VELZ_INDX) + phi(i,j,k,VELZ_INDX));

#ifdef IMM_BOUNDARY
    Real vfrac_L=phi(i,j,k-1,VFRAC_INDX);
    Real vfrac_R=phi(i,j,k,VFRAC_INDX);
#endif

#ifdef IMM_BOUNDARY
    //vfrac is zeroval on left or right, then face vel is zeroval
    if(vfrac_L == zeroval or vfrac_R == zeroval)
    {
        uface=zeroval;
        vface=zeroval;
        wface=zeroval;
    }
#endif

    if(nsflag)
    {
        for(int c=0;c<NCVARS;c++)
        {
            ul[c+RHO_IND]=phi(i,j,k-1,c+RHO_INDX);
            ur[c+RHO_IND]=phi(i,j,  k,c+RHO_INDX);
        }
        for(int c=0;c<NUM_SPECIES;c++)
        {
            ul[c+NCVARS]=phi(i,j,k-1,c+FLO_NVARS);
            ur[c+NCVARS]=phi(i,j,  k,c+FLO_NVARS);
        }
#ifdef IMM_BOUNDARY
        if(vfrac_L > zeroval and vfrac_R > zeroval)
#endif
        {
#ifdef SEC_ORDER_FLUX
#ifdef IMM_BOUNDARY
            if(phi(i,j,k-2,VFRAC_INDX) > zeroval and phi(i,j,k+1,VFRAC_INDX) > zeroval)
#endif
            {
                Real ulm1[NCVARS+NUM_SPECIES],urp1[NCVARS+NUM_SPECIES];
                for(int c=0;c<NCVARS;c++)
                {
                    ulm1[c+RHO_IND] = phi(i,j,k-2,c+RHO_INDX);
                    urp1[c+RHO_IND] = phi(i,j,k+1,c+RHO_INDX);
                }
                for(int c=0;c<NUM_SPECIES;c++)
                {
                    ulm1[c+NCVARS] = phi(i,j,k-2,c+FLO_NVARS);
                    urp1[c+NCVARS] = phi(i,j,k+1,c+FLO_NVARS);
                }
                get_higherorder_states(ul,ur,ulm1,urp1);
            }
#endif
        }
#ifdef IMM_BOUNDARY
        else if(vfrac_R == zeroval)
        {
            cons_to_prim(ul,pl);

            pr[DENS_IND] =  pl[DENS_IND];
            pr[VELX_IND] = -pl[VELX_IND];
            pr[VELY_IND] = -pl[VELY_IND];
            pr[VELZ_IND] = -pl[VELZ_IND];
            pr[PRES_IND] =  pl[PRES_IND];

            //zeroval gradient for species
            //note: this is only to make sure
            //rho,rhou,rhov,rhow,rhoe fluxes are correct
            for(int c=0;c<NUM_SPECIES;c++)
            {
                ur[c+NCVARS] = ul[c+NCVARS];
            }

            prim_to_cons(ur,pr);
        }
        else
        {
            cons_to_prim(ur,pr);

            pl[DENS_IND] =  pr[DENS_IND];
            pl[VELX_IND] = -pr[VELX_IND];
            pl[VELY_IND] = -pr[VELY_IND];
            pl[VELZ_IND] = -pr[VELZ_IND];
            pl[PRES_IND] =  pr[PRES_IND];

            //zeroval gradient for species
            //note: this is only to make sure
            //rho,rhou,rhov,rhow,rhoe fluxes are correct
            for(int c=0;c<NUM_SPECIES;c++)
            {
                ul[c+NCVARS] = ur[c+NCVARS];
            }

            prim_to_cons(ul,pl);
        }
#endif

#ifdef IMM_BOUNDARY
        //if both vfracs are 0 then fhalf should be 0
        if(vfrac_L > zeroval or vfrac_R > zeroval)
#endif
        {
            ausmp_up_flux(ul,ur,fhalf,normal);
            //laxf_flux(ul,ur,fhalf,normal);
        }

        for(int c=0;c<NCVARS;c++)
        {
            fz(i,j,k,c+RHO_INDX) = fhalf[c+RHO_IND];
        } 

#ifdef IMM_BOUNDARY
        if(vfrac_L > zeroval or vfrac_R > zeroval)
#endif
        {
            viscous_stencil_z(i,j,k,phi,VELX_INDX,dx,dudx,dudy,dudz);
            viscous_stencil_z(i,j,k,phi,VELY_INDX,dx,dvdx,dvdy,dvdz);
            viscous_stencil_z(i,j,k,phi,VELZ_INDX,dx,dwdx,dwdy,dwdz);

            dTdz = (phi(i,j,k,TEMP_INDX) - phi(i,j,k-1,TEMP_INDX))/dx[ZDIR];

            delu = dudx + dvdy + dwdz;

            tau[XDIR][XDIR]= visc*(two*dudx - two3rd*delu);
            tau[XDIR][YDIR]= visc*(dudy + dvdx);
            tau[XDIR][ZDIR]= visc*(dudz + dwdx);

            tau[YDIR][XDIR]= visc*(dvdx + dudy);
            tau[YDIR][YDIR]= visc*(two*dvdy - two3rd*delu);
            tau[YDIR][ZDIR]= visc*(dvdz + dwdy);

            tau[ZDIR][XDIR]= visc*(dwdx + dudz);
            tau[ZDIR][YDIR]= visc*(dwdy + dvdz);
            tau[ZDIR][ZDIR]= visc*(two*dwdz - two3rd*delu);

            fz(i,j,k,RHOU_INDX) -= ( tau[XDIR][XDIR]*normal[XDIR] 
                    + tau[YDIR][XDIR]*normal[YDIR] 
                    + tau[ZDIR][XDIR]*normal[ZDIR]);

            fz(i,j,k,RHOV_INDX) -= ( tau[XDIR][YDIR]*normal[XDIR] 
                    + tau[YDIR][YDIR]*normal[YDIR] 
                    + tau[ZDIR][YDIR]*normal[ZDIR]);

            fz(i,j,k,RHOW_INDX) -= ( tau[XDIR][ZDIR]*normal[XDIR] 
                    + tau[YDIR][ZDIR]*normal[YDIR] 
                    + tau[ZDIR][ZDIR]*normal[ZDIR]);

            fz(i,j,k,RHOE_INDX) -= thcond*dTdz;

            //viscous dissipation term

            //this term is the work done on the control volume by viscous forces - tau_{ji} n_j v_i
            fz(i,j,k,RHOE_INDX) -= (tau[XDIR][XDIR]*normal[XDIR] + 
                    tau[YDIR][XDIR]*normal[YDIR] +
                    tau[ZDIR][XDIR]*normal[ZDIR] )*uface;

            fz(i,j,k,RHOE_INDX) -= (tau[XDIR][YDIR]*normal[XDIR] + 
                    tau[YDIR][YDIR]*normal[YDIR] +
                    tau[ZDIR][YDIR]*normal[ZDIR] )*vface;

            fz(i,j,k,RHOE_INDX) -= (tau[XDIR][ZDIR]*normal[XDIR] + 
                    tau[YDIR][ZDIR]*normal[YDIR] +
                    tau[ZDIR][ZDIR]*normal[ZDIR] )*wface;
        }
    } 
    //Species convection and diffusion
    //convection
    for(int sp=0;sp<NUM_SPECIES;sp++)
    {
        Real phi_im2   = phi(i,j,k-2,FLO_NVARS+sp);
        Real phi_im1   = phi(i,j,k-1,FLO_NVARS+sp);
        Real phi_i     = phi(i,j,k,FLO_NVARS+sp);
        Real phi_ip1   = phi(i,j,k+1,FLO_NVARS+sp);
        Real rL = (phi_im1-phi_im2)/(phi_i-phi_im1+SMALLVAL);
        Real rR = (phi_i-phi_im1)/(phi_ip1-phi_i+SMALLVAL);
        Real phiL = phi_im1 + minmod_flux_limiter(rL)*half*(phi_im1-phi_im2);
        Real phiR = phi_i   - minmod_flux_limiter(rR)*half*(phi_ip1-phi_i);

        fz(i,j,k,FLO_NVARS+sp) = mflo_species::advect_flags[sp]* 
            (half*(wface+fabs(wface))*phiL +
             half*(wface-fabs(wface))*phiR);
    }

    //diffusion
    for(int sp=0;sp<NUM_SPECIES;sp++)
    {
        fz(i,j,k,FLO_NVARS+sp) -= specdiff(i,j,k,sp)*
            (phi(i,j,k,FLO_NVARS+sp)-phi(i,j,k-1,FLO_NVARS+sp))/dx[ZDIR];
    }
}
#endif

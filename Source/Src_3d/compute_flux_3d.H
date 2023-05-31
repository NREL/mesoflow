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
#include <ib_utils.H>

using namespace amrex;

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
Real harmonic_mean(Real d1,Real d2)
{
   Real heps=1e-30;
   return(2.0 * d1 * d2 / (d1 + d2 + heps));
}

AMREX_GPU_DEVICE
    AMREX_FORCE_INLINE
Real minmod_flux_limiter(Real r)
{
    return(std::max(zeroval,std::min(one,r)));
}

AMREX_GPU_DEVICE
    AMREX_FORCE_INLINE
void get_higherorder_states(Real uL[NCVARS+NUM_SPECIES],Real uR[NCVARS+NUM_SPECIES],
                            Real uim2[NCVARS+NUM_SPECIES],Real uip1[NCVARS+NUM_SPECIES])
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
void get_higherorder_states(Real uL,Real uR,
                            Real uim2,Real uip1)
{
    Real delt1,delt2;
    Real L,psi;
    Real sigma=two;

    delt1 = uip1   -   uR;
    delt2 = uL   - uim2;
    psi   = one - pow(fabs((delt1-delt2)/(fabs(delt1)+fabs(delt2)+TINYVAL)),sigma);
    L  = half*psi*(delt1+delt2); 

    uL = uL + half*L;    
    uR = uR - half*L;
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void compute_flux_weno(int i, int j, int k,int sweepdir,
                  Array4<Real> const& phi, //state variable
                  Array4<Real> const& fluid_transpcoeffs,
                  Array4<Real> const& specdiff,
                  Array4<Real> const& fx,
                  const GpuArray<Real, AMREX_SPACEDIM>& dx,
                  int hyperbolics_order,Real dissfactor,bool nsflag,
                  int spec_in_solid)
{
    Real ur[NCVARS+NUM_SPECIES],urp1[NCVARS+NUM_SPECIES],urp2[NCVARS+NUM_SPECIES];
    Real ul[NCVARS+NUM_SPECIES],ulm1[NCVARS+NUM_SPECIES],ulm2[NCVARS+NUM_SPECIES];
    
    Real fr[NCVARS+NUM_SPECIES],frp1[NCVARS+NUM_SPECIES],frp2[NCVARS+NUM_SPECIES];
    Real fl[NCVARS+NUM_SPECIES],flm1[NCVARS+NUM_SPECIES],flm2[NCVARS+NUM_SPECIES];

    Real fphalf[NCVARS+NUM_SPECIES]={zeroval};
    Real fmhalf[NCVARS+NUM_SPECIES]={zeroval};
    Real uphalf[NCVARS+NUM_SPECIES]={zeroval};
    Real umhalf[NCVARS+NUM_SPECIES]={zeroval};
    
    Real normal[AMREX_SPACEDIM]={zeroval};
    normal[sweepdir] = one;
    Real lambda_max;

    lambda_max=get_max_eigenvalue(ul,ur,normal);

    IntVect ivm2(i, j, k);
    IntVect ivm1(i, j, k);
    IntVect iv(i, j, k);
    IntVect ivp1(i, j, k);
    IntVect ivp2(i, j, k);
    IntVect ivp3(i, j, k);

    ivm2[sweepdir]          -= 2;
    ivm1[sweepdir]          -= 1;
    ivp1[sweepdir]          += 1;
    ivp2[sweepdir]          += 2;
    ivp3[sweepdir]          += 3;
    
    for(int c=0;c<NCVARS;c++)
    {
        ul[c+RHO_IND]   = phi(iv,  c+RHO_INDX);
        ulm1[c+RHO_IND] = phi(ivm1,c+RHO_INDX);
        ulm2[c+RHO_IND] = phi(ivm2,c+RHO_INDX);

        ur[c+RHO_IND]   = phi(ivp1, c+RHO_INDX);
        urp1[c+RHO_IND] = phi(ivp2, c+RHO_INDX);
        urp2[c+RHO_IND] = phi(ivp3, c+RHO_INDX);
    }
    for(int c=0;c<NUM_SPECIES;c++)
    {
        ul[c+NCVARS]   = phi(iv,  c+FLO_NVARS);
        ulm1[c+NCVARS] = phi(ivm1,c+FLO_NVARS);
        ulm2[c+NCVARS] = phi(ivm2,c+FLO_NVARS);

        ur[c+NCVARS]   = phi(ivp1, c+FLO_NVARS);
        urp1[c+NCVARS] = phi(ivp2, c+FLO_NVARS);
        urp2[c+NCVARS] = phi(ivp3, c+FLO_NVARS);
    }
    
    geteulerflux(ul,fl,normal);
    geteulerflux(ulm1,flm1,normal);
    geteulerflux(ulm2,flm2,normal);
    
    geteulerflux(ur,fr,normal);
    geteulerflux(urp1,frp1,normal);
    geteulerflux(urp2,frp2,normal);

    weno_reconstruct(flm2,flm1,fl,fr,frp1,frp2,fmhalf,fphalf);
    weno_reconstruct(ulm2,ulm1,ul,ur,urp1,urp2,umhalf,uphalf);

    for(int c=0;c<NCVARS;c++)
    {
        fx(iv,c+RHO_INDX)=half*(fphalf[c]+lambda_max*uphalf[c]) +
                          half*(fmhalf[c]-lambda_max*umhalf[c]);
    }
    
    for(int sp=0;sp<NUM_SPECIES;sp++)
    {
        fx(iv,FLO_NVARS+sp) = mflo_species::advect_flags[sp]*
                             (half*(fphalf[NCVARS+sp] + lambda_max*uphalf[NCVARS+sp]) +
                              half*(fmhalf[NCVARS+sp] - lambda_max*umhalf[NCVARS+sp]) );
    }
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void compute_flux(int i, int j, int k,int sweepdir,
                  Array4<Real> const& phi, //state variable
                  Array4<Real> const& fluid_transpcoeffs,
                  Array4<Real> const& specdiff,
                  Array4<Real> const& fx,
                  const GpuArray<Real, AMREX_SPACEDIM>& dx,
                  int hyperbolics_order,Real dissfactor,bool nsflag,
                  int spec_in_solid)
{
    //|--left1--|--left--|--right--|--right1--|
    //|--(i-2)--|--(i-1)--|--i--|--(i+1)--|
    //remember: face i separates (i-1) and i cells

    Real ul[NCVARS+NUM_SPECIES],ur[NCVARS+NUM_SPECIES];
    Real fhalf[NCVARS+NUM_SPECIES]={zeroval};
    Real pl[NCVARS],pr[NCVARS];
    Real ulm1[NCVARS+NUM_SPECIES],urp1[NCVARS+NUM_SPECIES];
    Real normal[AMREX_SPACEDIM]={zeroval};
    Real gradu[3],gradv[3],gradw[3];
    Real dudx,dudy,dudz;
    Real dvdx,dvdy,dvdz;
    Real dwdx,dwdy,dwdz;
    Real delu;
    Real dTdn;
    Real tau[AMREX_SPACEDIM][AMREX_SPACEDIM];
    Real visc,thcond;
    Real spec_dcoeff[NUM_SPECIES];

    normal[sweepdir] = one;

    IntVect left(i, j, k);
    IntVect right(i, j, k);
    IntVect left1(i, j, k);
    IntVect right1(i, j, k);

    IntVect top_left(i, j, k);
    IntVect bot_left(i, j, k);
    IntVect top_right(i, j, k);
    IntVect bot_right(i, j, k);

    IntVect frt_left(i, j, k);
    IntVect bck_left(i, j, k);
    IntVect frt_right(i, j, k);
    IntVect bck_right(i, j, k);

    int trans1dir = (sweepdir + 1) % AMREX_SPACEDIM;
    int trans2dir = (sweepdir + 2) % AMREX_SPACEDIM;

    left[sweepdir]          -= 1;
    top_left[sweepdir]      -= 1;
    bot_left[sweepdir]      -= 1;
    frt_left[sweepdir]      -= 1;
    bck_left[sweepdir]      -= 1;

    top_left[trans1dir]     += 1;
    top_right[trans1dir]    += 1;
    bot_left[trans1dir]     -= 1;
    bot_right[trans1dir]    -= 1;

    frt_left[trans2dir]     += 1;
    frt_right[trans2dir]    += 1;
    bck_left[trans2dir]     -= 1;
    bck_right[trans2dir]    -= 1;

    left1[sweepdir]         -= 2;
    right1[sweepdir]        += 1;
    
    visc   = harmonic_mean(fluid_transpcoeffs(i,j,k,VISC_INDX),fluid_transpcoeffs(i-1,j,k,VISC_INDX)); 
    thcond = harmonic_mean(fluid_transpcoeffs(i,j,k,THCOND_INDX),fluid_transpcoeffs(i-1,j,k,THCOND_INDX));
    for(int sp=0;sp<NUM_SPECIES;sp++)
    {
        spec_dcoeff[sp]=harmonic_mean(specdiff(left,sp),specdiff(right,sp));
    } 

    //viscous dissipation term
    Real facevel[AMREX_SPACEDIM];
    facevel[XDIR] = half*(phi(left,VELX_INDX) + phi(right,VELX_INDX));
    facevel[YDIR] = half*(phi(left,VELY_INDX) + phi(right,VELY_INDX));
    facevel[ZDIR] = half*(phi(left,VELZ_INDX) + phi(right,VELZ_INDX));
    
    for(int c=0;c<NCVARS;c++)
    {
        ul[c+RHO_IND] = phi(left,c+RHO_INDX);
        ur[c+RHO_IND] = phi(right,c+RHO_INDX);
        ulm1[c+RHO_IND] = phi(left1,  c+RHO_INDX);
        urp1[c+RHO_IND] = phi(right1, c+RHO_INDX);
    }
    for(int c=0;c<NUM_SPECIES;c++)
    {
        ul[c+NCVARS] = phi(left,c+FLO_NVARS);
        ur[c+NCVARS] = phi(right,c+FLO_NVARS);
        ulm1[c+NCVARS] = phi(left1,c+FLO_NVARS);
        urp1[c+NCVARS] = phi(right1,c+FLO_NVARS);
    }

    if(nsflag)
    {    

        if(   phi(right1,VFRAC_INDX) < one || phi(right,VFRAC_INDX) < one
           || phi(left, VFRAC_INDX)  < one || phi(left1,VFRAC_INDX) < one)
        {
            //make things first order
            for(int c=0;c<(NCVARS+NUM_SPECIES);c++)
            {
                ulm1[c+RHO_IND]=ul[c+RHO_IND];
                urp1[c+RHO_IND]=ur[c+RHO_IND];
            }
        }

        if(hyperbolics_order == 2)
        {
            get_higherorder_states(ul,ur,ulm1,urp1);
        }
        
        if(hyperbolics_order < 4)
        {
            ausmp_up_flux(ul,ur,fhalf,normal);
        }
        else
        {
            central_with_diss_flux(ul,ur,ulm1,urp1,fhalf,normal,dissfactor);
        }

        for(int c=0;c<NCVARS;c++)
        {
            fx(i,j,k,c+RHO_INDX) = fhalf[c+RHO_IND];
        }

        gradu[0]  = (phi(right,VELX_INDX) - phi(left,VELX_INDX))/dx[sweepdir];

        gradu[1] = ( half*(phi(top_left,VELX_INDX) + phi(top_right,VELX_INDX)) 
                    - half*(phi(bot_left,VELX_INDX) + phi(bot_right,VELX_INDX)) )/(two*dx[trans1dir]);

        gradu[2] = ( half*(phi(frt_left,VELX_INDX) + phi(frt_right,VELX_INDX)) 
                    - half*(phi(bck_left,VELX_INDX) + phi(bck_right,VELX_INDX)) )/(two*dx[trans2dir]);

        gradv[0]  = (phi(right,VELY_INDX) - phi(left,VELY_INDX))/dx[sweepdir];

        gradv[1] = ( half*(phi(top_left,VELY_INDX) + phi(top_right,VELY_INDX)) 
                    - half*(phi(bot_left,VELY_INDX) + phi(bot_right,VELY_INDX)) )/(two*dx[trans1dir]);

        gradv[2] = ( half*(phi(frt_left,VELY_INDX) + phi(frt_right,VELY_INDX)) 
                    - half*(phi(bck_left,VELY_INDX) + phi(bck_right,VELY_INDX)) )/(two*dx[trans2dir]);

        gradw[0]  = (phi(right,VELZ_INDX) - phi(left,VELZ_INDX))/dx[sweepdir];

        gradw[1] = ( half*(phi(top_left,VELZ_INDX) + phi(top_right,VELZ_INDX)) 
                    - half*(phi(bot_left,VELZ_INDX) + phi(bot_right,VELZ_INDX)) )/(two*dx[trans1dir]);

        gradw[2] = ( half*(phi(frt_left,VELZ_INDX) + phi(frt_right,VELZ_INDX)) 
                    - half*(phi(bck_left,VELZ_INDX) + phi(bck_right,VELZ_INDX)) )/(two*dx[trans2dir]);

        dudx=gradu[GET_XDIR(sweepdir)];
        dudy=gradu[GET_YDIR(sweepdir)];
        dudz=gradu[GET_ZDIR(sweepdir)];

        dvdx=gradv[GET_XDIR(sweepdir)];
        dvdy=gradv[GET_YDIR(sweepdir)];
        dvdz=gradv[GET_ZDIR(sweepdir)];

        dwdx=gradw[GET_XDIR(sweepdir)];
        dwdy=gradw[GET_YDIR(sweepdir)];
        dwdz=gradw[GET_ZDIR(sweepdir)];

        dTdn = (phi(right,TEMP_INDX) - phi(left,TEMP_INDX))/dx[sweepdir];

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

        fx(i,j,k,RHOE_INDX) -= thcond*dTdn;


        //this term is the work done on the control volume by viscous forces - tau_{ji} n_j v_i
        fx(i,j,k,RHOE_INDX) -= (tau[XDIR][XDIR]*normal[XDIR] + 
                                tau[YDIR][XDIR]*normal[YDIR] +
                                tau[ZDIR][XDIR]*normal[ZDIR] )*facevel[XDIR];

        fx(i,j,k,RHOE_INDX) -= (tau[XDIR][YDIR]*normal[XDIR] + 
                                tau[YDIR][YDIR]*normal[YDIR] +
                                tau[ZDIR][YDIR]*normal[ZDIR] )*facevel[YDIR];

        fx(i,j,k,RHOE_INDX) -= (tau[XDIR][ZDIR]*normal[XDIR] + 
                                tau[YDIR][ZDIR]*normal[YDIR] +
                                tau[ZDIR][ZDIR]*normal[ZDIR] )*facevel[ZDIR];
    }

    //Species convection and diffusion
    //convection
    for(int sp=0;sp<NUM_SPECIES;sp++)
    {
        Real phi_im2   = phi(left1,FLO_NVARS+sp);
        Real phi_im1   = phi(left,FLO_NVARS+sp);
        Real phi_i     = phi(right,FLO_NVARS+sp);
        Real phi_ip1   = phi(right1,FLO_NVARS+sp);

        Real phiL=phi_im1;
        Real phiR=phi_i;

        if(!spec_in_solid)
        {
           if(   phi(right1,VFRAC_INDX) < one || phi(right,VFRAC_INDX) < one
             || phi(left, VFRAC_INDX)  < one || phi(left1,VFRAC_INDX) < one)
           {
                //make things first order
                phi_im2=phi_im1;
                phi_ip1=phi_i;
           }
        }
        
        get_higherorder_states(phiL,phiR,phi_im2,phi_ip1);
        Real velface=facevel[sweepdir];
     
        //if there is transport of species in solid
        //dont use the reflected velocities
        if(spec_in_solid)
        {
            Real weightL=(phi(left,VFRAC_INDX)<one)?zeroval:one;
            Real weightR=(phi(left,VFRAC_INDX)<one)?zeroval:one;

            velface=(weightL*phi(left,VELX_INDX+sweepdir)
            +weightR*phi(right,VELX_INDX+sweepdir))/(weightL+weightR+TINYVAL);
            
        }
        
        fx(i,j,k,FLO_NVARS+sp) = mflo_species::advect_flags[sp]*
        (half*(velface+fabs(velface))*phiL +
         half*(velface-fabs(velface))*phiR);
    }

    for(int sp=0;sp<NUM_SPECIES;sp++)
    {
        fx(i,j,k,FLO_NVARS+sp) -= spec_dcoeff[sp]*(ur[NCVARS+sp]-ul[NCVARS+sp])/dx[sweepdir];
    }

}
#endif

#ifndef _inviscid_flux_schemes_H_
#define _inviscid_flux_schemes_H_
#include <globalDefines.H>
#include <mflo_constants.H>
#include <consprim_funcs.H>
#include <thermo.H>
#include <species.H>

//==========================================================================================
AMREX_GPU_DEVICE
    AMREX_FORCE_INLINE
amrex::Real get_max_eigenvalue(Real ul[NCVARS+NUM_SPECIES],Real ur[NCVARS+NUM_SPECIES],Real n[AMREX_SPACEDIM])
{
    Real presL,presR;
    Real rhoL,rhoR;
    Real ahalf;
    Real veln_L,veln_R,veln;
    Real gamaL,gamaR;
    Real rhoeL,rhoeR;

    rhoL=ul[RHO_IND];
    rhoR=ur[RHO_IND];
    
    rhoeL=ul[RHOE_IND] - half*(ul[RHOU_IND]*ul[RHOU_IND] +
                               ul[RHOV_IND]*ul[RHOV_IND] +
                               ul[RHOW_IND]*ul[RHOW_IND])/ul[RHO_IND];
    
    rhoeR=ur[RHOE_IND] - half*(ur[RHOU_IND]*ur[RHOU_IND] +
                               ur[RHOV_IND]*ur[RHOV_IND] +
                               ur[RHOW_IND]*ur[RHOW_IND])/ur[RHO_IND];

    presL = mflo_thermo::get_p_from_rec(rhoeL,ul[RHO_IND],ul+NCVARS);
    presR = mflo_thermo::get_p_from_rec(rhoeR,ur[RHO_IND],ul+NCVARS);

    veln_L  = (ul[RHOU_IND]*n[0] + ul[RHOV_IND]*n[1] + ul[RHOW_IND]*n[2])/ul[RHO_IND];
    veln_R  = (ur[RHOU_IND]*n[0] + ur[RHOV_IND]*n[1] + ur[RHOW_IND]*n[2])/ur[RHO_IND];

    veln = std::max(fabs(veln_L),fabs(veln_R));

    gamaL=mflo_thermo::get_gama_from_rpc(rhoL,presL,ul+NCVARS);
    gamaR=mflo_thermo::get_gama_from_rpc(rhoR,presR,ur+NCVARS);

    ahalf = std::max(sqrt(gamaL*presL/rhoL), sqrt(gamaR*presR/rhoR));

    return(ahalf+veln);
}
//==========================================================================================
AMREX_GPU_DEVICE
    AMREX_FORCE_INLINE
void laxf_flux(Real ul[NCVARS+NUM_SPECIES],Real ur[NCVARS+NUM_SPECIES],Real fhalf[NCVARS+NUM_SPECIES],Real n[AMREX_SPACEDIM])
{
    Real fl[NCVARS],fr[NCVARS];
    Real presL,presR;
    Real rhoL,rhoR;
    Real ahalf;
    Real veln_L,veln_R,veln;
    Real gamaL,gamaR;
    Real rhoeL,rhoeR;

    geteulerflux(ul,fl,n);
    geteulerflux(ur,fr,n);

    rhoL=ul[RHO_IND];
    rhoR=ur[RHO_IND];
    
    rhoeL=ul[RHOE_IND] - half*(ul[RHOU_IND]*ul[RHOU_IND] +
                               ul[RHOV_IND]*ul[RHOV_IND] +
                               ul[RHOW_IND]*ul[RHOW_IND])/ul[RHO_IND];
    
    rhoeR=ur[RHOE_IND] - half*(ur[RHOU_IND]*ur[RHOU_IND] +
                               ur[RHOV_IND]*ur[RHOV_IND] +
                               ur[RHOW_IND]*ur[RHOW_IND])/ur[RHO_IND];

    presL = mflo_thermo::get_p_from_rec(rhoeL,ul[RHO_IND],ul+NCVARS);
    presR = mflo_thermo::get_p_from_rec(rhoeR,ur[RHO_IND],ul+NCVARS);

    /*presL = (GAMA-one)*(ul[RHOE_IND] 
            - half*(ul[RHOU_IND]*ul[RHOU_IND] 
                + ul[RHOV_IND]*ul[RHOV_IND]
                + ul[RHOW_IND]*ul[RHOW_IND])/ul[RHO_IND]);

    presR = (GAMA-one)*(ur[RHOE_IND] 
            - half*(ur[RHOU_IND]*ur[RHOU_IND] 
                + ur[RHOV_IND]*ur[RHOV_IND]
                + ur[RHOW_IND]*ur[RHOW_IND])/ur[RHO_IND]);*/


    veln_L  = (ul[RHOU_IND]*n[0] + ul[RHOV_IND]*n[1] + ul[RHOW_IND]*n[2])/ul[RHO_IND];
    veln_R  = (ur[RHOU_IND]*n[0] + ur[RHOV_IND]*n[1] + ur[RHOW_IND]*n[2])/ur[RHO_IND];

    veln = std::max(fabs(veln_L),fabs(veln_R));

    gamaL=mflo_thermo::get_gama_from_rpc(rhoL,presL,ul+NCVARS);
    gamaR=mflo_thermo::get_gama_from_rpc(rhoR,presR,ur+NCVARS);

    ahalf = std::max(sqrt(gamaL*presL/rhoL), sqrt(gamaR*presR/rhoR));

    //FIXME: include species flux for consistency
    for(int i=0;i<NCVARS;i++)
    {
        fhalf[i] = half*(fl[i] + fr[i]) - half*(veln+ahalf)*(ur[i] - ul[i]);
    }
}
//==========================================================================================
AMREX_GPU_DEVICE
    AMREX_FORCE_INLINE
Real split_Mach(Real M,Real beta,Real plus_or_minus)
{
    Real Msplit;

    //Msplit(+-) = (+-)0.25*(M (+-) 1)**2   !subsonic
    //Msplit(+-) = 0.5 * (M (+-) |M|)       !supersonic

    if(plus_or_minus < 0)
    {

        if(fabs(M) <= one)
        {
            Msplit = -fourth*(M-one)*(M-one)
                -beta*(M*M-one)*(M*M-one);
        }
        else
        {
            Msplit = half*(M-fabs(M));
        }
    }
    else
    {
        if(fabs(M) <= 1)
        {
            Msplit = fourth*(M+one)*(M+one)
                +beta*(M*M-one)*(M*M-one);
        }
        else
        {
            Msplit = half*(M+fabs(M));
        }
    }

    return(Msplit);
}
//==========================================================================================
AMREX_GPU_DEVICE
    AMREX_FORCE_INLINE
Real split_pres(Real M,Real alpha,Real plus_or_minus)
{

    Real Psplit;

    //Psplit(+-) = 0.25 (M (+-) 1)**2 (2 (-+) M)  !subsonic
    //Psplit(+-) = 0.5 * (M (+-) |M|)/M           !supersonic

    if(plus_or_minus < 0)
    {

        if(fabs(M) <= one)
        {
            Psplit = fourth * (M-one)*(M-one) * (two+M) 
                - alpha*M*(M*M-one)*(M*M-one);
        }
        else
        {
            Psplit =  half * (M-fabs(M))/M;
        }
    }
    else
    {
        if(fabs(M) <= one)
        {
            Psplit = fourth * (M+one)*(M+one)*(two-M) 
                + alpha*M*(M*M-one)*(M*M-one);
        }
        else
        {
            Psplit = half * (M+fabs(M))/M;
        }
    }

    return(Psplit);

}
//==========================================================================================
AMREX_GPU_DEVICE
    AMREX_FORCE_INLINE
void ausmp_up_flux(Real ul[NCVARS+NUM_SPECIES],Real ur[NCVARS+NUM_SPECIES],
        Real fhalf[NCVARS+NUM_SPECIES],Real n[AMREX_SPACEDIM],Real Minf=threetenth)
{
    Real fluxL[NCVARS+NUM_SPECIES],fluxR[NCVARS+NUM_SPECIES];
    Real pl[NCVARS],pr[NCVARS];
    Real nx,ny,nz;
    Real vL,vR,ML,MR,HL,HR;
    Real M_half,MLplus,MRminus;
    Real PLplus,PRminus;
    Real astar2,ahatL,ahatR,a_half;
    Real mod_M_half,rho_half;
    Real gamaL,gamaR;

    Real M_bar2,Mo2,Mo,fa_Mo,Minf2;
    Real alpha,p_half;

    Real delta   = threetenth;
    Real beta    = eighth;
    Real sigma   = one;

    //Kp,Ku=0.d0 will make this scheme ausm_plus
    Real Kp=fourth;
    Real Ku=threefourth;

    Minf2=Minf*Minf;

    nx = n[XDIR];
    ny = n[YDIR];
    nz = n[ZDIR];

    cons_to_prim(ul,pl);
    cons_to_prim(ur,pr);

    gamaL=mflo_thermo::get_gama_from_rpc(pl[DENS_IND],pl[PRES_IND],ul+NCVARS);
    gamaR=mflo_thermo::get_gama_from_rpc(pr[DENS_IND],pr[PRES_IND],ur+NCVARS);

    vL = pl[VELX_IND]*nx + pl[VELY_IND]*ny + pl[VELZ_IND]*nz;
    vR = pr[VELX_IND]*nx + pr[VELY_IND]*ny + pr[VELZ_IND]*nz;

    HL = gamaL/(gamaL-one)*(pl[PRES_IND]/pl[DENS_IND]) +
        half*(pl[VELX_IND]*pl[VELX_IND] + pl[VELY_IND]*pl[VELY_IND]
                + pl[VELZ_IND]*pl[VELZ_IND]);

    HR = gamaR/(gamaR-one)*(pr[PRES_IND]/pr[DENS_IND]) + 
        half*(pr[VELX_IND]*pr[VELX_IND] + pr[VELY_IND]*pr[VELY_IND]
                + pr[VELZ_IND]*pr[VELZ_IND]);

    rho_half = half*(pl[DENS_IND]+pr[DENS_IND]);

    //left side
    astar2 = two*(gamaL-one)/(gamaL+one)*HL;
    ahatL  = astar2/std::max(sqrt(astar2),fabs(vL));

    //right side
    astar2 = two*(gamaR-one)/(gamaR+one)*HR;
    ahatR  = astar2/max(sqrt(astar2),fabs(vR));

    a_half = min(ahatL,ahatR);

    ML = vL/a_half;
    MR = vR/a_half;

    MLplus  = split_Mach(ML, beta,  one);
    MRminus = split_Mach(MR, beta, -one);

    M_half = MLplus + MRminus;

    //ausm_plus_up dissipation term========================
    M_bar2 = half*(vL*vL+vR*vR)/(a_half*a_half);
    Mo2    = std::min(one,std::max(M_bar2,Minf2));
    Mo     = sqrt(Mo2);
    fa_Mo  = Mo*(two-Mo);

    M_half = M_half - (Kp/fa_Mo)*std::max(one-sigma*M_bar2,zeroval) *
        (pr[PRES_IND]-pl[PRES_IND])/(rho_half*a_half*a_half);
    //======================================================

    fluxL[RHO_IND]   = ul[RHO_IND];
    fluxL[RHOU_IND]  = ul[RHOU_IND];
    fluxL[RHOV_IND]  = ul[RHOV_IND];
    fluxL[RHOW_IND]  = ul[RHOW_IND];
    fluxL[RHOE_IND]  = ul[RHOE_IND]+pl[PRES_IND];

    fluxR[RHO_IND]   = ur[RHO_IND];
    fluxR[RHOU_IND]  = ur[RHOU_IND];
    fluxR[RHOV_IND]  = ur[RHOV_IND];
    fluxR[RHOW_IND]  = ur[RHOW_IND];
    fluxR[RHOE_IND]  = ur[RHOE_IND]+pr[PRES_IND];
    
    for(int sp=0;sp<NUM_SPECIES;sp++)
    {
        fluxL[NCVARS+sp]=ul[NCVARS+sp];
        fluxR[NCVARS+sp]=ur[NCVARS+sp];
    }

    mod_M_half=fabs(M_half);

    //fix if M_half -> 0
    if(mod_M_half <= delta)
    {
        mod_M_half=(M_half*M_half+delta*delta)/(two*delta);
    }

    for(int i=0;i<(NCVARS+NUM_SPECIES);i++)
    {
        fhalf[i] = (half*M_half*(fluxL[i]+fluxR[i]) 
                - half*mod_M_half*(fluxR[i]-fluxL[i]))*a_half;
    }

    //pressure term
    //alpha=3.d0/16.d0 (original ausm_plus)
    alpha   = three*sixteenth*(-four + five*fa_Mo*fa_Mo);
    PLplus  = split_pres(ML, alpha,  one);
    PRminus = split_pres(MR, alpha, -one);
    p_half  = (PLplus*pl[PRES_IND]+PRminus*pr[PRES_IND]);

    //ausm_plus_up dissipation term========================
    p_half = p_half - Ku*PLplus*PRminus*(pl[DENS_IND]+pr[DENS_IND])
        *(fa_Mo*a_half)*(vR-vL);
    //======================================================

    fhalf[RHOU_IND] +=  p_half*nx;
    fhalf[RHOV_IND] +=  p_half*ny;
    fhalf[RHOW_IND] +=  p_half*nz;

}
//==========================================================================================
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void central_with_diss_flux(Real uim1[NCVARS+NUM_SPECIES],Real ui[NCVARS+NUM_SPECIES],
                            Real uim2[NCVARS+NUM_SPECIES],Real uip1[NCVARS+NUM_SPECIES],
        Real fhalf[NCVARS+NUM_SPECIES],Real n[AMREX_SPACEDIM],Real dissfactor)
{
    Real fim1[NCVARS+NUM_SPECIES],fi[NCVARS+NUM_SPECIES],
         fim2[NCVARS+NUM_SPECIES],fip1[NCVARS+NUM_SPECIES];
    Real presL,presR;
    Real rhoL,rhoR;
    Real ahalf;
    Real veln_L,veln_R,veln;
    Real fourth_order_flux,fourth_order_diss;
    Real rhoeL,rhoeR,gamaL,gamaR;

    geteulerflux(uim2,fim2,n);
    geteulerflux(uim1,fim1,n);
    geteulerflux(ui,fi,n);
    geteulerflux(uip1,fip1,n);

    rhoL=uim1[RHO_IND];
    rhoR=ui[RHO_IND];
    
    rhoeL=uim1[RHOE_IND] - half*(uim1[RHOU_IND]*uim1[RHOU_IND] +
                                 uim1[RHOV_IND]*uim1[RHOV_IND] +
                                 uim1[RHOW_IND]*uim1[RHOW_IND])/uim1[RHO_IND];
    
    rhoeR=ui[RHOE_IND] - half*(ui[RHOU_IND]*ui[RHOU_IND] +
                               ui[RHOV_IND]*ui[RHOV_IND] +
                               ui[RHOW_IND]*ui[RHOW_IND])/ui[RHO_IND];
    
    presL = mflo_thermo::get_p_from_rec(rhoeL,uim1[RHO_IND],uim1+NCVARS);
    presR = mflo_thermo::get_p_from_rec(rhoeR,ui[RHO_IND],ui+NCVARS);

    veln_L  = (uim1[RHOU_IND]*n[0] + uim1[RHOV_IND]*n[1] + uim1[RHOW_IND]*n[2])/uim1[RHO_IND];
    veln_R  = (ui[RHOU_IND]*n[0] + ui[RHOV_IND]*n[1] + ui[RHOW_IND]*n[2])/ui[RHO_IND];

    veln = std::max(fabs(veln_L),fabs(veln_R));
    
    gamaL=mflo_thermo::get_gama_from_rpc(rhoL,presL,uim1+NCVARS);
    gamaR=mflo_thermo::get_gama_from_rpc(rhoR,presR,ui+NCVARS);

    ahalf = std::max(sqrt(gamaL*presL/rhoL), sqrt(gamaR*presR/rhoR));

    for(int i=0;i<(NCVARS+NUM_SPECIES);i++)
    {
        fhalf[i] = twelveth*(7*fi[i]+7*fim1[i]-fim2[i]-fip1[i]);
        fhalf[i] -= twelveth*(3*ui[i]-3*uim1[i]-uip1[i]+ui[i])*(veln+ahalf)*dissfactor;
    }

}
//==========================================================================================
// WENO 5 JS (Original Jiang GS, Shu CW. Efficient Implementation of Weighted
// ENO Schemes.
//             J. Comput. Phys. 1996; 126(1): 202–228.
//             doi: 10.1006/jcph.1996.0130)
// @param s      Stencils i-2, i-1, i, i+1, i+2, i+3
// @param sm     The value of the interpolation on the left edge of the i cell
// @param sp     The value of the interpolation on the right edge of the i cell
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
weno_reconstruct(Real uim2[NCVARS+NUM_SPECIES],
                 Real uim1[NCVARS+NUM_SPECIES],
                 Real ui[NCVARS+NUM_SPECIES],
                 Real uip1[NCVARS+NUM_SPECIES],
                 Real uip2[NCVARS+NUM_SPECIES],
                 Real uip3[NCVARS+NUM_SPECIES],
                 Real umhalf[NCVARS+NUM_SPECIES],
                 Real uphalf[NCVARS+NUM_SPECIES])
{
    constexpr amrex::Real eps = 1e-6;
    amrex::Real vl[3];
    amrex::Real vr[3];
    amrex::Real beta[3];
    amrex::Real alpha[3];
    amrex::Real s[6];

    for(int c=0;c<(NCVARS+NUM_SPECIES);c++)
    {

        s[0] = uim2[c];
        s[1] = uim1[c];
        s[2] = ui[c];

        s[3] = uip1[c];
        s[4] = uip2[c];
        s[5] = uip3[c];

        beta[2] = (13.0 / 12.0) * pow(s[0] - 2.0 * s[1] + s[2], 2) +
        0.25 * pow(s[0] - 4.0 * s[1] + 3.0 * s[2], 2);
        beta[1] = (13.0 / 12.0) * pow(s[1] - 2.0 * s[2] + s[3], 2) +
        0.25 * (s[1] - s[3]) * (s[1] - s[3]);
        beta[0] = (13.0 / 12.0) * pow(s[2] - 2.0 * s[3] + s[4], 2) +
        0.25 * pow(3.0 * s[2] - 4.0 * s[3] + s[4], 2);

        beta[2] = 1.0 / ((eps + beta[2]) * (eps + beta[2]));
        beta[1] = 1.0 / ((eps + beta[1]) * (eps + beta[1]));
        beta[0] = 1.0 / ((eps + beta[0]) * (eps + beta[0]));

        alpha[2] = beta[2];
        alpha[1] = 6.0 * beta[1];
        alpha[0] = 3.0 * beta[0];
        amrex::Real alpha1 = 1.0 / (alpha[2] + alpha[1] + alpha[0]);

        vl[2] = 2.0 * s[0] - 7.0 * s[1] + 11.0 * s[2];
        vl[1] = -s[1] + 5.0 * s[2] + 2.0 * s[3];
        vl[0] = 2.0 * s[2] + 5.0 * s[3] - s[4];

        uphalf[c] = (1.0 / 6.0) * alpha1 *
        (alpha[2] * vl[2] + alpha[1] * vl[1] + alpha[0] * vl[0]);

        beta[2] = (13.0 / 12.0) * pow(s[5] - 2.0 * s[4] + s[3], 2) +
        0.25 * pow(s[5] - 4.0 * s[4] + 3.0 * s[3], 2);
        beta[1] = (13.0 / 12.0) * pow(s[4] - 2.0 * s[3] + s[2], 2) +
        0.25 * (s[4] - s[2]) * (s[4] - s[2]);
        beta[0] = (13.0 / 12.0) * pow(s[3] - 2.0 * s[2] + s[1], 2) +
        0.25 * pow(3.0 * s[3] - 4.0 * s[2] + s[1], 2);

        beta[2] = 1.0 / ((eps + beta[2]) * (eps + beta[2]));
        beta[1] = 1.0 / ((eps + beta[1]) * (eps + beta[1]));
        beta[0] = 1.0 / ((eps + beta[0]) * (eps + beta[0]));

        alpha[2] = beta[2];
        alpha[1] = 6.0 * beta[1];
        alpha[0] = 3.0 * beta[0];
        alpha1 = 1.0 / (alpha[2] + alpha[1] + alpha[0]);

        vr[2] = 11.0 * s[3] - 7.0 * s[4] + 2.0 * s[5];
        vr[1] = -s[4] + 5.0 * s[3] + 2.0 * s[2];
        vr[0] = 2.0 * s[3] + 5.0 * s[2] - s[1];

        umhalf[c] = (1.0 / 6.0) * alpha1 *
        (alpha[2] * vr[2] + alpha[1] * vr[1] + alpha[0] * vr[0]);
    }
}
//==========================================================================================
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void weno_flux(Real ulm2[NCVARS+NUM_SPECIES],
               Real ulm1[NCVARS+NUM_SPECIES],
               Real ul[NCVARS+NUM_SPECIES],
               Real ur[NCVARS+NUM_SPECIES],
               Real urp1[NCVARS+NUM_SPECIES],
               Real urp2[NCVARS+NUM_SPECIES],
               Real fhalf[NCVARS+NUM_SPECIES],
               Real normal[AMREX_SPACEDIM])
{
    //inspired from Desai et al., 
    //Direct numerical simulations of turbulent reacting flows with 
    //shock waves and stiff chemistry using many-core/GPU accelerationR
    //Computers and Fluids 215 (2021) 104787

    Real fr[NCVARS+NUM_SPECIES]={zeroval};
    Real frp1[NCVARS+NUM_SPECIES]={zeroval};
    Real frp2[NCVARS+NUM_SPECIES]={zeroval};
    Real fl[NCVARS+NUM_SPECIES]={zeroval};
    Real flm1[NCVARS+NUM_SPECIES]={zeroval};
    Real flm2[NCVARS+NUM_SPECIES]={zeroval};

    Real fphalf[NCVARS+NUM_SPECIES]={zeroval};
    Real fmhalf[NCVARS+NUM_SPECIES]={zeroval};
    Real uphalf[NCVARS+NUM_SPECIES]={zeroval};
    Real umhalf[NCVARS+NUM_SPECIES]={zeroval};

    Real lambda_max;

    lambda_max=get_max_eigenvalue(ul,ur,normal);
    geteulerflux(ul,fl,normal);
    geteulerflux(ulm1,flm1,normal);
    geteulerflux(ulm2,flm2,normal);

    geteulerflux(ur,fr,normal);
    geteulerflux(urp1,frp1,normal);
    geteulerflux(urp2,frp2,normal);

    weno_reconstruct(flm2,flm1,fl,fr,frp1,frp2,fmhalf,fphalf);
    weno_reconstruct(ulm2,ulm1,ul,ur,urp1,urp2,umhalf,uphalf);

    for(int c=0;c<(NCVARS+NUM_SPECIES);c++)
    {
        fhalf[c+RHO_IND]=half*(fphalf[c]+lambda_max*uphalf[c]) +
        half*(fmhalf[c]-lambda_max*umhalf[c]);
    }
}
//==========================================================================================
#endif

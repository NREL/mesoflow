#ifndef _consprim_funcs_h
#define _consprim_funcs_h
#include <globalDefines.H>
#include <mflo_constants.H>
#include <thermo.H>

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void cons_to_prim(Real u[NCVARS+NUM_SPECIES], Real p[NCVARS])
{
    p[DENS_IND] = u[RHO_IND];
    p[VELX_IND] = u[RHOU_IND] / u[RHO_IND];
    p[VELY_IND] = u[RHOV_IND] / u[RHO_IND];
    p[VELZ_IND] = u[RHOW_IND] / u[RHO_IND];

    //get only internal energy
    Real rhoe=u[RHOE_IND]-half*u[RHO_IND]*(p[VELX_IND]*p[VELX_IND] +
                                           p[VELY_IND]*p[VELY_IND] +
                                           p[VELZ_IND]*p[VELZ_IND]);

    p[PRES_IND] = mflo_thermo::get_p_from_rec(rhoe,u[RHO_IND],u+NCVARS);

}
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void prim_to_cons(Real u[NCVARS+NUM_SPECIES], Real p[NCVARS])
{
    u[RHO_IND] = p[DENS_IND];
    u[RHOU_IND] = p[DENS_IND] * p[VELX_IND];
    u[RHOV_IND] = p[DENS_IND] * p[VELY_IND];
    u[RHOW_IND] = p[DENS_IND] * p[VELZ_IND];

    u[RHOE_IND]=u[RHO_IND]*mflo_thermo::get_e_from_rpc(p[DENS_IND],p[PRES_IND],u+NCVARS);

    //add kinetic energy
    u[RHOE_IND] += half*u[RHO_IND]*(p[VELX_IND]*p[VELX_IND] +
                                    p[VELY_IND]*p[VELY_IND] +
                                    p[VELZ_IND]*p[VELZ_IND]);

    /*u[RHOE_IND] = p[PRES_IND] / (GAMA - one) +
                  half * p[DENS_IND] *
                      (p[VELX_IND] * p[VELX_IND] + p[VELY_IND] * p[VELY_IND] +
                       p[VELZ_IND] * p[VELZ_IND]);*/
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void geteulerflux(Real u[NCVARS+NUM_SPECIES], Real f[NCVARS], Real n[AMREX_SPACEDIM])
{
    Real nx, ny, nz;
    Real vx, vy, vz, rho, p;

    nx = n[XDIR];
    ny = n[YDIR];
    nz = n[ZDIR];

    rho = u[RHO_IND];
    vx = u[RHOU_IND] / rho;
    vy = u[RHOV_IND] / rho;
    vz = u[RHOW_IND] / rho;

    Real rhoe=u[RHOE_IND]-half*u[RHO_IND]*(vx*vx + vy*vy +vz*vz);
    //p = (GAMA - one) *
    //    (u[RHOE_IND] - half * rho * (vx * vx + vy * vy + vz * vz));

    p=mflo_thermo::get_p_from_rec(rhoe,u[RHO_IND],u+NCVARS);

    f[RHO_IND] = rho * vx * nx + rho * vy * ny + rho * vz * nz;
    f[RHOU_IND] =
        (rho * vx * vx + p) * nx + rho * vx * vy * ny + rho * vx * vz * nz;
    f[RHOV_IND] =
        rho * vy * vx * nx + (rho * vy * vy + p) * ny + rho * vy * vz * nz;
    f[RHOW_IND] =
        rho * vz * vx * nx + rho * vz * vy * ny + (rho * vz * vz + p) * nz;
    f[RHOE_IND] = (u[RHOE_IND] + p) * (vx * nx + vy * ny + vz * nz);
}

#endif

#include <AMReX_Box.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_Geometry.H>
#include <userfuncs.H>
#include <globalDefines.H>
#include <cflo_constants.H>

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void set_ib_state(int i, int j, int k, Array4<Real> const& phi)
{
    Real zerospec[NUM_SPECIES]={zero}; 

    if (phi(i, j, k, VFRAC_INDX) < one) 
    {
        phi(i, j, k, VFRAC_INDX) = zero;

        phi(i,j,k,RHO_INDX) = cflo_user_funcs::fs_rho;
        phi(i, j, k, RHOU_INDX) = zero;
        phi(i, j, k, RHOV_INDX) = zero;
        phi(i, j, k, RHOW_INDX) = zero;

        phi(i,j,k,RHOE_INDX)  = cflo_user_funcs::fs_rho*cflo_thermo
        ::get_e_from_rpc(cflo_user_funcs::fs_rho,cflo_user_funcs::fs_p,
                         zerospec);

        phi(i, j, k, DENS_INDX) = cflo_user_funcs::fs_rho;
        phi(i, j, k, VELX_INDX) = zero;
        phi(i, j, k, VELY_INDX) = zero;
        phi(i, j, k, VELZ_INDX) = zero;
        phi(i, j, k, PRES_INDX) = cflo_user_funcs::fs_p;
        phi(i, j, k, TEMP_INDX) = cflo_thermo::get_t_from_rpc(cflo_user_funcs::fs_rho,
                                                              cflo_user_funcs::fs_p,
                                                              zerospec);
        /*for(int sp=0;sp<NUM_SPECIES;sp++)
        {
            phi(i,j,k,CFLO_NVARS+sp) = zerospec[sp];
        }*/
    }
}

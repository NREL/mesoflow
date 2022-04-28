#include <AMReX_Box.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_Geometry.H>
#include <userfuncs.H>
#include <globalDefines.H>
#include <mflo_constants.H>

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void set_ib_state(int i, int j, int k, Array4<Real> const& phi)
{

    if (phi(i, j, k, VFRAC_INDX) < one) 
    {
        phi(i, j, k, VFRAC_INDX) = zeroval;

        phi(i,j,k,RHO_INDX) = mflo_user_funcs::dens_fs;
        phi(i, j, k, RHOU_INDX) = zeroval;
        phi(i, j, k, RHOV_INDX) = zeroval;
        phi(i, j, k, RHOW_INDX) = zeroval;

        //dont need to do this as it will remain the same as beginning
        phi(i,j,k,RHOE_INDX)  = phi(i,j,k,DENS_INDX)*mflo_thermo
            ::get_e_from_rpc(mflo_user_funcs::dens_fs,mflo_user_funcs::pres_fs,
                    mflo_user_funcs::spec_fs);

        phi(i, j, k, DENS_INDX) = mflo_user_funcs::dens_fs;
        phi(i, j, k, VELX_INDX) = zeroval;
        phi(i, j, k, VELY_INDX) = zeroval;
        phi(i, j, k, VELZ_INDX) = zeroval;
        phi(i, j, k, PRES_INDX) = mflo_user_funcs::pres_fs;
        phi(i, j, k, TEMP_INDX) = mflo_user_funcs::temp_fs;

        for(int sp=0;sp<NUM_SPECIES;sp++)
        {
            phi(i,j,k,FLO_NVARS+sp) = zeroval;
        }
    }
}

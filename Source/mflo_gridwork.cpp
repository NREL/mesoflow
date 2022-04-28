#include <AMReX_ParallelDescriptor.H>
#include <AMReX_ParmParse.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_FillPatchUtil.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_VisMF.H>
#include <AMReX_PhysBCFunct.H>
#include <kernels_3d.H>

#include <Prob.H>

#ifdef AMREX_MEM_PROFILING
#include <AMReX_MemProfiler.H>
#endif

#include <mflo.H>
#include <userfuncs.H>

// Make a new level using provided BoxArray and DistributionMapping and
// fill with interpolated coarse level data.
// overrides the pure virtual function in AmrCore
void mflo::MakeNewLevelFromCoarse(
    int lev, Real time, const BoxArray& ba, const DistributionMapping& dm)
{
    const int ncomp = phi_new[lev - 1].nComp();
    const int nghost = phi_new[lev - 1].nGrow();

    phi_new[lev].define(ba, dm, ncomp, nghost);
    phi_old[lev].define(ba, dm, ncomp, nghost);

    t_new[lev] = time;
    t_old[lev] = time - 1.e200;

    if (lev > 0 && do_reflux) {
        flux_reg[lev].reset(
            new FluxRegister(ba, dm, refRatio(lev - 1), lev, ncomp));
    }

    FillCoarsePatch(lev, time, phi_new[lev], 0, ncomp);
}

// Remake an existing level using provided BoxArray and DistributionMapping and
// fill with existing fine and coarse data.
// overrides the pure virtual function in AmrCore
void mflo::RemakeLevel(
    int lev, Real time, const BoxArray& ba, const DistributionMapping& dm)
{
    const int ncomp = phi_new[lev].nComp();
    const int nghost = phi_new[lev].nGrow();

    MultiFab new_state(ba, dm, ncomp, nghost);
    MultiFab old_state(ba, dm, ncomp, nghost);

    FillPatch(lev, time, new_state, 0, ncomp);

    std::swap(new_state, phi_new[lev]);
    std::swap(old_state, phi_old[lev]);

    t_new[lev] = time;
    t_old[lev] = time - 1.e200;

    if (lev > 0 && do_reflux) {
        flux_reg[lev].reset(
            new FluxRegister(ba, dm, refRatio(lev - 1), lev, ncomp));
    }
}

// Delete level data
// overrides the pure virtual function in AmrCore
void mflo::ClearLevel(int lev)
{
    phi_new[lev].clear();
    phi_old[lev].clear();
    flux_reg[lev].reset(nullptr);
}

// Make a new level from scratch using provided BoxArray and
// DistributionMapping. Only used during initialization. overrides the pure
// virtual function in AmrCore
void mflo::MakeNewLevelFromScratch(
    int lev, Real time, const BoxArray& ba, const DistributionMapping& dm)
{
    const int nghost = 0;
    int ncomp=TOTAL_NVARS;

    phi_new[lev].define(ba, dm, ncomp, nghost);
    phi_old[lev].define(ba, dm, ncomp, nghost);

    t_new[lev] = time;
    t_old[lev] = time - 1.e200;

    if (lev > 0 && do_reflux) {
        flux_reg[lev].reset(
            new FluxRegister(ba, dm, refRatio(lev - 1), lev, ncomp));
    }

    Real cur_time = t_new[lev];
    MultiFab& state = phi_new[lev];

    for (MFIter mfi(state); mfi.isValid(); ++mfi) {
        Array4<Real> fab = state[mfi].array();
        GeometryData geomData = geom[lev].data();
        const Box& box = mfi.validbox();

        amrex::launch(box,
        [=] AMREX_GPU_DEVICE (Box const& tbx)
        {
            initdata(tbx, fab, geomData);
        });
    }
}

// set covered coarse cells to be the average of overlying fine cells
void mflo::AverageDown()
{
    for (int lev = finest_level - 1; lev >= 0; --lev) {
        amrex::average_down(
            phi_new[lev + 1], phi_new[lev], geom[lev + 1], geom[lev], 0,
            phi_new[lev].nComp(), refRatio(lev));
    }
}

// more flexible version of AverageDown() that lets you average down across
// multiple levels
void mflo::AverageDownTo(int crse_lev)
{
    amrex::average_down(
        phi_new[crse_lev + 1], phi_new[crse_lev], geom[crse_lev + 1],
        geom[crse_lev], 0, phi_new[crse_lev].nComp(), refRatio(crse_lev));
}

// compute a new multifab by coping in phi from valid region and filling ghost
// cells works for single level and 2-level cases (fill fine grid ghost by
// interpolating from coarse)
void mflo::FillPatch(int lev, Real time, MultiFab& mf, int icomp, int ncomp)
{
    if (lev == 0) {
        Vector<MultiFab*> smf;
        Vector<Real> stime;
        GetData(0, time, smf, stime);

        GpuBndryFuncFab<AmrCoreFill> gpu_bndry_func(amrcore_fill_func);
        PhysBCFunct<GpuBndryFuncFab<AmrCoreFill>> physbc(
            geom[lev], bcs, gpu_bndry_func);
        amrex::FillPatchSingleLevel(
            mf, time, smf, stime, 0, icomp, ncomp, geom[lev], physbc, 0);
    } else {
        Vector<MultiFab*> cmf, fmf;
        Vector<Real> ctime, ftime;
        GetData(lev - 1, time, cmf, ctime);
        GetData(lev, time, fmf, ftime);

        Interpolater* mapper = &cell_cons_interp;

        GpuBndryFuncFab<AmrCoreFill> gpu_bndry_func(amrcore_fill_func);
        PhysBCFunct<GpuBndryFuncFab<AmrCoreFill>> cphysbc(
            geom[lev - 1], bcs, gpu_bndry_func);
        PhysBCFunct<GpuBndryFuncFab<AmrCoreFill>> fphysbc(
            geom[lev], bcs, gpu_bndry_func);

        amrex::FillPatchTwoLevels(
            mf, time, cmf, ctime, fmf, ftime, 0, icomp, ncomp, geom[lev - 1],
            geom[lev], cphysbc, 0, fphysbc, 0, refRatio(lev - 1), mapper, bcs,
            0);
    }
}

// fill an entire multifab by interpolating from the coarser level
// this comes into play when a new level of refinement appears
void mflo::FillCoarsePatch(
    int lev, Real time, MultiFab& mf, int icomp, int ncomp)
{
    BL_ASSERT(lev > 0);

    Vector<MultiFab*> cmf;
    Vector<Real> ctime;
    GetData(lev - 1, time, cmf, ctime);
    Interpolater* mapper = &cell_cons_interp;

    if (cmf.size() != 1) {
        amrex::Abort("FillCoarsePatch: how did this happen?");
    }

    GpuBndryFuncFab<AmrCoreFill> gpu_bndry_func(amrcore_fill_func);
    PhysBCFunct<GpuBndryFuncFab<AmrCoreFill>> cphysbc(
        geom[lev - 1], bcs, gpu_bndry_func);
    PhysBCFunct<GpuBndryFuncFab<AmrCoreFill>> fphysbc(
        geom[lev], bcs, gpu_bndry_func);

    amrex::InterpFromCoarseLevel(
        mf, time, *cmf[0], 0, icomp, ncomp, geom[lev - 1], geom[lev], cphysbc,
        0, fphysbc, 0, refRatio(lev - 1), mapper, bcs, 0);
}

#ifndef BCFILL_H
#define BCFILL_H

#include <AMReX_FArrayBox.H>
#include <AMReX_Geometry.H>
#include <AMReX_PhysBCFunct.H>
#include <userfuncs.H>
#include <globalDefines.H>
#include <externbc.H>

using namespace amrex;

struct AmrCoreFill
{
    AMREX_GPU_DEVICE
        void operator()(
                const IntVect& iv,
                Array4<Real> const& data,
                const int dcomp,
                const int numcomp,
                GeometryData const& geom,
                const Real time,
                const BCRec* bcr,
                const int bcomp,
                const int orig_comp) const
        {
            // do something for external Dirichlet (BCType::ext_dir)
            const int* domlo = geom.Domain().loVect();
            const int* domhi = geom.Domain().hiVect();
            const amrex::Real* prob_lo = geom.ProbLo();
            const amrex::Real* prob_hi = geom.ProbHi();
            const amrex::Real* dx = geom.CellSize();
            const amrex::Real x[AMREX_SPACEDIM] = {AMREX_D_DECL(
                    prob_lo[0] + (iv[0] + 0.5) * dx[0],
                    prob_lo[1] + (iv[1] + 0.5) * dx[1],
                    prob_lo[2] + (iv[2] + 0.5) * dx[2])};

            const int* bc = bcr->data();

            amrex::Real s_int[TOTAL_NVARS] = {0.0};
            amrex::Real s_ext[TOTAL_NVARS] = {0.0};

            // xlo and xhi
            int idir = 0;
            if ((bc[idir] == amrex::BCType::ext_dir) and (iv[idir] < domlo[idir])) 
            {
                amrex::IntVect loc(AMREX_D_DECL(domlo[idir], iv[1], iv[2]));
                for (int n = 0; n < TOTAL_NVARS; n++) 
                {
                    s_int[n] = data(loc, n);
                }
                externalbc(x, s_int, s_ext, idir, +1, time, geom);
                for (int n = 0; n < TOTAL_NVARS; n++) 
                {
                    data(iv, n) = s_ext[n];
                }
            } 
            else if ((bc[idir + AMREX_SPACEDIM] == amrex::BCType::ext_dir) and
                    (iv[idir] > domhi[idir])) 
            {
                amrex::IntVect loc(AMREX_D_DECL(domhi[idir], iv[1], iv[2]));
                for (int n = 0; n < TOTAL_NVARS; n++) 
                {
                    s_int[n] = data(loc, n);
                }
                externalbc(x, s_int, s_ext, idir, -1, time, geom);
                for (int n = 0; n < TOTAL_NVARS; n++) 
                {
                    data(iv, n) = s_ext[n];
                }
            }
            // ylo and yhi
            idir = 1;
            if ((bc[idir] == amrex::BCType::ext_dir) and (iv[idir] < domlo[idir])) 
            {
                amrex::IntVect loc(AMREX_D_DECL(iv[0], domlo[idir], iv[2]));
                for (int n = 0; n < TOTAL_NVARS; n++) 
                {
                    s_int[n] = data(loc, n);
                }
                externalbc(x, s_int, s_ext, idir, +1, time, geom);
                for (int n = 0; n < TOTAL_NVARS; n++) 
                {
                    data(iv, n) = s_ext[n];
                }
            } 
            else if ((bc[idir + AMREX_SPACEDIM] == amrex::BCType::ext_dir) and
                    (iv[idir] > domhi[idir])) 
            {
                amrex::IntVect loc(AMREX_D_DECL(iv[0], domhi[idir], iv[2]));
                for (int n = 0; n < TOTAL_NVARS; n++) 
                {
                    s_int[n] = data(loc, n);
                }
                externalbc(x, s_int, s_ext, idir, -1, time, geom);
                for (int n = 0; n < TOTAL_NVARS; n++) 
                {
                    data(iv, n) = s_ext[n];
                }
            }
            // zlo and zhi
            idir = 2;
            if ((bc[idir] == amrex::BCType::ext_dir) and (iv[idir] < domlo[idir])) 
            {
                for (int n = 0; n < TOTAL_NVARS; n++) 
                {
                    s_int[n] = data(iv[0], iv[1], domlo[idir], n);
                }
                externalbc(x, s_int, s_ext, idir, +1, time, geom);
                for (int n = 0; n < TOTAL_NVARS; n++) 
                {
                    data(iv, n) = s_ext[n];
                }
            } 
            else if ((bc[idir + AMREX_SPACEDIM] == amrex::BCType::ext_dir) and
                    (iv[idir] > domhi[idir])) 
            {
                for (int n = 0; n < TOTAL_NVARS; n++) 
                {
                    s_int[n] = data(iv[0], iv[1], domhi[idir], n);
                }
                externalbc(x, s_int, s_ext, idir, -1, time, geom);
                for (int n = 0; n < TOTAL_NVARS; n++) 
                {
                    data(iv, n) = s_ext[n];
                }
            }
        }
};

namespace {
    static AmrCoreFill amrcore_fill_func;
}

#endif

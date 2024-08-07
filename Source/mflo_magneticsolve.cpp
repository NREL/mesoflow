#include <AMReX_ParallelDescriptor.H>
#include <AMReX_ParmParse.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_FillPatchUtil.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_VisMF.H>
#include <AMReX_PhysBCFunct.H>
#include <AMReX_MLTensorOp.H>
#include <mflo.H>
#include <AMReX_MLABecLaplacian.H>
#include <AMReX_MLMG.H>
#ifdef AMREX_USE_HYPRE
#include <AMReX_Hypre.H>
#endif

void mflo::solve_magnetic_vecpot(Real current_time,int Adir)
{
    BL_PROFILE("mflo::solve_magnetic_vecpot()");
    
    // FIXME: add these as inputs
    int max_coarsening_level=10;
    int max_iter=200;
    int use_hypre=0;
    Real ascalar = 1.0;
    Real bscalar = 1.0;
    Real reltol=1e-10;
    Real abstol=1e-10;
    
    int linsolve_verbose=1;
    int dircn=Adir; //gpu capture
    
    ParmParse pp("mflo.potsolve");
    pp.query("max_coarsening_level",max_coarsening_level);
    pp.query("max_iter",max_iter);
    pp.query("use_hypre",use_hypre);
    pp.query("reltol",reltol);
    pp.query("abstol",abstol);

    //==================================================
    // amrex solves
    // read small a as alpha, b as beta

    //(A a - B del.(b del)) phi = f
    //
    // A and B are scalar constants
    // a and b are scalar fields
    // f is rhs
    // in this case: A=0,a=0,B=1,b=conductivity
    // note also the negative sign
    //====================================================

    const Real tol_rel = reltol;
    const Real tol_abs = abstol;
#ifdef AMREX_USE_HYPRE
    if(use_hypre)
    {
        amrex::Print()<<"using hypre\n";
    }
#endif

    std::array<LinOpBCType, AMREX_SPACEDIM> bc_potsolve_lo 
    = {AMREX_D_DECL(LinOpBCType::Robin, LinOpBCType::Robin, LinOpBCType::Robin)}; 

    std::array<LinOpBCType, AMREX_SPACEDIM> bc_potsolve_hi 
    = {AMREX_D_DECL(LinOpBCType::Robin, LinOpBCType::Robin, LinOpBCType::Robin)}; 
    
    int mixedbc=0;
    for (int idim = 0; idim < AMREX_SPACEDIM; idim++)
    {
        //lower side bcs
        if (bc_lo[idim] == PERBC)
        {
            bc_potsolve_lo[idim] = LinOpBCType::Periodic;
        }
        else
        {
            bc_potsolve_lo[idim] = LinOpBCType::Robin;
            mixedbc=1;
        }

        //higher side bcs
        if (bc_hi[idim] == PERBC)
        {
            bc_potsolve_hi[idim] = LinOpBCType::Periodic;
        }
        else
        {
            bc_potsolve_hi[idim] = LinOpBCType::Robin;
            mixedbc=1;
        }
    }

    Vector<MultiFab> potential(finest_level+1);
    Vector<MultiFab> acoeff(finest_level+1);
    Vector<MultiFab> solution(finest_level+1);
    Vector<MultiFab> rhs(finest_level+1);

    Vector<MultiFab> robin_a(finest_level+1);
    Vector<MultiFab> robin_b(finest_level+1);
    Vector<MultiFab> robin_f(finest_level+1);

    const int num_grow = 1;

    for (int ilev = 0; ilev <= finest_level; ilev++)
    {
        potential[ilev].define(grids[ilev], dmap[ilev], 1, num_grow);
        acoeff[ilev].define(grids[ilev], dmap[ilev], 1, 0);
        solution[ilev].define(grids[ilev], dmap[ilev], 1, num_grow);
        rhs[ilev].define(grids[ilev], dmap[ilev], 1, 0);

        robin_a[ilev].define(grids[ilev], dmap[ilev], 1, 1);
        robin_b[ilev].define(grids[ilev], dmap[ilev], 1, 1);
        robin_f[ilev].define(grids[ilev], dmap[ilev], 1, 1);
    }

    LPInfo info;
    info.setAgglomeration(true);
    info.setConsolidation(true);
    info.setMaxCoarseningLevel(max_coarsening_level);
    MLABecLaplacian linsolve_obj(Geom(0,finest_level), 
                               boxArray(0,finest_level), 
                 DistributionMap(0,finest_level), info);

    linsolve_obj.setMaxOrder(2);
    linsolve_obj.setDomainBC(bc_potsolve_lo, bc_potsolve_hi);
    linsolve_obj.setScalars(ascalar, bscalar);

    for (int ilev = 0; ilev <= finest_level; ilev++)
    {
        MultiFab Sborder(grids[ilev], dmap[ilev], phi_new[ilev].nComp(), num_grow);
        FillPatch(ilev, current_time, Sborder, 0, Sborder.nComp());
        potential[ilev].setVal(0.0);

        // Copy (FabArray<FAB>& dst, FabArray<FAB> const& src, int srccomp, 
        // int dstcomp, int numcomp, const IntVect& nghost)
        amrex::Copy(potential[ilev], Sborder, AX_INDX+Adir, 0, 1, num_grow);

        solution[ilev].setVal(0.0);
        amrex::MultiFab::Copy(solution[ilev], potential[ilev], 0, 0, 1, 1);
        rhs[ilev].setVal(0.0);
        acoeff[ilev].setVal(0.0);

        //default to homogenous Neumann
        robin_a[ilev].setVal(0.0);
        robin_b[ilev].setVal(1.0);
        robin_f[ilev].setVal(0.0);

        // Get the boundary ids
        const int* domlo_arr = geom[ilev].Domain().loVect();
        const int* domhi_arr = geom[ilev].Domain().hiVect();

        GpuArray<int,AMREX_SPACEDIM> domlo={AMREX_D_DECL(domlo_arr[0], domlo_arr[1], domlo_arr[2])};
        GpuArray<int,AMREX_SPACEDIM> domhi={AMREX_D_DECL(domhi_arr[0], domhi_arr[1], domhi_arr[2])};

        // fill dcoeff and rhs
        for (MFIter mfi(phi_new[ilev], TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const Box& bx = mfi.tilebox();
            const Box& gbx = amrex::grow(bx, 1);
            const auto dx = geom[ilev].CellSizeArray();
            auto prob_lo = geom[ilev].ProbLoArray();
            auto prob_hi = geom[ilev].ProbHiArray();
            const Box& domain = geom[ilev].Domain();

            Real time = current_time; // for GPU capture

            Array4<Real> phi_arr = Sborder.array(mfi);
            Array4<Real> rhs_arr = rhs[ilev].array(mfi);

            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {

                rhs_arr(i,j,k)=0.0;
                //add mu0 J here
               amrex::Real sigma=mflo_thermo::conductivity(phi_arr(i,j,k,TEMP_INDX),
                                                         phi_arr(i,j,k,VFRAC_INDX),phi_arr(i,j,k,FLO_NVARS+AR_ID),phi_arr(i,j,k,FLO_NVARS+H2_ID));
                amrex::Real J=sigma*phi_arr(i,j,k,EFLDX_INDX+dircn);
                //negative sign already accounted for in the diff term
                rhs_arr(i,j,k)+=MU0*J; 

            });
        }

        // average cell coefficients to faces, this includes boundary faces
        Array<MultiFab, AMREX_SPACEDIM> face_bcoeff;
        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim)
        {
            const BoxArray& ba = amrex::convert(acoeff[ilev].boxArray(), 
                                                IntVect::TheDimensionVector(idim));
            face_bcoeff[idim].define(ba, acoeff[ilev].DistributionMap(), 1, 0);
            face_bcoeff[idim].setVal(1.0);
        }
        // set boundary conditions
        for (MFIter mfi(phi_new[ilev], TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const Box& bx = mfi.tilebox();
            const auto dx = geom[ilev].CellSizeArray();
            auto prob_lo = geom[ilev].ProbLoArray();
            auto prob_hi = geom[ilev].ProbHiArray();
            const Box& domain = geom[ilev].Domain();

            Array4<Real> phi_arr = Sborder.array(mfi);

            Array4<Real> robin_a_arr = robin_a[ilev].array(mfi);
            Array4<Real> robin_b_arr = robin_b[ilev].array(mfi);
            Array4<Real> robin_f_arr = robin_f[ilev].array(mfi);

            Real time = current_time; // for GPU capture

            for (int idim = 0; idim < AMREX_SPACEDIM; ++idim)
            {
                //note: bdryLo/bdryHi grabs the face indices from bx that are the boundary
                //since they are face indices, the bdry normal index is 0/n+1, n is number of cells
                //so the ghost cell index at left side is i-1 while it is i on the right
                if (bx.smallEnd(idim) == domain.smallEnd(idim))
                {
                    amrex::ParallelFor(amrex::bdryLo(bx, idim), [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {

                        mflo_user_funcs::vec_potential_bc(i, j, k, dircn, idim, -1, 
                                                      phi_arr, robin_a_arr, 
                                                      robin_b_arr, robin_f_arr, 
                                                      prob_lo, prob_hi, dx, time);
                    });
                }
                if (bx.bigEnd(idim) == domain.bigEnd(idim))
                {
                    amrex::ParallelFor(amrex::bdryHi(bx, idim), [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                        mflo_user_funcs::vec_potential_bc(i, j, k, dircn, idim, +1, 
                                                      phi_arr, robin_a_arr, 
                                                      robin_b_arr, robin_f_arr, 
                                                      prob_lo, prob_hi, dx, time);
                    });
                }
            }
        }

        linsolve_obj.setACoeffs(ilev, acoeff[ilev]);

        // set b with diffusivities
        linsolve_obj.setBCoeffs(ilev, amrex::GetArrOfConstPtrs(face_bcoeff));

        // bc's are stored in the ghost cells of potential
        if(mixedbc)
        {
            linsolve_obj.setLevelBC(ilev, &potential[ilev], &(robin_a[ilev]), 
                                     &(robin_b[ilev]), &(robin_f[ilev]));
        }
        else
        {
            linsolve_obj.setLevelBC(ilev, &potential[ilev]);
        }

    }

    MLMG mlmg(linsolve_obj);
    mlmg.setMaxIter(max_iter);
    mlmg.setVerbose(linsolve_verbose);

#ifdef AMREX_USE_HYPRE
    if (use_hypre)
    {
        mlmg.setHypreOptionsNamespace("mflo.hypre");
        mlmg.setBottomSolver(MLMG::BottomSolver::hypre);
    }
#endif

    mlmg.solve(GetVecOfPtrs(solution), GetVecOfConstPtrs(rhs), tol_rel, tol_abs);
    amrex::Print()<<"Solved Magnetic vector potential along "<<Adir<<"\n";
    
    Vector<Array<MultiFab, AMREX_SPACEDIM>> gradsoln(finest_level+1);
    for (int ilev = 0; ilev <= finest_level; ilev++)
    {
        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim)
        {
            const BoxArray& faceba = amrex::convert(grids[ilev], 
                    IntVect::TheDimensionVector(idim));
            gradsoln[ilev][idim].define(faceba, dmap[ilev], 1, 0);
        }
    }
    mlmg.getGradSolution(GetVecOfArrOfPtrs(gradsoln));

    for (int ilev = 0; ilev <= finest_level; ilev++)
    {
        amrex::MultiFab::Copy(phi_new[ilev], solution[ilev], 0, AX_INDX+Adir, 1, 0);
        
        const Array<const MultiFab*, AMREX_SPACEDIM> allgrad = {&gradsoln[ilev][0], 
            &gradsoln[ilev][1], &gradsoln[ilev][2]};
        average_face_to_cellcenter(phi_new[ilev], GAXX_INDX+3*Adir, allgrad);
    }


    //clean-up
    potential.clear();
    acoeff.clear();
    solution.clear();
    rhs.clear();
    robin_a.clear();
    robin_b.clear();
    robin_f.clear();
}

void mflo::update_Bfields()
{
    for (int ilev = 0; ilev <= finest_level; ilev++)
    {
        for (MFIter mfi(phi_new[ilev], TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const Box& bx = mfi.tilebox();
            Array4<Real> phi_arr = phi_new[ilev].array(mfi);

            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                
                //Bx=dAz/dy-dAy/dz
                phi_arr(i,j,k,BX_INDX)=phi_arr(i,j,k,GAZY_INDX)-phi_arr(i,j,k,GAYZ_INDX);
                //By=dAx/dz-dAz/dx
                phi_arr(i,j,k,BY_INDX)=phi_arr(i,j,k,GAXZ_INDX)-phi_arr(i,j,k,GAZX_INDX);
                //Bz=dAy/dx-dAx/dy
                phi_arr(i,j,k,BZ_INDX)=phi_arr(i,j,k,GAYX_INDX)-phi_arr(i,j,k,GAXY_INDX);
            });
        }
    }
}

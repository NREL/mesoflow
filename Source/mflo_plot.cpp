#include <AMReX_ParallelDescriptor.H>
#include <AMReX_ParmParse.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_FillPatchUtil.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_VisMF.H>
#include <AMReX_PhysBCFunct.H>
#include <kernels_3d.H>

#ifdef AMREX_MEM_PROFILING
#include <AMReX_MemProfiler.H>
#endif

#include <mflo.H>
#include <userfuncs.H>

// utility to skip to next line in Header
void mflo::GotoNextLine(std::istream& is)
{
    constexpr std::streamsize bl_ignore_max{100000};
    is.ignore(bl_ignore_max, '\n');
}

// get plotfile name
std::string mflo::PlotFileName(int lev) const
{
    return amrex::Concatenate(plot_file, lev, 5);
}

// put together an array of multifabs for writing
Vector<const MultiFab*> mflo::PlotFileMF() const
{
    Vector<const MultiFab*> r;
    for (int i = 0; i <= finest_level; ++i) {
        r.push_back(&phi_new[i]);
    }
    return r;
}

// set plotfile variable names
Vector<std::string> mflo::PlotFileVarNames() const { return mflo_varnames; }

// write plotfile to disk
void mflo::WritePlotFile() const
{
    const std::string& plotfilename = PlotFileName(istep[0]);
    const auto& mf = PlotFileMF();
    const auto& varnames = PlotFileVarNames();

    amrex::Print() << "Writing plotfile " << plotfilename << "\n";

    amrex::WriteMultiLevelPlotfile(
        plotfilename, finest_level + 1, mf, varnames, Geom(), t_new[0], istep,
        refRatio());
}

void mflo::WriteCheckpointFile() const
{
    const auto& varnames = PlotFileVarNames();

    // chk00010            write a checkpoint file with this root directory
    // chk00010/Header     this contains information you need to save (e.g.,
    // finest_level, t_new, etc.) and also
    //                     the BoxArrays at each level
    // chk00010/Level_0/
    // chk00010/Level_1/
    // etc.                these subdirectories will hold the MultiFab data at
    // each level of refinement

    // checkpoint file name, e.g., chk00010
    const std::string& checkpointname = amrex::Concatenate(chk_file, istep[0]);

    amrex::Print() << "Writing checkpoint " << checkpointname << "\n";

    const int nlevels = finest_level + 1;

    // ---- prebuild a hierarchy of directories
    // ---- dirName is built first.  if dirName exists, it is renamed.  then
    // build
    // ---- dirName/subDirPrefix_0 .. dirName/subDirPrefix_nlevels-1
    // ---- if callBarrier is true, call ParallelDescriptor::Barrier()
    // ---- after all directories are built
    // ---- ParallelDescriptor::IOProcessor() creates the directories
    amrex::PreBuildDirectorHierarchy(checkpointname, "Level_", nlevels, true);

    // write Header file
    if (ParallelDescriptor::IOProcessor()) 
    {

        std::string HeaderFileName(checkpointname + "/Header");
        VisMF::IO_Buffer io_buffer(VisMF::IO_Buffer_Size);
        std::ofstream HeaderFile;
        HeaderFile.rdbuf()->pubsetbuf(io_buffer.dataPtr(), io_buffer.size());
        HeaderFile.open(HeaderFileName.c_str(),
            std::ofstream::out | std::ofstream::trunc | std::ofstream::binary);

        if (!HeaderFile.good()) 
        {
            amrex::FileOpenFailed(HeaderFileName);
        }

        HeaderFile.precision(17);

        // write out title line
        HeaderFile << "Checkpoint file for mflo\n";

        // write out finest_level
        HeaderFile << finest_level << "\n";

		//write current
		HeaderFile << current << "\n";
		
        // write out array of istep
        for (int i = 0; i < istep.size(); ++i) 
        {
            HeaderFile << istep[i] << " ";
        }
        HeaderFile << "\n";

        // write out array of dt
        for (int i = 0; i < dt.size(); ++i) 
        {
            HeaderFile << dt[i] << " ";
        }
        HeaderFile << "\n";

        // write out array of t_new
        for (int i = 0; i < t_new.size(); ++i) 
        {
            HeaderFile << t_new[i] << " ";
        }
        HeaderFile << "\n";

        // write the BoxArray at each level
        for (int lev = 0; lev <= finest_level; ++lev) 
        {
            boxArray(lev).writeOn(HeaderFile);
            HeaderFile << '\n';
        }
    }

    // write the MultiFab data to, e.g., chk00010/Level_0/
    for (int lev = 0; lev <= finest_level; ++lev) 
    {
        VisMF::Write(phi_new[lev], amrex::MultiFabFileFullPrefix(
                           lev, checkpointname, "Level_", varnames[0]));
    }
}

void mflo::ReadCheckpointFile()
{
    const auto& varnames = PlotFileVarNames();

    amrex::Print() << "Restart from checkpoint " << restart_chkfile << "\n";

    // Header
    std::string File(restart_chkfile + "/Header");

    VisMF::IO_Buffer io_buffer(VisMF::GetIOBufferSize());

    Vector<char> fileCharPtr;
    ParallelDescriptor::ReadAndBcastFile(File, fileCharPtr);
    std::string fileCharPtrString(fileCharPtr.dataPtr());
    std::istringstream is(fileCharPtrString, std::istringstream::in);

    std::string line, word;

    // read in title line
    std::getline(is, line);

    // read in finest_level
    is >> finest_level;
    GotoNextLine(is);

	//read in current
	is >> current;
	GotoNextLine(is);
	
    // read in array of istep
    std::getline(is, line);
    {
        std::istringstream lis(line);
        int i = 0;
        while (lis >> word) 
        {
            istep[i++] = std::stoi(word);
        }
    }

    // read in array of dt
    std::getline(is, line);
    {
        std::istringstream lis(line);
        int i = 0;
        while (lis >> word) 
        {
            dt[i++] = std::stod(word);
        }
    }

    // read in array of t_new
    std::getline(is, line);
    {
        std::istringstream lis(line);
        int i = 0;
        while (lis >> word) 
        {
            t_new[i++] = std::stod(word);
        }
    }

    for (int lev = 0; lev <= finest_level; ++lev) 
    {

        // read in level 'lev' BoxArray from Header
        BoxArray ba;
        ba.readFrom(is);
        GotoNextLine(is);

        // create a distribution mapping
        DistributionMapping dm{ba, ParallelDescriptor::NProcs()};

        // set BoxArray grids and DistributionMapping dmap in AMReX_AmrMesh.H
        // class
        SetBoxArray(lev, ba);
        SetDistributionMap(lev, dm);

        // build MultiFab and FluxRegister data
        int ncomp = varnames.size();
        int nghost = 0;
        phi_old[lev].define(grids[lev], dmap[lev], ncomp, nghost);
        phi_new[lev].define(grids[lev], dmap[lev], ncomp, nghost);
        if (lev > 0 && do_reflux) 
        {
            flux_reg[lev].reset(new FluxRegister(
                grids[lev], dmap[lev], refRatio(lev - 1), lev, ncomp));
        }
    }

    // read in the MultiFab data
    for (int lev = 0; lev <= finest_level; ++lev) 
    {
        VisMF::Read(phi_new[lev], amrex::MultiFabFileFullPrefix(
                    lev, restart_chkfile, "Level_", varnames[0]));
    }
}

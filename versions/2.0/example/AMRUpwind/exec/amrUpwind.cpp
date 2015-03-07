#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

#include <iostream>
#include <fstream> 
using std::ifstream;
using std::ios;

#include "FABView.H"

#include "ParmParse.H"
#include "CH_HDF5.H"
#include "parstream.H"

#include "AMR.H"
#include "AMRLevel.H"
#include "AMRLevelUpwindFactory.H"
#include "AMRLevelUpwind.H"

#ifdef CH_MPI
#include "CH_Attach.H"
#endif

// amrGodunov is a function (as opposed to inline in main()) to get
// around MPI scoping problems
void amrUpwind();

// setupFixedGrids allows fixed grids to be read in and used in this AMR
// computation example
void setupFixedGrids(Vector<Vector<Box> >& a_amrGrids,
                     const ProblemDomain&  a_domain,
                     int                   a_maxLevel,
                     int                   a_maxGridSize,
                     int                   a_blockFactor,
                     int                   a_verbosity,
                     std::string           a_gridFile);

// One more function for MPI
void dumpmemoryatexit();

int main(int a_argc, char* a_argv[])
{
#ifdef CH_MPI
  // Start MPI
  MPI_Init(&a_argc,&a_argv);
  setChomboMPIErrorHandler();
#endif

  // Check for an input file
  char* inFile = NULL;

  if (a_argc > 1)
  {
    inFile = a_argv[1];
  }
  else
  {
    pout() << "Usage:  amrUpwind...ex <inputfile>" << endl;
    pout() << "No input file specified" << endl;
    return -1;
  }

  // Parse the command line and the input file (if any)
  ParmParse pp(a_argc-2,a_argv+2,NULL,inFile);

  // Run amrUpwind, i.e., do the computation
  amrUpwind();

#ifdef CH_MPI
  // Exit MPI
  dumpmemoryatexit();
  MPI_Finalize();
#endif
}

void amrUpwind()
{
  // Read inputs that are prefixed with "main."
  ParmParse ppMain("main");

  // This determines the amount of diagnositic output generated
  int verbosity = 0;
  ppMain.query("verbosity",verbosity);
  CH_assert(verbosity >= 0);

  // Stop after this number of steps
  int nstop = 0;
  ppMain.get("max_step",nstop);

  // Stop when the simulation time get here
  Real stopTime = 0.0;
  ppMain.get("max_time",stopTime);

  // Set the physical size of the longest dimension of the domain
  Real domainLength = 1.0;
  ppMain.get("domain_length",domainLength);

  // Set the resolution of the coarsest level
  vector<int> numCells(SpaceDim);
  for (int i = 0; i < SpaceDim; ++i)
  {
    numCells[i] = 0;
  }
  ppMain.getarr("num_cells",numCells,0,SpaceDim);

  CH_assert(D_TERM(   (numCells[0] > 0),
                && (numCells[1] > 0),
                && (numCells[2] > 0)));
  CH_assert(D_TERM(   (numCells[0] % 2 == 0),
                && (numCells[1] % 2 == 0),
                && (numCells[2] % 2 == 0)));

  // Determine which spatial directions are periodic
  vector<int> isPeriodica(SpaceDim,1);
  bool isPeriodic[SpaceDim];

  ppMain.queryarr("is_periodic",isPeriodica,0,SpaceDim);
  // convert periodic from int->bool
  for (int dim = 0; dim < SpaceDim; dim++)
  {
    isPeriodic[dim] = (isPeriodica[dim] == 1);
    if (isPeriodic[dim] && verbosity >= 2 && procID() == 0)
    {
      pout() << "Using Periodic BCs in direction: " << dim << endl;
    }
  }

  // Maximum AMR level limit
  int maxLevel = 0;
  ppMain.get("max_level",maxLevel);
  int numReadLevels = Max(maxLevel,1);

  // Refinement ratios between levels
  std::vector<int> refRatios;
  // Note: this requires a refRatio to be defined for the finest level
  // (even though it will never be used)
  ppMain.getarr("ref_ratio",refRatios,0,numReadLevels+1);

  // Number of coarse time steps from one regridding to the next
  std::vector<int> regridIntervals;
  ppMain.getarr("regrid_interval",regridIntervals,0,numReadLevels);

  // How far to extend refinement from cells newly tagged for refinement
  int tagBufferSize = 3;
  ppMain.get("tag_buffer_size",tagBufferSize);

  // Threshold that triggers refinement
  Real refineThresh = 0.3;
  ppMain.get ("refine_thresh",refineThresh);

  // Minimum dimension of a grid
  int blockFactor = 1;
  ppMain.get("block_factor",blockFactor);

  // Maximum dimension of a grid
  int maxGridSize = 32;
  ppMain.get("max_grid_size",maxGridSize);

  Real fillRatio = 0.75;
  ppMain.get("fill_ratio",fillRatio);

  // Set up checkpointing
  int checkpointInterval = 0;
  ppMain.query("checkpoint_interval",checkpointInterval);

  // Set up plot file writing
  int plotInterval = 0;
  ppMain.query("plot_interval",plotInterval);

  // CFL number
  Real cfl = 0.8;
  ppMain.get("cfl",cfl);

  // Initial CFL multiplier
  Real initialDtMult = 0.1;
  ppMain.query("initial_dt_mult",initialDtMult);

  // Refinement ratios between levels
  std::vector<Real> advectVelVect;
  // Note: this requires a refRatio to be defined for the finest level
  // (even though it will never be used)
  ppMain.getarr("advection_velocity",advectVelVect,0,SpaceDim);
  RealVect advectionVel;
  for (int dir=0; dir<SpaceDim; dir++)
    {
      advectionVel[dir] = advectVelVect[dir];
    }


  // Determine if a fixed or variable time step will be used
  Real fixedDt = -1;
  ppMain.query("fixed_dt",fixedDt);

  // Limit the time step growth
  Real maxDtGrowth = 1.1;
  ppMain.get("max_dt_growth",maxDtGrowth);

  // Let the time step grow by this factor above the "maximum" before
  // reducing it
  Real dtToleranceFactor = 1.1;
  ppMain.get("dt_tolerance_factor",dtToleranceFactor);


  if (verbosity >= 2)
  {
    pout() << "maximum step = " << nstop << endl;
    pout() << "maximum time = " << stopTime << endl;

    pout() << "number of cells = " << D_TERM(numCells[0] << "  " <<,
                                             numCells[1] << "  " <<,
                                             numCells[2] << ) endl;

    pout() << "maximum level = " << maxLevel << endl;

    pout() << "refinement ratio = ";
    for (int i = 0; i < refRatios.size(); ++i)
    {
      pout() << refRatios[i] << " ";
    }
    pout() << endl;

    pout() << "regrid interval = ";
    for (int i = 0; i < regridIntervals.size(); ++i)
    {
      pout() << regridIntervals[i] << " ";
    }
    pout() << endl;

    pout() << "refinement threshold = " << refineThresh << endl;

    pout() << "blocking factor = " << blockFactor << endl;
    pout() << "max grid size = " << maxGridSize << endl;
    pout() << "fill ratio = " << fillRatio << endl;

    pout() << "checkpoint interval = " << checkpointInterval << endl;
    pout() << "plot interval = " << plotInterval << endl;
    pout() << "CFL = " << cfl << endl;
    pout() << "initial dt multiplier = " << initialDtMult << endl;
    if (fixedDt > 0)
    {
      pout() << "fixed dt = " << fixedDt << endl;
    }
    pout() << "maximum dt growth = " << maxDtGrowth << endl;
    pout() << "dt tolerance factor = " << dtToleranceFactor << endl;
  }

  ProblemDomain probDomain (IntVect::Zero,
                            IntVect(D_DECL(numCells[0]-1,
                                           numCells[1]-1,
                                           numCells[2]-1)),
                            isPeriodic);


  // Set up the AMRLevelUpwind factory
  AMRLevelUpwindFactory AMRLevelFact;
  AMRLevelFact.CFL(cfl);
  AMRLevelFact.domainLength(domainLength);
  AMRLevelFact.refinementThreshold(refineThresh);
  AMRLevelFact.tagBufferSize(tagBufferSize);
  AMRLevelFact.verbosity(verbosity);
  AMRLevelFact.initialDtMultiplier(initialDtMult);
  AMRLevelFact.advectionVel(advectionVel);

  AMR amr;

  // Set up the AMR object
  amr.define(maxLevel,refRatios,probDomain,&AMRLevelFact);

  if (fixedDt > 0)
  {
    amr.fixedDt(fixedDt);
  }

  // Set grid generation parameters
  amr.maxGridSize(maxGridSize);
  amr.blockFactor(blockFactor);
  amr.fillRatio(fillRatio);

  // The hyperbolic codes use a grid buffer of 1
  amr.gridBufferSize(1);

  // Set output parameters
  amr.checkpointInterval(checkpointInterval);
  amr.plotInterval(plotInterval);
  amr.regridIntervals(regridIntervals);
  amr.maxDtGrow(maxDtGrowth);
  amr.dtToleranceFactor(dtToleranceFactor);

  // Set up output files
  if (ppMain.contains("plot_prefix"))
  {
    std::string prefix;
    ppMain.query("plot_prefix",prefix);
    amr.plotPrefix(prefix);
  }

  if (ppMain.contains("chk_prefix"))
  {
    std::string prefix;
    ppMain.query("chk_prefix",prefix);
    amr.checkpointPrefix(prefix);
  }

  amr.verbosity(verbosity);

  // Set up input files
  if (!ppMain.contains("restart_file"))
  {
    if (!ppMain.contains("fixed_hierarchy"))
    {
      // initialize from scratch for AMR run
      // initialize hierarchy of levels
      amr.setupForNewAMRRun();
    }
    else
    {
      std::string gridFile;
      ppMain.query("fixed_hierarchy",gridFile);

      // initialize from a list of grids in "gridFile"
      Vector<Vector<Box> > amrGrids(maxLevel+1);
      setupFixedGrids(amrGrids,
                      probDomain,
                      maxLevel,
                      maxGridSize,
                      blockFactor,
                      verbosity,
                      gridFile);
      amr.setupForFixedHierarchyRun(amrGrids,1);
    }
  }
  else
  {
    std::string restartFile;
    ppMain.query("restart_file",restartFile);

#ifdef CH_USE_HDF5
    HDF5Handle handle(restartFile,HDF5Handle::OPEN_RDONLY);
    // read from checkpoint file
    amr.setupForRestart(handle);
    handle.close();
#else
    MayDay::Error("amrGodunov restart only defined with hdf5");
#endif
  }

  // Run the computation
  amr.run(stopTime,nstop);

  // Output the last plot file and statistics
  amr.conclude();
}

// setupFixedGrids allows fixed grids to be read in and used in this AMR
// computation example
void setupFixedGrids(Vector<Vector<Box> >& a_amrGrids,
                     const ProblemDomain&  a_domain,
                     int                   a_maxLevel,
                     int                   a_maxGridSize,
                     int                   a_blockFactor,
                     int                   a_verbosity,
                     std::string           a_gridFile)
{

  // Run this task on one processor
  if (procID() == uniqueProc(SerialTask::compute))
    {
      a_amrGrids.push_back(Vector<Box>(1,a_domain.domainBox()));
      
      // Read in predefined grids
      ifstream is(a_gridFile.c_str(), ios::in);
      
      if (is.fail())
        {
          MayDay::Error("Cannot open grids file");
        }
      
      // Format of file:
      //   number of levels, then for each level (starting with level 1):
      //   number of grids on level, list of boxes
      
      int inNumLevels;
      is >> inNumLevels;
      
      CH_assert (inNumLevels <= a_maxLevel+1);
      
      if (a_verbosity >= 3)
        {
          pout() << "numLevels = " << inNumLevels << endl;
        }
      
      while (is.get() != '\n');
      
      a_amrGrids.resize(inNumLevels);
      
      // Check to see if coarsest level needs to be broken up
      domainSplit(a_domain,a_amrGrids[0],a_maxGridSize,a_blockFactor);
      
      if (a_verbosity >= 3)
        {
          pout() << "level 0: ";
          for (int n = 0; n < a_amrGrids[0].size(); n++)
            {
              pout() << a_amrGrids[0][0] << endl;
            }
        }
      
      // Now loop over levels, starting with level 1
      int ngrid;
      for (int lev = 1; lev < inNumLevels; lev++)
        {
          is >> ngrid;
          
          if (a_verbosity >= 3)
            {
              pout() << "level " << lev << " numGrids = " << ngrid << endl;
              pout() << "Grids: ";
            }
          
          while (is.get() != '\n');
          
          a_amrGrids[lev].resize(ngrid);
          
          for (int i = 0; i < ngrid; i++)
            {
              Box bx;
              is >> bx;
              
              while (is.get() != '\n');
              
              // Quick check on box size
              Box bxRef(bx);
              
              if (bxRef.longside() > a_maxGridSize)
                {
                  pout() << "Grid " << bx << " too large" << endl;
                  MayDay::Error();
                }
              
              if (a_verbosity >= 3)
                {
                  pout() << bx << endl;
                }
              
              a_amrGrids[lev][i] = bx;
            } // End loop over boxes on this level
        } // End loop over levels
    }
  
  // Broadcast results to all the processors
  broadcast(a_amrGrids,uniqueProc(SerialTask::compute));

}





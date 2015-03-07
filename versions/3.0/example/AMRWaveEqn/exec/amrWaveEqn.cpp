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
using std::ifstream;
using std::ios;
#include <string>

#include "FABView.H"

#include "ParmParse.H"
#include "CH_HDF5.H"
#include "parstream.H"

#include "CH_Timer.H"

#include "AMR.H"
#include "AMRLevel.H"
#include "AMRLevelWaveEqnFactory.H"

#include "WaveIBC.H"

#ifdef CH_MPI
#include "CH_Attach.H"
#endif

#include "UsingNamespace.H"

// #define TRAP_FPE
#ifdef TRAP_FPE
extern "C" {
#include <fpu_control.h>
}
/* IM: Invalid operation mask
 * DM: Denormalized operand mask
 * ZM: Zero-divide mask
 * OM: Overflow mask
 * UM: Underflow mask
 * PM: Precision (inexact result) mask
  ---(pm is kinda stupid)
*/
static void __attribute__ ((constructor)) trapfpe(void)
{
  fpu_control_t cw =
    _FPU_DEFAULT & ~(_FPU_MASK_IM | _FPU_MASK_ZM | _FPU_MASK_OM | _FPU_MASK_UM);
  _FPU_SETCW(cw);
}
#endif

void amrWaveEqn();

// amrWaveEqn is a function (as opposed to inline in main()) to get
// around MPI scoping problems
// setupFixedGrids allows fixed grids to be read in and used in this AMR
// computation example
void setupFixedGrids(Vector<Vector<Box> >& a_amrGrids,
                     const ProblemDomain&  a_domain,
                     int                   a_maxLevel,
                     int                   a_maxGridSize,
                     int                   a_blockFactor,
                     int                   a_verbosity,
                     std::string           a_gridFile);

#ifndef CH_NTIMER
OldTimer *all_timer;
OldTimer *setup_timer;
OldTimer *solve_timer;
OldTimer *timestep_timer;
OldTimer *shutdown_timer;
#endif

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

int main(int a_argc, char* a_argv[])
{
#ifdef CH_MPI
  // Start MPI
  MPI_Init(&a_argc,&a_argv);
  // setChomboMPIErrorHandler();
#endif

#ifndef CH_NTIMER
  // timers
  all_timer      = new OldTimer("All",0) ;
  setup_timer    = new OldTimer("Setup",*all_timer) ;
  solve_timer    = new OldTimer("Solve",*all_timer,1) ;
  timestep_timer = new OldTimer("TimeStep",*solve_timer) ;
  shutdown_timer = new OldTimer("Shutdown",*all_timer) ;

  OldTimer::TimerInit(0);
  all_timer->start();
  setup_timer->start();
#endif

  // Check for an input file
  char* inFile = NULL;

#ifdef TRAP_FPE
  trapfpe();
#endif

  if (a_argc > 1)
  {
    inFile = a_argv[1];
  }
  else
  {
    pout() << "Usage:  amrWaveEqn...ex <inputfile>" << endl;
    pout() << "No input file specified" << endl;
    return -1;
  }

  // Parse the command line and the input file (if any)
  ParmParse pp(a_argc-2,a_argv+2,NULL,inFile);

  // Run amrWaveEqn, i.e., do the computation
  amrWaveEqn();

#ifdef CH_MPI
  // Exit MPI
  dumpmemoryatexit();
#endif

#ifndef CH_NTIMER
  shutdown_timer->stop() ;
  all_timer->stop() ;
  OldTimer::TimerSummary();
#endif

#ifdef CH_MPI
  MPI_Finalize();
#endif
}

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

void amrWaveEqn()
{

  // Read inputs that are prefixed with "wave."
  ParmParse ppwave("wave");

  // Determine the sample problem specified
  //XXX -- not used
  //XXXint problem = -1;
  std::string problemString;


  // This determines the amount of diagnositic output generated
  int verbosity = 0;
  ppwave.query("verbosity",verbosity);
 CH_assert(verbosity >= 0);

  // Parameters specific to different sample problems

  // Stop after this number of steps
  int nstop = 0;
  ppwave.query("max_step",nstop);

  // Stop when the simulation time get here
  Real stopTime = 0.0;
  ppwave.query("max_time",stopTime);

  // Set the physical size of the longest dimension of the domain
  Real domainLength = 1.0;
  ppwave.query("domain_length",domainLength);

  // Set the location of the lower left corner
  std::vector<Real> x0a(SpaceDim,0.0);
  ppwave.queryarr("x0",x0a,0,SpaceDim);
  RealVect x0;
  for( int d=0 ; d<SpaceDim ; ++d ) x0[d] = x0a[d] ;

  // Set the resolution of the coarsest level
  vector<int> numCells(SpaceDim);
  for (int i = 0; i < SpaceDim; ++i) numCells[i]=0;
  ppwave.queryarr("num_cells",numCells,0,SpaceDim);

 CH_assert(D_TERM(   (numCells[0] > 0),
                && (numCells[1] > 0),
                && (numCells[2] > 0)));
 CH_assert(D_TERM(   (numCells[0] % 2 == 0),
                && (numCells[1] % 2 == 0),
                && (numCells[2] % 2 == 0)));

  // Determine which spatial directions are periodic
  vector<int> isPeriodica(SpaceDim,0);
  bool isPeriodic[SpaceDim];

  ppwave.queryarr("is_periodic",isPeriodica,0,SpaceDim);
  // convert periodic from int->bool
  for (int dim=0; dim<SpaceDim; dim++)
    {
      isPeriodic[dim] = (isPeriodica[dim] == 1);
      if (isPeriodic[dim] && verbosity >= 2 && procID() == 0)
        pout() << "Using Periodic BCs in direction: " << dim << endl;
    }

  // Maximum AMR level limit
  int maxLevel = 0;
  ppwave.query("max_level",maxLevel);
  int numReadLevels = Max(maxLevel,1);

  // Refinement ratios between levels
  std::vector<int> refRatios;
  // Note: this requires a refRatio to be defined for the finest level
  // (even though it will never be used)
  ppwave.queryarr("ref_ratio",refRatios,0,numReadLevels+1);

  // Number of coarse time steps from one regridding to the next
  std::vector<int> regridIntervals;
  ppwave.queryarr("regrid_interval",regridIntervals,0,numReadLevels);

  // How far to extend refinement from cells newly tagged for refinement
  int tagBufferSize = 3;
  ppwave.query("tag_buffer_size",tagBufferSize);

  // Threshold that triggers refinement
  Real refineThresh = 0.3;
  ppwave.query ("refine_thresh",refineThresh);

  // Minimum dimension of a grid
  int blockFactor = 1;
  ppwave.query("block_factor",blockFactor);

  // Maximum dimension of a grid
  int maxGridSize = 32;
  ppwave.query("max_grid_size",maxGridSize);
  int maxBaseGridSize = 0;
  ppwave.query("max_base_grid_size",maxBaseGridSize);

  Real fillRatio = 0.75;
  ppwave.query("fill_ratio",fillRatio);

  // Set up checkpointing
  int checkpointInterval = 0;
  ppwave.query("checkpoint_interval",checkpointInterval);

  // Set up plot file writing
  int plotInterval = 0;
  ppwave.query("plot_interval",plotInterval);

  // CFL multiplier
  Real cfl = 0.8;
  ppwave.query("cfl",cfl);

  // Initial CFL multiplier
  Real initialCFL = 0.1;
  ppwave.query("initial_cfl",initialCFL);

  // Determine if a fixed or variable time step will be used
  Real fixedDt = -1;
  ppwave.query("fixed_dt",fixedDt);

  // Limit the time step growth
  Real maxDtGrowth = 1.1;
  ppwave.query("max_dt_growth",maxDtGrowth);

  // Let the time step grow by this factor above the "maximum" before
  // reducing it
  Real dtToleranceFactor = 1.1;
  ppwave.query("dt_tolerance_factor",dtToleranceFactor);

  // Print the parameters

  pout() << "maximum step = " << nstop << endl;
  pout() << "maximum time = " << stopTime << endl;

  pout() << "number of cells = " << D_TERM(numCells[0] << "  " <<,
                                               numCells[1] << "  " <<,
                                               numCells[2] << ) endl;

  pout() << "maximum level = " << maxLevel << endl;

  pout() << "refinement ratio = ";
  for (int i = 0; i < refRatios.size(); ++i) pout() << refRatios[i] << " ";
  pout() << endl;

  pout() << "regrid interval = ";
  for (int i = 0; i < regridIntervals.size(); ++i) pout() << regridIntervals[i] << " ";
  pout() << endl;

  pout() << "refinement threshold = " << refineThresh << endl;

  pout() << "blocking factor = " << blockFactor << endl;
  pout() << "max grid size = " << maxGridSize << endl;
  pout() << "max base grid size = " << maxBaseGridSize << endl ;
  pout() << "fill ratio = " << fillRatio << endl;

  pout() << "checkpoint interval = " << checkpointInterval << endl;
  pout() << "plot interval = " << plotInterval << endl;
  pout() << "CFL = " << cfl << endl;
  pout() << "initial CFL = " << initialCFL << endl;
  if (fixedDt > 0)
    {
      pout() << "fixed dt = " << fixedDt << endl;
    }
  pout() << "maximum dt growth = " << maxDtGrowth << endl;
  pout() << "dt tolerance factor = " << dtToleranceFactor << endl;


  // Create initial and boundary condition (IBC) object and initialize

  Real r0;
  ppwave.query("radius",r0);
  WaveIBC* waveibc = new WaveIBC;
  waveibc->setParams(r0,x0);
  PhysIBC* ibc  = waveibc;

  ProblemDomain probDomain (IntVect::Zero,
                            IntVect(D_DECL(numCells[0]-1,
                                           numCells[1]-1,
                                           numCells[2]-1)),
                            isPeriodic);

  // Set up the AMRLevel... factory
  //
  AMRLevelWaveEqnFactory amrWaveEqnFact;
  amrWaveEqnFact.CFL(cfl);
  amrWaveEqnFact.domainLength(domainLength);
  amrWaveEqnFact.refinementThreshold(refineThresh);
  amrWaveEqnFact.tagBufferSize(tagBufferSize);
  amrWaveEqnFact.verbosity(verbosity);
  amrWaveEqnFact.initialDtMultiplier(initialCFL);
  amrWaveEqnFact.IBC(ibc);
  AMR amr;

  // Set up the AMR object
  amr.define(maxLevel,refRatios,probDomain,&amrWaveEqnFact);

  if (fixedDt > 0)
    {
      amr.fixedDt(fixedDt);
    }

  // Set grid generation parameters
  amr.maxGridSize(maxGridSize);
  if( maxBaseGridSize != 0){ amr.maxBaseGridSize(maxBaseGridSize); }//defaults to maxGridSize
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
  if (ppwave.contains("plot_prefix"))
    {
      std::string prefix;
      ppwave.query("plot_prefix",prefix);
      amr.plotPrefix(prefix);
    }

  if (ppwave.contains("chk_prefix"))
    {
      std::string prefix;
      ppwave.query("chk_prefix",prefix);
      amr.checkpointPrefix(prefix);
    }

  amr.verbosity(verbosity);

  // Set up input files
  if (!ppwave.contains("restart_file"))
    {
      if (!ppwave.contains("fixed_hierarchy"))
        {
          // initialize from scratch for AMR run
          // initialize hierarchy of levels
          amr.setupForNewAMRRun();
        }
      else
        {
          std::string gridFile;
          ppwave.query("fixed_hierarchy",gridFile);

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
      ppwave.query("restart_file",restartFile);

#ifdef CH_USE_HDF5
      HDF5Handle handle(restartFile,HDF5Handle::OPEN_RDONLY);
      // read from checkpoint file
      amr.setupForRestart(handle);
      handle.close();
#else
      MayDay::Error("amrWaveEqn restart only defined with hdf5");
#endif
    }

#ifndef CH_NTIMER
  setup_timer->stop();
#endif

  // Run the computation
#ifndef CH_NTIMER
  solve_timer->start();
#endif

  amr.run(stopTime,nstop);

#ifndef CH_NTIMER
  solve_timer->stop();
#endif

  // Output the last plot file and statistics
#ifndef CH_NTIMER
  shutdown_timer->start() ;
#endif

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
//  // Run this task on one processor
//  if (procID() == uniqueProc(SerialTask::compute))
//  {
//    a_amrGrids.push_back(Vector<Box>(1,a_domain.domainBox()));
//
//    // Read in predefined grids
//    ifstream is(a_gridFile.c_str(), ios::in);
//
//    if (is.fail())
//    {
//      MayDay::Error("Cannot open grids file");
//    }
//
//    // Format of file:
//    //   number of levels, then for each level (starting with level 1):
//    //   number of grids on level, list of boxes
//
//    int inNumLevels;
//    is >> inNumLevels;
//
//   CH_assert (inNumLevels <= a_maxLevel+1);
//
//    if (a_verbosity >= 3)
//    {
//      pout() << "numLevels = " << inNumLevels << endl;
//    }
//
//    while (is.get() != '\n');
//
//    a_amrGrids.resize(inNumLevels);
//
//    // Check to see if coarsest level needs to be broken up
//    domainSplit(a_domain,a_amrGrids[0],a_maxGridSize,a_blockFactor);
//
//    if (a_verbosity >= 3)
//    {
//      pout() << "level 0: ";
//      for (int n = 0; n < a_amrGrids[0].size(); n++)
//      {
//        pout() << a_amrGrids[0][0] << endl;
//      }
//    }
//
//    // Now loop over levels, starting with level 1
//    int ngrid;
//    for (int lev = 1; lev < inNumLevels; lev++)
//    {
//      is >> ngrid;
//
//      if (a_verbosity >= 3)
//      {
//        pout() << "level " << lev << " numGrids = " << ngrid << endl;
//        pout() << "Grids: ";
//      }
//
//      while (is.get() != '\n');
//
//      a_amrGrids[lev].resize(ngrid);
//
//      for (int i = 0; i < ngrid; i++)
//      {
//        Box bx;
//        is >> bx;
//
//        while (is.get() != '\n');
//
//        // Quick check on box size
//        Box bxRef(bx);
//
//        if (bxRef.longside() > a_maxGridSize)
//        {
//          pout() << "Grid " << bx << " too large" << endl;
//          MayDay::Error();
//        }
//
//        if (a_verbosity >= 3)
//        {
//          pout() << bx << endl;
//        }
//
//        a_amrGrids[lev][i] = bx;
//      } // End loop over boxes on this level
//    } // End loop over levels
//  }
//
//  // Broadcast results to all the processors
//  broadcast(a_amrGrids,uniqueProc(SerialTask::compute));
}

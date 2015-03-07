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

#include "FABView.H"

#include "ParmParse.H"
#include "CH_HDF5.H"
#include "parstream.H"

#include "AMR.H"
#include "AMRLevel.H"
#include "AMRLevelMHDFactory.H"
#include "AMRLevelMHD.H"

#include "MHDPhysics.H"

#include "RotorMHDIBC.H"
#include "WaveMHDIBC.H"
// #include "RMMMHDIBC.H"
#include "ReconMHDIBC.H"

#include "BackwardEulerSolver.H"
#include "TGASolver.H"
#include "PoissonBC.H"

#ifdef CH_MPI
#include "CH_Attach.H"
#endif

#ifdef CH_LINUX
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
    _FPU_DEFAULT & ~(_FPU_MASK_IM | _FPU_MASK_ZM | _FPU_MASK_OM );  _FPU_SETCW(cw);
    // _FPU_DEFAULT & ~(_FPU_MASK_IM | _FPU_MASK_ZM | _FPU_MASK_OM | _FPU_MASK_UM);  _FPU_SETCW(cw);
}
#endif

// #define HALEM_PROC_SPEED
#ifdef HALEM_PROC_SPEED
#include <cstdio>
#include <sys/sysinfo.h>
#include <machine/hal_sysinfo.h>
#endif

#include "memusage.H"
#include "Timer.H"
Timer Everything    ("gov Everything", 0);
Timer TimeReadInput ("gov Read Input",   Everything);
Timer TimeSetupAMR  ("gov Setup AMR",    Everything);
Timer TimeRun       ("gov Run",          Everything);
Timer TimeConclude  ("gov Conclude",     Everything);

// Possible pressure relationships for the initial condition
#define PRESSURE_ISENTROPIC 0
#define PRESSURE_CONSTANT   1

// amrGodunov is a function (as opposed to inline in main()) to get
// around MPI scoping problems
void amrGodunov();

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
#ifdef CH_AIX
  H5dont_atexit();
#endif
  setChomboMPIErrorHandler();
#endif

  int rank, number_procs;

#ifdef CH_MPI
  MPI_Comm_rank(Chombo_MPI::comm, &rank);
  MPI_Comm_size(Chombo_MPI::comm, &number_procs);
#else
  rank = 0;
  number_procs = 1;
#endif

  if (rank == 0)
  {
    pout() << " number_procs = " << number_procs << endl;
  }

  Timer::TimerInit(rank);

  Everything.start();

#ifdef CH_LINUX
  // Abort when floating point exceptions occur
  trapfpe();
#endif

  // Check for an input file
  char* inFile = NULL;

  if (a_argc > 1)
  {
    inFile = a_argv[1];
  }
  else
  {
    pout() << "Usage:  amrGodunov...ex <inputfile>" << endl;
    pout() << "No input file specified" << endl;
    return -1;
  }

  // Parse the command line and the input file (if any)
  ParmParse pp(a_argc-2,a_argv+2,NULL,inFile);

  // Run amrGodunov, i.e., do the computation
  amrGodunov();

  Everything.stop();

#ifdef TIMER
  Real end_memory = get_memory_usage();

  pout() << endl
         << "Everything completed --- "
         << "mem: "
         << setw(8) << setprecision(3)
         << setiosflags(ios::fixed)
         << end_memory
         << " MB, time: "
         << setw(8) << setprecision(3)
         << setiosflags(ios::fixed)
         << Everything.wc_time()
         << " sec (wall-clock)" << endl << endl;
#endif

#if defined(TIMER) && defined(CH_MPI)
  Real avg_memory, min_memory, max_memory;
  gather_memory_from_procs(end_memory, avg_memory, min_memory, max_memory);
#endif

  Timer::TimerSummary();

#ifdef CH_MPI
  // Exit MPI
  dumpmemoryatexit();
  MPI_Finalize();
#endif
}

void amrGodunov()
{
  // Start timing the reading of the input file
  TimeReadInput.start();

  // Read inputs that are prefixed with "godonouv."
  ParmParse ppgodunov("godunov");

  // This determines the amount of diagnositic output generated
  int verbosity = 0;
  ppgodunov.query("verbosity",verbosity);
  CH_assert(verbosity >= 0);

  // For all gas dynamics
  // For all gas dynamics
  Real gamma = 1.667;
  ppgodunov.query("gamma",gamma);
  Real eta = 0.005;
  ppgodunov.query("eta",eta);
  Real mu = 0.01;
  ppgodunov.query("mu",mu);
  Real kappa = 0.01;
  ppgodunov.query("kappa",kappa);

  // Stop after this number of steps
  int nstop = 0;
  ppgodunov.get("max_step",nstop);

  // Stop when the simulation time get here
  Real stopTime = 0.0;
  ppgodunov.get("max_time",stopTime);

  // Set the physical size of the longest dimension of the domain
  Real domainLength = 1.0;
  ppgodunov.get("domain_length",domainLength);

  // Set the resolution of the coarsest level
  vector<int> numCells(SpaceDim);
  for (int i = 0; i < SpaceDim; ++i)
  {
    numCells[i] = 0;
  }
  ppgodunov.getarr("num_cells",numCells,0,SpaceDim);

  CH_assert(D_TERM(   (numCells[0] > 0),
                && (numCells[1] > 0),
                && (numCells[2] > 0)));
  CH_assert(D_TERM(   (numCells[0] % 2 == 0),
                && (numCells[1] % 2 == 0),
                && (numCells[2] % 2 == 0)));

  // Determine which spatial directions are periodic
  vector<int> isPeriodica(SpaceDim,0);
  bool isPeriodic[SpaceDim];

  ppgodunov.getarr("is_periodic",isPeriodica,0,SpaceDim);
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
  ppgodunov.get("max_level",maxLevel);
  int numReadLevels = Max(maxLevel,1);

  // Refinement ratios between levels
  std::vector<int> refRatios;
  // Note: this requires a refRatio to be defined for the finest level
  // (even though it will never be used)
  ppgodunov.getarr("ref_ratio",refRatios,0,numReadLevels+1);

  // Number of coarse time steps from one regridding to the next
  std::vector<int> regridIntervals;
  ppgodunov.getarr("regrid_interval",regridIntervals,0,numReadLevels);

  // How far to extend refinement from cells newly tagged for refinement
  int tagBufferSize = 3;
  ppgodunov.get("tag_buffer_size",tagBufferSize);

  // Threshold that triggers refinement
  Real refineThresh = 0.3;
  ppgodunov.get ("refine_thresh",refineThresh);

  // Minimum dimension of a grid
  int blockFactor = 1;
  ppgodunov.get("block_factor",blockFactor);

  // Maximum dimension of a grid
  int maxGridSize = 32;
  ppgodunov.get("max_grid_size",maxGridSize);

  Real fillRatio = 0.75;
  ppgodunov.get("fill_ratio",fillRatio);

  // Order of the normal predictor (PLM -> 1, PPM -> 2)
  std::string normalPred;
  int normalPredOrder;
  ppgodunov.get("normal_predictor",normalPred);
  if (normalPred == "PLM" || normalPred == "plm")
  {
    normalPredOrder = 1;
  }
  else if (normalPred == "PPM" || normalPred == "ppm")
  {
    normalPredOrder = 2;
  }
  else
  {
    MayDay::Error("Normal precdictor must by PLM or PPM");
  }

  // Use fourth order slopes
  int inFourthOrderSlopes = 1;
  bool useFourthOrderSlopes;
  ppgodunov.get("use_fourth_order_slopes",inFourthOrderSlopes);
  useFourthOrderSlopes = (inFourthOrderSlopes == 1);

  // Do slope limiting
  int inPrimLimiting = 1;
  bool usePrimLimiting;
  ppgodunov.get("use_prim_limiting",inPrimLimiting);
  usePrimLimiting = (inPrimLimiting == 1);

  // Do slope limiting using characteristics
  int inCharLimiting = 0;
  bool useCharLimiting;
  ppgodunov.get("use_char_limiting",inCharLimiting);
  useCharLimiting = (inCharLimiting == 1);

  // Do slope flattening
  int inFlattening = 1;
  bool useFlattening;
  ppgodunov.get("use_flattening",inFlattening);
  useFlattening = (inFlattening == 1);

  // Apply artificial viscosity
  int inArtificialViscosity = 1;
  bool useArtificialViscosity;
  ppgodunov.get("use_artificial_viscosity",inArtificialViscosity);
  useArtificialViscosity = (inArtificialViscosity == 1);

  // Artificial viscosity coefficient/multiplier
  Real artificialViscosity = 0.1;
  ppgodunov.get("artificial_viscosity",artificialViscosity);

  // Apply filter to magnetic field
  int inFilterBField=1;
  bool doFilterBField;
  ppgodunov.get("filter_BField",inFilterBField);
  doFilterBField = (inFilterBField==1);

  // Set up checkpointing
  int checkpointInterval = 0;
  ppgodunov.query("checkpoint_interval",checkpointInterval);

  // Set up plot file writing
  int plotInterval = 0;
  ppgodunov.query("plot_interval",plotInterval);

  // CFL multiplier
  Real cfl = 0.8;
  ppgodunov.get("cfl",cfl);

  // Initial CFL multiplier
  Real initialCFL = 0.1;
  ppgodunov.get("initial_cfl",initialCFL);

  // Determine if a fixed or variable time step will be used
  Real fixedDt = -1;
  ppgodunov.query("fixed_dt",fixedDt);

  // Limit the time step growth
  Real maxDtGrowth = 1.1;
  ppgodunov.get("max_dt_growth",maxDtGrowth);

  // Let the time step grow by this factor above the "maximum" before
  // reducing it
  Real dtToleranceFactor = 1.1;
  ppgodunov.get("dt_tolerance_factor",dtToleranceFactor);

  // End timing the reading of the input file
  TimeReadInput.stop();

#ifdef TIMER
  pout() << "Input Read completed --- "
         << "mem: "
         << setw(8) << setprecision(3)
         << setiosflags(ios::fixed)
         << get_memory_usage()
         << " MB, time: "
         << setw(8) << setprecision(3)
         << setiosflags(ios::fixed)
         << TimeReadInput.wc_time()
         << " sec (wall-clock)" << endl;
#endif

  // Start timing AMR solver setup
  TimeSetupAMR.start();

  // Create and define IBC (initial and boundary condition) object
  PhysIBC* ibc;

  // A minimum pressure needed to construct MHDPhysics - used in slope
  // flattening
  Real smallPressure;

  // Domain ghost BCs for specific problems - used for diffusion solves
  DomainGhostBC dombcB;
  DomainGhostBC dombcV;
  DomainGhostBC dombcT;
  //Default domain ghost BCs
  for (int dir=0; dir<SpaceDim; dir++)
    {
      SideIterator sit;
      for (sit.reset(); sit.ok(); ++sit)
        {
	  //	           scalarDirichletBC thisBC(dir,sit());
	  NeumannBC thisBC(dir,sit(),Interval(0,2));
          dombcB.setBoxGhostBC(thisBC);
          dombcV.setBoxGhostBC(thisBC);
	  NeumannBC thisBCT(dir,sit());
          dombcT.setBoxGhostBC(thisBCT);
        }
    }
  // Determine the sample problem specified
  std::string problemString;
  if (ppgodunov.contains("problem"))
  {
    ppgodunov.query("problem",problemString);

    // Print some parameters
    if (verbosity >= 2)
    {
      pout() << "problem = " << problemString << endl;
      pout() << "gamma = " << gamma << endl;
    }

    if (problemString == "rotor")
    {
      // Rotor problem

      // Define IBC for rotor problem
      RotorMHDIBC* rotorMHDibc = new RotorMHDIBC;
      rotorMHDibc->setFortranCommon(smallPressure,
                                    gamma,mu,eta,kappa);
      ibc = rotorMHDibc;
     
      //Setup domain ghost BCs for diffusion solves
      // Zero gradient for all for the Rotor problem
      for (int dir=0; dir<SpaceDim; dir++)
	{
	  SideIterator sit;
	  for (sit.reset(); sit.ok(); ++sit)
	    {
	      //	           scalarDirichletBC thisBC(dir,sit());
	      NeumannBC thisBC(dir,sit(),Interval(0,2));
	      dombcB.setBoxGhostBC(thisBC);
	      dombcV.setBoxGhostBC(thisBC);
	      NeumannBC thisBCT(dir,sit());
	      dombcT.setBoxGhostBC(thisBCT);
	    }
	}

    }
    else if (problemString == "wave")
    {
      // Plane wave problem
      Real alpha = 45;
      ppgodunov.get("alpha",alpha);

      int  pdir = 1;
      ppgodunov.get("pdir",pdir);

      int  kratio = 1;
      ppgodunov.get("kratio",kratio);

      int  waveNumber = 3;
      ppgodunov.get("waveNumber",waveNumber);

      Real pertAmplitude = 0.001;
      ppgodunov.get("pertAmplitude",pertAmplitude);

      if (verbosity >= 2)
      {
        pout() << "alpha = " << alpha << endl;
        pout() << "pdir = " << pdir << endl;
        pout() << "kratio = " << kratio << endl;
        pout() << "waveNumber = " << waveNumber << endl;
        pout() << "pertAmplitude = " << pertAmplitude << endl;
      }

      // Define IBC for plane wave problem
      WaveMHDIBC* waveMHDibc = new WaveMHDIBC;
      waveMHDibc->setFortranCommon(smallPressure,
                                   gamma,
				   mu,
				   eta,
				   kappa,
                                   alpha,
                                   pdir,
                                   kratio,
                                   waveNumber,
                                   pertAmplitude);
      ibc = waveMHDibc;
      // Domain ghost BCs can be default - these are not
      // used anyway because the wave problem is computed with
      // periodic BCs
    }
#if 0
    else if (problemString == "rm")
    {
      // RM problem
      Real alpha=45;
      ppgodunov.get("alpha",alpha);

      Real mach=1.2;
      ppgodunov.get("mach",mach);

      Real dratio = 3;
      ppgodunov.get("dratio",dratio);

      Real betainv=1;
      ppgodunov.get("betainv",betainv);

      if (verbosity >= 2)
      {
        pout() << "alpha = " << alpha << endl;
        pout() << "mach = " << mach << endl;
        pout() << "dratio = " << dratio << endl;
        pout() << "betainv = " << betainv << endl;
      }

      // Define IBC for RM problem
      RMMHDIBC* rmMHDibc = new RMMHDIBC;
      rmMHDibc->setFortranCommon(smallPressure,
                                 gamma,
                                 alpha,
                                 mach,
                                 dratio,
                                 betainv);
      ibc = rmMHDibc;
    }
#endif
    else if (problemString == "recon")
    {
      // Reconnection problem

      // Define IBC for reconnection problem
      ReconMHDIBC* reconMHDibc = new ReconMHDIBC;
      //      reconMHDibc->setFortranCommon(smallPressure,
      //                                    gamma,);
      
      reconMHDibc->setParameters(gamma,mu,eta,kappa);
      ibc = reconMHDibc;
      // Zero gradient for all except 
      // no slip for velocity
      // specified wall temperature 
      // zero normal B
      // Periodic in x
      //BCs for Y
      
      {
      int dir=1;
      SideIterator sit;
      for (sit.reset(); sit.ok(); ++sit)
	{
	  pout() << "Recon BCs" <<sit() << "dir=" << dir<<endl;
	  // No slip on all components of velocity
	  DirichletBC velBC(dir,sit(),Interval(0,2));
	  dombcV.setBoxGhostBC(velBC);
	  NeumannBC BxBC(dir,sit(),Interval(0,0));
	  dombcB.setBoxGhostBC(BxBC);
	  DirichletBC ByBC(dir,sit(),Interval(1,1));
	  dombcB.setBoxGhostBC(ByBC);
	  NeumannBC BzBC(dir,sit(),Interval(2,2));
	  dombcB.setBoxGhostBC(BzBC);
	  NeumannBC thisBCT(dir,sit());
	  dombcT.setBoxGhostBC(thisBCT);
	}
      }
    }
    else
    {
      // The sample problem name given isn't valid
      pout() << "Invalid problem, \"" << problemString << "\", specified in input file" << endl << endl;
      return;
    }
  }
  else
  {
    // A sample problem must be specified
    pout() << "\"godunov.problem\" not specified in input file" << endl << endl;
    return;
  }

  if (verbosity >= 2)
  {
    pout() << "verbosity = " << verbosity << endl;

    pout() << "maximum step = " << nstop << endl;
    pout() << "maximum time = " << stopTime << endl;
    if (fixedDt > 0)
    {
      pout() << "fixed dt = " << fixedDt << endl;
    }

    pout() << "number of cells = " << D_TERM(numCells[0] << "  " <<,
                                             numCells[1] << "  " <<,
                                             numCells[2] << ) endl;
    pout() << "is period = " << D_TERM(isPeriodic[0] << "  " <<,
                                       isPeriodic[1] << "  " <<,
                                       isPeriodic[2] << ) endl;

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
    pout() << "tag buffer size = " << tagBufferSize << endl;

    pout() << "refinement threshold = " << refineThresh << endl;

    pout() << "blocking factor = " << blockFactor << endl;
    pout() << "max grid size = " << maxGridSize << endl;
    pout() << "fill ratio = " << fillRatio << endl;

    pout() << "normal predictor = ";
    if (normalPredOrder == 1)
    {
      pout() << "PLM" << endl;
    }
    else if (normalPredOrder == 2)
    {
      pout() << "PPM" << endl;
    }
    else
    {
      pout() << "Unknown (" << normalPredOrder << ")" << endl;
    }

    pout() << "slope order = "
           << (useFourthOrderSlopes ? "2nd" : "4th") << endl;
    pout() << "use primitive slope limiting = "
           << (usePrimLimiting ? "yes" : "no") << endl;
    pout() << "use characteristic slope limiting = "
           << (useCharLimiting ? "yes" : "no") << endl;
    pout() << "use slope flattening = "
           << (useFlattening ? "yes" : "no") << endl;

    pout() << "use artificial viscosity = "
           << (useArtificialViscosity ? "yes" : "no") << endl;
    if (useArtificialViscosity)
    {
      pout() << "artificial viscosity = " << artificialViscosity << endl;
    }

    pout() << "checkpoint interval = " << checkpointInterval << endl;
    pout() << "plot interval = " << plotInterval << endl;

    pout() << "CFL = " << cfl << endl;
    pout() << "initial CFL = " << initialCFL << endl;

    pout() << "maximum dt growth = " << maxDtGrowth << endl;
    pout() << "dt tolerance factor = " << dtToleranceFactor << endl;
  }

  ProblemDomain probDomain (IntVect::Zero,
                            IntVect(D_DECL(numCells[0]-1,
                                           numCells[1]-1,
                                           numCells[2]-1)),
                            isPeriodic);

  // Set up the physics for magneto-hydrodynamics
  MHDPhysics mhdPhysics(smallPressure);
  mhdPhysics.setPhysIBC(ibc);

  // Cast to physics base class pointer for technical reasons
  GodunovPhysics* godunovPhysics = static_cast<GodunovPhysics*> (&mhdPhysics);

  // Set up the AMRLevel... factory
  AMRLevelMHDFactory amrGodFact;

  BaseHeatSolver* diffusionSolverPtr;

  diffusionSolverPtr= (BaseHeatSolver *) new TGASolver;

//   DomainGhostBC dombcB;
//   DomainGhostBC dombcV;
//   DomainGhostBC dombcT;
//   for (int dir=0; dir<SpaceDim; dir++)
//     {
//       SideIterator sit;
//       for (sit.reset(); sit.ok(); ++sit)
//         {
// 	  //	           scalarDirichletBC thisBC(dir,sit());
// 	  NeumannBC thisBC(dir,sit(),Interval(0,2));
//           dombcB.setBoxGhostBC(thisBC);
//           dombcV.setBoxGhostBC(thisBC);
// 	  NeumannBC thisBCT(dir,sit());
//           dombcT.setBoxGhostBC(thisBCT);
//         }
//     }
  amrGodFact.define(cfl,
                    domainLength,
                    verbosity,
                    refineThresh,
                    tagBufferSize,
                    initialCFL,
                    godunovPhysics,
                    normalPredOrder,
                    useFourthOrderSlopes,
                    usePrimLimiting,
                    useCharLimiting,
                    useFlattening,
                    useArtificialViscosity,
                    artificialViscosity,
                    doFilterBField,
		    diffusionSolverPtr,
		    dombcB,
		    dombcV,
		    dombcT,
		    gamma,
		    mu,
		    eta,
		    kappa);
  
  AMR amr;

  // Set up the AMR object
  amr.define(maxLevel,refRatios,probDomain,&amrGodFact);

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
  if (ppgodunov.contains("plot_prefix"))
  {
    std::string prefix;
    ppgodunov.query("plot_prefix",prefix);
    amr.plotPrefix(prefix);
  }

  if (ppgodunov.contains("chk_prefix"))
  {
    std::string prefix;
    ppgodunov.query("chk_prefix",prefix);
    amr.checkpointPrefix(prefix);
  }

  amr.verbosity(verbosity);

  // Set up input files
  if (!ppgodunov.contains("restart_file"))
  {
    if (!ppgodunov.contains("fixed_hierarchy"))
    {
      // initialize from scratch for AMR run
      // initialize hierarchy of levels
      amr.setupForNewAMRRun();
    }
    else
    {
      std::string gridFile;
      ppgodunov.query("fixed_hierarchy",gridFile);

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
    ppgodunov.query("restart_file",restartFile);

#ifdef CH_USE_HDF5
    HDF5Handle handle(restartFile,HDF5Handle::OPEN_RDONLY);
    // read from checkpoint file
    amr.setupForRestart(handle);
    handle.close();
#else
    MayDay::Error("amrGodunov restart only defined with hdf5");
#endif
  }

  // End timing AMR solver setup
  TimeSetupAMR.stop();

#ifdef TIMER
  pout() << "AMR Setup completed ---- "
         << "mem: "
         << setw(8) << setprecision(3)
         << setiosflags(ios::fixed)
         << get_memory_usage()
         << " MB, time: "
         << setw(8) << setprecision(3)
         << setiosflags(ios::fixed)
         << TimeSetupAMR.wc_time()
         << " sec (wall-clock)" << endl;

  if (verbosity >= 1)
  {
    pout() << endl;
  }
#endif

  // Run and time the computation
  TimeRun.start();
  amr.run(stopTime,nstop);
  TimeRun.stop();

#ifdef TIMER
  if (verbosity >= 1)
  {
    pout() << endl;
  }

  pout() << "AMR Run completed ------ "
         << "mem: "
         << setw(8) << setprecision(3)
         << setiosflags(ios::fixed)
         << get_memory_usage()
         << " MB, time: "
         << setw(8) << setprecision(3)
         << setiosflags(ios::fixed)
         << TimeRun.wc_time()
         << " sec (wall-clock)" << endl;
#endif

  // Output the last plot file and statistics - time the process
  TimeConclude.start();
  amr.conclude();
  TimeConclude.stop();

#ifdef TIMER
  pout() << "AMR Conclude completed - "
         << "mem: "
         << setw(8) << setprecision(3)
         << setiosflags(ios::fixed)
         << get_memory_usage()
         << " MB, time: "
         << setw(8) << setprecision(3)
         << setiosflags(ios::fixed)
         << TimeConclude.wc_time()
         << " sec (wall-clock)" << endl;
#endif
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

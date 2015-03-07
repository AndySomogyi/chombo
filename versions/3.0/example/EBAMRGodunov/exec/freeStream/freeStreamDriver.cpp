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

#include "ParmParse.H"
#include "CH_HDF5.H"
#include "parstream.H"

#include "EBPlanarShockIBCFactory.H"
#include "ModianoIBCFactory.H"
#include "EBPatchPolytropicFactory.H"

#include "EBPatchPolytropic.H"

#include "EBAMRGodunovFactory.H"
#include "EBAMRGodunov.H"
#include "AMRLevel.H"
#include "AMR.H"
#include "BaseIVFactory.H"
#include "EBIndexSpace.H"
#include "EBISLayout.H"
#include "BoxIterator.H"
#include "ParmParse.H"
#include "BRMeshRefine.H"
#include "LoadBalance.H"
#include "GeometryShop.H"
#include "PolyGeom.H"
#include "LevelData.H"
#include "EBCellFAB.H"
#include "EBCellFactory.H"
#include "VoFIterator.H"
#include "EBArith.H"
#include "AllRegularService.H"
#include "EBLevelRedist.H"
#include "RedistStencil.H"
#include "SlabService.H"
#include "TiltedCylinderIF.H"

#include "DebugDump.H"
#include "EBDebugDump.H"

#ifdef CH_MPI
#include "CH_Attach.H"
#endif

#if CH_SPACEDIM==2
IntVect ivdebug(3, 5);
#elif CH_SPACEDIM==3
IntVect ivdebug(1, 1, 0);
#endif

void makeGeometry(Box& a_coarsestDomain,
                  Real& a_dx);

void amrGodunov(const Box&  coarsestDomain,
                const Real& dx);

void setupFixedGrids(Vector<Vector<Box> >& a_amrGrids,
                     const ProblemDomain&  a_coarsestDomain,
                     int                   a_max_level,
                     int                   a_max_grid_size,
                     int                   a_block_factor,
                     int                   a_verbosity,
                     std::string           a_grid_file);

/***************/
/***************/
int
main(int a_argc, char* a_argv[])
{
#ifdef CH_MPI
  MPI_Init(&a_argc,&a_argv);
  // setChomboMPIErrorHandler();
#endif

  // Check for an input file
  char* inFile = NULL;

  if (a_argc > 1)
    {
      inFile = a_argv[1];
    }
  else
    {
      pout() << "Usage: <executable name> <inputfile>" << endl;
      pout() << "No input file specified" << endl;
      return -1;
    }
  // Parse the command line and the input file (if any)
  ParmParse pp(a_argc-2,a_argv+2,NULL,inFile);

  Box coarsestDomain;
  Real dx;
  // run amrGodunov
  makeGeometry(coarsestDomain, dx);
  amrGodunov(coarsestDomain, dx);

#ifdef CH_MPI
  dumpmemoryatexit();
  MPI_Finalize();
#endif
}

/***************/
/***************/
void
makeGeometry(Box& a_coarsestDomain,
             Real& a_dx)
{
  ParmParse ppgodunov;
  //parse input file
  int max_level = 0;
  ppgodunov.get("max_level",max_level);

  int num_read_levels = Max(max_level,1);
  std::vector<int> refRatios; // (num_read_levels,1);
  // note this requires a refRatio to be defined for the
  // finest level (even though it will never be used)

  ppgodunov.getarr("ref_ratio",refRatios,0,num_read_levels+1);
  ParmParse pp;
  RealVect origin = RealVect::Zero;
  Vector<int> n_cell(SpaceDim);
  pp.getarr("n_cell",n_cell,0,SpaceDim);

  CH_assert(n_cell.size() == SpaceDim);
  IntVect lo = IntVect::Zero;
  IntVect hi;
  for(int ivec = 0; ivec < SpaceDim; ivec++)
    {
      if(n_cell[ivec] <= 0)
        {
          pout() << " bogus number of cells input = " << n_cell[ivec];
          MayDay::Error();
        }
      hi[ivec] = n_cell[ivec] - 1;
    }

  a_coarsestDomain = Box(lo, hi);
  Box finestDomain = a_coarsestDomain;
  Real prob_hi;
  int numOpen = n_cell[0];
  pp.get("domain_length",prob_hi);
  a_dx = prob_hi/numOpen;
  Real fineDx = a_dx;
  for(int ilev = 0; ilev < max_level; ilev++)
    {
      finestDomain.refine(refRatios[ilev]);
      fineDx /= refRatios[ilev];
    }

  RealVect cylinderAxis;
  vector<Real>  cylinderAxisVect(SpaceDim);
  pp.getarr("cylinder_axis",cylinderAxisVect, 0, SpaceDim);
  for(int idir = 0; idir < SpaceDim; idir++)
    {
      cylinderAxis[idir] = cylinderAxisVect[idir];
    }
  Real cylinderRadius;
  pp.get("cylinder_radius",cylinderRadius);
  pout() << "using a tilted cylinder implicit function" << endl;
  RealVect corner = RealVect::Zero;
  bool negativeInside = true;
  TiltedCylinderIF tunnel(cylinderRadius, cylinderAxis, corner, negativeInside);

  RealVect vectDx = RealVect::Unit;
  vectDx *= fineDx;

  GeometryShop workshop(tunnel,0,vectDx);
  EBIndexSpace* ebisPtr = Chombo_EBIS::instance();
  int ebMax = 1024;
  ebisPtr->define(finestDomain, origin, fineDx, workshop, ebMax, ebMax);
}
/***************/
/***************/

void amrGodunov(const Box& a_domain,
                const Real& a_dx)
{
  // read inputs
  ParmParse pp;

  int verbosity;
  pp.get("verbosity",verbosity);
  CH_assert(verbosity >= 0);

  Real gamma = 1.4;
  pp.get("gamma",gamma);


  int nstop = 0;
  pp.get("max_step",nstop);

  int redistRad = 0;
  pp.get("redist_radius",redistRad);

  Real stopTime = 0.0;
  pp.get("max_time",stopTime);

  bool is_periodic[SpaceDim];
  for (int i = 0; i < SpaceDim; ++i)
    {
      is_periodic[i] = false;
    }

  int max_level = 0;
  pp.get("max_level",max_level);

  int num_read_levels = Max(max_level,1);
  std::vector<int> ref_ratios; // (num_read_levels,1);
  // note this requires a ref_ratio to be defined for the
  // finest level (even though it will never be used)
  pp.getarr("ref_ratio",ref_ratios,0,num_read_levels+1);

  std::vector<int> regrid_intervals; // (num_read_levels,1);
  pp.getarr("regrid_interval",regrid_intervals,0,num_read_levels);

  Real refineThresh = 0.3;
  pp.get ("refine_thresh",refineThresh);

  int block_factor = 1;
  pp.get("block_factor",block_factor);

  int tagBufferSize = 3;
  pp.get("tag_buffer_size",tagBufferSize);

  int max_grid_size = 32;
  pp.get("max_grid_size",max_grid_size);

  Real fill_ratio = 0.75;
  pp.get("fill_ratio",fill_ratio);

  int checkpoint_interval = 0;
  pp.get("checkpoint_interval",checkpoint_interval);

  int plot_interval = 0;
  pp.get("plot_interval",plot_interval);

  Real cfl = 0.8;
  pp.get("cfl",cfl);

  Real initialCFL = 0.1;
  pp.get("initial_cfl",initialCFL);

  Real fixed_dt = -1;
  pp.get("fixed_dt",fixed_dt);

  Real max_dt_growth = 1.1;
  pp.get("max_dt_growth",max_dt_growth);

  Real dt_tolerance_factor = 1.1;
  pp.get("dt_tolerance_factor",dt_tolerance_factor);

  vector<Real> length(SpaceDim);
  pp.getarr("domain_length",length,0,SpaceDim);
  RealVect domainLength;
  for(int idir=0;idir<SpaceDim;idir++)
    {
      domainLength[idir] = length[idir];
    }

  int ifourth, iflatten, iartvisc;
  pp.get("use_fourth_order_slopes", ifourth);
  pp.get("use_flattening"         , iflatten);
  pp.get("use_art_visc"           , iartvisc);
  bool useFourthOrderSlopes = (ifourth  ==1);
  bool useFlattening        = (iflatten ==1);
  bool useArtificialVisc    = (iartvisc ==1);

  ProblemDomain prob_domain(a_domain.smallEnd(),
                            a_domain.bigEnd(),
                            is_periodic);

  // create initial/boundary condition object

  // create initial/boundary condition object
  pout() << "Modiano initial and boundary conditions" << endl;

  RealVect modianoAxis;

  Real waveAmp, waveWidth;
  pp.get("wave_amplitude", waveAmp);
  pp.get("wave_width", waveWidth);

  pout() << "uses cylinder axis for direction" << endl;
  vector<Real>  cylinderAxisVect(SpaceDim);
  pp.getarr("cylinder_axis",cylinderAxisVect, 0, SpaceDim);
  RealVect cylinderAxis;
  for(int idir = 0; idir < SpaceDim; idir++)
    {
      cylinderAxis[idir] = cylinderAxisVect[idir];
    }
  Real sum;
  PolyGeom::unifyVector(cylinderAxis, sum);
  modianoAxis = cylinderAxis;

  vector<Real> centerpp(SpaceDim,0.5);
  RealVect center;
  pp.getarr("wave_center",centerpp,0,SpaceDim);
  for (int i = 0; i < SpaceDim; i++)
    center[i] = centerpp[i];


  int idoNegWave;
  pp.get("use_negative_wave",idoNegWave);
  bool useNegativeWave = (idoNegWave==1);
  if(useNegativeWave)
    {
      pout() << "u - c wave " << endl;
    }
  else
    {
      pout() << "u + c wave " << endl;
    }

  int idoFreeStream;
  pp.get("free_stream_prob",idoFreeStream);
  bool doFreeStream = (idoFreeStream==1);
  if(doFreeStream)
    {
      pout() << "free stream modiano bump " << endl;
    }
  else
    {
      pout() << "simple wave prob " << endl;
    }

  bool doSmushing = true;
  bool doRZCoords = false;
  bool useLimiting = true;
  bool hasSourceTerm = false;
  bool tagAll = false;

  ModianoIBCFactory bcfactory(gamma, waveAmp, waveWidth,
                              center, modianoAxis, doFreeStream, useNegativeWave);

  int iusemassredist;
  pp.get("use_mass_redist", iusemassredist);
  bool useMassRedist    = (iusemassredist ==1);
  //create patch integrator
  EBPatchPolytropicFactory patchGamma(&bcfactory,
                                      gamma,
                                      useFourthOrderSlopes,
                                      useFlattening,
                                      useArtificialVisc,
                                      useLimiting,
                                      doRZCoords);

  EBAMRGodunovFactory amrg_fact(initialCFL,
                                cfl,
                                redistRad,
                                domainLength,
                                refineThresh,
                                tagBufferSize,
                                verbosity,
                                useMassRedist,
                                doSmushing,
                                doRZCoords,
                                hasSourceTerm,
                                &patchGamma,
                                tagAll);

  AMR amr;
  amr.define(max_level,ref_ratios,prob_domain,&amrg_fact);

  if (fixed_dt > 0)
    {
      amr.fixedDt(fixed_dt);
    }

  // set grid generation parameters
  amr.maxGridSize(max_grid_size);
  amr.blockFactor(block_factor);
  amr.fillRatio(fill_ratio);

  // the hyperbolic codes use a grid buffer of 1
  amr.gridBufferSize(1);

  // set output parameters
  amr.checkpointInterval(checkpoint_interval);
  amr.plotInterval(plot_interval);
  amr.regridIntervals(regrid_intervals);
  amr.maxDtGrow(max_dt_growth);
  amr.dtToleranceFactor(dt_tolerance_factor);

  if (pp.contains("plot_prefix"))
    {
      std::string prefix;
      pp.get("plot_prefix",prefix);
      amr.plotPrefix(prefix);
    }

  if (pp.contains("chk_prefix"))
    {
      std::string prefix;
      pp.get("chk_prefix",prefix);
      amr.checkpointPrefix(prefix);
    }

  amr.verbosity(verbosity);

  std::string base_string;
  pp.get("fixed_hierarchy",base_string);

  //append .spacedim d

  const char* base_cstr  = base_string.c_str();
  char filename[100];
  sprintf(filename,"%s.%dd",base_cstr, SpaceDim );
  string grid_file(filename);

  // initialize from a list of grids in "grid_file"
  Vector<Vector<Box> > amrGrids(max_level+1);
  setupFixedGrids(amrGrids,
                  prob_domain,
                  max_level,
                  max_grid_size,
                  block_factor,
                  verbosity,
                  grid_file);
  amr.setupForFixedHierarchyRun(amrGrids,1);

  // run
  amr.run(stopTime,nstop);

  // output last pltfile and statistics
  amr.conclude();
  //cleanup
}

void setupFixedGrids(Vector<Vector<Box> >& a_amrGrids,
                     const ProblemDomain&  a_domain,
                     int                   a_max_level,
                     int                   a_max_grid_size,
                     int                   a_block_factor,
                     int                   a_verbosity,
                     std::string           a_grid_file)
{
  if (procID() == uniqueProc(SerialTask::compute))
    {
      a_amrGrids.push_back(Vector<Box>(1,a_domain.domainBox()));

      // read in predefined grids
      ifstream is(a_grid_file.c_str(), ios::in);

      if (is.fail())
        {
          MayDay::Error("Cannot open grids file");
        }

      // format of file -- number of levels, then for each level starting
      // with level 1, number of grids on level, list of boxes

      int in_numLevels;
      is >> in_numLevels;

      CH_assert (in_numLevels <= a_max_level+1);

      pout() << "numLevels = " << in_numLevels << endl;

      while (is.get() != '\n');

      a_amrGrids.resize(in_numLevels);

      // check to see if coarsest level needs to be broken up
      domainSplit(a_domain,a_amrGrids[0],a_max_grid_size,a_block_factor);

      if (a_verbosity >= 3)
        {
          pout() << "level 0: ";
          for (int n = 0; n < a_amrGrids[0].size(); n++)
            {
              pout() << a_amrGrids[0][0] << endl;
            }
        }

      // now loop over levels, starting with level 1
      int ngrid;
      for (int lev = 1; lev < in_numLevels; lev++)
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

              // quick check on box size
              Box bxRef(bx);

              if (bxRef.longside() > a_max_grid_size)
                {
                  pout() << "Grid " << bx << " too large" << endl;
                  MayDay::Error();
                }

              if (a_verbosity >= 3) {
                pout() << bx << endl;
              }

              a_amrGrids[lev][i] = bx;
            } // end loop over boxes on this level
        } // end loop over levels
    }

  broadcast(a_amrGrids,uniqueProc(SerialTask::compute));
}


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
#include <stdio.h>

#include "ParmParse.H"
#include "CH_HDF5.H"
#include "parstream.H"

#include "EBExplosionIBCFactory.H"
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
#include "EBAMRIO.H"
#include "DebugDump.H"
#include "EBDebugDump.H"
#include "AllRegularService.H"
#include "TiltedCylinderIF.H"
#include "GodunovGeom.H"

#ifdef CH_MPI
#include "CH_Attach.H"
#endif

#if CH_SPACEDIM==2
IntVect ivdebug(3, 5);
#elif CH_SPACEDIM==3
IntVect ivdebug(1, 1, 0);
#endif

/***************/
Vector<string> getNames()
{
  Vector<string> retval;
  EBPatchPolytropic patchInt;
  Vector<string> consNames  = patchInt.stateNames();
  Vector<string> primNames  = patchInt.primNamesNoLog();
  retval = consNames;
  retval.append(primNames);
  return retval;
}
/***************/
void
generateExactData(LevelData<EBCellFAB>& a_exactData,
                  const DisjointBoxLayout& a_grids,
                  const EBISLayout&        a_ebisl,
                  const Box& a_domain,
                  const Real& a_dx,
                  const Real& a_finalTime,
                  const Real& a_dt,
                  const int&  a_nvar)
{
  EBCellFactory ebcf(a_ebisl);
  a_exactData.define(a_grids, a_nvar, IntVect::Zero, ebcf);

  RealVect modianoAxis, center;
  Real waveAmp, waveWidth, gamma;

  ParmParse pp;
  pp.get("wave_amplitude", waveAmp);
  pp.get("wave_width", waveWidth);
  pp.get("gamma",gamma);
  int idoNegWave;
  pp.get("use_negative_wave",idoNegWave);
  bool useNegativeWave = (idoNegWave==1);

  vector<Real> centerpp(SpaceDim,0.5);
  pp.getarr("wave_center",centerpp,0,SpaceDim);

  for(int idir = 0; idir < SpaceDim; idir++)
    {
      center[idir] = centerpp[idir];
    }

  int testverbosity;
  pp.get("testverbosity", testverbosity);
  int whichGeom;
  pp.get("which_geom", whichGeom);
  if(whichGeom == 3)
    {
      if(testverbosity > 1)
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
    }
  else
    {
      if(testverbosity > 1)
        pout() << "uses channel slope for direction" << endl;
      RealVect channelNormal;
      vector<Real>  channelNormalVect(SpaceDim);
      pp.getarr("channel_normal",channelNormalVect, 0, SpaceDim);
      for(int idir = 0; idir < SpaceDim; idir++)
        {
          channelNormal[idir] = channelNormalVect[idir];
        }
      modianoAxis = mibcGetNormalVector(channelNormal);;
    }


  int idoFreeStream;
  pp.get("free_stream_prob",idoFreeStream);
  bool doFreeStream = (idoFreeStream==1);

  ModianoIBC mibc(gamma, waveAmp, waveWidth,
                  center, modianoAxis, doFreeStream, useNegativeWave);

  bool is_periodic[SpaceDim];
  for (int i = 0; i < SpaceDim; ++i)
    {
      is_periodic[i] = false;
    }
  ProblemDomain probdomain(a_domain.smallEnd(), a_domain.bigEnd(), is_periodic);
  mibc.define(probdomain, a_dx);
  mibc.setToExactConsAndPrim(a_exactData, a_ebisl, a_finalTime);
}
/************/
/************/
void
generateCalculatedData(LevelData<EBCellFAB>& a_calcData,
                       Real& a_dt,
                       Real& a_finalTime,
                       DisjointBoxLayout& a_grids,
                       EBISLayout&        a_ebisl,
                       const Box&  a_domain,
                       const Real& a_dx,
                       const int&  a_nstop)
{
  AMR amr;
  // read inputs
  ParmParse pp;

  int testverbosity;
  pp.get("testverbosity", testverbosity);
  int verbosity;
  pp.get("verbosity",verbosity);
 CH_assert(verbosity >= 0);

  Real gamma = 1.4;
  pp.get("gamma",gamma);


  int redistRad = 0;
  pp.get("redist_radius",redistRad);

  int nstop = a_nstop;

  //do not want this stopping calc
  Real stopTime = 1.0e6;

  bool is_periodic[SpaceDim];
  for (int i = 0; i < SpaceDim; ++i)
    {
      is_periodic[i] = false;
    }

  //single level calc
  int max_level = 0;
  // note this requires a ref_ratio to be defined for the
  // finest level (even though it will never be used)
  std::vector<int> ref_ratios(1,1);

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
  Real domainLength;
  pp.get("domain_length",domainLength);

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
  if(testverbosity > 1)
    pout() << "Modiano initial and boundary conditions" << endl;

  RealVect modianoAxis;

  Real waveAmp, waveWidth;
  pp.get("wave_amplitude", waveAmp);
  pp.get("wave_width", waveWidth);
  int whichGeom;
  pp.get("which_geom", whichGeom);
  if(whichGeom == 3)
    {
      if(testverbosity > 1)
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
    }
  else
    {
      if(testverbosity > 1)
        pout() << "uses channel slope for direction" << endl;
      RealVect channelNormal;
      vector<Real>  channelNormalVect(SpaceDim);
      pp.getarr("channel_normal",channelNormalVect, 0, SpaceDim);
      for(int idir = 0; idir < SpaceDim; idir++)
        {
          channelNormal[idir] = channelNormalVect[idir];
        }
      modianoAxis = mibcGetNormalVector(channelNormal);;
    }

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
      if(testverbosity > 1)
        pout() << "u - c wave " << endl;
    }
  else
    {
      if(testverbosity > 1)
        pout() << "u + c wave " << endl;
    }

  int idoFreeStream;
  pp.get("free_stream_prob",idoFreeStream);
  bool doFreeStream = (idoFreeStream==1);
  if(doFreeStream)
    {
      if(testverbosity > 1)
        pout() << "free stream modiano bump " << endl;
    }
  else
    {
      if(testverbosity > 1)
        pout() << "simple wave prob " << endl;
    }

  ModianoIBCFactory bcfactory(gamma, waveAmp, waveWidth,
                              center, modianoAxis, doFreeStream, useNegativeWave);

  int uselim;
  pp.get("use_limiting", uselim);
  bool useLimiting = (uselim==1);
  bool doRZCoords = false;
  bool hasSourceTerm = false;
  bool doSmushing = true;
  //create patch integrator
  EBPatchPolytropicFactory patchGamma(&bcfactory,
                                      gamma,
                                      useFourthOrderSlopes,
                                      useFlattening,
                                      useArtificialVisc,
                                      useLimiting,
                                      doRZCoords);

  //dummies
  Real refineThresh = 0.7;
  int tagBufferSize = 1;
  int iusemassredist;
  pp.get("use_mass_redist", iusemassredist);
  bool useMassRedist    = (iusemassredist ==1);

  EBAMRGodunovFactory amrg_fact(initialCFL,
                                cfl,
                                redistRad,
                                domainLength*RealVect::Unit,
                                refineThresh,
                                tagBufferSize,
                                verbosity,
                                useMassRedist,
                                doSmushing,
                                doRZCoords,
                                hasSourceTerm,
                                &patchGamma);

  int max_grid_size = 32;
  pp.get("max_grid_size",max_grid_size);

  amr.define(max_level,ref_ratios,prob_domain,&amrg_fact);
  amr.maxGridSize(max_grid_size);

  if (fixed_dt > 0)
    {
      amr.fixedDt(fixed_dt);
    }

  // the hyperbolic codes use a grid buffer of 1
  std::string prefix("plt");
  amr.gridBufferSize(1);
  int plotInterval;
  pp.get("plot_interval", plotInterval);
  amr.checkpointInterval(-1);
  amr.plotInterval(plotInterval);
  amr.maxDtGrow(max_dt_growth);
  amr.dtToleranceFactor(dt_tolerance_factor);
  amr.verbosity(verbosity);

  // initialize hierarchy of levels
  amr.setupForNewAMRRun();
  // run
  amr.run(stopTime,nstop);
  // output last pltfile and stapptistics
  amr.conclude();

  Vector<AMRLevel*> amrLevelData = amr.getAMRLevels();
 CH_assert(amrLevelData.size() == 1);

  EBAMRGodunov* amrGodunovData = dynamic_cast<EBAMRGodunov*>(amrLevelData[0]);
  LevelData<EBCellFAB>& amrData = amrGodunovData->getStateNew();

  a_grids = amrData.disjointBoxLayout();
  a_finalTime = amr.getCurrentTime();
  a_dt = amrGodunovData->getDt();
  a_ebisl = amrGodunovData->getEBISLayout();
  amrGodunovData->fillConsAndPrim(a_calcData);
}
/***************/
/***************/
void
generateErrorOverTime(LevelData<EBCellFAB>& a_error,
                      DisjointBoxLayout&    a_grids,
                      EBISLayout&           a_ebisl,
                      const Box&            a_domain,
                      const Real&           a_dx,
                      const int&            a_nstop)
{
  LevelData<EBCellFAB> calc;
  LevelData<EBCellFAB> exac;

  Real finalTime, dt;
  generateCalculatedData(calc, dt,
                         finalTime, a_grids, a_ebisl,
                         a_domain, a_dx, a_nstop);


  EBPatchPolytropic patchInt;
  bool is_periodic[SpaceDim];
  for (int i = 0; i < SpaceDim; ++i)
    {
      is_periodic[i] = false;
    }
  ProblemDomain probdomain(a_domain.smallEnd(), a_domain.bigEnd(), is_periodic);
  patchInt.define(probdomain, a_dx);

  int nCons = patchInt.numConserved();
  int nPrim = patchInt.numPrimitives();
  int consAndPrim = nCons+nPrim;
  generateExactData( exac, a_grids, a_ebisl,
                    a_domain, a_dx, finalTime, dt, consAndPrim);


  EBCellFactory ebcf(a_ebisl);

  a_error.define(a_grids, consAndPrim, IntVect::Zero, ebcf);

  Interval interv(0, consAndPrim-1);
  for(DataIterator dit = a_grids.dataIterator(); dit.ok(); ++dit)
    {
      EBCellFAB& errorFAB = a_error[dit()];
      const EBCellFAB& calcFAB = calc[dit()];
      const EBCellFAB& exacFAB = exac[dit()];
      //error is exact-calculated data
      errorFAB.copy(exacFAB);
      errorFAB -= calcFAB;

      //finally divide by the final time because that is
      //what we want
      if(finalTime > 0)
        {
          errorFAB /= finalTime;
        }
      else
        {
          MayDay::Error("final time = 0");
        }

      for(int icomp = 0; icomp < a_error[dit()].nComp(); icomp++)
        {
          a_error[dit()].setCoveredCellVal(0.0, icomp);
        }
    }
}

/***************/
/***************/
void
outputError(const LevelData<EBCellFAB>& a_errorFine,
            const EBISLayout& a_ebislFine,
            const DisjointBoxLayout& a_gridsFine,
            const Box& a_domainFine,
            const int& a_nstepmax)
{
#ifdef CH_USE_HDF5
  char charstr[100];
  sprintf(charstr, "pltError%d.%dd.hdf5", a_nstepmax, SpaceDim);
  string fileFine(charstr);

  EBPatchPolytropic patchInt;
  Vector<string> names  = getNames();
  bool rCov = true;
  int nCons = patchInt.numConserved();
  int nPrim = patchInt.numPrimitives();
  int consAndPrim = nCons+nPrim;
  Vector<Real> cVal(consAndPrim,0.0);
  //values that don't matter in output file
  Real dxFine = 1.0;
  Real time   = 1.0;
  Real dt     = 1.0;

  writeEBHDF5(fileFine, a_gridsFine, a_errorFine, names,
              a_domainFine, dxFine, dt, time,
              rCov, cVal);
#endif
}
/***************/
/***************/
int iscript(int icomp, int inorm, int ncomp)
{
  return icomp + inorm*ncomp;
}
/***************/
/***************/
void
makeFinestDomain(Box& a_domain,
                 Real& a_dx)
{
  //parse input file.  single level
  ParmParse pp;
  Vector<int> n_cell(SpaceDim);
  pp.getarr("n_cell",n_cell,0,SpaceDim);

  int testverbosity;
  pp.get("testverbosity", testverbosity);
 CH_assert(n_cell.size() == SpaceDim);
  IntVect lo = IntVect::Zero;
  IntVect hi;
  for(int ivec = 0; ivec < SpaceDim; ivec++)
    {
      if(n_cell[ivec] <= 0)
        {
          if(testverbosity > 1)
            pout() << " bogus number of cells input = " << n_cell[ivec];
          MayDay::Error();
        }
      hi[ivec] = n_cell[ivec] - 1;
    }

  a_domain = Box(lo, hi);

  Real prob_hi;
  int numOpen = n_cell[0];
  pp.get("domain_length",prob_hi);
  a_dx = prob_hi/numOpen;
}
/***/
void
makeGeometry(const Box& a_domain,
             const Real& a_dx)
{

  RealVect origin = RealVect::Zero;
  ParmParse pp;
  Real channelRadius;
  pp.get("channel_radius", channelRadius);

  EBIndexSpace* ebisPtr = Chombo_EBIS::instance();
  int biggridsize = 2048;
  char ebfilename[80];
  int  fineLength = a_domain.size(0);
  pp.get("max_grid_size", biggridsize);
  pout() << "using a tilted cylinder" << endl;
  vector<Real>  cylinderAxisVect(SpaceDim);
  pp.getarr("cylinder_axis",cylinderAxisVect, 0, SpaceDim);
  RealVect cylinderAxis;
  for(int idir = 0; idir < SpaceDim; idir++)
    {
      cylinderAxis[idir] = cylinderAxisVect[idir];
    }
  Real sum;
  PolyGeom::unifyVector(cylinderAxis, sum);

  RealVect corner = RealVect::Zero;
  bool negativeInside = true;
  TiltedCylinderIF tunnel(channelRadius, cylinderAxis, corner, negativeInside);

  RealVect vectDx = RealVect::Unit;
  vectDx *= a_dx;

  GeometryShop workshop(tunnel,0,vectDx);
  ebisPtr->define(a_domain, origin, a_dx, workshop, biggridsize);
  sprintf(ebfilename, "cyl%d.%dd.hdf5",fineLength, SpaceDim);
}
/***************/
/***************/
void
multiplyByKappa(LevelData<EBCellFAB>& a_error,
                const EBISLayout& a_ebisl,
                const DisjointBoxLayout& a_grids)
{
 for(DataIterator dit = a_grids.dataIterator(); dit.ok(); ++dit)
    {
      const Box& box = a_grids.get(dit());
      const EBISBox& ebisBox = a_ebisl[dit()];
      IntVectSet ivs(box);
      for(VoFIterator vofit(ivs, ebisBox.getEBGraph()); vofit.ok(); ++vofit)
        {
          const VolIndex& vof = vofit();
          Real volFrac = ebisBox.volFrac(vof);
          for(int ivar = 0; ivar < a_error.nComp(); ivar++)
            {
              a_error[dit()](vof, ivar) *= volFrac;
            }
        }
    }
}
/***************/
int
main(int a_argc, char* a_argv[])
{
#ifdef CH_MPI
  MPI_Init(&a_argc,&a_argv);
  {
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
    ParmParse pp(a_argc-2,a_argv+2,NULL,inFile);
    int nstop;
    pp.get("max_step",nstop);
    Box domainFine;
    Real dxFine;
    makeFinestDomain(domainFine, dxFine);

    LevelData<EBCellFAB>     errorFine;
    DisjointBoxLayout gridsFine;
    EBISLayout        ebislFine;
    int testverbosity;

    pp.get("testverbosity", testverbosity);
    if(testverbosity >= 1)
      pout() << "generating fine domain" << endl;

    makeGeometry(domainFine, dxFine);


    if(testverbosity >= 1)
      pout() << "starting to generate error" << endl;

    for(int istepmax = 1; istepmax < nstop; istepmax++)
      {
        if(testverbosity >= 1)
          pout() << "generating error/time for stepmax = " << istepmax << endl;

        generateErrorOverTime(errorFine, gridsFine, ebislFine, domainFine, dxFine, istepmax);

        if(testverbosity > 1)
          pout() << "Outputting error to file" << endl;

        outputError(errorFine, ebislFine, gridsFine, domainFine, istepmax);
      }

#ifdef CH_MPI
    dumpmemoryatexit();
  }
  MPI_Finalize();
#endif
}


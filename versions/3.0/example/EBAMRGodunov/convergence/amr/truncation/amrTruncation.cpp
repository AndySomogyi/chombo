#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

#include <cmath>
#include <cstdio>
#include <iostream>

#include "ParmParse.H"
#include "CH_HDF5.H"
#include "parstream.H"
#include "ModianoIBCFactory.H"
#include "EBPatchPolytropicFactory.H"
#include "EBPatchPolytropic.H"
#include "EBAMRGodunovFactory.H"
#include "EBAMRGodunov.H"
#include "AMRLevel.H"
#include "AMR.H"
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
#include "EBAMRDataOps.H"
#include "RedistStencil.H"
#include "EBAMRIO.H"
#include "DebugDump.H"
#include "EBDebugDump.H"
#include "TiltedCylinderIF.H"
#include "EBAMRGodunovFactory.H"
#include "EBAMRIO.H"
#include "GodunovGeom.H"

#ifdef CH_MPI
#include "CH_Attach.H"
#endif

/************/
void
generateData(Vector< LevelData<EBCellFAB>* >& a_datum,
             Vector< DisjointBoxLayout     >& a_grids,
             Vector< EBISLayout            >& a_ebisl,
             Real                           & a_dt,
             const Vector<Vector<Box>      >& a_boxes,
             const Vector<int              >& a_refRat,
             const Box                      & a_domain,
             const int                      & a_refToFinestCalc)
{
  //get IBC factory
  RefCountedPtr<ModianoIBCFactory> ibc;
  getModianoIBCFactory(ibc);

  //get patch integrator
  RefCountedPtr<EBPatchPolytropicFactory> patch;
  getEBPPFactoryXY(patch, &(*ibc));

  //get level integrator
  RefCountedPtr<EBAMRGodunovFactory>  amrg;
  getEBAMRGFactory(amrg, &(*patch));

  // read inputs
  ParmParse pp;
  int max_level = a_boxes.size() - 1;


  bool is_periodic[SpaceDim];
  for (int i = 0; i < SpaceDim; ++i)
    {
      is_periodic[i] = false;
    }
  ProblemDomain prob_domain(a_domain.smallEnd(),
                            a_domain.bigEnd(),
                            is_periodic);


  AMR amr;
  amr.define(max_level,a_refRat,prob_domain,&(*amrg));

  if(a_refToFinestCalc > 1)
    {
      amr.fixedDt(a_dt);
    }

  int max_grid_size = 32;
  pp.get("max_grid_size",max_grid_size);
  char prefixchar[100];
  sprintf(prefixchar, "pltnx%dSoln", a_domain.size(0));
  std::string prefix(prefixchar);
  Real max_dt_growth = 1.1;
  pp.get("max_dt_growth",max_dt_growth);
  Real dt_tolerance_factor = 1.1;
  pp.get("dt_tolerance_factor",dt_tolerance_factor);
  amr.maxGridSize(max_grid_size);
  amr.checkpointInterval(-1);
  amr.plotInterval(1);
  amr.maxDtGrow(max_dt_growth);
  amr.dtToleranceFactor(dt_tolerance_factor);
  int verbosity;
  pp.get("verbosity", verbosity);
  amr.verbosity(verbosity);

  // the hyperbolic codes use a grid buffer of 1
  amr.gridBufferSize(1);

  // initialize hierarchy of levels
  amr.setupForFixedHierarchyRun(a_boxes);
  // run
  Real stopTime = 10000000;
  int nstop  = 1;

  amr.run(stopTime,nstop);
  // output last pltfile and stapptistics
  amr.conclude();

  Vector<AMRLevel*> amrLevelData = amr.getAMRLevels();
  CH_assert(amrLevelData.size() == a_boxes.size());

  a_grids.resize(amrLevelData.size());
  a_ebisl.resize(amrLevelData.size());
  a_datum.resize(amrLevelData.size(), NULL);
  if(a_refToFinestCalc == 1)
    {
      a_dt = amrLevelData[0]->dt();
    }
  for(int ilev = 0; ilev < amrLevelData.size(); ilev++)
    {
      EBAMRGodunov* amrGodunovData =  (EBAMRGodunov*) (amrLevelData[ilev]);
      LevelData<EBCellFAB>& amrData = amrGodunovData->getStateNew();

      a_grids[ilev] = amrData.disjointBoxLayout();
      a_ebisl[ilev] = amrGodunovData->getEBISLayout();

      const LevelData<EBCellFAB>& stateNew = amrGodunovData->getStateNew();
      const LevelData<EBCellFAB>& stateOld = amrGodunovData->getStateOld();
      EBCellFactory fact(a_ebisl[ilev]);
      a_datum[ilev] = new LevelData<EBCellFAB>(a_grids[ilev], stateNew.nComp(), stateNew.ghostVect(), fact);
      CH_assert(a_dt > 0);
      for(DataIterator dit = a_grids[ilev].dataIterator(); dit.ok(); ++dit)
        {
          EBCellFAB& datFAB = (*a_datum[ilev])[dit()];
          //returning an estimate to du/dt
          datFAB.copy(stateNew[dit()]);
          datFAB -=   stateOld[dit()];
          datFAB /= a_dt;
        }
    }
}
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

/*****/
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
  RealVect cylinderOrigin = RealVect::Zero;

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

  bool negativeInside = true;
  TiltedCylinderIF tunnel(channelRadius, cylinderAxis, cylinderOrigin, negativeInside);

  RealVect vectDx = RealVect::Unit;
  vectDx *= a_dx;

  GeometryShop workshop(tunnel,0,vectDx);
  ebisPtr->define(a_domain, origin, a_dx, workshop, biggridsize);
}
/*****/
void
coarsenBoxes(Vector< Vector<Box>      >&    a_boxesCoar,
             const Vector<Vector<Box> >&    a_boxesFine,
             int a_refToCoar)
{
  a_boxesCoar.resize(a_boxesFine.size());
  for(int ilev = 0; ilev < a_boxesFine.size(); ilev++)
    {
      a_boxesCoar[ilev].resize(a_boxesFine[ilev].size());
      for(int ibox = 0; ibox < a_boxesFine[ilev].size(); ibox++)
        {
          a_boxesCoar[ilev][ibox] = coarsen(a_boxesFine[ilev][ibox], a_refToCoar);
        }
    }
}
/*****/
void
getBoxes(Vector<Vector<Box> >&   a_boxes,
         Vector<int>&            a_refRat,
         const Box&              a_domain)
{
  int amrRef = 2;
  a_refRat.resize(2, amrRef);
  a_boxes.resize(2);
  Box fineBox = refine(a_domain, amrRef);
  int ishrink = fineBox.size(0);
  //this leaves 1/4 refined.
  ishrink *= 3;
  ishrink /= 8;
  fineBox.grow(-ishrink);
  a_boxes[0] = Vector<Box>(1, a_domain);
  a_boxes[1] = Vector<Box>(1, fineBox);
}
/*****/
void
solutionErrorTest(int testverbosity, int fileout)
{
  Box domaiFine, domaiMedi, domaiCoar;
  Real dxFine, dxMedi, dxCoar;
  makeFinestDomain(domaiFine, dxFine);
  domaiMedi = coarsen(domaiFine, 2);
  domaiCoar = coarsen(domaiMedi, 2);
  dxMedi = 2.0*dxFine;
  dxCoar = 2.0*dxMedi;

  Vector< LevelData<EBCellFAB> *>      solutFine, solutMedi, solutCoar;
  Vector< LevelData<EBCellFAB> *>                 errorMedi, errorCoar;
  Vector<DisjointBoxLayout      >      gridsFine, gridsMedi, gridsCoar;
  Vector<EBISLayout             >      ebislFine, ebislMedi, ebislCoar;
  Vector<Vector<Box>            >      boxesFine, boxesMedi, boxesCoar;
  Vector<int>  refRat;

  getBoxes(    boxesFine, refRat, domaiFine);
  coarsenBoxes(boxesMedi, boxesFine, 2);
  coarsenBoxes(boxesCoar, boxesMedi, 2);

  if(testverbosity >= 1)
    pout() << "generating geometry" << endl;

  //because we are in amr land
  Box domainRFF = domaiFine;
  for(int ilev = 1; ilev < refRat.size(); ilev++)
    {
      domainRFF.refine(refRat[ilev-1]);
    }

  makeGeometry(domainRFF, dxFine);

  int refFine = 1;  int refMedi = 2;   int refCoar = 4;
  //everyone runs the same dt
  Real dt;
  if(testverbosity >= 1) pout() << "generating fine solution" << endl;
  generateData(solutFine, gridsFine, ebislFine, dt, boxesFine, refRat, domaiFine, refFine);

  if(testverbosity >= 1) pout() << "generating medi solution" << endl;
  generateData(solutMedi, gridsMedi, ebislMedi, dt, boxesMedi, refRat, domaiMedi, refMedi);

  if(testverbosity >= 1) pout() << "generating coar solution" << endl;
  generateData(solutCoar, gridsCoar, ebislCoar, dt, boxesCoar, refRat, domaiCoar, refCoar);


  if(testverbosity >= 1) pout() << "generating medi error from medi and fine solutions" << endl;
  EBAMRDataOps::getErrorFromCoarseAndFine(errorMedi,
                                          solutMedi, gridsMedi, ebislMedi, domaiMedi,
                                          solutFine, gridsFine, ebislFine, domaiFine,
                                          refRat);

  if(testverbosity >= 1) pout() << "generating coar error from medi and coar solutions" << endl;
  EBAMRDataOps::getErrorFromCoarseAndFine(errorCoar,
                                          solutCoar, gridsCoar, ebislCoar, domaiCoar,
                                          solutMedi, gridsMedi, ebislMedi, domaiMedi,
                                          refRat);

  Vector<Real> order;
  EBArith::compareError(order,
                        errorMedi, errorCoar,
                        gridsMedi, gridsCoar,
                        ebislMedi, ebislCoar,
                        refRat,    domaiCoar, testverbosity);

#ifdef CH_USE_HDF5
  if(fileout == 1)
    {
      if(testverbosity > 1)  pout() << "Outputting error to file" << endl;
      const char fileCoar[] = "plotCoarError.hdf5";
      const char fileMedi[] = "plotMediError.hdf5";
      writeEBAMRname(&errorCoar, fileCoar);
      writeEBAMRname(&errorMedi, fileMedi);
    }
#else
  pout() << "hdf5 not compiled in so I cannot output error to the file" << endl;
#endif

  for(int ilev = 0; ilev < refRat.size(); ilev++)
    {
      delete errorMedi[ilev];
      delete errorCoar[ilev];
      delete solutFine[ilev];
      delete solutMedi[ilev];
      delete solutCoar[ilev];
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
    int testverbosity;
    int fileout;
    pp.get("do_erroroutput", fileout);
    pp.get("testverbosity", testverbosity);

    solutionErrorTest(testverbosity, fileout);

#ifdef CH_MPI
    dumpmemoryatexit();
  }
  MPI_Finalize();
#endif
}


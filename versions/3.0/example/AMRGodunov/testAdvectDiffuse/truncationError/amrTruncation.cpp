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
#include "AMR.H"
#include "AMRLevelAdvectDiffuseFactory.H"
#include "AdvectTestIBC.H"
#include "AdvectionFunctions.H"
#include "DebugDump.H"
#include "memtrack.H"
#include "CH_Attach.H"
#include "FABView.H"
#include "AdvectDiffuseUtils.H"
#include "parstream.H"
/************/
void
generateData(Vector< LevelData<FArrayBox>* >& a_datum,
             Vector< DisjointBoxLayout     >& a_grids,
             Real                           & a_dt,
             const Vector<Vector<Box>      >& a_boxes,
             const Vector<int              >& a_refRat,
             const ProblemDomain            & a_domain,
             const int                      & a_refToFinestCalc)
{
  RefCountedPtr<AdvectTestIBC> ibc;
  getAdvectTestIBC(ibc);

  AdvectPhysics advPhys;
  advPhys.setPhysIBC(&(*ibc));

  AdvectionVelocityFunction velFunc;
  getAdvectionVelocityFunction(velFunc);

  RefCountedPtr<AMRLevelAdvectDiffuseFactory>  amrg_fact;
  getAMRLADFactory(amrg_fact, velFunc, advPhys);

  AMR amr;
  defineAMR(amr, amrg_fact, a_domain, a_refRat);

  //see if dt is fixed (not so if you are trying to figure out what the stable
  // time step is
  if(a_refToFinestCalc > 1)
    {
      amr.fixedDt(a_dt);
    }

  // initialize hierarchy of levels
  amr.setupForFixedHierarchyRun(a_boxes);

  // run fixed number of steps
  Real stopTime = 1.0e10;
  int nstop = 1;
  Vector<AMRLevel*> amrLevelData = amr.getAMRLevels();
  if(a_refToFinestCalc == 1)
    {
      a_dt = amrLevelData[0]->dt();
    }
  amr.run(stopTime,nstop);
  // output last pltfile and stapptistics
  amr.conclude();

  CH_assert(amrLevelData.size() == a_boxes.size());

  a_grids.resize(amrLevelData.size());
  a_datum.resize(amrLevelData.size(), NULL);
  for(int ilev = 0; ilev < amrLevelData.size(); ilev++)
    {
      AMRLevelAdvectDiffuse* amrGodunovData =  (AMRLevelAdvectDiffuse*) (amrLevelData[ilev]);
      LevelData<FArrayBox>& stateNew = amrGodunovData->getStateNew();
      LevelData<FArrayBox>& stateOld = amrGodunovData->getStateOld();
      a_grids[ilev] = stateNew.disjointBoxLayout();
      a_datum[ilev] = new LevelData<FArrayBox>(a_grids[ilev], stateNew.nComp(), IntVect::Zero);
      for(DataIterator dit = a_grids[ilev].dataIterator(); dit.ok(); ++dit)
        {
          FArrayBox& datFAB = (*a_datum[ilev])[dit()];
          //returning an estimate to du/dt
          datFAB.copy(stateNew[dit()]);
          datFAB -=   stateOld[dit()];
          datFAB /= a_dt;
        }
    }
}

/*****/
void
truncationErrorTest(int testverbosity, int fileout)
{
  ProblemDomain domaiFine, domaiMedi, domaiCoar;
  Real dxFine, dxMedi, dxCoar;
  makeFinestDomain(domaiFine, dxFine);
  domaiMedi = coarsen(domaiFine, 2);
  domaiCoar = coarsen(domaiMedi, 2);
  dxMedi = 2.0*dxFine;
  dxCoar = 2.0*dxMedi;

  Vector< LevelData<FArrayBox> *>      solutFine, solutMedi, solutCoar;
  Vector< LevelData<FArrayBox> *>                 errorMedi, errorCoar;
  Vector<DisjointBoxLayout      >      gridsFine, gridsMedi, gridsCoar;
  Vector<Vector<Box>            >      boxesFine, boxesMedi, boxesCoar;
  Vector<int>  refRat;

  getBoxes(    boxesFine, refRat, domaiFine.domainBox());
  coarsenBoxes(boxesMedi, boxesFine, 2);
  coarsenBoxes(boxesCoar, boxesMedi, 2);

  int refFine = 1;  int refMedi = 2;   int refCoar = 4;
  Real dt; //uses finest dt calculation
  if(testverbosity >= 1) pout() << "generating fine solution" << endl;
  generateData(solutFine, gridsFine, dt, boxesFine, refRat, domaiFine,  refFine);

  if(testverbosity >= 1) pout() << "generating medi solution" << endl;
  generateData(solutMedi, gridsMedi, dt, boxesMedi, refRat, domaiMedi,  refMedi);

  if(testverbosity >= 1) pout() << "generating coar solution" << endl;
  generateData(solutCoar, gridsCoar, dt, boxesCoar, refRat, domaiCoar,  refCoar);


  if(testverbosity >= 1) pout() << "generating medi error from medi and fine solutions" << endl;
  getErrorFromCoarseAndFine(errorMedi,
                            solutMedi, gridsMedi,  domaiMedi,
                            solutFine, gridsFine,  domaiFine,
                            refRat);

  if(testverbosity >= 1) pout() << "generating coar error from medi and coar solutions" << endl;
  getErrorFromCoarseAndFine(errorCoar,
                            solutCoar, gridsCoar, domaiCoar,
                            solutMedi, gridsMedi, domaiMedi, refRat);

  Vector<Real> order;
  compareError(order,
               errorMedi, errorCoar,
               gridsMedi, gridsCoar,
               refRat, domaiCoar, testverbosity);

#ifdef CH_USE_HDF5
  if(fileout == 1)
    {
      if(testverbosity > 1)  pout() << "Outputting error to file" << endl;
      const char fileCoar[] = "plotCoarError.hdf5";
      const char fileMedi[] = "plotMediError.hdf5";
      writeVectorLevelName(&errorCoar, &refRat, fileCoar);
      writeVectorLevelName(&errorMedi, &refRat, fileMedi);
    }
#else
  MayDay::Warning("hdf5 not compiled in so I cannot output error to the file");
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

    truncationErrorTest(testverbosity, fileout);

#ifdef CH_MPI
    dumpmemoryatexit();
  }
  MPI_Finalize();
#endif
}


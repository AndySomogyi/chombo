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
#include "EBBackwardEuler.H"
#include "AllRegularService.H"
#include "EBLevelRedist.H"
#include "EBAMRDataOps.H"
#include "RedistStencil.H"
#include "EBAMRIO.H"
#include "DebugDump.H"
#include "EBDebugDump.H"
#include "TiltedCylinderIF.H"
#include "EBAMRIO.H"
#include "PoissonUtilities.H"
#include "EBAMRPoissonOp.H"
#include "EBEllipticLoadBalance.H"

#ifdef CH_MPI
#include "CH_Attach.H"
#endif

/************/
void
coarsenEBBEB(Vector< Vector<Box>      >&    a_boxesCoar,
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

/************/
void
getSimpleEBBEB(Vector<Vector<Box> >&   a_boxes,
               Vector<int>&            a_refRat,
               const Box&              a_domain)
{
  ParmParse pp;
  int maxlev;
  pp.get("max_level", maxlev);
  int amrRef = 2;
  
  a_refRat.resize(maxlev+1, amrRef);
  a_boxes.resize( maxlev+1);
  Box domlev = a_domain;
  int ishrink = domlev.size(0);
  ishrink /= 2;
  a_boxes[0] = Vector<Box>(1, a_domain);
  for(int ilev = 1 ; ilev <= maxlev; ilev++)
    {
      Box fineBox = refine(domlev, amrRef);
      fineBox.grow(-ishrink);

      a_boxes[ilev] = Vector<Box>(1, fineBox);

      domlev.refine(amrRef);
      ishrink  += 8;
    }
}
/************/
void
defineEBBEG(Vector< DisjointBoxLayout     >& a_grids,
            Vector< EBISLayout            >& a_ebisl,
            const Vector<Vector<Box>      >& a_boxes,
            const Vector<int              >& a_refRat,
            const PoissonParameters        & a_params)
{
  a_grids.resize(a_boxes.size());
  a_ebisl.resize(a_boxes.size());
  ProblemDomain domLev = a_params.coarsestDomain;
  const EBIndexSpace* const  ebisPtr = Chombo_EBIS::instance();
  for(int ilev = 0; ilev < a_boxes.size(); ilev++)
    {
      Vector<int> procs;
      EBEllipticLoadBalance(procs, a_boxes[ilev], domLev);
      a_grids[ilev] = DisjointBoxLayout(a_boxes[ilev], procs);
      int numGhost = 4;
      ebisPtr->fillEBISLayout(a_ebisl[ilev],
                              a_grids[ilev],
                              domLev,
                              numGhost);
      domLev.refine(a_refRat[ilev]);
    }
}
/************/
void
defineEBBackS(AMRMultiGrid<LevelData<EBCellFAB> >&         a_solver,
             const Vector<DisjointBoxLayout>&             a_grids,
             const Vector<EBISLayout>&                    a_ebisl,
             LinearSolver<LevelData<EBCellFAB > >&        a_bottomSolver,
             const PoissonParameters&                     a_params,
             const Real&                                  a_time,
             Real                                         a_alpha,
             Real                                         a_beta,
             RefCountedPtr<AMRLevelOpFactory<LevelData<EBCellFAB> > >& operatorFactory)
{
  CH_TIME("PoissonUtilities::defineSolver");
  ParmParse pp;
  int numPreCondIters = 40;
  pp.get("num_pre_cond_iters",numPreCondIters);
  int relaxtype;
  pp.get("mg_relax_type", relaxtype);
  if(relaxtype == 0)
    {
      pout() << "Using levelJacobi" << endl;
    }
  else if(relaxtype == 1)
    {
      pout() << "Using levelMultiColorGS" << endl;
    }
  else if(relaxtype == 2)
    {
      pout() << "Using levelGSRB" << endl;
    }
  else if(relaxtype == 3)
    {
      pout() << "Using slow multicolor relax" << endl;
    }
  else if(relaxtype == 4)
    {
      pout() << "Using levelMultiColorGS with clone" << endl;
    }
  else
    {
      MayDay::Error("bogus relax type in input file");
    }

  int numSmooths;
  pp.get("mg_num_smooths", numSmooths);
  Real hang, eps;
  pp.get("mg_hang", hang);
  pp.get("mg_eps",  eps);
  int numMG, iterMax;
  pp.get("mg_num_cycles", numMG);

  pp.get("mg_iter_max", iterMax);

  getEBAMRPFactory(operatorFactory,
                   a_grids,  a_ebisl,
                   a_params, numPreCondIters, relaxtype,
                   a_time,a_alpha,  a_beta);

  ProblemDomain coarsestDomain(a_params.coarsestDomain);
  a_solver.define(coarsestDomain,
                  *operatorFactory,
                  &a_bottomSolver,
                  a_params.numLevels);
  Real normThresh = 1.0e-30;
  a_solver.setSolverParameters(numSmooths, numSmooths, numSmooths,
                               numMG, iterMax, eps, hang, normThresh);
  a_solver.m_verbosity = 4;
  if(pp.contains("multigrid_verbosity"))
    {
      pp.get("multigrid_verbosity", a_solver.m_verbosity);
    }
}
/************/
void
newAndFillSmoothCellData(Vector<LevelData<EBCellFAB>* >&        a_scalar,
                         const Vector<DisjointBoxLayout>&       a_grids,
                         const Vector<EBISLayout>&              a_ebisl,
                         const PoissonParameters&               a_params,
                         int a_ncomp,
                         bool a_usePhiGhost)
{
  int nlevels = a_grids.size();
  a_scalar.resize(nlevels, NULL);

  RealVect dxLev = RealVect::Unit;
  dxLev *=a_params.coarsestDx;
  IntVect ivghost = a_params.ghostPhi;
  if(!a_usePhiGhost)
    ivghost = a_params.ghostRHS;
  ProblemDomain domLev = a_params.coarsestDomain;
  for(int ilev = 0; ilev < nlevels; ilev++)
    {
      EBCellFactory ebcellfact(a_ebisl[ilev]);
      a_scalar[ilev] = new LevelData<EBCellFAB>(a_grids[ilev], a_ncomp, ivghost, ebcellfact);

      EBLevelDataOps::setToZero(*(a_scalar[ilev]));
      for(int icomp= 0; icomp < a_ncomp; icomp++)
        {
          setTrigPhi(*a_scalar[ilev], dxLev, a_params, 0.0, icomp);
        }

      //prepare dx, domain for next level
      dxLev /=      a_params.refRatio[ilev];
      domLev.refine(a_params.refRatio[ilev]);
    }
}
/************/
void
generateData(Vector< LevelData<EBCellFAB>* >& a_datum,
             Vector< DisjointBoxLayout     >& a_grids,
             Vector< EBISLayout            >& a_ebisl,
             const Vector<Vector<Box>      >& a_boxes,
             const PoissonParameters        & a_params,
             const int                      & a_refToFinestCalc)
{
  defineEBBEG(a_grids, a_ebisl, a_boxes, a_params.refRatio, a_params);
  pout() << "filling data" << endl;
  Vector<LevelData<EBCellFAB>* > phiOld, source;
  Vector<LevelData<EBCellFAB>* >& phiNew = a_datum;
  newAndFillSmoothCellData(phiOld, a_grids, a_ebisl, a_params, 1, true);
  newAndFillSmoothCellData(phiNew, a_grids, a_ebisl, a_params, 1, true);
  newAndFillSmoothCellData(source, a_grids, a_ebisl, a_params, 1, true);


  pout() << "defining  backward Euler  solver" << endl;
  RefCountedPtr< AMRMultiGrid<LevelData<EBCellFAB> > >  
    solver = RefCountedPtr< AMRMultiGrid<LevelData<EBCellFAB> > > (new AMRMultiGrid<LevelData<EBCellFAB> >);

  Real beta;
  ParmParse pp;
  pp.get("diffusion_coeff", beta);

  RefCountedPtr<AMRLevelOpFactory<LevelData<EBCellFAB> > > factory;
  Real alpha = 1;  
  BiCGStabSolver<LevelData<EBCellFAB> > bottomSolver;
  bottomSolver.m_verbosity = 0;
  defineEBBackS(*solver, a_grids, a_ebisl, bottomSolver, a_params, 0.0,alpha, beta, factory);
  EBBackwardEuler               integrator(solver, *factory, a_params.coarsestDomain, a_params.refRatio, -1, 3);

  Real dt = a_params.coarsestDx[0];
  int nstop;
  pp.get("max_steps", nstop);
  nstop /= a_refToFinestCalc;

  Real time = 0;
  int maxlev = a_datum.size() - 1;

  solver->init(phiOld, source, maxlev, 0);
  for(int istep = 0; istep < nstop; istep++)
    {
      integrator.oneStep(phiNew, phiOld, source, dt, 0, maxlev);
      time += dt;
      EBAMRDataOps::assign(phiOld,phiNew);
      pout() << "step = " << istep << ", time = "  << time << endl;
    }

  //phinew = datum lives on.
  for(int ilev = 0; ilev <= maxlev; ilev++)
    {
      delete phiOld[ilev];
      delete source[ilev];
    }
}

/*****/
void
solutionErrorTest(int testverbosity, int fileout)
{
  
  if(testverbosity >= 1)
    pout() << "getting parameters" << endl;
  PoissonParameters         paramFine, paramMedi, paramCoar;
  getPoissonParameters(paramFine);
  paramCoar = paramFine;
  paramMedi = paramFine;
  paramMedi.coarsen(2);
  paramCoar.coarsen(4);
  
  if(testverbosity >= 1)
    pout() << "generating geometry" << endl;

  definePoissonGeometry(paramFine);

  Box domaiFine = paramFine.coarsestDomain.domainBox();
  Box domaiMedi = paramMedi.coarsestDomain.domainBox();
  Box domaiCoar = paramCoar.coarsestDomain.domainBox();

  if(testverbosity >= 1)
    pout() << "getting grids" << endl;

  Vector<DisjointBoxLayout> gridsFine, gridsMedi, gridsCoar;
  Vector<EBISLayout>        ebislFine, ebislMedi, ebislCoar;
  Vector<Vector<Box>            >      boxesFine, boxesMedi, boxesCoar;

  getSimpleEBBEB(boxesFine, paramFine.refRatio, domaiFine);
  coarsenEBBEB(boxesMedi, boxesFine, 2);
  coarsenEBBEB(boxesCoar, boxesMedi, 2);

  Vector< LevelData<EBCellFAB> *>      solutFine, solutMedi, solutCoar;
  Vector< LevelData<EBCellFAB> *>                 errorMedi, errorCoar;



  int refFine = 1;  int refMedi = 2;   int refCoar = 4;

  if(testverbosity >= 1) pout() << "generating fine solution" << endl;
  generateData(solutFine, gridsFine, ebislFine, boxesFine, paramFine, refFine);

  if(testverbosity >= 1) pout() << "generating medi solution" << endl;
  generateData(solutMedi, gridsMedi, ebislMedi, boxesMedi, paramMedi, refMedi);

  if(testverbosity >= 1) pout() << "generating coar solution" << endl;
  generateData(solutCoar, gridsCoar, ebislCoar, boxesCoar,  paramCoar, refCoar);


  if(testverbosity >= 1) pout() << "generating medi error from medi and fine solutions" << endl;
  EBAMRDataOps::getErrorFromCoarseAndFine(errorMedi,
                                          solutMedi, gridsMedi, ebislMedi, domaiMedi,
                                          solutFine, gridsFine, ebislFine, domaiFine,
                                          paramFine.refRatio);

  if(testverbosity >= 1) pout() << "generating coar error from medi and coar solutions" << endl;
  EBAMRDataOps::getErrorFromCoarseAndFine(errorCoar,
                                          solutCoar, gridsCoar, ebislCoar, domaiCoar,
                                          solutMedi, gridsMedi, ebislMedi, domaiMedi,
                                          paramMedi.refRatio);

  if(testverbosity >= 1)
    pout() << "comparing error " << endl;
  Vector<Real> order;
  EBArith::compareError(order,
                        errorMedi, errorCoar,
                        gridsMedi, gridsCoar,
                        ebislMedi, ebislCoar,
                        paramMedi.refRatio,    domaiCoar, testverbosity);

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

  for(int ilev = 0; ilev < paramFine.refRatio.size(); ilev++)
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
    pp.get("do_error_output", fileout);
    pp.get("testverbosity", testverbosity);

    solutionErrorTest(testverbosity, fileout);

#ifdef CH_MPI
    dumpmemoryatexit();
  }
  MPI_Finalize();
#endif
}


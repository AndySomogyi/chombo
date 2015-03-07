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
#include "LoadBalance.H"
#include "BRMeshRefine.H"

#include "EBIndexSpace.H"
#include "EBISLayout.H"
#include "EBCellFactory.H"

#include "AMRINSUtils.H"
#include "EBAMRPoissonOp.H"

#include "EBFABView.H"
#include "EBDebugDump.H"
#include "DebugDump.H"
#include "EBAMRNoSubcycle.H"
#include "EBAMRDataOps.H"
#include "EBCoarseAverage.H"
#include "NoFlowVortex.H"
#include "UsingNamespace.H"

using std::cerr;

#define CHECK_FLOATING_PT 0

#if CHECK_FLOATING_PT==1
#include <fenv.h>
#endif

/******/
void getError(Vector< LevelData<EBCellFAB>* >&           a_errorCoar,
              const Vector< LevelData<EBCellFAB>* >&     a_solnCoar,
              const Vector< DisjointBoxLayout >&         a_gridsCoar,
              const Vector< EBISLayout >&                a_ebislCoar,
              const ProblemDomain&                       a_level0DomainCoar,
              const Vector< LevelData<EBCellFAB>* >&     a_solnFine,
              const Vector< DisjointBoxLayout >&         a_gridsFine,
              const Vector< EBISLayout >&                a_ebislFine,
              const ProblemDomain&                       a_level0DomainFine,
              Vector<string>&                            a_names,
              const AMRParameters&                       a_params)
{

  EBAMRDataOps::getErrorFromCoarseAndFine(a_errorCoar,
                                          a_solnCoar,
                                          a_gridsCoar,
                                          a_ebislCoar,
                                          a_level0DomainCoar,
                                          a_solnFine,
                                          a_gridsFine,
                                          a_ebislFine,
                                          a_level0DomainFine,
                                          a_params.m_refRatio);

}
/****/
void getSolution(Vector< LevelData<EBCellFAB>* >&           a_soln,
                 Vector< DisjointBoxLayout >&               a_grids,
                 Vector< EBISLayout >&                      a_ebisl,
                 const ProblemDomain&                       a_level0Domain,
                 Vector<string>&                            a_names,
                 const AMRParameters&                       a_params,
                 const Real&                                a_fixedDt,
                 int                                        a_gphiIterations,
                 int                                        a_maxStep,
                 const Vector<Vector<Box> > &               a_vectBoxes)
{
  ParmParse pp;
  Vector<Real> centerVect(SpaceDim);
  IntVect hiSide;
  Vector<int> nCells;
  pp.getarr("n_cell",  nCells,0,SpaceDim);
  pp.getarr("vortex_center",centerVect,0,SpaceDim);
  RealVect center;
  for(int idir = 0; idir < SpaceDim; idir++)
    {
      center[idir] = centerVect[idir];
      hiSide[idir] = nCells[idir] -1;
    }

  Real coreRadius, coreStrength;
  pp.get("vortex_strength",coreStrength);
  pp.get("core_radius",coreRadius);

  Real viscosity = 0.0;
  pp.get("viscosity", viscosity);

  NoFlowVortexFactory ibc(center, coreRadius, coreStrength);

  //time controlled by maxstep, fixedDt
  Real stopTime = 1.0e10;

  AMRParameters params =  a_params;
  params.m_initIterations = a_gphiIterations;

  EBAMRNoSubcycle amr(params, ibc, a_level0Domain, viscosity);

  amr.setupForFixedHierarchyRun(a_vectBoxes);

  if(a_fixedDt > 0.0)
    {
      amr.useFixedDt(a_fixedDt);
    }
  amr.run(stopTime, a_maxStep);


  Vector<string> presnames(1, string("pressure"));
  Vector<string> velonames(SpaceDim);
  Vector<string> gphinames(SpaceDim);
  int nlevels = a_params.m_maxLevel + 1;
  for(int idir = 0; idir < SpaceDim; idir++)
    {
      char velochar[100];
      char gphichar[100];
      sprintf(velochar, "velocity%d", idir);
      sprintf(gphichar, "gradPres%d", idir);
      velonames[idir] = string(velochar);
      gphinames[idir] = string(gphichar);
    }

  a_names = velonames;
  a_names.append(gphinames);
  a_names.append(presnames);

  a_grids = amr.getGrids();
  a_ebisl = amr.getEBISLayouts();
  a_soln.resize(nlevels);
  for(int ilev = 0; ilev < nlevels; ilev++)
    {
      EBCellFactory ebcellfact(a_ebisl[ilev]);
      a_soln[ilev] = new LevelData<EBCellFAB>(a_grids[ilev], 2*SpaceDim + 1,  IntVect::Zero, ebcellfact);
    }

  const Vector<LevelData<EBCellFAB>* >& velo = amr.getVeloNew();
  const Vector<LevelData<EBCellFAB>* >& gphi = amr.getGphiNew();
  const Vector<LevelData<EBCellFAB>* >& pres = amr.getPresNew();
  for(int ilev = 0; ilev < nlevels; ilev++)
    {
      Interval srcInterv, dstInterv;
      srcInterv = Interval(0, SpaceDim-1);
      dstInterv = Interval(0, SpaceDim-1);
      velo[ilev]->copyTo(srcInterv, *a_soln[ilev], dstInterv);

      srcInterv = Interval(0, SpaceDim-1);
      dstInterv = Interval(SpaceDim, 2*SpaceDim-1);
      gphi[ilev]->copyTo(srcInterv, *a_soln[ilev], dstInterv);

      srcInterv = Interval(0, 0);
      dstInterv = Interval(2*SpaceDim , 2*SpaceDim );
      pres[ilev]->copyTo(srcInterv, *a_soln[ilev], dstInterv);
    }
}
/***************/
int main(int argc, char* argv[])
{
#ifdef CH_MPI
  MPI_Init(&argc, &argv);
#endif
  // Scoping trick
  {
#ifdef CH_MPI
    MPI_Barrier(Chombo_MPI::comm);
#endif

#if CHECK_FLOATING_PT==1
    //    int except =  FE_DIVBYZERO | FE_UNDERFLOW | FE_OVERFLOW |  FE_INVALID ;
    int except =  FE_DIVBYZERO | FE_OVERFLOW |  FE_INVALID ;
    feenableexcept(except);
#endif

    if (argc < 2)
      {
        cerr << " usage " << argv[0] << " <input_file_name> " << endl;
        exit(0);
      }

    char* inFile = argv[1];
    ParmParse pp(argc-2,argv+2,NULL,inFile);

    ProblemDomain domainCoar, domainFine, domainMedi;
    AMRParameters params;
    getAMRINSParameters(params, domainCoar);
    int numFiltrations;
    pp.get("num_filtrations", numFiltrations);
    params.m_numFilterIterations = numFiltrations;
    params.m_subtractOffMean = true;

    domainMedi = refine(domainCoar, 2);
    domainFine = refine(domainMedi, 2);

    Vector<DisjointBoxLayout> gridsFine, gridsCoar, gridsMedi;
    Vector<EBISLayout>        ebislFine, ebislCoar, ebislMedi;

    int nlevels = params.m_maxLevel + 1;
    Vector<LevelData<EBCellFAB>* > solnFine( nlevels, NULL);
    Vector<LevelData<EBCellFAB>* > solnMedi( nlevels, NULL);
    Vector<LevelData<EBCellFAB>* > solnCoar( nlevels, NULL);
    Vector<LevelData<EBCellFAB>* > errorMedi(nlevels, NULL);
    Vector<LevelData<EBCellFAB>* > errorCoar(nlevels, NULL);


    pout() << "generating geometry only on finest level" << endl;
    AMRINSGeometry(params, domainFine);

    Vector<Vector<Box> > fineBoxes, mediBoxes, coarBoxes;
    getFixedGrids(fineBoxes, mediBoxes, coarBoxes, params, domainCoar);

    Real fixedDt;
    pp.get("fixed_dt", fixedDt);
    CH_assert(fixedDt > 0.0);
    int maxStep;
    pp.get("max_step", maxStep);



    pout() << "generating fine solution " << endl;
    Vector<string>   names;
    int gphiIter = 0;
    pp.get("init_iterations", gphiIter);
    getSolution(solnFine, gridsFine, ebislFine, domainFine, names, params, fixedDt, gphiIter, maxStep, fineBoxes);

    pout() << "generating medi solution and error" << endl;
    //AMRINSGeometry(params, domainMedi);
    //    gphiIter = 0;
    fixedDt *= 2.;
    maxStep /= 2;
    getSolution(solnMedi, gridsMedi, ebislMedi, domainMedi, names, params, fixedDt, gphiIter, maxStep, mediBoxes);

    getError(errorMedi,
             solnMedi, gridsMedi, ebislMedi, domainMedi,
             solnFine, gridsFine, ebislFine, domainFine,
             names, params);

    pout() << "generating coar solution and error" << endl;
    //AMRINSGeometry(params, domainCoar);
    //    gphiIter = 0;
    fixedDt *= 2;
    maxStep /= 2;
    getSolution(solnCoar, gridsCoar, ebislCoar, domainCoar, names, params, fixedDt, gphiIter, maxStep, coarBoxes);

    getError(errorCoar,
             solnCoar, gridsCoar, ebislCoar, domainCoar,
             solnMedi, gridsMedi, ebislMedi, domainMedi,
             names, params);

    string testName("Solution error");
    compareError(errorMedi,   errorCoar,
                 gridsMedi,   gridsCoar,
                 ebislMedi,   ebislCoar,
                 domainMedi,  domainCoar,
                 names, testName, params);

#if CH_SPACEDIM==2
    string fileFine("pltFineError.2d.hdf5");
    string fileCoar("pltCoarError.2d.hdf5");
#else
    string fileFine("pltFineError.3d.hdf5");
    string fileCoar("pltCoarError.3d.hdf5");
#endif
    int dofileout;
    pp.get("do_error_output", dofileout);
    if(dofileout == 1)
      {
        outputError(errorMedi,   errorCoar,
                    gridsMedi,   gridsCoar,
                    domainMedi,  domainCoar,
                    fileFine,     fileCoar,
                    names, params);
      }

    for(int ilev = 0; ilev < nlevels; ilev++)
      {
        delete solnFine[ilev];
        delete solnMedi[ilev];
        delete solnCoar[ilev];
        solnFine[ilev] = NULL;
        solnMedi[ilev] = NULL;
        solnCoar[ilev] = NULL;

        delete errorMedi[ilev];
        delete errorCoar[ilev];
        errorMedi[ilev] = NULL;
        errorCoar[ilev] = NULL;
      }

  }// End scoping trick

#ifdef CH_MPI
  MPI_Finalize();
#endif
}


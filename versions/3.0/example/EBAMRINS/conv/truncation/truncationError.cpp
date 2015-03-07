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
using std::cerr;

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
#include "EBCoarseAverage.H"
#include "InflowOutflowIBC.H"
#include "UsingNamespace.H"


#define CHECK_FLOATING_PT 0

#if CHECK_FLOATING_PT==1
#include <fenv.h>
#endif

/****/
void getSoln(Vector< LevelData<EBCellFAB>* >&           a_soln,
             Real&                                      a_dt,
             Vector< DisjointBoxLayout >&               a_grids,
             Vector< EBISLayout >&                      a_ebisl,
             const ProblemDomain&                       a_level0Domain,
             Vector<string>&                            a_names,
             const AMRParameters&                       a_params,
             int                                        a_initIterations,
             bool                                       a_atFinestLevel,
             const Vector<Vector<Box> > &               a_vectBoxes)
{
  int flowDir;
  ParmParse pp;
  pp.get("flow_dir", flowDir);

  Real inflowVel;
  pp.get("inflow_vel",  inflowVel);

  Real viscosity = 0.0;
  pp.get("viscosity", viscosity);

  int idoSlipWalls;
  pp.get("do_slip_walls", idoSlipWalls);
  bool doSlip = (idoSlipWalls == 1);
  pout() << "uniform inflow condition" << endl;
  int orderEBBC=2;
  pp.query("order_ebbc", orderEBBC);
  InflowOutflowIBCFactory ibc(flowDir, inflowVel, orderEBBC, doSlip);

  int maxStep = 1;

  Real stopTime = 1.0e10;


  AMRParameters params =  a_params;
  params.m_initIterations = a_initIterations;

  EBAMRNoSubcycle amr(params, ibc, a_level0Domain, viscosity);

  amr.setupForFixedHierarchyRun(a_vectBoxes);
  ///klugy mechanism to make sure everyone is using the same time step
  if(!a_atFinestLevel)
    {
      amr.useFixedDt(a_dt);
    }

  amr.run(stopTime, maxStep);


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
      a_soln[ilev]  = new LevelData<EBCellFAB>(a_grids[ilev], 2*SpaceDim + 1,  IntVect::Zero, ebcellfact);
    }

  const Vector<LevelData<EBCellFAB>* >& veloNew = amr.getVeloNew();
  const Vector<LevelData<EBCellFAB>* >& gphiNew = amr.getGphiNew();
  const Vector<LevelData<EBCellFAB>* >& presNew = amr.getPresNew();

  ///klugy mechanism to make sure everyone is using the same time step
  if(a_atFinestLevel)
    {
      a_dt = amr.getDt();
    }
  for(int ilev = 0; ilev < nlevels; ilev++)
    {
      Interval srcInterv, dstInterv;

      srcInterv = Interval(0, SpaceDim-1);
      dstInterv = Interval(0, SpaceDim-1);
      veloNew[ilev]->copyTo(srcInterv, *a_soln[ilev], dstInterv);

      srcInterv = Interval(0, SpaceDim-1);
      dstInterv = Interval(SpaceDim, 2*SpaceDim-1);
      gphiNew[ilev]->copyTo(srcInterv, *a_soln[ilev], dstInterv);

      srcInterv = Interval(0, 0);
      dstInterv = Interval(2*SpaceDim, 2*SpaceDim);
      presNew[ilev]->copyTo(srcInterv, *a_soln[ilev], dstInterv);
    }
}
/***************/
void getError(Vector< LevelData<EBCellFAB>* >& a_errorCoar,
              Vector< DisjointBoxLayout >&     a_gridsCoar,
              Vector< EBISLayout >&            a_ebislCoar,
              Vector<string>&                  a_names,
              const Vector<Vector<Box> > &     a_vectBoxesCoar,
              const ProblemDomain&             a_domainCoar,
              const Vector<Vector<Box> > &     a_vectBoxesFine,
              const ProblemDomain&             a_domainFine,
              const AMRParameters&             a_params,
              Real& a_dt,
              bool a_atFinestLevel)
{
  //this get generated and thrown away
  Vector< DisjointBoxLayout >  gridsFine;
  Vector< EBISLayout >         ebislFine;

  int nlevels = a_params.m_maxLevel + 1;
  a_errorCoar.resize(nlevels);
  int nref = 2;  //nothing to do with param refinement ratio. this is the refinement between the two solutions

  Vector<LevelData<EBCellFAB>* > solnFine( nlevels, NULL);
  Vector<LevelData<EBCellFAB>* > solnCoar( nlevels, NULL);

  pout() << "generating fine soln " << endl;
  int initIter = 0;
  getSoln(solnFine, a_dt, gridsFine, ebislFine, a_domainFine, a_names, a_params,  initIter, a_atFinestLevel, a_vectBoxesFine);

  pout() << "generating coarse soln " << endl;
  initIter = 0;
  getSoln(solnCoar, a_dt, a_gridsCoar, a_ebislCoar, a_domainCoar, a_names, a_params,  initIter, false, a_vectBoxesCoar);

  pout() << "generating truncation error = (Ave(fine) - coar)/dt" << endl;

  ProblemDomain domLevCoar = a_domainCoar;
  for(int ilev = 0; ilev < nlevels; ilev++)
    {
      EBCellFactory ebcellfact(a_ebislCoar[ilev]);
      int nvar = solnCoar[ilev]->nComp();
      Interval interv(0, nvar-1);

      a_errorCoar[ilev] = new LevelData<EBCellFAB>(a_gridsCoar[ilev], nvar,  IntVect::Zero, ebcellfact);

      EBCoarseAverage averageOp(gridsFine[ilev], a_gridsCoar[ilev], ebislFine[ilev], a_ebislCoar[ilev], domLevCoar, nref, nvar, Chombo_EBIS::instance());
      //here make error = Ave(fine)
      averageOp.average(*a_errorCoar[ilev], *solnFine[ilev], interv);
      //now subtract off coarse so error= Ave(Fine) - coar
      for(DataIterator dit = a_gridsCoar[ilev].dataIterator(); dit.ok(); ++dit)
        {
          EBCellFAB&       errorFAB = (*a_errorCoar[ilev])[dit()];
          const EBCellFAB& solnFAB  = (   *solnCoar[ilev])[dit()];

          errorFAB -= solnFAB;

          errorFAB /= a_dt;
        }
      domLevCoar.refine(a_params.m_refRatio[ilev]);
    }

    for(int ilev = 0; ilev < nlevels; ilev++)
      {
        delete solnFine[ilev];
        delete solnCoar[ilev];
        solnFine[ilev] = NULL;
        solnCoar[ilev] = NULL;
      }
}
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
    pp.get("num_filter_iterations", numFiltrations);
    params.m_numFilterIterations  = numFiltrations;
    bool doStokesFlow;
    pp.get("do_stokes_flow", doStokesFlow);
    params.m_stokesFlow = doStokesFlow;

    //need this so we can say du/dt = new-old/dt
    params.m_copyOverOld = false;
    domainMedi = refine(domainCoar, 2);
    domainFine = refine(domainMedi, 2);

    Vector<DisjointBoxLayout> gridsCoar, gridsMedi;
    Vector<EBISLayout>        ebislCoar, ebislMedi;

    int nlevels = params.m_maxLevel + 1;
    Vector<LevelData<EBCellFAB>* > errorMedi(nlevels, NULL);
    Vector<LevelData<EBCellFAB>* > errorCoar(nlevels, NULL);


    pout() << "generating geometry only on finest level" << endl;
    AMRINSGeometry(params, domainFine);

    Vector<Vector<Box> > boxesFine, boxesMedi, boxesCoar;
    getFixedGrids(boxesFine, boxesMedi, boxesCoar, params, domainCoar);

    Vector<string> names;
    Real dt = 1.2345678e9;
    getError(errorMedi, gridsMedi, ebislMedi, names,
             boxesMedi, domainMedi,
             boxesFine, domainFine,
             params, dt, true);


    getError(errorCoar, gridsCoar, ebislCoar, names,
             boxesCoar, domainCoar,
             boxesMedi, domainMedi,
             params, dt, false);

    string testName("Truncation error");
    compareError(errorMedi,   errorCoar,
                 gridsMedi,   gridsCoar,
                 ebislMedi,   ebislCoar,
                 domainMedi,  domainCoar,
                 names, testName, params);

#if CH_SPACEDIM==2
    string fileFine("plotFineError.2d.hdf5");
    string fileCoar("plotCoarError.2d.hdf5");
#else
    string fileFine("plotFineError.3d.hdf5");
    string fileCoar("plotCoarError.3d.hdf5");
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


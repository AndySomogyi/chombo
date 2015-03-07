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
#include "AMRPoissonOp.H"
#include "MACProjBC.H"
#include "FABView.H"
#include "DebugDump.H"
#include "PoissonUtilities.H"
#include "functionsF_F.H"
#include "PoissProbF_F.H"
#include "BiCGStabSolver.H"
#include "LevelMACProjector.H"
#include "LevelCCProjector.H"


void setExactVel(LevelData<FArrayBox>&      a_vel,
                 const Real&                a_dx,
                 const PoissonParameters&   a_params)
{
  CH_assert(a_vel.nComp() == SpaceDim);

  for (DataIterator dit = a_vel.dataIterator(); dit.ok(); ++dit)
    {
          FArrayBox& curPhiFAB = a_vel[dit()];
          Box curPhiBox = curPhiFAB.box();

          RealVect probLo = RealVect::Zero;
          RealVect probHi = RealVect::Unit; //never used.
          const RealVect&     trig = getTrigRV();

          for(int idir = 0; idir < SpaceDim; idir++)
            {
              RealVect dxVect  = a_dx*RealVect::Unit;
              FORT_GETPHI(CHF_FRA1(curPhiFAB,idir),
                          CHF_CONST_REALVECT(trig),
                          CHF_CONST_REALVECT(dxVect),
                          CHF_CONST_REALVECT(probLo),
                          CHF_CONST_REALVECT(probHi),
                          CHF_BOX(curPhiBox));
            }
    }
}
/******/
void getError(Vector< LevelData<FArrayBox>* >&     a_error,
              const Vector< DisjointBoxLayout >&   a_grids,
              const PoissonParameters&             a_params)
{
  ParmParse pp;

  int nlevels = a_params.numLevels;
  a_error.resize(nlevels);
  Vector<LevelData<FArrayBox>* >  velocity(nlevels, NULL);
  Vector<LevelData<FArrayBox>* >  gradient(nlevels, NULL);
  Vector<LevelData<FArrayBox>* >  pressure(nlevels, NULL);

  RefCountedPtr< AMRMultiGrid<LevelData<FArrayBox> > >  solver =
    RefCountedPtr< AMRMultiGrid<LevelData<FArrayBox> > >( new AMRMultiGrid<LevelData<FArrayBox> >());
  BiCGStabSolver<LevelData<FArrayBox> >   bottomSolver;
  bottomSolver.m_verbosity = 0;
  defineSolver(*solver, a_grids, bottomSolver, a_params);
  {
    Real dxLev = a_params.coarsestDx;
    ProblemDomain domLev(a_params.coarsestDomain);
    for(int ilev = 0; ilev < nlevels; ilev++)
      {
        velocity[ilev] = new LevelData<FArrayBox>(a_grids[ilev], SpaceDim, IntVect::Unit);
        gradient[ilev] = new LevelData<FArrayBox>(a_grids[ilev], SpaceDim, IntVect::Zero);
        a_error[ilev]  = new LevelData<FArrayBox>(a_grids[ilev],        1, IntVect::Zero);
        pressure[ilev] = new LevelData<FArrayBox>(a_grids[ilev],        1, IntVect::Unit);
        for(DataIterator dit = a_grids[ilev].dataIterator(); dit.ok(); ++dit)
          {
            (*a_error[ilev])[dit()].setVal(0.);
          }

        setExactVel(         *velocity[ilev], dxLev, a_params);

        //prepare dx, domain for next level
        dxLev /=      a_params.refRatio[ilev];
        domLev.refine(a_params.refRatio[ilev]);
      }
  }
  //create error = divergence of the velocity field after projection.
  {
    Real dxLev = a_params.coarsestDx;
    ProblemDomain domLev(a_params.coarsestDomain);
    for(int ilev = 0; ilev < nlevels; ilev++)
      {
        RefCountedPtr<LevelMACProjector> macProj = RefCountedPtr<LevelMACProjector>
          (new LevelMACProjector(a_grids[ilev], domLev, solver, dxLev, ilev, nlevels,
                                 noFlowVelocityBC, noFlowGradientBC));

        LevelCCProjector  ccProj(a_grids[ilev], domLev, macProj, dxLev, ilev, nlevels, a_params.refRatio[ilev]);
        LevelData<FArrayBox>* presCoar = NULL;
        LevelData<FArrayBox>* veloCoar = NULL;
        if(ilev > 0)
          {
            presCoar = pressure[ilev-1];
            veloCoar = velocity[ilev-1];
          }
        ccProj.project(*velocity[ilev],
                       *gradient[ilev],
                       *pressure[ilev],
                       veloCoar,
                       presCoar);

        ccProj.divergence(*a_error [ilev],
                          *velocity[ilev],
                          veloCoar);

        //prepare dx, domain for next level
        dxLev /=      a_params.refRatio[ilev];
        domLev.refine(a_params.refRatio[ilev]);
      }
  }

  //delete the local news.  error must survive outside this scope
  for(int ilev = 0; ilev < nlevels; ilev++)
    {
      delete velocity[ilev];
      delete pressure[ilev];
      delete gradient[ilev];
    }
}

/***************/
void outputError(const Vector< LevelData<FArrayBox>* >&   a_errorFine,
                 const Vector< LevelData<FArrayBox>* >&   a_errorCoar,
                 const Vector< DisjointBoxLayout >&       a_gridsFine,
                 const Vector< DisjointBoxLayout >&       a_gridsCoar,
                 const PoissonParameters&                 a_paramsFine,
                 const PoissonParameters&                 a_paramsCoar)
{
#if CH_SPACEDIM==2
  string fileFine("pltFineError.2d.hdf5");
  string fileCoar("pltCoarError.2d.hdf5");
#else
  string fileFine("pltFineError.3d.hdf5");
  string fileCoar("pltCoarError.3d.hdf5");
#endif
  string phiname("error");
  outputData(a_errorFine, a_gridsFine,
             a_paramsFine.coarsestDomain, a_paramsFine.refRatio,
             a_paramsFine.coarsestDx, a_paramsFine.numLevels,
             fileFine, phiname);
  outputData(a_errorCoar, a_gridsCoar,
             a_paramsCoar.coarsestDomain, a_paramsCoar.refRatio,
             a_paramsCoar.coarsestDx, a_paramsCoar.numLevels,
             fileCoar, phiname);
}
int main(int argc, char* argv[])
{
#ifdef CH_MPI
  MPI_Init(&argc, &argv);
#endif
  // Scoping trick
  {

    if (argc < 2)
      {
        cerr << " usage " << argv[0] << " <input_file_name> " << endl;
        exit(0);
      }

    char* inFile = argv[1];
    ParmParse pp(argc-2,argv+2,NULL,inFile);

    PoissonParameters paramFine, paramCoar;
    Vector<DisjointBoxLayout> gridsFine, gridsCoar;

    //read params from file
    getPoissonParameters(paramFine);
    paramCoar = paramFine;
    paramCoar.coarsen(2);
    int nlevels = paramCoar.numLevels;
    Vector<LevelData<FArrayBox>* > errorFine(nlevels, NULL);
    Vector<LevelData<FArrayBox>* > errorCoar(nlevels, NULL);

    setGrids(gridsFine,  paramFine);

    pout() << "generating fine error" << endl;
    getError(errorFine, gridsFine,  paramFine);

    getCoarseLayoutsFromFine(gridsCoar, gridsFine, paramCoar);

    pout() << "generating coarse error" << endl;
    getError(errorCoar, gridsCoar, paramCoar);

    int dofileout;
    pp.get("do_error_output", dofileout);
    if(dofileout == 1)
      {
        outputError(errorFine,   errorCoar,
                    gridsFine,   gridsCoar,
                    paramFine,   paramCoar);
      }

    string testname("Solution error");
    compareError(errorFine,   errorCoar,
                 gridsFine,   gridsCoar,
                 paramFine,   paramCoar,
                 testname);

    for(int ilev = 0; ilev < nlevels; ilev++)
      {
        delete errorFine[ilev];
        delete errorCoar[ilev];
        errorFine[ilev] = NULL;
        errorCoar[ilev] = NULL;
      }

  }// End scoping trick

#ifdef CH_MPI
  MPI_Finalize();
#endif
}

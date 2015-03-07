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
#include "FABView.H"
#include "DebugDump.H"
#include "PoissonUtilities.H"
#include "BiCGStabSolver.H"


/******/
void getError(Vector< LevelData<FArrayBox>* >&     a_error,
              const Vector< DisjointBoxLayout >&   a_grids, 
              const PoissonParameters&             a_params)
{
  ParmParse pp;

  int nlevels = a_params.numLevels;
  a_error.resize(nlevels);
  Vector<LevelData<FArrayBox>* > phiExac(nlevels, NULL);
  Vector<LevelData<FArrayBox>* > klpExac(nlevels, NULL);
                                         
  AMRMultiGrid<LevelData<FArrayBox> > solver;
  BiCGStabSolver<LevelData<FArrayBox> >   bottomSolver;
  bottomSolver.m_verbosity = 0;
  defineSolver(solver, a_grids, bottomSolver, a_params);

  RealVect dxLev = RealVect::Unit;
  dxLev *=a_params.coarsestDx;
  ProblemDomain domLev(a_params.coarsestDomain);
  for(int ilev = 0; ilev < nlevels; ilev++)
    {
      a_error[ilev] = new LevelData<FArrayBox>(a_grids[ilev], 1,IntVect::Unit);
      klpExac[ilev] = new LevelData<FArrayBox>(a_grids[ilev], 1,IntVect::Zero);
      phiExac[ilev] = new LevelData<FArrayBox>(a_grids[ilev], 1,IntVect::Unit);
      for(DataIterator dit = a_grids[ilev].dataIterator(); dit.ok(); ++dit)
        {
          (*a_error[ilev])[dit()].setVal(0.);
        }

      //set phi = phiExact, rhs=lphiexact  This makes AMRResidual return lphiexact-Lphi
      setTrigPhi(         *phiExac[ilev], dxLev, a_params);
      setTrigKappaLOfPhi (*klpExac[ilev], dxLev, a_params);

      //prepare dx, domain for next level
      dxLev /=      a_params.refRatio[ilev];
      domLev.refine(a_params.refRatio[ilev]);
    }

  //create calculated data and error
  int lbase = 0;
  solver.solve(a_error, klpExac, a_params.numLevels-1, lbase);
  for(int ilev = 0; ilev < nlevels; ilev++)
    {
      for(DataIterator dit = a_grids[ilev].dataIterator(); dit.ok(); ++dit)
        {
          (*a_error[ilev])[dit()] -= (*phiExac[ilev])[dit()];
        }
    }
  //delete the local news.  error must survive outside this scope
  for(int ilev = 0; ilev < nlevels; ilev++)
    {
      delete phiExac[ilev];
      delete klpExac[ilev];
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

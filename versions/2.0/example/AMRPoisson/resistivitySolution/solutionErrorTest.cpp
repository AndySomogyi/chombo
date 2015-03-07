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
#include "ResistivityOp.H"
#include "FABView.H"
#include "DebugDump.H"
#include "PoissonUtilities.H"
#include "BiCGStabSolver.H"


/******/
void getError(Vector< LevelData<FArrayBox>* >&     a_error,
              Vector< LevelData<FArrayBox>* >&     a_divB,
              const Vector< DisjointBoxLayout >&   a_grids, 
              const PoissonParameters&             a_params,
              bool a_isFine)
{
  ParmParse pp;

  int nlevels = a_params.numLevels;
  a_error.resize(nlevels);
  Vector<LevelData<FArrayBox>* > magExac(nlevels, NULL);
  Vector<LevelData<FArrayBox>* > klpExac(nlevels, NULL);
  Vector<ResistivityOp*>         ebamrpo(nlevels, NULL);
  Vector<RefCountedPtr<LevelData<FluxBox> > >  etaExac(nlevels);
                                         
  AMRMultiGrid<LevelData<FArrayBox> > solver;
  BiCGStabSolver<LevelData<FArrayBox> >   bottomSolver;
  bottomSolver.m_verbosity = 0;
  defineResistivitySolver(solver, a_grids, bottomSolver, a_params);

  RealVect dxLev = RealVect::Unit;
  dxLev *=a_params.coarsestDx;
  ProblemDomain domLev(a_params.coarsestDomain);
  for(int ilev = 0; ilev < nlevels; ilev++)
    {
      a_error[ilev] = new LevelData<FArrayBox>(a_grids[ilev], SpaceDim, IntVect::Unit);
      a_divB [ilev] = new LevelData<FArrayBox>(a_grids[ilev],        1, IntVect::Zero);
      klpExac[ilev] = new LevelData<FArrayBox>(a_grids[ilev], SpaceDim, IntVect::Zero);
      magExac[ilev] = new LevelData<FArrayBox>(a_grids[ilev], SpaceDim, IntVect::Unit);
      etaExac[ilev] = RefCountedPtr<LevelData<FluxBox> >(new LevelData<FluxBox>  (a_grids[ilev],        1,  IntVect::Zero));
      for(DataIterator dit = a_grids[ilev].dataIterator(); dit.ok(); ++dit)
        {
          (*a_error[ilev])[dit()].setVal(0.);
        }

      //set mag = magExact, rhs=lmagexact  This makes AMRResidual return lmagexact-Lmag
      setMagResistive(*magExac[ilev], dxLev, a_params);
      setKLBResistive(*klpExac[ilev], dxLev, a_params);
      setEtaResistive(*etaExac[ilev], dxLev, a_params);

      //prepare dx, domain for next level
      dxLev /=      a_params.refRatio[ilev];
      domLev.refine(a_params.refRatio[ilev]);

    }

  ResistivityOpFactory opFactory(a_grids, etaExac,
                                 a_params.alpha, a_params.beta,
                                 a_params.refRatio,
                                 a_params.coarsestDomain,
                                 a_params.coarsestDx, 
                                 &ParseBC);

  domLev = ProblemDomain(a_params.coarsestDomain);
  for(int ilev = 0; ilev < nlevels; ilev++)
    {
      ebamrpo[ilev] = opFactory.AMRnewOp(domLev);
      domLev.refine(a_params.refRatio[ilev]);
    }

  //create calculated data and error
  int lbase = 0;
  solver.m_verbosity = 4;
  solver.solve(a_error, klpExac, a_params.numLevels-1, lbase);

  //compute div B
  Vector< LevelData<FArrayBox> *>& mag = a_error;
  for(int ilev = 0; ilev < nlevels; ilev++)
    {
      const LevelData<FArrayBox>* magC = NULL;
      if(ilev > 0)
        {
          magC = mag[ilev-1];
        }
      ebamrpo[ilev]->divergenceCC((*a_divB[ilev]), (*mag[ilev]), magC);
    }

  string fileKLBExac;
  string fileMagCalc;
  string fileMagExac;
  if(a_isFine)
    {
#if CH_SPACEDIM==2
      fileKLBExac =   string("pltKLBExacFine.2d.hdf5");
      fileMagCalc =   string("pltMagCalcFine.2d.hdf5");
      fileMagExac =   string("pltMagExacFine.2d.hdf5");
#else                       
      fileKLBExac =   string("pltKLBExacFine.3d.hdf5");
      fileMagCalc =   string("pltMagCalcFine.3d.hdf5");
      fileMagExac =   string("pltMagExacFine.3d.hdf5");
#endif
    }
  else
    {
#if CH_SPACEDIM==2
      fileKLBExac =   string("pltKLBExacCoar.2d.hdf5");
      fileMagCalc =   string("pltMagCalcCoar.2d.hdf5");
      fileMagExac =   string("pltMagExacCoar.2d.hdf5");
#else                                 
      fileKLBExac =   string("pltKLBExacCoar.3d.hdf5");
      fileMagCalc =   string("pltMagCalcCoar.3d.hdf5");
      fileMagExac =   string("pltMagExacCoar.3d.hdf5");
#endif
    }

  outputData(klpExac, a_grids, 
             a_params.coarsestDomain, a_params.refRatio,
             a_params.coarsestDx,     a_params.numLevels, 
             fileKLBExac, string("klbExac"));
  outputData(a_error, a_grids, 
             a_params.coarsestDomain, a_params.refRatio,
             a_params.coarsestDx,     a_params.numLevels, 
             fileMagCalc, string("magCalc"));
  outputData(magExac, a_grids, 
             a_params.coarsestDomain, a_params.refRatio,
             a_params.coarsestDx,     a_params.numLevels, 
             fileMagExac, string("magExac"));

  for(int ilev = 0; ilev < nlevels; ilev++)
    {
      for(DataIterator dit = a_grids[ilev].dataIterator(); dit.ok(); ++dit)
        {
          (*a_error[ilev])[dit()] -= (*magExac[ilev])[dit()];
        }
    }

  //delete the local news.  error must survive outside this scope
  for(int ilev = 0; ilev < nlevels; ilev++)
    {
      delete magExac[ilev];
      delete klpExac[ilev];
      delete ebamrpo[ilev];
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
    Vector<LevelData<FArrayBox>* >  divBFine(nlevels, NULL);
    Vector<LevelData<FArrayBox>* >  divBCoar(nlevels, NULL);

    setGrids(gridsFine,  paramFine);

    pout() << "generating fine error" << endl;
    getError(errorFine, divBFine, gridsFine,  paramFine, true);

    getCoarseLayoutsFromFine(gridsCoar, gridsFine, paramCoar);
    
    pout() << "generating coarse error" << endl;
    getError(errorCoar, divBCoar, gridsCoar, paramCoar, false);
    
    int dofileout;
    pp.get("do_error_output", dofileout);
    if(dofileout == 1)
      {
        outputError(errorFine,   errorCoar,
                    gridsFine,   gridsCoar, 
                    paramFine,   paramCoar);
      }

    string testname("Solution error");
    pout() << "convergence test for error in solution of L(B) = rho" << endl;
    compareError(errorFine,   errorCoar,
                 gridsFine,   gridsCoar, 
                 paramFine,   paramCoar,
                 testname);

    
    pout() << "convergence test for divergence(B) in solution of L(B) = rho" << endl;
    pout() << "Warning: this will only look second order if divExact(Bexact) = 0" << endl;
    string divbname("Divergence of B");
    compareError(divBFine,    divBCoar,
                 gridsFine,   gridsCoar, 
                 paramFine,   paramCoar,
                 divbname);

    for(int ilev = 0; ilev < nlevels; ilev++)
      {
        delete errorFine[ilev];
        delete errorCoar[ilev];
        delete divBFine [ilev];
        delete divBCoar [ilev];

        errorFine[ilev] = NULL;
        errorCoar[ilev] = NULL;
        divBFine[ilev]  = NULL;
        divBCoar[ilev]  = NULL;
      }

  }// End scoping trick

#ifdef CH_MPI
  MPI_Finalize();
#endif
}

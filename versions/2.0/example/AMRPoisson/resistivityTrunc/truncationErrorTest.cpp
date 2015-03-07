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


/******/
void getError(Vector< LevelData<FArrayBox>* >&     a_error,
              const Vector< DisjointBoxLayout >&   a_grids, 
              const PoissonParameters&             a_params)
{
  ParmParse pp;

  int nlevels = a_params.numLevels;
  a_error.resize(nlevels);
  Vector<LevelData<FArrayBox>* > magExac(nlevels, NULL);
  Vector<LevelData<FArrayBox>* > klbExac(nlevels, NULL);
  Vector<ResistivityOp*>         ebamrpo(nlevels, NULL);
  Vector<RefCountedPtr<LevelData<FluxBox> > >  etaExac(nlevels);
  RealVect dxLev = RealVect::Unit;
  dxLev *=a_params.coarsestDx;
  for(int ilev = 0; ilev < nlevels; ilev++)
    {
      a_error[ilev] = new LevelData<FArrayBox>(a_grids[ilev], SpaceDim,  IntVect::Zero);
      klbExac[ilev] = new LevelData<FArrayBox>(a_grids[ilev], SpaceDim,  IntVect::Zero);
      magExac[ilev] = new LevelData<FArrayBox>(a_grids[ilev], SpaceDim,3*IntVect::Unit);
      etaExac[ilev] = RefCountedPtr<LevelData<FluxBox> >(new LevelData<FluxBox>  (a_grids[ilev],        1,  IntVect::Zero));
      //set mag = magExact, rhs=lmagexact  This makes AMRResidual return lmagexact-Lmag
      setEtaResistive(*etaExac[ilev], dxLev, a_params);
      setMagResistive(*magExac[ilev], dxLev, a_params);
      setKLBResistive(*klbExac[ilev], dxLev, a_params);
      dxLev /=      a_params.refRatio[ilev];
    }
                                         
  ResistivityOpFactory opFactory(a_grids, etaExac,
                                 a_params.alpha, a_params.beta,
                                 a_params.refRatio,
                                 a_params.coarsestDomain,
                                 a_params.coarsestDx, 
                                 &ParseBC);

  ProblemDomain domLev(a_params.coarsestDomain);
  for(int ilev = 0; ilev < nlevels; ilev++)
    {
      ebamrpo[ilev] = opFactory.AMRnewOp(domLev);
      domLev.refine(a_params.refRatio[ilev]);
    }

    
  //create calculated data and error
  for(int ilev = 0; ilev < nlevels; ilev++)
    {
      bool homogeneous = false;
      if(nlevels == 1)
        {
          //single level calc
          ebamrpo[ilev]->residual(*a_error[ilev], *magExac[ilev], *klbExac[ilev], homogeneous);
        }
      else if(ilev == nlevels-1)  //at finest level  of multi level calc
        {
          ebamrpo[ilev]->AMRResidualNF(*a_error[ilev], *magExac[ilev], *magExac[ilev-1],  
                                       *klbExac[ilev], homogeneous);
        }
      else if(ilev == 0)  //at lowest level of multi level calc
        {
          ebamrpo[ilev]->AMRResidualNC(*a_error[ilev], *magExac[ilev+1], *magExac[ilev],
                                        *klbExac[ilev], homogeneous, ebamrpo[ilev+1]);
        }
      else  //mid level of multi level calc 
        {
          ebamrpo[ilev]->AMRResidual  (*a_error[ilev], *magExac[ilev+1], *magExac[ilev], *magExac[ilev-1],
                                       *klbExac[ilev],   homogeneous, ebamrpo[ilev+1]);
        
        }
    }  

  //delete the local news.  error must survive outside this scope and eta
  // is refcounted 
  for(int ilev = 0; ilev < nlevels; ilev++)
    {
      delete klbExac[ilev];
      delete magExac[ilev];
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

    string testname("Truncation error");
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

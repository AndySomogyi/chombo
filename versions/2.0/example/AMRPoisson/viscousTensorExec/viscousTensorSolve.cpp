#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

#include  <iostream>
#include "LevelData.H"
#include "FArrayBox.H"
#include "ParmParse.H"
#include "Vector.H"
#include "AMRIO.H"
#include "SPMD.H"
#include "BRMeshRefine.H"
#include "LoadBalance.H"
#include "computeSum.H"
#include "PoissProbF_F.H"
#include "ProblemDomain.H"
#include "BCFunc.H"
#include "ViscousTensorOp.H"
#include "AMRMultiGrid.H"
#include "BoxIterator.H"
#include "BiCGStabSolver.H"
#include "DebugDump.H"
#include "PoissonUtilities.H"
#include <fstream>

void
outputData(const Vector<LevelData<FArrayBox>* >& vectPhi,
           const Vector<LevelData<FArrayBox>* >& vectRhs,
           const Vector<DisjointBoxLayout>& vectGrids,
           const PoissonParameters&   params)
{
  string phiname("phi");
  string phifile("phi.hdf5");
  string rhsname("rhs");
  string rhsfile("rhs.hdf5");
  outputData(vectPhi, vectGrids, params.coarsestDomain, params.refRatio,
             params.coarsestDx, params.numLevels, phifile, phiname);
  outputData(vectRhs, vectGrids, params.coarsestDomain, params.refRatio,
             params.coarsestDx, params.numLevels, rhsfile, rhsname);
}

/*****/
void runSolver()
{
  PoissonParameters params;
  getPoissonParameters(params);

  ParmParse pp2;

  Vector<DisjointBoxLayout> vectGrids;
  setGrids(vectGrids,  params);
  
  int numlevels = params.numLevels;
  for(int ilev = 0; ilev < vectGrids.size(); ilev++)
    {
      pout() << "grid at level " << ilev << " = " << vectGrids[ilev] << endl;
    }

  Vector<LevelData<FArrayBox>* > phi(numlevels, NULL);
  Vector<LevelData<FArrayBox>* > rhs(numlevels, NULL);

  int ncomp  = SpaceDim;
  for(int ilev = 0; ilev < numlevels; ilev++)
    {
      phi[ilev] = new LevelData<FArrayBox>(vectGrids[ilev], ncomp,  2*IntVect::Unit);
      rhs[ilev] = new LevelData<FArrayBox>(vectGrids[ilev], ncomp,    IntVect::Zero);
      LevelData<FArrayBox>& phifabs= *phi[ilev];
      LevelData<FArrayBox>& rhsfabs= *rhs[ilev];
      for(DataIterator dit = phifabs.dataIterator(); dit.ok(); ++dit)
        {
          phifabs[dit()].setVal(0.);
          rhsfabs[dit()].setVal(1.);
        }
    }

  //  setRHS(rhs, params);
  AMRMultiGrid<LevelData<FArrayBox> > solver;
  BiCGStabSolver<LevelData<FArrayBox> >   bottomSolver;
  bottomSolver.m_verbosity = 0;

  defineViscousTensorSolver(solver, vectGrids, bottomSolver, params);
  int lbase;
  pp2.get("lbase", lbase);
  solver.m_verbosity =5;
  //  solver.solve(phi, rhs, numlevels-1, lbase);
  solver.solve(phi, rhs, numlevels-1, lbase);
  //solver.relaxOnlyHomogeneous(phi, rhs, numlevels-1, lbase);

  outputData(phi, rhs, vectGrids, params);

  for(int ilev = 0; ilev < numlevels; ilev++)
    {
      delete phi[ilev];
      delete rhs[ilev];
    }
}
/*****/
int main(int argc, char* argv[])
{
#ifdef CH_MPI
  MPI_Init(&argc, &argv);
#endif
  //scoping trick
  {
    if (argc < 2) 
      {
        cerr<< " usage " << argv[0] << " <input_file_name> " << endl;
        exit(0);
      }
    char* in_file = argv[1];
    ParmParse  pp(argc-2,argv+2,NULL,in_file);

    int iterations;
    pp.get("iterations", iterations);
    for(int iiter = 0; iiter < iterations; iiter++)
      {
        runSolver();
      }

  }//end scoping trick

#ifdef CH_MPI
  MPI_Finalize();
#endif

  return(0);
}

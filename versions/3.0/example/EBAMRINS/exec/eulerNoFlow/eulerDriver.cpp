#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

#include "ParmParse.H"
#include "CH_HDF5.H"
#include "parstream.H"


#include "BaseIVFactory.H"
#include "EBIndexSpace.H"
#include "EBISLayout.H"
#include "BoxIterator.H"
#include "ParmParse.H"
#include "BRMeshRefine.H"
#include "LoadBalance.H"
#include "GeometryShop.H"
#include "LevelData.H"
#include "EBCellFAB.H"
#include "EBCellFactory.H"
#include "VoFIterator.H"
#include "EBArith.H"
#include "DebugDump.H"
#include "EBDebugDump.H"
#include "EBFABView.H"
#include "memtrack.H"
#include "AMRINSUtils.H"
#include "CH_Attach.H"
#include "EBAMRNoSubcycle.H"
#include "NoFlowVortex.H"
#include <iostream>
#include "UsingNamespace.H"

/***************/
void ebamrieuler(const AMRParameters& a_params,
                 const ProblemDomain& a_coarsestDomain)
{

  // read inputs
  ParmParse ppebamrieuler;

  int flowDir;
  ppebamrieuler.get("flow_dir", flowDir);
  Vector<Real> centerVect(SpaceDim);
  IntVect hiSide;
  Vector<int> nCells;
  ppebamrieuler.getarr("n_cell",  nCells,0,SpaceDim);
  ppebamrieuler.getarr("spot_center",centerVect,0,SpaceDim);
  RealVect center;
  for(int idir = 0; idir < SpaceDim; idir++)
    {
      center[idir] = centerVect[idir];
      hiSide[idir] = nCells[idir] -1;
    }

  Real inflowVel = 100.0;
  Real coreRadius;
  ppebamrieuler.get("spot_radius", coreRadius);

  Real viscosity = 0.0;
  ppebamrieuler.get("viscosity", viscosity);


  NoFlowVortexFactory ibc(center, coreRadius, inflowVel);


  EBAMRNoSubcycle  kahuna(a_params, ibc, a_coarsestDomain, viscosity);
  if (!ppebamrieuler.contains("restart_file"))
    {
      pout() << "starting fresh AMR run" << endl;
      kahuna.setupForAMRRun();
    }
  else
    {
      std::string restart_file;
      ppebamrieuler.get("restart_file",restart_file);
      pout() << " restarting from file " << restart_file << endl;
      kahuna.setupForRestart(restart_file);
    }

  Real fixedDt;
  ppebamrieuler.get("fixed_dt", fixedDt);
  if(fixedDt > 0)
    {
      kahuna.useFixedDt(fixedDt);
    }
  int maxStep;
  ppebamrieuler.get("max_step", maxStep);

  Real stopTime = 0.0;
  ppebamrieuler.get("max_time",stopTime);

  kahuna.run(stopTime, maxStep);
}

/***************/
int
main(int a_argc, char* a_argv[])
{
#ifdef CH_MPI
  MPI_Init(&a_argc,&a_argv);
#endif
  { //scoping trick
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
    // Parse the command line and the input file (if any)
    ParmParse pp(a_argc-2,a_argv+2,NULL,inFile);

    ProblemDomain coarsestDomain;
    AMRParameters params;
    getAMRINSParameters(params, coarsestDomain);
    //define geometry
    AMRINSGeometry(params, coarsestDomain);

    ebamrieuler(params, coarsestDomain);
  }//end scoping trick
#ifdef CH_MPI
  MPI_Finalize();
#endif
}

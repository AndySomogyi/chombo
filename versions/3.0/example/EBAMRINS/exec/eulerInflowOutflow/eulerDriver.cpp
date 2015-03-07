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
#include "InflowOutflowIBC.H"
#include "EBFABView.H"
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

  Real coreRadius, inflowVel;
  ppebamrieuler.get("spot_radius", coreRadius);
  ppebamrieuler.get("inflow_vel", inflowVel);

  Real viscosity = 0.0;
  ppebamrieuler.get("viscosity", viscosity);

  int idoSlipWalls;
  ppebamrieuler.get("do_slip_walls", idoSlipWalls);
  IntVect doSlipWallsLo = idoSlipWalls*IntVect::Unit;
  IntVect doSlipWallsHi = idoSlipWalls*IntVect::Unit;


  int orderEBBC =1;
  ParmParse pp;
  pp.query("order_ebbc", orderEBBC);
  bool doSlip = idoSlipWalls==1;
  InflowOutflowIBCFactory ibc(flowDir, inflowVel, orderEBBC, doSlip);


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

  int maxStep;
  ppebamrieuler.get("max_step", maxStep);

  Real stopTime = 0.0;
  ppebamrieuler.get("max_time",stopTime);
  //debug
  //kahuna.useFixedDt(1.0e-4);
  //end debug
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
    int tagOnScalar;
    pp.get("tag_on_scalar", tagOnScalar);
    params.m_tagOnScalarToo = (tagOnScalar==1);
    int numFilt;
    pp.get("num_filter_iterations", numFilt);
    params.m_numFilterIterations = numFilt;

    int gphiIterations;
    pp.get("num_gphi_iterations", gphiIterations);
    params.m_gphiIterations = gphiIterations;

    //hack to keep outflow from getting refined
    params.m_chopXTagsHi = true;
    //define geometry
    AMRINSGeometry(params, coarsestDomain);

    ebamrieuler(params, coarsestDomain);
  }//end scoping trick
#ifdef CH_MPI
  MPI_Finalize();
#endif
}

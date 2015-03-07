#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

// Purpose:
//  Test the infrastucture for periodic boundary conditions
//
// Usage:
//  <program-name> [-q|-v] ...
//
//  where:
//    -q means run quietly (only pass/fail messages printed)
//    -v means run verbosely (all messages printed)
//    -write means write out the boxes; otherwise will read in
//           boxes from a file and compare with boxes generated
//           in tests.
//    ... all non-option arguments are ignored (Chombo convention)
//
//  Default is `-v'
//
//  Unknown options are treated as errors.
//

// Include files:
#include <cstdio>

#include "parstream.H"
#include "LoadBalance.H"
#include "FArrayBox.H"
#include "LevelData.H"
#include "ProblemDomain.H"
#include "BoxIterator.H"
#include "BRMeshRefine.H"
#include "AMRIO.H"

#include "FineInterp.H"
#include "PiecewiseLinearFillPatch.H"
#include "GhostBC.H"
#include "PoissonBC.H"

#ifdef CH_MPI
#include "mpi.h"
#endif
#include "UsingNamespace.H"

//////////////////////////////////////////////////////////////
using std::endl;

void parseTestOptions(int argc, char* argv[]);

void initData(LevelData<FArrayBox>& a_data, const Real a_dx);

/// Global variables for handling output

static const char *pgmname = "testPeriodic";
static const char *indent = "   " ,*indent2 = "      " ;

static bool verbose = true ;

/// Code:

int
main(int argc ,char *argv[] )
{
#if 0   // Problem with header consistency, branch vs trunk.  Fix on merge.


#ifdef CH_MPI
  MPI_Init(&argc, &argv);
#endif

  parseTestOptions(argc, argv);

  // establish periodic domain -- first multiply periodic in all directions
  int status = 0;
  int baseDomainSize = 8;
  int numGhost = 2;
  IntVect ghostVect(numGhost*IntVect::Unit);

  Box baseDomBox(IntVect::Zero, (baseDomainSize-1)*IntVect::Unit);
  ProblemDomain baseDomain(baseDomBox);

  // set periodic in all directions
  for (int dir=0; dir<SpaceDim; dir++) {
    baseDomain.setPeriodic(dir, true);
  }

  // set up physical BC class for partially periodic tests
  DomainGhostBC physBC;
  for (int dir=0; dir<SpaceDim; dir++) {
    SideIterator sit;
    for (sit.begin(); sit.ok(); ++sit)
      {
        // it should be easy to see if Neumann bcs are doing the right thing
        NeumannBC thisBC(dir,sit());
        physBC.setBoxGhostBC(thisBC);
      }
  }

  {
    // FineInterp test
    int nRef = 2;
    ProblemDomain fineDomain(baseDomain);
    fineDomain.refine(nRef);

    const Box domainBox = baseDomain.domainBox();
    Vector<Box> crseBoxes(1, domainBox);
    Vector<int> crseProcAssign(1,0);
    DisjointBoxLayout crseGrids(crseBoxes, crseProcAssign, baseDomain);

    // note the lack of ghost cells on the base level
    LevelData<FArrayBox> crseData(crseGrids, 1, IntVect::Zero);

    Real dx = 1.0/(baseDomainSize);
    initData(crseData, dx);

    // try multiple boxes for fine level
    Vector<Box> fineBoxes(3);
    Vector<int> fineProcAssign(3);

    fineBoxes[0] = Box(IntVect::Zero, 9*IntVect::Unit);
    fineBoxes[1] = Box(IntVect(D_DECL(6,10,10)), 15*IntVect::Unit);
    // this box should fail the disjointness test
    //fineBoxes[1] = Box(IntVect(D_DECL(7,10,10)), 17*IntVect::Unit);
    fineBoxes[2] = Box(IntVect(D_DECL(12,0,0)), IntVect(D_DECL(15,7,7)));

    int loadbalancestatus = LoadBalance(fineProcAssign, fineBoxes);
    CH_assert (loadbalancestatus == 0);

    DisjointBoxLayout fineGrids;
    fineGrids.define(fineBoxes, fineProcAssign, fineDomain);
    LevelData<FArrayBox> fineData(fineGrids, 1, ghostVect);

    // set these to a bogus value to start with
    DataIterator fineDit = fineData.dataIterator();
    for (fineDit.begin(); fineDit.ok(); ++fineDit)
      {
        fineData[fineDit()].setVal(1.0e9);
      }

    FineInterp interpolator(fineGrids, 1, nRef);

    interpolator.interpToFine(fineData, crseData);

    // now fill in ghost cells
    PiecewiseLinearFillPatch filler(fineGrids, crseGrids, 1, nRef,
                                    numGhost);

    filler.fillInterp(fineData, crseData, crseData, 1.0, 0, 0, 1);
  }

  {
    // now test partially periodic case (also coarsen back to original domain
    baseDomain.setPeriodic(0,false);
    // try similar tests to the previous case
    int nRef = 2;
    ProblemDomain fineDomain(baseDomain);
    fineDomain.refine(nRef);

    const Box domainBox = baseDomain.domainBox();
    Vector<Box> crseBoxes(1, domainBox);
    Vector<int> crseProcAssign(1,0);
    DisjointBoxLayout crseGrids(crseBoxes, crseProcAssign, baseDomain);

    // note the lack of ghost cells on the base level
    LevelData<FArrayBox> crseData(crseGrids, 1, IntVect::Zero);

    Real dx = 1.0/(baseDomainSize);
    initData(crseData, dx);

    // try multiple boxes for fine level
    Vector<Box> fineBoxes(3);
    Vector<int> fineProcAssign(3);

    fineBoxes[0] = Box(IntVect::Zero, 9*IntVect::Unit);
    fineBoxes[1] = Box(IntVect(D_DECL(6,10,10)), 15*IntVect::Unit);
    // this box should fail the disjointness test
    //fineBoxes[1] = Box(IntVect(D_DECL(7,10,10)), 17*IntVect::Unit);
    fineBoxes[2] = Box(IntVect(D_DECL(12,0,0)), IntVect(D_DECL(15,7,7)));

    int loadbalancestatus = LoadBalance(fineProcAssign, fineBoxes);
    CH_assert (loadbalancestatus == 0);

    DisjointBoxLayout fineGrids;
    fineGrids.define(fineBoxes, fineProcAssign, fineDomain);
    LevelData<FArrayBox> fineData(fineGrids, 1, ghostVect);

    // set these to a bogus value to start with
    DataIterator fineDit = fineData.dataIterator();
    for (fineDit.begin(); fineDit.ok(); ++fineDit)
      {
        fineData[fineDit()].setVal(1.0e9);
      }

    FineInterp interpolator(fineGrids, 1, nRef);

    interpolator.interpToFine(fineData, crseData);

    // now fill in ghost cells
    PiecewiseLinearFillPatch filler(fineGrids, crseGrids, 1, nRef,
                                    numGhost);

    filler.fillInterp(fineData, crseData, crseData, 1.0, 0, 0, 1);

    fineData.exchange(fineData.interval());

  }

#ifdef CH_MPI
  MPI_Finalize();
#endif

  pout() << indent << pgmname << ": "
     << ( (status == 0) ? "passed all tests" : "failed at least one test,")
     << endl;

  return status ;
#endif // 0
}

///
// Parse the standard test options (-v -q) out of the command line.
// Stop parsing when a non-option argument is found.
///
void
parseTestOptions( int argc ,char* argv[] )
{
  for( int i = 1 ; i < argc ; ++i )
    {
      if( argv[i][0] == '-' ) //if it is an option
        {
          // compare 3 chars to differentiate -x from -xx
          if( strncmp( argv[i] ,"-v" ,3 ) == 0 )
            {
              verbose = true ;
              // argv[i] = "" ;
            }
          else
            {
              break ;
            }
        }
    }
  return ;
}

void
initData(LevelData<FArrayBox>& a_data, const Real a_dx)
{

  DataIterator dit = a_data.dataIterator();
  const DisjointBoxLayout& interiorBoxes = a_data.getBoxes();

  for (dit.begin(); dit.ok(); ++dit)
    {
      // first set to a bogus value which will persist in ghost cells
      // after initialization
      a_data[dit()].setVal(1.0e9);

      // this will be slow, but who cares?
      FArrayBox& localData = a_data[dit()];
      BoxIterator boxIt(interiorBoxes[dit()]);
      Real localVal;
      for (boxIt.begin(); boxIt.ok(); ++boxIt)
        {
          const IntVect& loc = boxIt();
          localVal = a_dx*(D_TERM(loc[0]+0.5,
                                + loc[1]+0.5,
                                + loc[2]+0.5));
          localData(loc) = localVal;
        }

    }
}

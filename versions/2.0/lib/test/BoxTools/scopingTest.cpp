#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

#include <cmath>
#include "Misc.H"
#include "FArrayBox.H"
#include "DisjointBoxLayout.H"
#include "LevelData.H"
#include "BoxIterator.H"
#include "LayoutIterator.H"
#include "LoadBalance.H"
#include "UsingNamespace.H"

/// Global variables for test output:

static const char *pgmname = "scopingTest" ;
static const char *indent = "   " ,*indent2 = "      " ;
static bool verbose = true ;

/// Prototypes:

void
parseTestOptions( int argc ,char* argv[] ) ;
void
setGrids( Vector<Box>& coarboxes);
void
makeLevelData(const Vector<Box>& boxes,
              LevelData<FArrayBox>& levelFab,
              Vector<DataIndex>& vectInd);

/**
   scopingTest returns:
    0: all tests passed
 */
int scopingTest(void);

int main(int argc, char* argv[])
{
#ifdef CH_MPI
  MPI_Init(&argc, &argv);
#endif
  bool passed = true;
  int icode = 0;

  parseTestOptions( argc ,argv ) ;
  if( verbose )
    cout << indent2 << "Beginning " << pgmname << " ..." << endl ;

#ifndef CH_MPI
  icode = scopingTest();
#endif

  if(icode != 0)
    {
      cout << indent << pgmname
           << ": failed with error code " << icode
           << endl;
      passed = false;
    }

  if(passed)
    cout << indent << pgmname
           << ": passed all tests"
           << endl;
#ifdef CH_MPI
  MPI_Finalize();
#endif

  return 0;
}

int scopingTest()
{
  Vector<Box> boxes;
  int retflag = 0;
  setGrids(boxes);

  Vector<DataIndex> vectInd;
  LevelData<FArrayBox> levelFab;
  makeLevelData(boxes, levelFab, vectInd);
  const DisjointBoxLayout& dbl = levelFab.getBoxes();
  for(int ivec = 0; ivec < vectInd.size(); ivec++)
    {
      if(!dbl.check(vectInd[ivec]))
        {
          if(verbose)
            cout << indent2 << pgmname
                 << ": failed at box " << ivec << endl;
          retflag = 1;
        }
      else
        {
          levelFab[vectInd[ivec]].setVal(0.);
        }
      if(retflag > 0) break;
    }
  return retflag;
}
/*****************/
/*****************/
void
makeLevelData(const Vector<Box>& boxes,
              LevelData<FArrayBox>& levelFab,
              Vector<DataIndex>& vectInd)
{
  vectInd.resize(boxes.size());
  Vector<int> procIDs(boxes.size());
  LoadBalance(procIDs, boxes);

  DisjointBoxLayout dbl1;
  for(int ibox = 0; ibox < boxes.size(); ++ibox)
    {
     vectInd[ibox] = dbl1.addBox(boxes[ibox], procIDs[ibox]);
    }
  dbl1.close();
  levelFab.define(dbl1, 1);
}
/*****************/
/*****************/
void
setGrids(Vector<Box>& boxes)
{
  int nc = 16;
  IntVect ivclo = IntVect::Zero;
  IntVect ivchi = (nc-1)*IntVect::Unit;

  Box domc = Box(ivclo, ivchi);
  int ichop = domc.size(0)/2;
  Box bc1 = domc;
  Box bc2 = bc1.chop(0, ichop);
#if CH_SPACEDIM != 1
  Box bc3 = bc2.chop(1, ichop);
#endif
  boxes.resize(0);
  boxes.push_back(bc1);
  boxes.push_back(bc2);
#if CH_SPACEDIM != 1
  boxes.push_back(bc3);
#endif
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
              argv[i] = "" ;
            }
          else if( strncmp( argv[i] ,"-q" ,3 ) == 0 )
            {
              verbose = false ;
              argv[i] = "" ;
            }
          else
            {
              break ;
            }
        }
    }
  return ;
}

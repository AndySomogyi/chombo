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

#include "DebugDump.H"
#include "QuadCFInterp.H"
#include "LayoutIterator.H"
#include "LoadBalance.H"
#include "UsingNamespace.H"

/*****************/
/*****************/
void setHierarchy( int& nref,
                   Vector<Box>& coarboxes,
                   Vector<Box>& fineboxes,
                   Box& domc, Box& domf);
int quadCFInterpTest();
Real getQuadValue(const IntVect& index,Real dx, int idir);
void parseTestOptions(int argc ,char *argv[]) ;
/*****************/
/*****************/

/// Global variables for handling output:
static const char *pgmname = "quadCFInterpTest" ;
static const char *indent = "   " ,*indent2 = "      " ;
static bool verbose = true ;
/*****************/
/*****************/

int
main(int argc, char *argv[])
{
#ifdef CH_MPI
  MPI_Init(&argc, &argv);
#endif
  parseTestOptions( argc ,argv ) ;
  if( verbose )
    cout << indent2 << "Beginning " << pgmname << " ..." << endl ;

  int eekflag = quadCFInterpTest();
  if(eekflag != 0)
    {
      cout << indent << pgmname << ": failed with error code: "
           << eekflag << endl;
    }
  else
    {
      cout << indent << pgmname << ": passed. " << endl;
    }
#ifdef CH_MPI
  MPI_Finalize();
#endif
  return eekflag;
}

/*****************/
/*****************/

int
quadCFInterpTest()
{
  int nref;
  Vector<Box> fineboxes;
  Vector<Box> coarboxes;
  Box domf, domc;

  //set the grids and refinement ratio
  setHierarchy(nref, coarboxes,  fineboxes, domc, domf);
  int retflag = 0;

  Interval interv(0,0);
  Vector<int> procAssignCoar(coarboxes.size(), 0);
  Vector<int> procAssignFine(fineboxes.size(), 0);
  LoadBalance(procAssignCoar, coarboxes);
  LoadBalance(procAssignFine, fineboxes);
  DisjointBoxLayout dblFine,dblCoar;
  dblCoar.define(coarboxes, procAssignCoar);
  dblFine.define(fineboxes, procAssignFine);

  dblCoar.close();
  dblFine.close();
  if(verbose)
    {
      cout << indent2 << "dblCoar:" << endl;
      dumpDBL(&dblCoar);
      cout << indent2 << "dblFine:" << endl;
      dumpDBL(&dblFine);
    }

  Real dxf = 3.14159;
  Real dxc = nref*dxf;

  LevelData<FArrayBox>
    coarData(dblCoar, 1, IntVect::Unit),
    fineData(dblFine, 1, IntVect::Unit);

  DataIterator coarIt = coarData.dataIterator();
  DataIterator fineIt = fineData.dataIterator();

  QuadCFInterp QCFI(dblFine, &dblCoar, dxf, nref, 1, domf);

  for(int idir = 0; idir < SpaceDim; idir++)
    {
      //make coarse data = quadratic everywhere
      for(coarIt.reset(); coarIt.ok(); ++coarIt)
        {
          FArrayBox& coarFab=coarData[coarIt()];
          BoxIterator bit(coarFab.box());
          for(bit.reset(); bit.ok(); ++bit)
            {
              coarFab(bit(), 0) = getQuadValue(bit(), dxc, idir);
            }
        }

      //make fine data = quadratic
      //everywhere EXCEPT at ghost cells
      for(fineIt.reset(); fineIt.ok(); ++fineIt)
        {
          FArrayBox&  fineFab = fineData[fineIt()];
          fineFab.setVal(0.);
          Box fabIntBox= dblFine.get(fineIt());
          BoxIterator bit(fabIntBox);
          for(bit.reset(); bit.ok(); ++bit)
            {
              fineFab(bit(), 0) = getQuadValue(bit(), dxf, idir);
            }
        }

      //exchange data where ghost data overlaps valid data
      fineData.exchange(interv);

      //here is the big coarse-fine interpolation call
      QCFI.coarseFineInterp(fineData, coarData);

      for(fineIt.reset(); fineIt.ok(); ++fineIt)
        {
          FArrayBox&  fineFab = fineData[fineIt()];
          Box iterBox= dblFine.get(fineIt());
          iterBox.grow(idir,1);
          iterBox &= domf;
          BoxIterator bit(iterBox);
          for(bit.reset(); bit.ok(); ++bit)
            {
              Real rightAns = getQuadValue(bit(), dxf, idir);
              Real fabVal = fineFab(bit(), 0);
              Real diff =  Abs(fabVal - rightAns) / rightAns ;
              Real eps = 1.0e-12;
#ifdef CH_USE_FLOAT
              eps = sqrt(eps);
#endif
              if(diff > eps)
                {
                  cout << indent
                       << "QCFI failed  at"
                       << "   idir = " << idir
                       << ",  iv   = " << bit() << endl;
                  if(verbose)
                    {
                      cout << indent2
                           << "RightAns = " << rightAns
                           << ",  fabval = " << fabVal
                           << ",  diff   = " << diff
                           << endl;
                    }
                  retflag = 1;
                }
              //if(retflag > 0) break;
            }
        }
    }
  return retflag;
}

Real
getQuadValue(const IntVect& index, Real dx, int idir)
{
  Real Acoeff = 2.23;
  Real Bcoeff = -7.389;
  Real Ccoeff = 3.98;

  Real x = dx*(index[idir] + 0.5);
  Real val = Acoeff*x*x + Bcoeff*x + Ccoeff;
  return val;
}

void
setHierarchy( int& nref,
              Vector<Box>& coarboxes,
              Vector<Box>& fineboxes,
              Box& domc, Box& domf)
{
  nref = 2;
  int nc = 16;
  int nc16 = nc/2;
  int nc8 = nc/4;
  int nc24 = 3*nc/4;
  int nc12 = 3*nc/8;
  int nc20 = 5*nc/8;

  IntVect ivclo = IntVect::Zero;

  fineboxes.resize(0);

#if(CH_SPACEDIM ==1)

  IntVect ivchi(nc-1);

  fineboxes.resize(3);
  fineboxes[0].define(IntVect(0), IntVect(nc8-1));
  fineboxes[1].define(IntVect(nc12), IntVect(nc20-1));
  fineboxes[2].define(IntVect(nc24), ivchi);

  // this is to prevent unused variable warning in 1d
  int temp =nc16;
  nc16 = temp;

#elif(CH_SPACEDIM ==2)

  IntVect ivchi(nc-1, nc-1);
  Box boxf1(IntVect(0 ,nc8), IntVect(nc16-1,nc24));
  Box boxf2(IntVect(nc16,nc12), IntVect(nc24-1,nc20-1));
  Box boxf3(IntVect(nc8 ,0 ), IntVect(nc20-1,nc8-1));
  Box boxf4(IntVect(nc20,nc8 ), IntVect(nc-1,nc12-1));

  fineboxes.resize(4);
  fineboxes[0] = boxf1;
  fineboxes[1] = boxf2;
  fineboxes[2] = boxf3;
  fineboxes[3] = boxf4;  

#elif(CH_SPACEDIM==3)
  IntVect ivchi(nc-1, nc-1, nc-1);
  Box boxf1(IntVect(0 ,nc8,0), IntVect(nc16-1,nc24,nc-1));
  Box boxf2(IntVect(nc16,nc12,0), IntVect(nc24-1,nc20-1,nc-1));
  Box boxf3(IntVect(nc8 ,0 ,0), IntVect(nc20-1,nc8-1,nc-1));
  Box boxf4(IntVect(nc20,nc8,0 ), IntVect(nc-1,nc12-1,nc-1));

  fineboxes.resize(4);
  fineboxes[0] = boxf1;
  fineboxes[1] = boxf2;
  fineboxes[2] = boxf3;
  fineboxes[3] = boxf4;  

#else
  error bogus spacedim;
#endif
  
  for (int i=0; i<fineboxes.size(); ++i)
    {
      fineboxes[i].refine(nref);
    }

  domc = Box(ivclo, ivchi);
  int ichop = domc.size(0)/2;
  Box bc1 = domc;
  Box bc2 = bc1.chop(0, ichop);
#if (CH_SPACEDIM > 1)
  Box bc3 = bc2.chop(1, ichop);
  coarboxes.push_back(bc3);
#endif
  coarboxes.push_back(bc1);
  coarboxes.push_back(bc2);


  domf = refine(domc,nref);
  return ;
}

///
// Parse the standard test options (-v -q) out of the command line.
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
        }
    }
  return ;
}

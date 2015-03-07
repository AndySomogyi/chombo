#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

// Test program for AMR/LoadBalance
// The program simulates different processor counts so
// it should not be run under MPI.

// Test 1: single level, small number of loads.
// Test 2: single level, large number of loads.
// Test 3: single level, loads from Brian's HDF test code

#include <iostream>
using std::cout;
#include "SPMD.H"  // to get num_procs global variable
#include "LoadBalance.H"
#include "Misc.H"
#include "parstream.H"
#include "UsingNamespace.H"

/// Prototypes:
int
testLB1(void);
int
testLB2(void);
int
testLB3(void);

using std::endl;

void
parseTestOptions(int argc ,char* argv[]) ;

/// Global variables for handling output:
static const char* pgmname = "testLoadBalance" ;
static const char* indent = "   ";
static const char* indent2 = "      " ;
static bool verbose = true ;

int
main(int argc ,char* argv[])
{
  int stat_all = 0 ;
#ifdef CH_MPI
  MPI_Init(&argc, &argv);
#endif

  parseTestOptions( argc ,argv ) ;
  if( verbose )
    pout() << indent2 << "Beginning " << pgmname << endl ;

  int status = testLB1();

  if( status == 0 ){
    if( verbose ) pout() << indent << pgmname << " passed test 1." << endl ;
  }else{
    pout() << indent << pgmname << " failed test 1 with return code " << status << endl ;
    stat_all = status ;
  }

  status = testLB2();

  if( status == 0 ){
    if( verbose ) pout() << indent << pgmname << " passed test 2." << endl ;
  }else{
    pout() << indent << pgmname << " failed test 2 with return code " << status << endl ;
    stat_all = status ;
  }

  status = testLB3();

  if( status == 0 ){
    if( verbose ) pout() << indent << pgmname << " passed test 3." << endl ;
  }else{
    pout() << indent << pgmname << " failed test 3 with return code " << status << endl ;
    stat_all = status ;
  }

  if( stat_all == 0 )
    pout() << indent << pgmname << ": passed all tests." << endl ;
  else
    pout() << indent << pgmname << ": failed one or more tests." << endl ;

#ifdef CH_MPI
  MPI_Finalize();
#endif

  return stat_all ;
}

int
testLB1()
{
  const int numLoads = 9 ;
  Vector<Vector<long> > loads( 1 ) ;
  Real total_load = 0.0 ;
  for( int i=0 ; i<numLoads ; ++i )
    {
      loads[0].push_back( (i+1) ) ;
      total_load += (i+1) ;
    }
  Vector<Vector<Box> > grids( 1 ) ;
  grids[0].resize( numLoads ) ;
  Vector<int> refratios( 1,2 ) ;
  Vector<Vector<int> > assignments( 1 ) ;
  Real eff_ratio ;
  int status ;

  status = LoadBalance( assignments ,eff_ratio ,grids ,loads ,refratios,1 ) ;
  if( status != 0 ) return status ;
  if( verbose ) {
    pout() << indent2 << "test1: 1 processors, total load " << total_load
         << ", efficiency " << (int)eff_ratio*100.0 << "%, assignments are" << endl ;
    pout() << indent2 ;
    for( int i=0 ; i<numLoads ; ++i )
      {
        pout() << i << "{" << loads[0][i] << "}->" << assignments[0][i] << " " ;
      }
    pout() << endl ;
  }
  if( eff_ratio != 1.0 )
    {
      if( verbose )
        pout() << indent2 << "test1: 1 processor, expected efficiency 100%, got "
             << eff_ratio*100 << endl ;
      return -11 ;
    }
  for( int i=0 ; i<numLoads ; ++i )
    {
      if( assignments[0][i] != 0 )
        {
          if( verbose )
            pout() << indent2 << "test1: 1 processor, expected assignment 0, got "
                 << assignments[0][i] << ", for load " << i << endl ;
          return -12 ;
        }
    }

  status = LoadBalance( assignments ,eff_ratio ,grids ,loads ,refratios, 2 ) ;
  if( status != 0 ) return status ;
  if( verbose ) {
    pout() << indent2 << "test1: 2 processors, total load " << total_load
         << ", efficiency " << (int)(eff_ratio*100.0) << "%, assignments are" << endl ;
    pout() << indent2 ;
    for( int i=0 ; i<numLoads ; ++i )
      {
        pout() << i << "{" << loads[0][i] << "}->" << assignments[0][i] << " " ;
      }
    pout() << endl ;
  }
  if( eff_ratio < 0.95 )
    {
      if( verbose )
        pout() << indent2 << "test1: 2 processor, expected efficiency >95%, got "
             << eff_ratio*100 << endl ;
      return -21 ;
    }
  for( int i=0 ; i<numLoads ; ++i )
    {
      if( assignments[0][i] < 0 || assignments[0][i] > 1 )
        {
          if( verbose )
            pout() << indent2 << "test1: 2 processor, got assignment " << assignments[0][i] << ", for load " << i << endl ;
          return -22 ;
        }
    }

  status = LoadBalance( assignments ,eff_ratio ,grids ,loads ,refratios, 4 ) ;
  if( verbose ) {
    pout() << indent2 << "test1: 4 processors, total load " << total_load
         << ", efficiency " << (int)(eff_ratio*100.0) << "%, assignments are" << endl ;
    pout() << indent2 ;
    for( int i=0 ; i<numLoads ; ++i )
      {
        pout() << i << "{" << loads[0][i] << "}->" << assignments[0][i] << " " ;
      }
    pout() << endl ;
  }
  if( eff_ratio < 0.80 )
    {
      if( verbose )
        pout() << indent2 << "test1: 4 processor, expected efficiency >80%, got "
             << eff_ratio*100 << endl ;
      return -41 ;
    }
  for( int i=0 ; i<numLoads ; ++i )
    {
      if( assignments[0][i] < 0 || assignments[0][i] > 3 )
        {
          if( verbose )
            pout() << indent2 << "test1: 4 processor: got assignment " << assignments[0][i] << ", for load " << i << endl ;
          return -42 ;
        }
    }

  return status ;
}

int
testLB2()
{
  const int numLoads = 8 ;
  Vector<Vector<long> > loads( 1 ) ;
  Real total_load = 0.0 ;
  loads[0].push_back( 7 ) ; total_load += 7 ;
  loads[0].push_back( 6 ) ; total_load += 6 ;
  loads[0].push_back( 5 ) ; total_load += 5 ;
  loads[0].push_back( 4 ) ; total_load += 4 ;
  loads[0].push_back( 4 ) ; total_load += 4 ;
  loads[0].push_back( 3 ) ; total_load += 3 ;
  loads[0].push_back( 3 ) ; total_load += 3 ;
  loads[0].push_back( 2 ) ; total_load += 2 ;

  Vector<Vector<Box> > grids( 1 ) ;
  grids[0].resize( numLoads ) ;
  Vector<int> refratios( 1,2 ) ;
  Vector<Vector<int> > assignments( 1 ) ;
  Real eff_ratio ;
  int status ;

  status = LoadBalance( assignments ,eff_ratio ,grids ,loads ,refratios, 3 ) ;
  if( status != 0 ) return status ;
  if( verbose ) {
    pout() << indent2 << "test2: 3 processors, total load " << total_load
         << ", efficiency " << (int)(eff_ratio*100.0) << "%, assignments are" << endl ;
    pout() << indent2 ;
    for( int i=0 ; i<numLoads ; ++i )
      {
        pout() << i << "{" << loads[0][i] << "}->" << assignments[0][i] << " " ;
      }
    pout() << endl ;
  }
  Real correct_eff = (Real)11 / (Real)12 ;
  if( Abs(eff_ratio - correct_eff) > (Real)1e-5 )
    {
      if( verbose )
        pout() << indent2 << "test2: 3 processor, expected efficiency " << correct_eff
             << ", got " << eff_ratio*100 << endl ;
      return -221 ;
    }
  for( int i=0 ; i<numLoads ; ++i )
    {
      if( assignments[0][i] < 0 || assignments[0][i] > 2 )
        {
          if( verbose )
            pout() << indent2 << "test2: 3 processor, assignment " << assignments[0][i] << ", for load " << i
                 << " is out of range" << endl ;
          return -222 ;
        }
    }

  status = LoadBalance( assignments ,eff_ratio ,grids ,loads ,refratios, 4 ) ;
  if( verbose ) {
    pout() << indent2 << "test2: 4 processors, total load " << total_load
         << ", efficiency " << (int)(eff_ratio*100.0) << "%, assignments are" << endl ;
    pout() << indent2 ;
    for( int i=0 ; i<numLoads ; ++i )
      {
        pout() << i << "{" << loads[0][i] << "}->" << assignments[0][i] << " " ;
      }
    pout() << endl ;
  }
  if( eff_ratio < 0.80 )
    {
      if( verbose )
        pout() << indent2 << "test2: 4 processor, expected efficiency >80%, got "
             << eff_ratio*100 << endl ;
      return -241 ;
    }
  for( int i=0 ; i<numLoads ; ++i )
    {
      if( assignments[0][i] < 0 || assignments[0][i] > 3 )
        {
          if( verbose )
            pout() << indent2 << "test2: 4 processor: got assignment " << assignments[0][i] << ", for load " << i << endl ;
          return -242 ;
        }
    }

  return status ;
}

int
testLB3()
{
  const int numLoads = 16 ;
  Vector<Vector<long> > loads( 1 ) ;
  Real total_load = 0.0 ;
  loads[0].push_back( 1 ) ; total_load += 1 ; //box ((4,10) (4,10) (0,0))
  loads[0].push_back( 1 ) ; total_load += 1 ; //box ((5,7) (5,7) (0,0))
  loads[0].push_back( 1 ) ; total_load += 1 ; //box ((5,13) (5,13) (0,0))
  loads[0].push_back( 1 ) ; total_load += 1 ; //box ((6,6) (6,6) (0,0))
  loads[0].push_back( 1 ) ; total_load += 1 ; //box ((6,14) (6,14) (0,0))
  loads[0].push_back( 3 ) ; total_load += 3 ; //box ((7,5) (9,5) (0,0))
  loads[0].push_back( 3 ) ; total_load += 3 ; //box ((7,15) (9,15) (0,0))
  loads[0].push_back( 1 ) ; total_load += 1 ; //box ((10,4) (10,4) (0,0))
  loads[0].push_back( 8 ) ; total_load += 8 ; //box ((10,15) (13,16) (0,0))
  loads[0].push_back( 2 ) ; total_load += 2 ; //box ((11,5) (12,5) (0,0))
  loads[0].push_back( 1 ) ; total_load += 1 ; //box ((13,5) (13,5) (0,0))
  loads[0].push_back( 1 ) ; total_load += 1 ; //box ((14,6) (14,6) (0,0))
  loads[0].push_back( 1 ) ; total_load += 1 ; //box ((14,14) (14,14) (0,0))
  loads[0].push_back( 1 ) ; total_load += 1 ; //box ((15,7) (15,7) (0,0))
  loads[0].push_back( 1 ) ; total_load += 1 ; //box ((15,13) (15,13) (0,0))
  loads[0].push_back( 1 ) ; total_load += 1 ; //box ((16,10) (16,10) (0,0))

  Vector<Vector<Box> > grids( 1 ) ;
  grids[0].resize( numLoads ) ;
  Vector<int> refratios( 1,2 ) ;
  Vector<Vector<int> > assignments( 1 ) ;
  Real eff_ratio ;
  int status ;

  status = LoadBalance( assignments ,eff_ratio ,grids ,loads ,refratios, 3 ) ;
  if( status != 0 ) return status ;
  if( verbose ) {
    pout() << indent2 << "test3: 3 processors, total load " << total_load
         << ", efficiency " << (int)(eff_ratio*100.0) << "%, assignments are" << endl ;
    pout() << indent2 ;
    for( int i=0 ; i<numLoads ; ++i )
      {
        pout() << i << "{" << loads[0][i] << "}->" << assignments[0][i] << " " ;
      }
    pout() << endl ;
  }
  Real correct_eff = (Real)9 / (Real)10 ;
  if( Abs(eff_ratio - correct_eff) > (Real)1e-5 )
    {
      if( verbose )
        pout() << indent2 << "test3: 3 processor, expected efficiency " << correct_eff
             << ", got " << eff_ratio*100 << endl ;
      return -321 ;
    }
  for( int i=0 ; i<numLoads ; ++i )
    {
      if( assignments[0][i] < 0 || assignments[0][i] > 3 )
        {
          if( verbose )
            pout() << indent2 << "test3: 3 processor, assignment " << assignments[0][i] << ", for load " << i
                 << " is out of range" << endl ;
          return -322 ;
        }
    }

  return 0 ;
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

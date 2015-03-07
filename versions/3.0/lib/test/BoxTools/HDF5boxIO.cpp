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
using std::sqrt;
#include <iostream>
using std::endl;

#include "DisjointBoxLayout.H"
#include "LayoutIterator.H"
#include "CH_HDF5.H"
#include "BRMeshRefine.H"
#include "IntVectSet.H"
#include "LoadBalance.H"
#include "UsingNamespace.H"

/// Prototypes:
void
setCircleTags(IntVectSet& ivs, int circleR, int thickness = 1,
              const IntVect& center= IntVect::Zero);
void
buildDisjointBoxLayout(DisjointBoxLayout& plan,
                       const IntVectSet& tags, const Box& domain);

void
parseTestOptions( int argc ,char* argv[] ) ;

int
test();

/// Global variables for handling output:
static const char *pgmname = "HDF5boxIO" ;
static const char *indent = "   ", *indent2 = "      " ;
static bool verbose = true ;

/// Code:

int
main(int argc, char* argv[])
{
#ifdef CH_MPI
  MPI_Init(&argc, &argv);
#endif
  parseTestOptions( argc ,argv ) ;

  if( verbose )
    pout() << indent2 << "Beginning " << pgmname << " ..." << endl ;

  ///
  // Run the tests
  ///
  int icode = test();
  if(icode != 0)
    {
      pout() << indent << pgmname <<" failed"<<endl;
    }
  else
    {
      pout() << indent << pgmname <<" passed"<<endl;
    }
#ifdef CH_MPI
  MPI_Finalize();
#endif
  return icode;
}

// returns 0 on all tests passed.

int test()
{
#ifdef CH_USE_HDF5
  int error;
  HDF5Handle testFile;

  CH_assert(!testFile.isOpen());

  error = testFile.open("boxIO.h5", HDF5Handle::CREATE);
  if(error != 0)
    {
      if( verbose )
        pout() << indent2 << "File creation failed "<<error<<endl;
      return error;
    }

 CH_assert(testFile.isOpen());

  DisjointBoxLayout plan1, plan2;
  {
    Box domain(IntVect::Zero, 20*IntVect::Unit);

    IntVectSet tags;

    IntVect center = 10*IntVect::Unit;

    setCircleTags(tags, 6, 3, center); //circle/sphere

    buildDisjointBoxLayout(plan1, tags, domain);

    tags.makeEmpty();

    setCircleTags(tags, 5, 2, center);

    buildDisjointBoxLayout(plan2, tags, domain);
  }

  testFile.setGroupToLevel(0);
  error = write(testFile, plan1);
   if(error != 0)
    {
      if( verbose )
        pout() << indent2 << "box write failed "<<error<<endl;
      return error;
    }

  testFile.setGroupToLevel(1);

  error = write(testFile, plan2);
    if(error != 0)
    {
      if( verbose )
        pout() << indent2 << "box write failed "<<error<<endl;
      return error;
    }

  testFile.close();

 CH_assert(!testFile.isOpen());

  //====================================================================

  DisjointBoxLayout readplan1, readplan2;

  testFile.open("boxIO.h5", HDF5Handle::OPEN_RDONLY);

  testFile.setGroupToLevel(1);
  Vector<Box> boxes;
  Vector<int> assignments;
  error = read(testFile, boxes);
  if(error != 0)
    {
      if( verbose )
        pout() << indent2 << "box read failed "<<error<<endl;
      return error;
    }
  LoadBalance(assignments, boxes);
  readplan2.define(boxes, assignments);

  testFile.setGroupToLevel(0);

  error = read(testFile, boxes);
   if(error != 0)
    {
      if( verbose )
        pout() << indent2 << "box read failed "<<error<<endl;
      return error;
    }
   LoadBalance(assignments, boxes);
   readplan1.define(boxes, assignments);

   // now the test:
   LayoutIterator p1 = plan1.layoutIterator();
   LayoutIterator rp1 = readplan1.layoutIterator();
   LayoutIterator p2 = plan2.layoutIterator();
   LayoutIterator rp2 = readplan2.layoutIterator();

   if(plan1.size() != readplan1.size())
     {
       if( verbose )
         pout() << indent2 << "plan1 size different on read "<<endl;
       return 1;
     }

   for(;p1.ok(); ++p1, ++rp1)
     {
       if(plan1.get(p1()) != readplan1.get(rp1()))
         {
           if( verbose )
             {
               pout() << indent2 << "plan1 != readplan1  on read \n";
               pout()<<plan1<<"\n\n"<<readplan1<<endl;
             }
           return 2;
         }

     }

   if(plan2.size() != readplan2.size())
     {
       if( verbose )
         pout() << indent2 << "plan2 size different on read "<<endl;
       return 1;
     }

   for(;p2.ok(); ++p2, ++rp2)
     {
       if(plan2.get(p2()) != readplan2.get(rp2()))
         {
           if( verbose )
             {
               pout() << indent2 << "plan2 != readplan2  on read \n";
               pout()<<plan2<<"\n\n"<<readplan2<<endl;
             }
           return 2;

         }
     }

   testFile.close();
#endif
   return 0;

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
      else if( strncmp( argv[i] ,"-q" ,3 ) == 0 )
      {
        verbose = false ;
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
setCircleTags(IntVectSet& ivs, int circleR, int thickness,
              const IntVect& center)
{
 CH_assert(circleR > 0);
 CH_assert(thickness < circleR && thickness > 0);
  for(int t = 0; t<thickness; ++t){
    for(int x=-circleR-t; x<=circleR+t; ++x)
    {
      // not perfect, but I'm not fussy here
#if CH_SPACEDIM == 1
      //interval
      ivs |= center + IntVect(x);
#elif CH_SPACEDIM == 2
      // circle
      int Y = (int)sqrt((Real)abs(circleR*circleR - x*x));
      ivs |= center+IntVect(x,Y);
      ivs |= center+IntVect(x,-Y);
#elif CH_SPACEDIM == 3
      //sphere
      int Y = (int)sqrt((Real)abs(circleR*circleR - x*x));
      for(int y=-Y; y<=Y; ++y)
      {
        int Z = (int)sqrt((Real)abs(circleR*circleR - x*x - y*y));
        ivs |= center+IntVect(x,y,Z);
        ivs |= center+IntVect(x,y,-Z);
      }
#elif CH_SPACEDIM == 4
      //hypersphere
      int Y = (int)sqrt((Real)abs(circleR*circleR - x*x));
      for(int y=-Y; y<=Y; ++y)
      {
        int Z = (int)sqrt((Real)abs(circleR*circleR - x*x - y*y));
        for (int z=-Z; z<=Z; ++z)
          {
            int U = (int)sqrt((Real)abs(circleR*circleR -x*x -y*y -z*z));
            ivs |= center+IntVect(x,y,z,U);
            ivs |= center+IntVect(x,y,z,-U);
          }
      }
#elif CH_SPACEDIM == 5
      //hypersphere
      int Y = (int)sqrt((Real)abs(circleR*circleR - x*x));
      for(int y=-Y; y<=Y; ++y)
      {
        int Z = (int)sqrt((Real)abs(circleR*circleR - x*x - y*y));
        for (int z=-Z; z<=Z; ++z)
          {
            int U = (int)sqrt((Real)abs(circleR*circleR -x*x -y*y -z*z));
            for (int u=-U; u<=U; ++u)
              {
                int V = (int)sqrt((Real)abs(circleR*circleR -x*x -y*y -z*z -u*u));
                ivs |= center+IntVect(x,y,z,u,V);
                ivs |= center+IntVect(x,y,z,u,-V);
              }
          }
      }
#elif CH_SPACEDIM == 6
      //hypersphere
      int Y = (int)sqrt((Real)abs(circleR*circleR - x*x));
      for(int y=-Y; y<=Y; ++y)
      {
        int Z = (int)sqrt((Real)abs(circleR*circleR - x*x - y*y));
        for (int z=-Z; z<=Z; ++z)
          {
            int U = (int)sqrt((Real)abs(circleR*circleR -x*x -y*y -z*z));
            for (int u=-U; u<=U; ++u)
              {
                int V = (int)sqrt((Real)abs(circleR*circleR -x*x -y*y -z*z -u*u));
                for (int v=-V; v<=V; ++v)
                  {
                    int W = (int)sqrt((Real)abs(circleR*circleR -x*x -y*y -z*z -u*u -v*v));
                    ivs |= center+IntVect(x,y,z,u,v,-W);
                    ivs |= center+IntVect(x,y,z,u,v,-W);
                  }
              }
          }
      }
#else
#error HDF5boxIO failed: implemented only for 1-6D
#endif
    }
  }
}

void
buildDisjointBoxLayout(DisjointBoxLayout& plan,
                       const IntVectSet& tags, const Box& domain)
{
  Vector<Vector<Box> > vectBox(1);
  Vector<Vector<int> > assignments(1);
  Vector<Vector<long> > computeLoads(1);

  IntVectSet domainivs(domain);
  int maxsize = 128;
  Vector<int> fakeNRef(2,2);
  Real fillRatio = 0.5;
  int blockFactor = 1;
  int bufferSize = 1;

  // this is constructed with bogus values, since we'll be calling
  // the makeBoxes function directly from here
  BRMeshRefine mr(domain, fakeNRef, fillRatio, blockFactor, bufferSize, maxsize);
  mr.makeBoxes(vectBox[0], tags, domainivs, domain, maxsize, 0);

  assignments[0].resize(vectBox[0].size());
  computeLoads[0].resize(vectBox[0].size());

  for(int i=0; i<vectBox[0].size(); ++i)
    {
      computeLoads[0][i] = vectBox[0][i].numPts();
    }
  Real eff; Vector<int> refRatio(1,1);
  int stat = LoadBalance(assignments, eff, vectBox, computeLoads, refRatio);
  if( stat != 0 ) {
    MayDay::Error("loadBalance() FAILED",stat);
  }

  plan.define(vectBox[0], assignments[0]);

  plan.close();
}

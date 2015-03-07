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
using std::endl;

#include "BRMeshRefine.H"
#include "LoadBalance.H"
#include "CH_HDF5.H"
#include "parstream.H"
#include "BoxIterator.H"
#include "FABView.H"

#include "NewPoissonOp.H"
#include "AMRPoissonOp.H"
#include "BCFunc.H"

#include "UsingNamespace.H"

/// Global variables for handling output:
static const char* pgmname = "testAMRPoissonOp" ;
static const char* indent = "   ";
static const char* indent2 = "      " ;
static bool verbose = true ;

/// Application-specific global variables:
static Box domain = Box(IntVect(D_DECL(0,0,0)), IntVect(D_DECL(63,63,63)));
static Real dx = 1.0/64;
static int blockingFactor = 8;

///
// Parse the standard test options (-v -q -h) and
// app-specific options (-S <domain_size>) out of the command line
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
          else if( strncmp( argv[i] ,"-S" ,2 ) == 0 )
            {
              // see if the option has the size value attached
              if( strlen( argv[i] ) > 2 )
                {
                  // skip the "-S" and use whatever is left
                  argv[i] += 2 ;
                }
              else
                {
                  // clear the option arg and move to the next one
                  // argv[i++] = "" ;
                }
              int size = atoi( argv[i] );
              domain = Box(IntVect(D_DECL(-size,-size,-size)), IntVect(D_DECL(size-1,size-1,size-1)));
              dx = (Real)1.6 / (Real)(size*2) ;
              // argv[i] = "" ;
            }
          else if( strncmp( argv[i] ,"-h" ,3 ) == 0 )
            {
              pout() << "usage: " << pgmname << " [-hqv] [-S domain_size]" << std::endl ;
              exit( 99 ) ;
            }
        }
    }
  return ;
}

int
testPoissonOp();

int
testAMRPoissonOp();

int
main(int argc ,char* argv[])
{
#ifdef CH_MPI
  MPI_Init (&argc, &argv);
#endif
  parseTestOptions( argc ,argv ) ;
  if( verbose )
    pout () << indent2 << "Beginning " << pgmname << " ..." << endl ;

  int overallStatus = 0;
  int status = testPoissonOp();

  if( status == 0 )
  {
    pout() << indent << "PoissonOp" << " passed." << endl ;
  }
  else
  {
    overallStatus = 1;
    pout() << indent << "PoissonOp" << " failed with return code " << status << endl ;
  }

  status = testAMRPoissonOp();

  if( status == 0 )
  {
    pout() << indent << pgmname << " passed." << endl ;
  }
  else
  {
    overallStatus = 1;
    pout() << indent << pgmname << " failed with return code " << status << endl ;
  }


#ifdef CH_MPI
  MPI_Finalize ();
#endif
  return overallStatus;
}


//   u = x*x + y*y + z*z
//   du/dx = 2*x
//   du/dy = 2*y
//   du/dz = 2*z
//   Laplace(u) = 2*CH_SPACEDIM


extern "C"
{

  void Parabola_neum(Real* pos,
                     int* dir,
                     Side::LoHiSide* side,
                     Real* a_values)
  {
    switch (*dir)
      {
      case 0:
        a_values[0]=2*pos[0];
        return;
      case 1:
        a_values[0]=2*pos[1];
        return;
      case 2:
        a_values[0]=2*pos[2];
        return;
      default:
        MayDay::Error("no such dimension");
      };
  }
  void Parabola_diri(Real* pos,
                     int* dir,
                     Side::LoHiSide* side,
                     Real* a_values)
  {
    a_values[0] = D_TERM(pos[0]*pos[0],+pos[1]*pos[1],+pos[2]*pos[2]);
  }

  //this sets ghost cells in domain or out.
  //this makes zero sense outside of test environment
  void DirParabolaBC(FArrayBox& a_state,
                     const Box& valid,
                     const ProblemDomain& a_domain,
                     Real a_dx,
                     bool a_homogeneous)
  {
    for(int i=0; i<CH_SPACEDIM; ++i)
      {
        for(SideIterator sit; sit.ok(); ++sit)
          {
            DiriBC(a_state,
                   valid,
                   dx,
                   a_homogeneous,
                   Parabola_diri,
                   i,
                   sit(),2);
          }
      }
  }

  //this sets ghost cells in domain or out.
  //this makes zero sense outside of test environment
  void NeumParabolaBC(FArrayBox& a_state,
                      const Box& valid,
                      const ProblemDomain& a_domain,
                      Real a_dx,
                      bool a_homogeneous)
  {
    for(int i=0; i<CH_SPACEDIM; ++i)
      {
        for(SideIterator sit; sit.ok(); ++sit)
          {

            NeumBC(a_state,
                   valid,
                   dx,
                   a_homogeneous,
                   Parabola_neum,
                   i,
                   sit());
          }
      }
  }

}

static BCValueFunc pointFunc = Parabola_diri;

void parabola(const Box& box, int comps, FArrayBox& t)
{
  RealVect pos;
  int dir;
  Side::LoHiSide side;
  int num = 1;
  ForAllXBNN(Real,t, box, 0, comps)
    {
      num=nR;
      D_TERM(pos[0]=dx*(iR+0.5);, pos[1]=dx*(jR+0.5);, pos[2]=dx*(kR+0.5));
      pointFunc(&(pos[0]), &dir, &side, &tR);
    }EndFor;
}

void makeGrids(DisjointBoxLayout& a_dbl, const Box& a_domain)
{

  BRMeshRefine br;
  Box domain = a_domain;
  domain.coarsen(blockingFactor);
  domain.refine(blockingFactor);
  CH_assert(domain == a_domain);
  domain.coarsen(blockingFactor);

  ProblemDomain junk(domain);
  IntVectSet pnd(domain);
  IntVectSet tags;
  for(BoxIterator bit(domain); bit.ok(); ++bit)
    {
      const IntVect& iv = bit();
      if(D_TERM(true, && iv[1]< 2*iv[0] && iv[1]>iv[0]/2, && iv[2] < domain.bigEnd(2)/2))
        {
          tags|=iv;
        }
    }
  Vector<Box> boxes;
  br.makeBoxes(boxes, tags, pnd, junk, 32/blockingFactor, 1);
  Vector<int> procs;
  LoadBalance(procs, boxes);
  for(int i=0; i<boxes.size(); ++i) boxes[i].refine(blockingFactor);
  a_dbl.define(boxes, procs);
}

int
testAMRPoissonOp()
{
#ifdef CH_USE_HDF5
  writeLevel(NULL);
#endif
  ProblemDomain regularDomain(domain);

  DisjointBoxLayout  dbl;

  makeGrids(dbl, domain);
  // make interesting grids here.

  dbl.close();

  DataIterator dit(dbl);
  LevelData<FArrayBox> phi(dbl, 1, IntVect::Unit);
  LevelData<FArrayBox> rhs(dbl, 1);
  LevelData<FArrayBox> lofphi(dbl, 1);

  struct setvalue{
    static void setFunc(const Box& box,
                        int comps, FArrayBox& t) {t.setVal(4.0);}
  };

  rhs.apply(setvalue::setFunc);
  phi.apply(parabola);

  RealVect pos(IntVect::Unit);
  pos*=dx;

  AMRPoissonOp amrop;

  amrop.define(dbl, pos[0], regularDomain, DirParabolaBC);

  amrop.applyOp(lofphi, phi);

  for(dit.begin(); dit.ok(); ++dit)
    {
      lofphi[dit]-=2*CH_SPACEDIM;
      Real norm0 = lofphi[dit].norm(0);
      Real norm1 = lofphi[dit].norm(1);
      Real norm2 = lofphi[dit].norm(2);
      if(norm0 > 1E-6)
        {
          pout()<< "applyOp for AMR parabolic dirichlet function seems busted on Box " << lofphi[dit].box()
                << ", max(L(phi))=" << norm0 <<std::endl;
          pout() << "applyOp: 1-norm, 2-norm = " << norm1 << " , " << norm2 << std::endl ;
          return 1;
        }
      else if(verbose)
      {
        pout() << indent << "applyOp for AMR parabolic dirichlet on " << lofphi[dit].box()
               << ": norms(0,1,2) = " << norm0 << " , " << norm1 << " , " << norm2 << std::endl ;
      }
    }

  amrop.define(dbl, pos[0], regularDomain, NeumParabolaBC);

  amrop.applyOp(lofphi, phi);

  for(dit.begin(); dit.ok(); ++dit)
    {
      lofphi[dit]-=2*CH_SPACEDIM;
      Real norm0 = lofphi[dit].norm(0);
      Real norm1 = lofphi[dit].norm(1);
      Real norm2 = lofphi[dit].norm(2);
      if(norm0 > 1E-6)
        {
          pout()<< "applyOp for AMR parabolic neumann function seems busted on Box " << lofphi[dit].box()
                << ", max(L(phi))=" << norm0 <<std::endl;
          pout() << "applyOp: 1-norm, 2-norm = " << norm1 << " , " << norm2 << std::endl ;
          return 1;
        }
      else if(verbose)
      {
        pout() << indent << "applyOp for AMR parabolic neumann on " << lofphi[dit].box()
               << ": norms(0,1,2) = " << norm0 << " , " << norm1 << " , " << norm2 << std::endl ;
      }
    }

  return 0;
}

int
testPoissonOp()
{
  if(verbose) pout() << indent << "solving on domain " << domain << " with dx = " << dx << std::endl ;

#ifdef CH_USE_HDF5
  writeLevel(NULL);
#endif
  ProblemDomain regularDomain(domain);

  Box phiBox = domain;
  phiBox.grow(1);

  FArrayBox phi(phiBox, 1);
  FArrayBox rhs(domain, 1);
  FArrayBox lofphi(domain, 1);

  rhs.setVal(2*CH_SPACEDIM);
  parabola(phiBox, 1, phi);

  RealVect pos(IntVect::Unit);
  pos*=dx;

  NewPoissonOp op;

  op.define(pos, regularDomain, DirParabolaBC);

  op.applyOp(lofphi, phi);

  lofphi-=2*CH_SPACEDIM;
  Real norm0 = lofphi.norm(0);
  Real norm1 = lofphi.norm(1);
  Real norm2 = lofphi.norm(2);

  if(norm0 > 1E-6)
    {
      pout() << "applyOp for parabolic dirichlet function seems busted, max(L(phi))=" << norm0
             <<std::endl;
      pout() << "applyOp: 1-norm, 2-norm = " << norm1 << " , " << norm2 << std::endl ;
      return 1;
    }
  else if(verbose)
  {
    pout() << indent << "applyOp for parabolic dirichlet: norms(0,1,2) = "
           << norm0 << " , " << norm1 << " , " << norm2 << std::endl ;
  }

  op.define(pos, regularDomain, NeumParabolaBC);

  op.applyOp(lofphi, phi);

  lofphi-=2*CH_SPACEDIM;
  norm0 = lofphi.norm(0);
  norm1 = lofphi.norm(1);
  norm2 = lofphi.norm(2);
  if(norm0 > 1E-6)
    {
      pout()<< "applyOp for parabolic neumann function seems busted, max(L(phi))=" << norm0
            <<std::endl;
      pout() << "applyOp: 1-norm, 2-norm = " << norm1 << " , " << norm2 << std::endl ;
      return 1;
    }
  else if(verbose)
  {
    pout() << indent << "applyOp for parabolic neumann: norms(0,1,2) = "
           << norm0 << " , " << norm1 << " , " << norm2 << std::endl ;
  }

  return 0;
}

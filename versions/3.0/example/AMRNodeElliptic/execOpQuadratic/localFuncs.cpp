#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

// localFuncs.cpp
// petermc, 24 March 2003

#include <string>
using std::string;

#include "localFuncs.H"
#include "NodePoissonBC.H"
#include "ParmParse.H"
#include "testProb_F.H"
using std::cerr;
using std::cout;
using std::endl;

// --------------------------------------------------------------
int
setDomainBC(DomainNodeBC& a_dombc,
            const bool a_verbose)
{
  /*
    Read input file and broadcast the vectors ibclo, ibchi.
  */
#ifdef CH_MPI
  MPI_Barrier(Chombo_MPI::comm);
#endif
  Vector<int> ibclo;
  Vector<int> ibchi;
  if (procID() == 0)
    {
      ParmParse pp("main");
      pp.getarr("bc_lo", ibclo, 0, SpaceDim);
      pp.getarr("bc_hi", ibchi, 0, SpaceDim);
    }
  broadcast(ibclo, 0);
  broadcast(ibchi, 0);

  /*
    On each side, set a_dombc to be an object of class
    DirichletNodeBC or NeumannNodeBC.
  */
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      SideIterator sit;
      for (sit.begin(); sit.ok(); sit.next())
        {
          // settings for lo or hi
          int ibc;
          string sideString;
          switch (sit())
            {
            case Side::Lo:
              {
                ibc = ibclo[idir];
                sideString = "lo";
                break;
              }
            case Side::Hi:
              {
                ibc = ibchi[idir];
                sideString = "hi";
                break;
              }
            default:
              {
                cerr << "localFuncs::setDomainBC(): bogus side" << endl;
                return(1);
              }
            }
          // now check bc_lo or bc_hi.
          string bcString;
          if (ibc == 0) // 0 for Dirichlet
            {
              bcString = "Dirichlet";
              DirichletNodeBC dirbc(idir, sit());
              a_dombc.setFaceNodeBC(dirbc);
            }
          else if (ibc == 1) // 1 for Neumann
            {
              bcString = "Neumann";
              NeumannNodeBC neumbc(idir, sit());
              a_dombc.setFaceNodeBC(neumbc);
            }
          else if (ibc == 2) // 2 for inhomogeneous Dirichlet
            {
              bcString = "inhomogeneous Dirichlet";
              DirichletNodeBC dirbc(idir, sit());
              dirbc.setInhomogeneous(true);
              a_dombc.setFaceNodeBC(dirbc);
            }
          else if (ibc == 3) // 3 for inhomogeneous Neumann
            {
              bcString = "inhomogeneous Neumann";
              NeumannNodeBC neumbc(idir, sit());
              neumbc.setInhomogeneous(true);
              a_dombc.setFaceNodeBC(neumbc);
            }
          else
            {
              if (a_verbose)
                cerr << "setDomainBC: bogus input bc_" << sideString
                     << " flag" << endl;
              return(1);
            }

          if (a_verbose && procID() == 0)
            cout << bcString << " boundary conditions in "
                 << sideString << " direction for side " << idir << endl;
        } // side iterator
    } // direction
  return(0);
}


// --------------------------------------------------------------
int
setPhiQuadratic(Vector<LevelData<NodeFArrayBox>* >& a_phi,
                const Vector<Real>&                 a_vectDx,
                const Vector<ProblemDomain>&        a_vectDomain,
                const int a_numlevels,
                const bool a_verbose)
{
  for (int ilev = 0; ilev < a_numlevels; ilev++)
    {
      Real dxlev = a_vectDx[ilev];
      // ProblemDomain domlev = a_vectDomain[ilev];
      // domlev.surroundingNodes();
      if (a_verbose)
        cout << "level " << ilev << " dx = " << dxlev << endl;
      LevelData<NodeFArrayBox>* phiThisLev = a_phi[ilev];
      for (DataIterator it(phiThisLev->dataIterator()); it.ok(); ++it)
        {
          FArrayBox& phiFab = phiThisLev->operator[](it()).getFab();
          // phiFab.box() is node-centered.

          FORT_GETPHINODEQUADRATIC(CHF_FRA(phiFab),
                                   CHF_BOX(phiFab.box()),
                                   CHF_CONST_REAL(dxlev));
        }
    }
  return(0);
}

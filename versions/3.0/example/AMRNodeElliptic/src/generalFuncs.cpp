#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

// generalFuncs.cpp
// petermc, 24 March 2003

#include <cmath>
#include <string>
using std::string;

#include "generalFuncs.H"
#include "BoxIterator.H"
#include "ParmParse.H"
#include "BRMeshRefine.H"
#include "SPMD.H"
#include "LoadBalance.H"
#include "NodeSetOperations.H"
using std::cout;
using std::cerr;
using std::endl;

// ---------------------------------------------------------
Box
boxFromVector(const Vector<int>&  a_ints,
              int                 a_start) // default a_start == 0
{
  IntVect lo(D_DECL(a_ints[a_start+0], a_ints[a_start+1], a_ints[a_start+2]));
  IntVect hi(D_DECL(a_ints[a_start+3], a_ints[a_start+4], a_ints[a_start+5]));
  Box bx(lo, hi);
  return bx;
}


// ---------------------------------------------------------
int
readGrids(Vector<DisjointBoxLayout>& a_vectGrids,
          Vector<Box>&               a_vectDomain,
          Vector<Real>&              a_vectDx,
          Vector<int>&               a_vectRefRatio,
          RealVect&                  a_probLo,
          RealVect&                  a_probHi,
          int& a_numlevels,
          const bool a_verbose,
          int a_refined)
{
  Vector<ProblemDomain> probdomain;
  int eekflag = readGrids(a_vectGrids, probdomain, a_vectDx, a_vectRefRatio,
                          a_probLo, a_probHi,
                          a_numlevels, a_verbose, a_refined);

  for (int ilev = 0; ilev < a_numlevels; ilev++)
    a_vectDomain[ilev] = probdomain[ilev].domainBox();

  return eekflag;
}


// ---------------------------------------------------------
int
readGrids(Vector<DisjointBoxLayout>& a_vectGrids,
          Vector<ProblemDomain>&     a_vectDomain,
          Vector<Real>&              a_vectDx,
          Vector<int>&               a_vectRefRatio,
          RealVect&                  a_probLo,
          RealVect&                  a_probHi,
          int& a_numlevels,
          const bool a_verbose,
          int a_refined)
{
#ifdef CH_MPI
  MPI_Barrier(Chombo_MPI::comm);
#endif
  Vector<int> ancells(SpaceDim);
  Vector<Real> prob_loa(SpaceDim);
  Vector<Real> prob_hia(SpaceDim);
  Vector<int> level1domainIndices, level1boxIndices;
  Vector<int> level2domainIndices, level2boxIndices;

  // variables with defaults

  int max_level = 1;
  // a_numlevels = max_level + 1;
  Real fillRat = 0.77;
  Vector<int> periodic(SpaceDim);
  for (int idir = 0; idir < SpaceDim; idir++) periodic[idir] = 0;
  int maxboxsize = 32;
  int level1boxcount = 0;
  int level2boxcount = 0;

  if(procID() == 0)
    {
      ParmParse pp("main");
      pp.query("max_level", max_level);
      a_numlevels = max_level + 1;
      pp.getarr("n_cell", ancells, 0, SpaceDim);
      pp.getarr("ref_ratio", a_vectRefRatio, 0, a_numlevels); // note a_numlevels here
      pp.getarr("prob_lo", prob_loa, 0, SpaceDim);
      pp.getarr("prob_hi", prob_hia, 0, SpaceDim);
      pp.query("fill_ratio", fillRat);
      pp.query("maxboxsize", maxboxsize);
      pp.queryarr("is_periodic", periodic, 0, SpaceDim);
      if (max_level >= 1)
        {
          if (pp.query("level_1_boxcount", level1boxcount))
            {
              pp.getarr("level_1_boxes", level1boxIndices,
                        0, 2*3*level1boxcount);
              pp.queryarr("level_1_domain", level1domainIndices,
                          0, 2*3);
            }
          if (max_level >= 2)
            {
              if (pp.query("level_2_boxcount", level2boxcount))
                {
                  pp.getarr("level_2_boxes", level2boxIndices,
                            0, 2*3*level2boxcount);
                  pp.queryarr("level_2_domain", level2domainIndices,
                              0, 2*3);
                }
            }
        }
    }
  broadcast(a_vectRefRatio, 0);
  broadcast(a_numlevels, 0);

  broadcast(ancells, 0);
  broadcast(prob_loa, 0);
  broadcast(prob_hia, 0);
  broadcast(level1domainIndices, 0);
  broadcast(level1boxIndices, 0);
  broadcast(level2domainIndices, 0);
  broadcast(level2boxIndices, 0);

  broadcast(max_level, 0);
  broadcast(fillRat, 0);
  broadcast(periodic, 0);
  broadcast(maxboxsize, 0);
  broadcast(level1boxcount, 0);
  broadcast(level2boxcount, 0);

  int eekflag = 0;

  a_probLo = RealVect(D_DECL(prob_loa[0], prob_loa[1], prob_loa[2]));
  a_probHi = RealVect(D_DECL(prob_hia[0], prob_hia[1], prob_hia[2]));

  // corners of basedom
  IntVect ivlo = IntVect::Zero;
  IntVect ivhi;
  for (int idir = 0; idir < SpaceDim; idir++)
    ivhi[idir] = a_refined * ancells[idir] - 1;

  // basedom is CELL-centered, [ivlo : ivhi]
  // in each dimension, [0 : ancells - 1]
  Box basedomBox = Box(ivlo, ivhi, IntVect::Zero);
  ProblemDomain basedom(basedomBox);
  for (int idir = 0; idir < SpaceDim; idir++)
    if (periodic[idir] == 1) basedom.setPeriodic(idir, true);

  // assign domains at all levels:  CELL-centered
  eekflag = assignDomains(a_vectDomain, a_vectRefRatio, basedom, a_numlevels);
  if (eekflag != 0)
    {
      cerr << "readGrids: assignDomains returned error code "
           << eekflag << endl;
      return(1);
    }

  int ncbase = a_refined * ancells[0]; // 32
  Vector< Vector<Box> > gridsAllLevels(a_numlevels);
  gridsAllLevels[0].push_back(basedomBox);

  int maxRefRatio = 2;
  for (int irat = 0; irat < a_vectRefRatio.size(); irat++)
    {
      maxRefRatio = Max(a_vectRefRatio[irat], maxRefRatio);
    }
  int blockFactor = maxRefRatio * 2;
  // If very small, then need smaller blockFactor.
  if (ncbase <= 16) blockFactor = maxRefRatio;

  // max_level == 0:  nothing to do.
  // max_level == 1:  the four familiar boxes.
  // setGridsFull:  at every level, cover the whole domain.
  // setGridsOne:  one grid in the center

  if (max_level >= 1)
    {
      eekflag = getGridsLevel(gridsAllLevels[1], a_vectDomain[1],
                              level1domainIndices, level1boxcount, level1boxIndices);
      if (eekflag != 0)
        {
          cerr << "getGridsLevel on level 1 returned error code "
               << eekflag << endl;
          return(1);
        }

      if (max_level >= 2)
        {
          eekflag = getGridsLevel(gridsAllLevels[2], a_vectDomain[2],
                                  level2domainIndices, level2boxcount, level2boxIndices);
          if (eekflag != 0)
            {
              cerr << "getGridsLevel on level 2 returned error code "
                   << eekflag << endl;
              return(1);
            }
        }
    }

  // set grids, splitting them up if necessary, and assign procs
  eekflag = assignAllGrids(a_vectGrids, a_vectDomain, gridsAllLevels,
                           maxboxsize, blockFactor);

  if (eekflag != 0)
    {
      cerr << "setGrids: assignAllGrids returned error code " << eekflag
           << endl;
      return(1);
    }

  // calculate grid spacing at each level
  Real dx0 = (prob_hia[0] - prob_loa[0]) / Real(ncbase);
  eekflag = assignDx(a_vectDx, a_vectRefRatio, dx0, a_numlevels);

  if (eekflag != 0)
    {
      cerr << "setGrids: assignDx returned error code " << eekflag
           << endl;
      return(1);
    }

  return(0);
}


// ---------------------------------------------------------
int
readDomain(Vector<ProblemDomain>&     a_vectDomain,
           Vector<Real>&              a_vectDx,
           Vector<int>&               a_vectRefRatio,
           RealVect&                  a_probLo,
           RealVect&                  a_probHi,
           int& a_numlevels,
           const bool a_verbose,
           int a_refined)
{
#ifdef CH_MPI
  MPI_Barrier(Chombo_MPI::comm);
#endif
  Vector<int> ancells(SpaceDim);
  Vector<Real> prob_loa(SpaceDim);
  Vector<Real> prob_hia(SpaceDim);

  // variables with defaults

  int max_level = 1;
  // a_numlevels = max_level + 1;
  Real fillRat = 0.77;
  Vector<int> periodic(SpaceDim);
  for (int idir = 0; idir < SpaceDim; idir++) periodic[idir] = 0;
  int maxboxsize = 32;

  if(procID() == 0)
    {
      ParmParse pp("main");
      pp.query("max_level", max_level);
      a_numlevels = max_level + 1;
      pp.getarr("n_cell", ancells, 0, SpaceDim);
      pp.getarr("ref_ratio", a_vectRefRatio, 0, a_numlevels); // note a_numlevels here
      pp.getarr("prob_lo", prob_loa, 0, SpaceDim);
      pp.getarr("prob_hi", prob_hia, 0, SpaceDim);
      pp.query("fill_ratio", fillRat);
      pp.query("maxboxsize", maxboxsize);
      pp.queryarr("is_periodic", periodic, 0, SpaceDim);
    }
  broadcast(a_vectRefRatio, 0);
  broadcast(a_numlevels, 0);

  broadcast(ancells, 0);
  broadcast(prob_loa, 0);
  broadcast(prob_hia, 0);

  broadcast(max_level, 0);
  broadcast(fillRat, 0);
  broadcast(periodic, 0);
  broadcast(maxboxsize, 0);

  int eekflag = 0;

  a_probLo = RealVect(D_DECL(prob_loa[0], prob_loa[1], prob_loa[2]));
  a_probHi = RealVect(D_DECL(prob_hia[0], prob_hia[1], prob_hia[2]));

  // corners of basedom
  IntVect ivlo = IntVect::Zero;
  IntVect ivhi;
  for (int idir = 0; idir < SpaceDim; idir++)
    ivhi[idir] = a_refined * ancells[idir] - 1;

  // basedom is CELL-centered, [ivlo : ivhi]
  // in each dimension, [0 : ancells - 1]
  Box basedomBox = Box(ivlo, ivhi, IntVect::Zero);
  ProblemDomain basedom(basedomBox);
  for (int idir = 0; idir < SpaceDim; idir++)
    if (periodic[idir] == 1) basedom.setPeriodic(idir, true);

  // assign domains at all levels:  CELL-centered
  eekflag = assignDomains(a_vectDomain, a_vectRefRatio, basedom, a_numlevels);
  if (eekflag != 0)
    {
      cerr << "readGrids: assignDomains returned error code "
           << eekflag << endl;
      return(1);
    }

  int ncbase = a_refined * ancells[0]; // 32

  // calculate grid spacing at each level
  Real dx0 = (prob_hia[0] - prob_loa[0]) / Real(ncbase);
  eekflag = assignDx(a_vectDx, a_vectRefRatio, dx0, a_numlevels);

  if (eekflag != 0)
    {
      cerr << "setGrids: assignDx returned error code " << eekflag
           << endl;
      return(1);
    }

  return(0);
}


// ---------------------------------------------------------
int
getGridsLevel(Vector<Box>&          a_grids,
              const ProblemDomain&  a_domain,
              const Vector<int>&    a_subdomainIndices,
              const int             a_boxCount,
              const Vector<int>&    a_boxIndices)
{
  Box domain = a_domain.domainBox();

  if (a_boxCount == 0)
    {
      // Make level 1 cover the whole domain.
      a_grids.push_back(domain);
    }
  else
    {
      // Figure out the coarsening factor.
      Box subdomain = boxFromVector(a_subdomainIndices);
      int multiplier = 1;
      Box expandedDomain(subdomain);

      Box refinedExpandedDomain = refine(expandedDomain, 2);
      // keep refining expandedDomain until it's just big enough to fit
      // inside domain.
      while (domain.contains(refinedExpandedDomain))
        {
          expandedDomain = refinedExpandedDomain;
          refinedExpandedDomain = refine(expandedDomain, 2);
          multiplier *= 2;
        }

      for (int indbox = 0; indbox < a_boxCount; indbox++)
        {
          Box baseBox = boxFromVector(a_boxIndices, 2*3*indbox);
          if (! subdomain.contains(baseBox))
            {
              cerr << "readGrids: level 1 box " << baseBox
                   << " not in " << subdomain << endl;
              return(1);
            }
          // Now scale thisBox up to size of the domain.
          Box thisBox = refine(baseBox, multiplier);

          a_grids.push_back(thisBox);
        }
    } // if (level1boxcount > 0)

  return(0);
}


// ---------------------------------------------------------
int
setGridsTwoLevel(const IntVect&             a_lengths,
                 const IntVect&             a_offsets,
                 const IntVect&             a_sublengths,
                 const int                  a_refRatio,
                 const Real                 a_dx,
                 Vector<DisjointBoxLayout>& a_vectGrids,
                 Vector<ProblemDomain>&     a_vectDomain,
                 Vector<Real>&              a_vectDx,
                 Vector<int>&               a_vectRefRatio,
                 const bool a_verbose)
{
#ifdef CH_MPI
  MPI_Barrier(Chombo_MPI::comm);
#endif
  int numLevels = 2;
  Vector<int> ancells(SpaceDim);
  Vector<Real>  prob_loa(SpaceDim);
  Vector<Real>  prob_hia(SpaceDim);

  int eekflag = 0;

  // corners of basedom
  IntVect ivlo = IntVect::Zero;
  IntVect ivhi = a_lengths - IntVect::Unit;

  // basedom is CELL-centered, [ivlo : ivhi]
  // in each dimension, [0 : ancells - 1]
  Box basedom = Box(ivlo, ivhi, IntVect::Zero);

  for (int ilev = 0; ilev < numLevels; ilev++)
    a_vectRefRatio.push_back(a_refRatio);

  // assign domains at all levels:  CELL-centered
  eekflag = assignDomains(a_vectDomain, a_vectRefRatio, basedom, numLevels);
  if (eekflag != 0)
    {
      cerr << "setGridsTwoLevel: assignDomains returned error code " << eekflag << endl;
      return(1);
    }

  // assign grids
  Vector<Vector<Box> > gridsAllLevels(numLevels);
  Box grid1(a_offsets, a_offsets + a_sublengths - IntVect::Unit);
  grid1.refine(a_refRatio);

  gridsAllLevels[0].push_back(a_vectDomain[0].domainBox());
  gridsAllLevels[1].push_back(grid1);

  // set grids and assign procs
  a_vectGrids.resize(numLevels);
  for (int ilev = 0; ilev < numLevels; ilev++)
    {
      eekflag = assignLevelGrids(a_vectGrids, a_vectDomain[ilev], gridsAllLevels[ilev], ilev);
      if (eekflag != 0)
        {
          cerr << "setGrids: assignLevelGrids returned error code " << eekflag
               << endl;
          return(1);
        }
    }

  // calculate grid spacing at each level
  eekflag = assignDx(a_vectDx, a_vectRefRatio, a_dx, numLevels);
  if (eekflag != 0)
    {
      cerr << "setGrids: assignDx returned error code " << eekflag
           << endl;
      return(1);
    }

  return(0);
}


// ---------------------------------------------------------
int
assignAllGrids(Vector<DisjointBoxLayout>&   a_vectGrids,
               const Vector<ProblemDomain>& a_vectDomain,
               const Vector<Vector<Box> >&  a_boxes,
               int a_maxboxsize, // if omitted then 0, meaning no maximum
               int a_blockFactor) // if omitted then 1
{
  int numlevels = a_vectDomain.size();
  CH_assert( a_boxes.size() == numlevels );
  // set grids, splitting them up if necessary, and assign procs
  a_vectGrids.resize(numlevels);
  for (int ilev = 0; ilev < numlevels; ilev++)
    {
      int eekflag = assignLevelGrids(a_vectGrids, a_vectDomain[ilev],
                                     a_boxes[ilev], ilev,
                                     a_maxboxsize, a_blockFactor);
      if (eekflag != 0)
        {
          cerr << "assignAllGrids: assignLevelGrids at level "
               << ilev << " returned error code " << eekflag
               << endl;
          return(1);
        }
    }
  return(0);
}


// ---------------------------------------------------------
int
assignLevelGrids(Vector<DisjointBoxLayout>&  a_vectGrids,
                 const ProblemDomain&        a_domain,
                 const Vector<Box>&          a_boxes,
                 int a_level,
                 int a_maxboxsize, // if omitted then 0, meaning no maximum
                 int a_blockFactor) // if omitted then 1
{
  Vector<Box> newboxes(0);
  if (a_maxboxsize == 0)
    {
      newboxes.append(a_boxes);
    }
  else
    {
      for (int gindex = 0; gindex < a_boxes.size(); gindex++)
        {
          Box biggerBox = a_boxes[gindex];
          Vector<Box> splitBoxes(0);
          domainSplit(biggerBox, splitBoxes, a_maxboxsize, a_blockFactor);
          newboxes.append(splitBoxes);
        }
    }

  Vector<int> procAssign;
  int eekflag = LoadBalance(procAssign, newboxes);

  if (eekflag != 0)
    {
      cerr << "setLevelGrids(): LoadBalance returned error code " << eekflag
           << endl;
      return(1);
    }

  // grids at level a_level
  a_vectGrids[a_level].define(newboxes, procAssign, a_domain);
  a_vectGrids[a_level].close();
  return(0);
}


// ---------------------------------------------------------
int
assignLevelGrids(DisjointBoxLayout&          a_vectGrids,
                 const ProblemDomain&        a_domain,
                 const Vector<Box>&          a_boxes,
                 int a_maxboxsize, // if omitted then 0, meaning no maximum
                 int a_blockFactor) // if omitted then 1
{
  Vector<Box> newboxes(0);
  if (a_maxboxsize == 0)
    {
      newboxes.append(a_boxes);
    }
  else
    {
      for (int gindex = 0; gindex < a_boxes.size(); gindex++)
        {
          Box biggerBox = a_boxes[gindex];
          Vector<Box> splitBoxes(0);
          domainSplit(biggerBox, splitBoxes, a_maxboxsize, a_blockFactor);
          newboxes.append(splitBoxes);
        }
    }

  Vector<int> procAssign;
  int eekflag = LoadBalance(procAssign, newboxes);

  if (eekflag != 0)
    {
      cerr << "setLevelGrids(): LoadBalance returned error code " << eekflag
           << endl;
      return(1);
    }

  // grids at level a_level
  a_vectGrids.define(newboxes, procAssign, a_domain);
  a_vectGrids.close();
  return(0);
}


// ---------------------------------------------------------
int
assignDx(Vector<Real>&              a_vectDx,
         const Vector<int>&         a_vectRefRatio,
         const Real                 a_dx0,
         const int a_numlevels)
{
  a_vectDx.resize(a_numlevels, 0.0);
  a_vectDx[0] = a_dx0; // (prob_hia[0]-prob_loa[0])/ancells[0];
  for (int ilev = 1; ilev < a_numlevels; ilev++)
    {
      a_vectDx[ilev] = a_vectDx[ilev-1]/a_vectRefRatio[ilev-1];
    }
  return(0);
}


// ---------------------------------------------------------
int
assignDomains(Vector<ProblemDomain>&     a_vectDomain,
              const Vector<int>&         a_vectRefRatio,
              const ProblemDomain&       a_baseDomain,
              const int a_numlevels)
{
  a_vectDomain.resize(a_numlevels);
  a_vectDomain[0] = a_baseDomain;
  for (int ilev = 1;ilev < a_numlevels; ilev++)
    {
      CH_assert(a_vectRefRatio[ilev-1] > 0);
      a_vectDomain[ilev] = refine(a_vectDomain[ilev-1], a_vectRefRatio[ilev-1]);
    }
  return(0);
}


// --------------------------------------------------------------
int
getDiff(Vector<LevelData<NodeFArrayBox>* >& a_vectDiff,
        const Vector<LevelData<NodeFArrayBox>* >& a_vectPhi1,
        const Vector<LevelData<NodeFArrayBox>* >& a_vectPhi2)
{
  int numlevels = a_vectPhi1.size();
  for (int ilev = 0; ilev < numlevels; ilev++)
    {
      const LevelData<NodeFArrayBox>& phi1 = *a_vectPhi1[ilev];
      const LevelData<NodeFArrayBox>& phi2 = *a_vectPhi2[ilev];
      const int nComp = phi1.nComp();
      CH_assert(nComp == phi2.nComp());
      const DisjointBoxLayout& grids1 = phi1.getBoxes();
      // const DisjointBoxLayout& grids2 = phi2.getBoxes();
      // CH_assert(grids1 == grids2);
      a_vectDiff[ilev] = new
        LevelData<NodeFArrayBox>(grids1, nComp, IntVect::Zero);
      LevelData<NodeFArrayBox>& diffLevel = *a_vectDiff[ilev];
      for (DataIterator dit(grids1.dataIterator()); dit.ok(); ++dit)
        {
          const FArrayBox& phi1Fab = phi1[dit()].getFab();
          const FArrayBox& phi2Fab = phi2[dit()].getFab();
          FArrayBox& diffFab = diffLevel[dit()].getFab();
          diffFab.copy(phi1Fab);
          diffFab -= phi2Fab;
        }
    }
  return(0);
}


// --------------------------------------------------------------
int
getMagnitude(Vector<LevelData<NodeFArrayBox>* >& a_vectMag,
             const Vector<LevelData<NodeFArrayBox>* >& a_vectField)
{
  int numlevels = a_vectField.size();
  for (int ilev = 0; ilev < numlevels; ilev++)
    {
      const LevelData<NodeFArrayBox>& field = *a_vectField[ilev];
      const int nComp = field.nComp();
      CH_assert(nComp == SpaceDim);
      const DisjointBoxLayout& grids = field.getBoxes();
      // const DisjointBoxLayout& grids2 = phi2.getBoxes();
      // CH_assert(grids1 == grids2);
      a_vectMag[ilev] = new
        LevelData<NodeFArrayBox>(grids, 1, IntVect::Zero);
      LevelData<NodeFArrayBox>& magLevel = *a_vectMag[ilev];
      for (DataIterator dit(grids.dataIterator()); dit.ok(); ++dit)
        {
          const FArrayBox& fieldFab = field[dit()].getFab();
          FArrayBox& magFab = magLevel[dit()].getFab();

          Box bx(magFab.box());
          for (BoxIterator bit(bx); bit.ok(); ++bit)
            {
              IntVect iv = bit();

              Real magPoint2 = 0.;
              for (int idim = 0; idim < SpaceDim; idim++)
                magPoint2 += fieldFab(iv, idim) * fieldFab(iv, idim);

              magFab(iv, 0) = sqrt(magPoint2);
            }
        }
    }
  return(0);
}


// --------------------------------------------------------------
int
getTruncError(Vector<LevelData<NodeFArrayBox>* >& a_vectErr,
              const Vector<LevelData<NodeFArrayBox>* >& a_vectLap,
              const Vector<LevelData<NodeFArrayBox>* >& a_vectRhs,
              const Vector<DisjointBoxLayout>& a_vectGrids,
              const Vector<Box>& a_vectDomain,
              const Vector<int>& a_vectRatio,
              const Vector<Real>& a_vectDx,
              const int a_numlevels,
              const bool a_verbose)
{
  Vector<ProblemDomain> probdomain(a_numlevels);
  for (int ilev = 0; ilev < a_numlevels; ilev++)
    probdomain[ilev] = ProblemDomain(a_vectDomain[ilev]);

  int eekflag = getTruncError(a_vectErr, a_vectLap, a_vectRhs,
                              a_vectGrids, probdomain, a_vectRatio, a_vectDx,
                              a_numlevels, a_verbose);
  return eekflag;
}


// --------------------------------------------------------------
int
getTruncError(Vector<LevelData<NodeFArrayBox>* >& a_vectErr,
              const Vector<LevelData<NodeFArrayBox>* >& a_vectLap,
              const Vector<LevelData<NodeFArrayBox>* >& a_vectRhs,
              const Vector<DisjointBoxLayout>& a_vectGrids,
              const Vector<ProblemDomain>& a_vectDomain,
              const Vector<int>& a_vectRatio,
              const Vector<Real>& a_vectDx,
              const int a_numlevels,
              const bool a_verbose)
{
  for (int ilev = 0; ilev < a_numlevels; ilev++)
    {
      const DisjointBoxLayout& grids = a_vectGrids[ilev];
      const ProblemDomain& domain = a_vectDomain[ilev];

      const LevelData<NodeFArrayBox>& lofPhiLevel = *a_vectLap[ilev];
      const LevelData<NodeFArrayBox>& rhsLevel = *a_vectRhs[ilev];
      int ncomp = rhsLevel.nComp();
      a_vectErr[ilev] = new
        LevelData<NodeFArrayBox>(grids, ncomp, IntVect::Zero);
      LevelData<NodeFArrayBox>& errLevel = *a_vectErr[ilev];

      // Need IVSVext for zeroing out boundary nodes.
      LayoutData< Vector<IntVectSet> > IVSV;
      LayoutData< Vector<IntVectSet> > IVSVext;
      interiorBoundaryNodes(IVSV, grids, domain);
      exteriorBoundaryNodes(IVSVext, IVSV, grids);

      for (DataIterator dit(lofPhiLevel.dataIterator()); dit.ok(); ++dit)
        {
          FArrayBox& errFab = errLevel[dit()].getFab();
          const FArrayBox& lapFab = lofPhiLevel[dit()].getFab();
          const FArrayBox& rhsFab = rhsLevel[dit()].getFab();

          errFab.copy(rhsFab);
          errFab -= lapFab;
        }
      zeroBoundaryNodes(errLevel, IVSVext);
    }
  return(0);
}


// --------------------------------------------------------------
int
project2(Vector<LevelData<NodeFArrayBox>* >& a_vectProj,
         const Vector<LevelData<NodeFArrayBox>* >& a_vectPhi,
         const Vector<DisjointBoxLayout>& a_vectGrids)
{
  int numlevels = a_vectPhi.size();
  for (int ilev = 0; ilev < numlevels; ilev++)
    {
      const LevelData<NodeFArrayBox>& phiLevel = *a_vectPhi[ilev];
      const int nComp = phiLevel.nComp();
      const DisjointBoxLayout& grids = phiLevel.getBoxes();
      DisjointBoxLayout gridsCoarse;
      coarsen(gridsCoarse, grids, 2);
      // a_vectProj[ilev] = new
      // LevelData<NodeFArrayBox>(a_vectGrids[ilev], nComp, IntVect::Zero);
      LevelData<NodeFArrayBox>* projPtr = new
        LevelData<NodeFArrayBox>(gridsCoarse, nComp, IntVect::Zero);
      LevelData<NodeFArrayBox>& proj = *projPtr;
      for (DataIterator dit(phiLevel.dataIterator()); dit.ok(); ++dit)
        {
          const FArrayBox& phiFab = phiLevel[dit()].getFab();
          FArrayBox& projFab = proj[dit()].getFab();
          Box projBox = projFab.box();

          for (BoxIterator bit(projBox); bit.ok(); ++bit)
            {
              IntVect ivc = bit();
              IntVect ivf = 2 * ivc;
              for (int ivar = 0; ivar < nComp; ivar++)
                projFab(ivc, ivar) = phiFab(ivf, ivar);
            }
        }
      // Unfortunately you have to do this copyTo business because
      // of some quirk that doesn't recognize the same indices
      // used in the two LevelDatas when DEBUG==TRUE.
      a_vectProj[ilev] = new
        LevelData<NodeFArrayBox>(a_vectGrids[ilev], nComp, IntVect::Zero);
      proj.copyTo(proj.interval(),
                  *a_vectProj[ilev],
                  a_vectProj[ilev]->interval());
      delete projPtr;
    }
  return(0);
}


// --------------------------------------------------------------
int
printErrorNorms(const string& a_prefix,
                Real a_dxFine,
                int a_deginterp,
                Real a_normMaxerrFine,
                Real a_norm1errFine,
                Real a_norm2errFine,
                Real a_normMaxerrCoarse,
                Real a_norm1errCoarse,
                Real a_norm2errCoarse)
{
  int dxinv = int((1.0 + 1.e-10)/a_dxFine);
  Real logTwo = log(2.);
  const char* cprefix = a_prefix.c_str();
  printf("%s %d & %d & $L_\\infty$ & 1/%-4d & %.5E & %.5E & %.2f \\\\ \\hline\n",
         cprefix, SpaceDim, a_deginterp, dxinv,
         a_normMaxerrFine, a_normMaxerrCoarse,
         log(a_normMaxerrCoarse/a_normMaxerrFine)/logTwo);

  printf("%s %d & %d & $L_1     $ & 1/%-4d & %.5E & %.5E & %.2f \\\\ \\hline\n",
         cprefix, SpaceDim, a_deginterp, dxinv,
         a_norm1errFine, a_norm1errCoarse,
         log(a_norm1errCoarse/a_norm1errFine)/logTwo);

  printf("%s %d & %d & $L_2     $ & 1/%-4d & %.5E & %.5E & %.2f \\\\ \\hline\n",
         cprefix, SpaceDim, a_deginterp, dxinv,
         a_norm2errFine, a_norm2errCoarse,
         log(a_norm2errCoarse/a_norm2errFine)/logTwo);
  return(0);
}


// --------------------------------------------------------------
int
printDiffNorms(const string& a_prefix,
               Real a_dxFine,
               int a_deginterp,
               Real a_normMaxdiff,
               Real a_norm1diff,
               Real a_norm2diff)
{
  int dxinv = int((1.0 + 1.e-10)/a_dxFine);
  int halfdxinv = int((1.0 + 1.e-10)/(0.5 * a_dxFine));
  const char* cprefix = a_prefix.c_str();
  printf("%s %d & %d & $L_\\infty$ & 1/%-4d & %.5E &             &   \\\\ \\hline\n",
         cprefix, SpaceDim, a_deginterp, dxinv, a_normMaxdiff);
  printf("%s %d & %d & $L_\\infty$ & 1/%-4d &             & %.5E &   \\\\ \\hline\n",
         cprefix, SpaceDim, a_deginterp, halfdxinv, a_normMaxdiff);
  printf("%s %d & %d & $L_1     $ & 1/%-4d & %.5E &             &   \\\\ \\hline\n",
         cprefix, SpaceDim, a_deginterp, dxinv, a_norm1diff);
  printf("%s %d & %d & $L_1     $ & 1/%-4d &             & %.5E &   \\\\ \\hline\n",
         cprefix, SpaceDim, a_deginterp, halfdxinv, a_norm1diff);
  printf("%s %d & %d & $L_2     $ & 1/%-4d & %.5E &             &   \\\\ \\hline\n",
         cprefix, SpaceDim, a_deginterp, dxinv, a_norm2diff);
  printf("%s %d & %d & $L_2     $ & 1/%-4d &             & %.5E &   \\\\ \\hline\n",
         cprefix, SpaceDim, a_deginterp, halfdxinv, a_norm2diff);
  return(0);
}

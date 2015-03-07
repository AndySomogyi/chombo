#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

//  ANAG, LBNL, DTG

#ifndef _EBINDEXSPACE_H_
#define _EBINDEXSPACE_H_

#include "REAL.H"
#include "IntVect.H"
#include "IntVectSet.H"
#include "VolIndex.H"
#include "FaceIndex.H"
#include "IntVectSet.H"
#include "DisjointBoxLayout.H"
#include "LevelData.H"
#include "EBISBox.H"
#include "GeometryService.H"
#include "EBISLayout.H"
#include "CH_HDF5.H"
#include <map>
#include "NamespaceHeader.H"

typedef std::map<DisjointBoxLayout, EBISLayout> dmap;

class EBIndexSpace;

///singleton version of EBIndexSpace.
class Chombo_EBIS
{
private:
  ///
  /**
   This is NULL until the first time
     instance() is called.
     On the first time it is called it
     is set to  new EBIndexSpace();
     After that it points to the same
     thing.
   */
  static EBIndexSpace* s_instance;

  Chombo_EBIS()
  {;}

  ~Chombo_EBIS()
  {;}
public:
  ///
  /**
     Return a pointer to the singleton
     EBIndexSpace.  Until this is called
     once, the singleton is null.  See
     Design patterns for why this is done
     this way.
   */
  static EBIndexSpace* instance();

};

///
/**
   EBISLevel represents the geometric information
   of the domain at a particular level of refinement.
*/
class EBISLevel
{
public:
  void coarsenVoFs(EBISLevel& a_fineEBIS);

  void fixFineToCoarse(EBISLevel& a_fineEBIS);

  void coarsenFaces(EBISLevel& a_fineEBIS);

  int numVoFsOnProc() const;

  void levelStitch(EBISLevel&       a_otherPhase,
                   const EBISLevel* a_fineThisPtr,
                   const EBISLevel* a_fineOtherPtr); // MF Addition  bvs

  void cellStitch(EBData&           a_ebdataCoarA,
                  EBData&           a_ebdataCoarB,
                  const EBGraph&    a_ebgrapCoarA,
                  const EBGraph&    a_ebgrapCoarB,
                  const EBISBox&    a_ebisbxFineA,
                  const EBISBox&    a_ebisbxFineB,
                  const IntVect&    a_iv,
                  const int&        a_aphase,
                  const int&        a_bphase);
  EBISLevel();
  ///
  /**
     Defines the level.  Calls the geoserver
     to fill the layout with geometric information.
   */
  EBISLevel(const ProblemDomain& a_domain,
            const RealVect& a_origin,
            const Real& a_dx,
            const GeometryService& a_geoserver,
            const EBIndexSpace* const ebisPtr = Chombo_EBIS::instance());

  ///
  /**
     Defines the level from a finer level
     (factor of two refinement ratio).
     Fills in all the coarse-fine information
     in both levels.
   */
  EBISLevel(EBISLevel& a_finerLevel, const GeometryService& a_geoserver,
            const EBIndexSpace* const ebisPtr = Chombo_EBIS::instance());

#ifdef CH_USE_HDF5
  /// read from file
  EBISLevel(HDF5Handle& handle);

  ///
  void
  write(HDF5Handle&   a_handle);
#endif

  ///
  /**
     Checks to see the vofs are in the correct cells.
     Checks to see that the faces are over the correct cells.
     Checks that volume fractions, area fractions are positive.
     Bail out with MayDay::Error if anything fails.
  */
  void sanityCheck(const EBIndexSpace* const ebisPtr = Chombo_EBIS::instance());

  ///
  void fillEBISLayout(EBISLayout& a_ebis,
                      const DisjointBoxLayout& a_grids,
                      const int& nghost) const;
  ~EBISLevel();

  const ProblemDomain& getDomain() const;

  const DisjointBoxLayout& getGrids() const;

  DisjointBoxLayout getIrregGrids() const;

  DisjointBoxLayout getFlowGrids() const;

  IntVectSet irregCells() const;

  const Real& getDX() const;

  const RealVect& getOrigin() const;

  ///
  void fixRegularNextToMultiValued();

  void clearMultiBoundaries();
  void setBoundaryPhase(int phase);

  int m_phase;



private:

  ///
  DisjointBoxLayout m_grids;

  ///
  ProblemDomain m_domain;

  ///
  RealVect m_origin;

  ///
  Real m_dx;

  ///
  LevelData<EBGraph> m_graph;

  static bool s_verbose;

  static Real s_tolerance;
  ///
  LevelData<EBData>  m_data;

  //make grids for this level.

  static void makeBoxes(Vector<Box>& a_boxes,
                        Vector<long>& a_irregCount,
                        const Box& a_region,
                        const ProblemDomain& a_domain,
                        const GeometryService& a_geoserver,
                        const RealVect&        a_origin,
                        const Real&            a_dx,
                        const int              a_ncellmax,
                        const EBIndexSpace* const ebisPtr);
  static void makeBoxes(std::list<Box>& a_boxes,
                        const Box& a_region,
                        const ProblemDomain& a_domain,
                        const GeometryService& a_geoserver,
                        const RealVect&        a_origin,
                        const Real&            a_dx,
                        const int              a_ncellmax);

  //these are disallowed
  void operator=(const EBISLevel& ebiin)
  {
    MayDay::Error("EBISLevel::operator= not allowed");
  }

  EBISLevel(const EBISLevel& ebiin)
  {
    MayDay::Error("EBISLevel::copy constructor not allowed");
  }

  Real m_tolerance; // used in Multifluid face matching algorithm.

  bool match(const Real& a, const Real& b) const;

  bool match(const RealVect& a, const RealVect& b) const;

  void refreshCache() const ;
  mutable std::map<int, dmap> cache;
  mutable int m_cacheMisses, m_cacheHits, m_cacheStale;

  friend class EBIndexSpace;
};

///
/**
   EBIndexSpace represents the geometric information
   of the domain.    It should be generated through
   the Chombo_EBIS class.  This follows the singleton
   pattern.
*/
class EBIndexSpace
{
public:

  EBIndexSpace();
  ///
  /**
      If a_ncellMax is set, that is the max width of
      an internal grid.  Otherwise use defaults
      of (16 in 3D, 64 in 2d)
   */
  void
  define(const ProblemDomain& a_domain,
         const RealVect& a_origin,
         const Real& a_dx,
         const GeometryService& a_geoserver,
         int a_nCellMax = -1,
         int maxCoarsenings = -1);

#ifdef CH_USE_HDF5
  ///
  void
  define(HDF5Handle& handle, int maxCoarsenings = -1);
#endif

  void
  define(EBISLevel* level0,
         int a_nCellMax = -1,
         int a_maxCoarsenings = -1);

  EBISLevel* buildFirstLevel(const ProblemDomain& a_domain,
                             const RealVect& a_origin,
                             const Real& a_dx,
                             const GeometryService& a_geoserver,
                             int a_nCellMax,
                             int a_maxCoarsenings);

  EBISLevel* buildNextLevel(const GeometryService& a_geoserver);

  ///
  ~EBIndexSpace();

  ///
  int numLevels() const;

  ///
  /**
     Get the number of vofs over the entire domain.
     This is blocking as a broadcast and gather are required.
   */
  int numVoFs(const ProblemDomain& a_domain) const;

  ///
  /**
     return level index of domain.
     return -1 if a_domain does not
     correspond to any refinement of EBIS.
   */
  int getLevel(const ProblemDomain& a_domain) const;

  RealVect getOrigin() const {return m_ebisLevel[0]->getOrigin();}
  ///
  /**
     returns the problem domain box at a level of refinement
  */
  const ProblemDomain& getBox(int level) const
  {
    return m_domainLevel[level];
  }

  ///
  void fillEBISLayout(EBISLayout& a_ebis,
                      const DisjointBoxLayout& a_grids,
                      const ProblemDomain& a_domain,
                      const int & nghost) const;

  ///
  /**
     Return true if the define function has been called.
   */
  bool isDefined() const;

  ///
  void clear();

  int getNCellMax() const;
  ///
  /**
     This makes Chombo_EBIS be the only function
     that can call the private constructor.
   */
  friend class Chombo_EBIS;

#ifdef CH_USE_HDF5
  ///
  /**
     Writes out finest level as a default.   Writes at the level you are domain
     you want if a_outputLevel is defined.
   */
  void
  write(HDF5Handle& handle, 
        ProblemDomain a_outputLevel = ProblemDomain());

  void
  writeInfo(HDF5Handle& handle) const ;

#endif

  DisjointBoxLayout getIrregGrids(int level) const
  { return m_ebisLevel[level]->getIrregGrids();}

  DisjointBoxLayout getFlowGrids(int level) const
  { return m_ebisLevel[level]->getFlowGrids();}

  DisjointBoxLayout levelGrids(int level) const
  {  return m_ebisLevel[level]->m_grids;}

  Real dx(int level) const
  {  return m_ebisLevel[level]->m_dx;}

  IntVectSet irregCells(int depth=0) const
  {
    return m_ebisLevel[depth]->irregCells();
  }

  void setCellMax(int max) {m_nCellMax = max;}

  void resetLevels(int nLevel);

  static bool s_MFSingleBox;
private:
  //this is private to force the singleton thing.
  //EBIndexSpace();

  int  m_nCellMax;
  bool m_isDefined;

  Vector<EBISLevel*> m_ebisLevel;
  Vector<ProblemDomain> m_domainLevel;
  int m_nlevels;
  static Real s_tolerance;
  static bool s_verbose;

  //disallowed for performance reasons
  void operator=(const EBIndexSpace& ebiin)
  {
    MayDay::Error("no assignment of EBIndexSpace allowed");
  }
  EBIndexSpace(const EBIndexSpace& ebiin)
  {
    MayDay::Error("no copy construction of EBIndexSpace allowed");
  }
};

inline
bool EBISLevel::match(const Real& a, const Real& b) const
{
  return a<b+m_tolerance && a>b-m_tolerance;
}

inline
bool EBISLevel::match(const RealVect& a, const RealVect& b) const
{
  return D_TERM( a[0]<b[0]+m_tolerance && a[0]>b[0]-m_tolerance, &&
                 a[1]<b[1]+m_tolerance && a[1]>b[1]-m_tolerance, &&
                 a[2]<b[2]+m_tolerance && a[2]>b[2]-m_tolerance) ;
}

#include "NamespaceFooter.H"
#endif

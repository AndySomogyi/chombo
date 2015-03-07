#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

#ifndef _CFIVS_H_
#define _CFIVS_H_

#include <iostream>
#include <cmath>
#include "SPACE.H"
#include <cstdlib>

#include "REAL.H"
#include "IntVect.H"
#include "Box.H"
#include "LoHiSide.H"
#include "DisjointBoxLayout.H"
#include "IntVectSet.H"
#include "ProblemDomain.H"
#include "NamespaceHeader.H"

/// Internal Chombo class to get fine intvectset at CF interface
/**
    Class to get fine intvectset at cf interface for one side of a given
    fine-level Box.
    This class should be considered internal to AMRSolver and should
    not be considered part of the Chombo API.
*/
class CFIVS
{
public:
  /// Null constructor
  CFIVS();

  /// Destructor
  ~CFIVS();

  /// Deprecated constructor (The ProblemDomain class has replaced the Box domain)
  /**
     Replaces Box domain with a ProblemDomain and calls define fn.
  */
  CFIVS(const Box&               a_domain,
        const Box&               a_boxIn,
        const DisjointBoxLayout& a_fineBoxes,
        int                      a_direction,
        Side::LoHiSide           a_hiorlo);

  /// Full constructor
  /**
     Calls corresponding define function.
   */
  CFIVS(const ProblemDomain&     a_domain,
        const Box&               a_boxIn,
        const DisjointBoxLayout& a_fineBoxes,
        int                      a_direction,
        Side::LoHiSide           a_hiorlo);

  /// Deprecated define function
  /**
      Replaces Box domain with a ProblemDomain and calls define fn.
  */
  void define(const Box&               a_domain,
              const Box&               a_boxIn,
              const DisjointBoxLayout& a_fineBoxes,
              int                      a_direction,
              Side::LoHiSide           a_hiorlo);

  /// Full define function
  /**
   */
  void define (const ProblemDomain&     a_domain,    /// Problem domain at the fine level
               const Box&               a_boxIn,     /// fine-level grid box
               const DisjointBoxLayout& a_fineBoxes, /// fine-level grids
               int                      a_direction, /// direction of the face of interest
               Side::LoHiSide           a_hiorlo);   /// is the face on the hi or lo side of the box?

  /// faster define function.  Requires the use of CFStencil::buildPeriodicVector static function.
  void define(const ProblemDomain&     a_domain,
              const Box&               a_boxIn,
              const Vector<Box>&       a_periodicfineBoxes,
              int                      a_direction,
              Side::LoHiSide           a_hiorlo);

  /// Returns "true" if the object has been defined
  bool isDefined() const;

  /// Returns "true" if this CF IntVectSet can be represented as a simple Box.
  /**
      Returns 'true' if this coarse-fine IntVectSet can be represented
      as just a Box, ie.  a user can perform a dense computation instead of
      a pointwise calculation using IVSIterator.
  */
  bool isPacked() const;

  /// Returns packed box for when isPacked() is true.
  /** if isPacked() returns 'true', then the box returned by this function
      is a suitable substitute for this IntVectSet
  */
  const Box& packedBox() const;

  /// Get fine intvects which need to be interpolated.
  /**
      This will be empty if isEmpty() returns true
  */
  const IntVectSet& getFineIVS() const ;

  static long long m_packCount, m_sparseCount;

  void coarsen(int a_ref);

protected:
  // For internal use
  void  setDefaultValues();

  // The following variables are mainly for self-identification

  /**
     Fine intvects which need to be interpolated
  */
  IntVectSet m_fiinterpIVS;

  /**
     Can the coarse-fine IntVectSet be represented as a Box?
  */
  bool m_packed;

  /**
      If m_packed is true, the Box with which the coarse-fine
      IntVectSet can be represented
  */
  Box  m_packedBox;

  /**
     Has this object been defined?
  */
  bool m_isDefined;


};

inline bool CFIVS::isPacked() const
{
  return m_packed;
}

inline const Box& CFIVS::packedBox() const
{
  return m_packedBox;
}

#include "NamespaceFooter.H"
#endif

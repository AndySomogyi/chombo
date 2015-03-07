#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

#ifndef _REFINEMENTCRITERION_H_
#define _REFINEMENTCRITERION_H_

#if defined(CH_Darwin) && defined(__GNUC__) && ( __GNUC__ == 3 )
// deal with the broken isnan()/isinf() in GCC on MacOS
#include <unistd.h>
#define _GLIBCPP_USE_C99 1
#endif

#include "MayDay.H"
#include "REAL.H"
#include "Vector.H"

#include "Notation.H"

#include "NamespaceHeader.H"

///
/**
   This is the base class for refinement criterion used for the subdivision of
   cells in geometry generation.  If the result if "doRefine()" is true, the
   cell is subdivided in the directions in which "a_refineInDir" is non-zero.
   The results of this subdivision are then combined to get results for the
   original cell.
 */
class RefinementCriterion
{
public:
  /// Null constructor
  RefinementCriterion()
  {
  }

  /// Destructor
  virtual ~RefinementCriterion()
  {
  }

  /// Should a cell be subdivided and in which directions
  /**
     This method returns true if the current cell should be subdivided.  The
     subdivsion should occur in all the directions where "a_refineInDir" is
     non-zero.
   */
  virtual bool doRefine(Vector<int>                 & a_refineInDir,
                        const int                   & a_dim,
                        const Vector<Real>          & a_dx,
                        const Vector<Vector<Real> > & a_residual) = 0;

protected:
};

#include "NamespaceFooter.H"

#endif

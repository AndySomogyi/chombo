#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

#ifndef _FACESREFINEMENT_H_
#define _FACESREFINEMENT_H_

#if defined(CH_Darwin) && defined(__GNUC__) && ( __GNUC__ == 3 )
// deal with the broken isnan()/isinf() in GCC on MacOS
#include <unistd.h>
#define _GLIBCPP_USE_C99 1
#endif

#include "RefinementCriterion.H"

#include "NamespaceHeader.H"

class FacesRefinement: public RefinementCriterion
{
public:
  // constructors
  FacesRefinement(const bool & a_refinedOnce=true,
                  const bool & a_refinedTwice=true,
                  const bool & a_refinedThree=true,
                  const bool & a_refinedFour=true);

  FacesRefinement(const FacesRefinement& a_facesRefinement);

  // destructor
  virtual ~FacesRefinement();

  virtual bool doRefine(Vector<int>                 & a_refineInDir,
                        const int                   & a_dim,
                        const Vector<Real>          & a_dx,
                        const Vector<Vector<Real> > & a_residual);

  virtual bool doRefine(const int & a_dim,
                        const Real & a_dxRatio);

  virtual bool doRefine(const Vector< Vector<Real> > & a_residual);

  void setRefinedOnce();

  void setRefinedTwice();

  void setRefinedThree();

  void setRefinedFour();

protected:
  // member data
  bool m_refinedOnce;
  bool m_refinedTwice;
  bool m_refinedThree;
  bool m_refinedFour;
};

#include "NamespaceFooter.H"

#endif


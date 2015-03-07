#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

#ifndef _REFERENCEHEIGHTIF_H_
#define _REFERENCEHEIGHTIF_H_

#include "Notation.H"
#include "RealVect.H"
#include "IndexTM.H"
#include "BaseIF.H"
#include "NamespaceHeader.H"

///
/**
   This adapter wraps any BaseIF subclass and gives it the
   the concept of a reference height. The adapter passes through
   all function calls to the adapted object, except value(RealVect)..
   The reference height is used to evaluate the value function when
   SpaceDim=GLOBALDIM-1,

 */
class ReferenceHeightIF : public BaseIF
{
public:

  /**
   Construct the ReferenceHeightIF using another IF plus
   a reference height. The delegate IF will be copied.
  */
  ReferenceHeightIF(const  BaseIF & a_implicitFunction,
                    const  Real   & a_referenceHeight):
    m_implicitFunction(a_implicitFunction.newImplicitFunction()),
    m_referenceHeight(a_referenceHeight){;}

  /// Destructor cleans up the member implicit function
  virtual ~ReferenceHeightIF()
  {
    delete m_implicitFunction;
  }

  ///
  /**
     Return the value of the function at a_point. value(RealVect)
     to value(IndexTM). However, in the case where GLOBALDIM == 3 and
     SpaceDim == 2 the missing dimension is filled in with the reference height.
     Subclasses need not implement this function, but if

  */
  virtual Real value(const RealVect& a_point) const
  {
    IndexTM<Real,GLOBALDIM> pt;

    for (int idir = 0; idir < SpaceDim; ++idir)
      {
        pt[idir] = a_point[idir];
      }

    if (GLOBALDIM == 3 && SpaceDim == 2)
    {
      pt[SpaceDim] = m_referenceHeight;
    }

    return value(pt);
  }


   ///
  /**
     Return the value of the function at a_point (of type INdexTM).
  */
  virtual Real value(const IndexTM<Real,GLOBALDIM>& a_point) const
  {
    return m_implicitFunction->value(a_point);
  };

  ///
  /**
     Return the gradient of the function at a_point.
  */
  virtual IndexTM<Real,GLOBALDIM> grad(const IndexTM<Real,GLOBALDIM>& a_point) const
  {
    return m_implicitFunction->grad(a_point);
  };

  ///
  /**
     Return the normalized gradient of the function at a_point.
  */
  virtual IndexTM<Real,GLOBALDIM> normal(const IndexTM<Real,GLOBALDIM>& a_point) const
  {
    return m_implicitFunction->normal(a_point);
  };

  ///
  /**
     Return the gradient of the gradient at a_point.
  */
  virtual Vector<IndexTM <Real,GLOBALDIM> > gradGrad(const IndexTM<Real,GLOBALDIM>& a_point) const
  {
    return m_implicitFunction->gradGrad(a_point);
  };

  ///
  /**
     Return the gradient of the normal at a_point.
  */
  virtual Vector<IndexTM <Real,GLOBALDIM> > gradNormal(const IndexTM<Real,GLOBALDIM>& a_point) const
  {
    return m_implicitFunction->gradNormal(a_point);
  };

  ///
  /**
     Return a newly allocated derived class.  The responsibility
     for deleting the memory is left to the calling function.
  */
  virtual ReferenceHeightIF* newImplicitFunction() const
  {
    return new ReferenceHeightIF(*m_implicitFunction, m_referenceHeight);
  }

  ///
  /**
    Return the reference height for this IF
  */
  Real getReferenceHeight() const
  {
    return m_referenceHeight;
  }

private:

  // Default constructor is disabled in order to prevent initialization
  // without reference height. Base classes should all use the constructor
  // version that takes a reference height
  ReferenceHeightIF()
  {
  }

  BaseIF* m_implicitFunction;
  Real    m_referenceHeight;
};

#include "NamespaceFooter.H"
#endif

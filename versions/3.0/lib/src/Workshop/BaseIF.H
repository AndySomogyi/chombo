#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

#ifndef _BASEIF_H_
#define _BASEIF_H_

#include "RealVect.H"
#include "ProblemDomain.H"
#include "IndexTM.H"

#include "Notation.H"
#include "GeometryService.H"

#include "NamespaceHeader.H"

///
/**
   This is the base class for an implicit function specification of geometry.
   All that is needed is a constructor/destructor, a method to give the value
   of the function at any point in space (in 2D or 3D), and a factory method
   (these will probably all be very similar).
 */
class BaseIF
{

public:

  /// Default constructor
  BaseIF()
  {
  }

  /// Default destructor
  virtual ~BaseIF()
  {
  }

  ///
  /**
   Return the value of the function at a_point.  When delineating a domain,
   the level set value=0 represents the boundary and value<0 is inside the
   fluid.
  */
  virtual Real value(const RealVect& a_point) const = 0;


  virtual bool fastIntersection(const Box&           a_region,
                                const ProblemDomain& a_domain,
                                const RealVect&      a_origin,
                                const Real&          a_dx) const
  {
    RealVect low, high;
    corners(a_region, a_origin, a_dx, low, high);
    return fastIntersection(low, high);
  }

  virtual bool fastIntersection(const RealVect& a_low,
                                const RealVect& a_high) const
  {return false;}

  virtual GeometryService::InOut InsideOutside(const Box&           a_region,
                                               const ProblemDomain& a_domain,
                                               const RealVect&      a_origin,
                                               const Real&          a_dx) const
  {
    RealVect low, high;
    corners(a_region, a_origin, a_dx, low, high);
    return InsideOutside(low, high);
  }

  virtual GeometryService::InOut InsideOutside(const RealVect& a_low,
                                               const RealVect& a_high) const

  {
    MayDay::Error("This class has not implemented a fastIntersection operation");
    return GeometryService::Irregular;
  }
   ///
  /**
     Return the value of the function at a_point (of type INdexTM).
  */
  virtual Real value(const IndexTM<Real,GLOBALDIM>& a_point) const
  {
    MayDay::Abort("'value' called with 'IndexTM' not implemented in BaseIF");
    Real retval = LARGEREALVAL;
    return retval;
  };

  ///
  /**
     Return the gradient of the function at a_point.
  */
  virtual IndexTM<Real,GLOBALDIM> grad(const IndexTM<Real,GLOBALDIM>& a_point) const
  {
    MayDay::Abort("'grad' not implemented in BaseIF");
    return IndexTM<Real,GLOBALDIM>();
  };

  ///
  /**
     Return the normalized gradient of the function at a_point.
     The gradients and normals are based on value, which yields
     outward-pointing normals (which is different the computational
     computation elsewhere in Chombo).
  */
  virtual IndexTM<Real,GLOBALDIM> normal(const IndexTM<Real,GLOBALDIM>& a_point) const
  {
    MayDay::Abort("'normal' not implemented in BaseIF");
    return IndexTM<Real,GLOBALDIM>();
  };

  ///
  /**
     Return the gradient of the gradient at a_point.
  */
  virtual Vector<IndexTM <Real,GLOBALDIM> > gradGrad(const IndexTM<Real,GLOBALDIM>& a_point) const
  {
    MayDay::Abort("'gradGrad' not implemented in BaseIF");
    return Vector<IndexTM <Real,GLOBALDIM> >();
  };

  ///
  /**
     Return the gradient of the normal at a_point.
  */
  virtual Vector<IndexTM <Real,GLOBALDIM> > gradNormal(const IndexTM<Real,GLOBALDIM>& a_point) const
  {
    MayDay::Abort("'gradNormal' not implemented in BaseIF");
    return Vector<IndexTM <Real,GLOBALDIM> >();
  };

  ///
  /**
     Return a newly allocated derived class.  The responsibility
     for deleting the memory is left to the calling function.
  */
  virtual BaseIF* newImplicitFunction() const = 0;


  static void corners(const Box& a_region, const RealVect& a_origin, const Real& a_dx,
                      RealVect& a_lo, RealVect& a_hi)
  {
    a_lo = a_origin + RealVect(a_region.smallEnd())*a_dx;
    a_hi = a_origin + RealVect(a_region.bigEnd()+IntVect::Unit)*a_dx;
  }
};

#include "NamespaceFooter.H"

#endif

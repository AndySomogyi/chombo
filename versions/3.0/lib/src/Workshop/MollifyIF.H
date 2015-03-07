#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

#ifndef _MOLLIFYIF_H_
#define _MOLLIFYIF_H_

#include "MayDay.H"
#include "RealVect.H"
#include "Vector.H"

#include "BaseIF.H"

#include "NamespaceHeader.H"

///
/**
    This implicit function mollifies on implicit function using another
    implicit function (e.g., "GaussianSphereIF").  Basically, the IF being
    mollified is sampled in a small (specified) neighborhood of a point and
    convolved with the mollifying function in this neighborhood.
 */
class MollifyIF: public BaseIF
{
public:
  ///
  /**
      Constructor specifying the implicit function to mollify (a_impFunc),
      the mollifier (a_mollifier), a neighorhood to use (a_min, a_max) is
      each dimension, and a sampling (a_numPts).
   */
  MollifyIF(const BaseIF& a_impFunc,
            const BaseIF& a_mollifier,
            const Real&   a_min,
            const Real&   a_max,
            const int&    a_numPts);

  /// Copy constructor
  MollifyIF(const MollifyIF& a_inputIF);

  /// Destructor
  virtual ~MollifyIF();

  ///
  /**
      Return the value of the function at a_point.
   */
  virtual Real value(const RealVect& a_point) const;

  ///
  /**
     Return the value of the function at a_point.
   */
  virtual Real value(const IndexTM<Real,GLOBALDIM>& a_point) const;

  ///
  /**
     Return the gradient of the function at a_point.
  */
  virtual IndexTM<Real,GLOBALDIM> grad(const IndexTM<Real,GLOBALDIM>& a_point) const;

  ///
  /**
     Return the normalized gradient of the function at a_point.
  */
  virtual IndexTM<Real,GLOBALDIM> normal(const IndexTM<Real,GLOBALDIM>& a_point) const;

  ///
  /**
     Return the gradient of the gradient at a_point.
  */
  virtual Vector<IndexTM <Real,GLOBALDIM> > gradGrad(const IndexTM<Real,GLOBALDIM>& a_point) const;

  ///
  /**
     Return the gradient of the normal at a_point.
  */
  virtual Vector<IndexTM <Real,GLOBALDIM> > gradNormal(const IndexTM<Real,GLOBALDIM>& a_point) const;

  virtual BaseIF* newImplicitFunction() const;

protected:
  BaseIF*   m_impFunc;          // function to mollify
  BaseIF*   m_mollifier;        // mollifier

  Real      m_min;              // minimum of the mollifying neighborhood
  Real      m_max;              // maximum of the mollifying neighborhood
  int       m_numPts;           // number of sample points in meighborhood

  Real      m_dx;               // spacing in the sampled neighborhood
  Box       m_sampleBox;        // the sample box

  // These are "mutable" so they can be modified in "const" methods - this is
  // done so that they can be cached but ONLY if the function which needs them
  // is called.
  mutable Real m_mollifierSum;                  // the total of the samples
  mutable FArrayBox m_sampledMollifier;         // sampled mollifier

  mutable FArrayBox m_sampledGradMollifier;     // sampled grad(mollifier)
  mutable bool m_cachedGradMollifier;           // flag it has been cached

  mutable FArrayBox m_sampledGradGradMollifier; // sampled gradGrad(mollifier)
  mutable bool m_cachedGradGradMollifier;       // flag it has been cached

private:
  MollifyIF()
  {
    MayDay::Error("MollifyIF uses strong construction");
  }

  void operator=(const MollifyIF& a_inputIF)
  {
    MayDay::Error("MollifyIF doesn't allow assignment");
  }
};

#include "NamespaceFooter.H"
#endif

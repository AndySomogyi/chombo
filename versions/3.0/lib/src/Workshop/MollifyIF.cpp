#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

#include "MollifyIF.H"

#include "NamespaceHeader.H"

MollifyIF::MollifyIF(const BaseIF& a_impFunc,
                     const BaseIF& a_mollifier,
                     const Real&   a_min,
                     const Real&   a_max,
                     const int&    a_numPts)
{
  // Copy the implicit functions
  m_impFunc   = a_impFunc.newImplicitFunction();
  m_mollifier = a_mollifier.newImplicitFunction();

  // Save the sampled neighborhood specifications
  m_min    = a_min;
  m_max    = a_max;
  m_numPts = a_numPts;

  // Spacing of the sample points
  if (m_numPts > 1)
  {
    m_dx = (m_max - m_min) / (m_numPts - 1);
  }
  else
  {
    m_dx = 0.0;
  }

  // Box for the sampled mollifier
  m_sampleBox.define(IntVect::Zero, (m_numPts-1) * IntVect::Unit);

  // Create the FArrayBox for the sampled mollifier
  int nComp = 1;

  m_sampledMollifier.define(m_sampleBox,nComp);

  // Sum the mollifier values so they can be normalized
  m_mollifierSum = 0.0;

  // Sample the mollifier on the necessary grid
  for (BoxIterator bit(m_sampleBox); bit.ok(); ++bit)
  {
    IntVect iv = bit();
    RealVect point;

    for (int idir = 0; idir < SpaceDim; idir++)
    {
      point[idir] = m_max - iv[idir] * m_dx;
    }

    Real value = m_mollifier->value(point);

    m_mollifierSum += value;
    m_sampledMollifier(iv,0) = value;
  }

  // Normalize the sampled mollifier to be a partition of unity (i.e. to sum
  // to 1.0)
  m_sampledMollifier *= 1.0 / m_mollifierSum;

  // grad and gradGrad haven't been computed/cached yet
  m_cachedGradMollifier     = false;
  m_cachedGradGradMollifier = false;
}

MollifyIF::MollifyIF(const MollifyIF& a_inputIF)
{
  // Copy things over
  if (a_inputIF.m_impFunc == NULL)
  {
    m_impFunc = NULL;
  }
  else
  {
    m_impFunc = a_inputIF.m_impFunc->newImplicitFunction();
  }

  if (a_inputIF.m_mollifier == NULL)
  {
    m_mollifier = NULL;
  }
  else
  {
    m_mollifier = a_inputIF.m_mollifier->newImplicitFunction();
  }

  m_min    = a_inputIF.m_min;
  m_max    = a_inputIF.m_max;
  m_numPts = a_inputIF.m_numPts;

  m_dx     = a_inputIF.m_dx;

  m_sampleBox = a_inputIF.m_sampleBox;

  m_sampledMollifier.define(a_inputIF.m_sampledMollifier.box(),
                            a_inputIF.m_sampledMollifier.nComp());
  m_sampledMollifier.copy  (a_inputIF.m_sampledMollifier);

  m_cachedGradMollifier = a_inputIF.m_cachedGradMollifier;
  m_sampledGradMollifier.define(a_inputIF.m_sampledGradMollifier.box(),
                                a_inputIF.m_sampledGradMollifier.nComp());
  m_sampledGradMollifier.copy  (a_inputIF.m_sampledGradMollifier);

  m_cachedGradGradMollifier = a_inputIF.m_cachedGradGradMollifier;
  m_sampledGradGradMollifier.define(a_inputIF.m_sampledGradGradMollifier.box(),
                                    a_inputIF.m_sampledGradGradMollifier.nComp());
  m_sampledGradGradMollifier.copy  (a_inputIF.m_sampledGradGradMollifier);
}

MollifyIF::~MollifyIF()
{
  // Delete the IFs (if they exist)
  if (m_impFunc != NULL)
  {
    delete m_impFunc;
  }

  if (m_mollifier != NULL)
  {
    delete m_mollifier;
  }
}

Real MollifyIF::value(const RealVect& a_point) const
{
  Real retval = 0.0;

  // Apply the mollifier to the function at the current point (a_point)
  for (BoxIterator bit(m_sampleBox); bit.ok(); ++bit)
  {
    IntVect iv = bit();
    RealVect point(a_point);

    for (int idir = 0; idir < SpaceDim; idir++)
    {
      point[idir] += m_min + iv[idir] * m_dx;
    }

    retval += m_sampledMollifier(iv,0) * m_impFunc->value(point);
  }

  return retval;
}

Real MollifyIF::value(const IndexTM<Real,GLOBALDIM>& a_point) const
{
  RealVect point;
  for (int idir = 0 ; idir < GLOBALDIM ; idir++)
    {
      point[idir] = a_point[idir];
    }

  return value(point);
}

IndexTM<Real,GLOBALDIM> MollifyIF::grad(const IndexTM<Real,GLOBALDIM>& a_point) const
{
  // Check if the sampled grad(mollifier) has been computed/cached
  if (!m_cachedGradMollifier)
  {
    // Create the FArrayBox for the sampled mollifier
    int nComp = GLOBALDIM;

    m_sampledGradMollifier.define(m_sampleBox,nComp);

    // Sample the mollifier on the necessary grid
    for (BoxIterator bit(m_sampleBox); bit.ok(); ++bit)
    {
      IntVect iv = bit();
      IndexTM<Real,GLOBALDIM> point;

      for (int idir = 0; idir < SpaceDim; idir++)
      {
        point[idir] = m_max - iv[idir] * m_dx;
      }

      IndexTM<Real,GLOBALDIM> value = m_mollifier->grad(point);

      for (int idir = 0; idir < GLOBALDIM; idir++)
      {
        m_sampledGradMollifier(iv,idir) = value[idir];
      }
    }

    // Normalize the sampled mollifier to be a partition of unity (i.e. to sum
    // to 1.0)
    m_sampledGradMollifier *= 1.0 / m_mollifierSum;

    m_cachedGradMollifier = true;
  }

  IndexTM<Real,GLOBALDIM> retval;

  for (int idir = 0; idir < GLOBALDIM; idir++)
  {
    retval[idir] = 0.0;
  }

  // Apply grad(mollifier) to the function at the current point (a_point)
  for (BoxIterator bit(m_sampleBox); bit.ok(); ++bit)
  {
    IntVect iv = bit();
    RealVect point;

    for (int idir = 0; idir < SpaceDim; idir++)
    {
      point[idir] = a_point[idir] + m_min + iv[idir] * m_dx;
    }

    Real value = m_impFunc->value(point);

    for (int idir = 0; idir < GLOBALDIM; idir++)
    {
      retval[idir] += m_sampledGradMollifier(iv,idir) * value;
    }
  }

  return retval;
}

IndexTM<Real,GLOBALDIM> MollifyIF::normal(const IndexTM<Real,GLOBALDIM>& a_point) const
{
  // Get the gradient
  IndexTM<Real,GLOBALDIM> retval = grad(a_point);

  // Get the length (norm) of the gradient
  Real norm = 0.0;
  for (int idir = 0; idir < GLOBALDIM; idir++)
  {
    norm += retval[idir] * retval[idir];
  }
  norm = sqrt(norm);

  // Normalize the gradient to get the normal
  if (norm != 0.0)
  {
    for (int idir = 0; idir < GLOBALDIM; idir++)
    {
      retval[idir] /= norm;
    }
  }

  return retval;
}

Vector<IndexTM<Real,GLOBALDIM> > MollifyIF::gradGrad(const IndexTM<Real,GLOBALDIM>& a_point) const
{
  // Check if the sampled grad(mollifier) has been computed/cached
  if (!m_cachedGradGradMollifier)
  {
    // Create the FArrayBox for the sampled mollifier
    int nComp = GLOBALDIM * GLOBALDIM;

    m_sampledGradGradMollifier.define(m_sampleBox,nComp);

    // Sample the mollifier on the necessary grid
    for (BoxIterator bit(m_sampleBox); bit.ok(); ++bit)
    {
      IntVect iv = bit();
      IndexTM<Real,GLOBALDIM> point;

      for (int idir = 0; idir < SpaceDim; idir++)
      {
        point[idir] = m_max - iv[idir] * m_dx;
      }

      Vector<IndexTM<Real,GLOBALDIM> > value = m_mollifier->gradGrad(point);

      int index = 0;
      for (int idir = 0; idir < GLOBALDIM; idir++)
      {
        for (int jdir = 0; jdir < GLOBALDIM; jdir++)
        {
          m_sampledGradGradMollifier(iv,index) = value[idir][jdir];
          index++;
        }
      }
    }

    // Normalize the sampled mollifier to be a partition of unity (i.e. to sum
    // to 1.0)
    m_sampledGradGradMollifier *= 1.0 / m_mollifierSum;

    m_cachedGradGradMollifier = true;
  }

  Vector<IndexTM<Real,GLOBALDIM> > retval(GLOBALDIM);

  for (int idir = 0; idir < GLOBALDIM; idir++)
  {
    for (int jdir = 0; jdir < GLOBALDIM; jdir++)
    {
      retval[idir][jdir] = 0.0;
    }
  }

  // Apply gradGrad(mollifier) to the function at the current point (a_point)
  for (BoxIterator bit(m_sampleBox); bit.ok(); ++bit)
  {
    IntVect iv = bit();
    RealVect point;

    for (int idir = 0; idir < SpaceDim; idir++)
    {
      point[idir] = a_point[idir] + m_min + iv[idir] * m_dx;
    }

    Real value = m_impFunc->value(point);

    int index = 0;
    for (int idir = 0; idir < GLOBALDIM; idir++)
    {
      for (int jdir = 0; jdir < GLOBALDIM; jdir++)
      {
        retval[idir][jdir] += m_sampledGradGradMollifier(iv,index) * value;
        index++;
      }
    }
  }

  return retval;
}

Vector<IndexTM<Real,GLOBALDIM> > MollifyIF::gradNormal(const IndexTM<Real,GLOBALDIM>& a_point) const
{
  Vector<IndexTM<Real,GLOBALDIM> > retval(GLOBALDIM);

  // Get the gradient
  IndexTM<Real,GLOBALDIM> curGrad = grad(a_point);

  // Get the gradient of the gradient
  Vector<IndexTM<Real,GLOBALDIM> > curGradGrad = gradGrad(a_point);

  // Get the square of the length of the gradient
  Real sum2 = 0.0;
  for (int idir = 0; idir < GLOBALDIM; idir++)
  {
    sum2 += curGrad[idir] * curGrad[idir];
  }

  // Get the cube of length of the gradient
  Real sum3 = pow(sum2,(Real)1.5);

  // Compute the deriatives of the normal
  for (int idir = 0; idir < GLOBALDIM; idir++)
  {
    for (int jdir = 0; jdir < GLOBALDIM; jdir++)
    {
      // A sum of first and second derivatives needed
      Real mixedSum = 0.0;

      for (int kdir = 0; kdir < GLOBALDIM; kdir++)
      {
        mixedSum += curGrad[kdir] * curGradGrad[kdir][jdir];
      }

      if (sum3 != 0) {
        retval[idir][jdir] = (curGradGrad[idir][jdir]*sum2 - curGrad[idir]*mixedSum) / sum3;
      }
      else
      {
        retval[idir][jdir] = 0.0;
      }
    }
  }

  return retval;
}

BaseIF* MollifyIF::newImplicitFunction() const
{
  MollifyIF* mollifyPtr = new MollifyIF(*m_impFunc, *m_mollifier,
                                        m_min, m_max, m_numPts);

  return static_cast<BaseIF*>(mollifyPtr);
}

#include "NamespaceFooter.H"

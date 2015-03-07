#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

#include "PolynomialIF.H"

#include "NamespaceHeader.H"

PolynomialIF::PolynomialIF(const Vector<PolyTerm>& a_polynomial,
                           const bool&             a_inside)
{
  // Copy polynomial
  int size = a_polynomial.size();

  m_polynomial.resize(size);

  for (int iterm = 0; iterm < size; iterm++)
  {
    m_polynomial[iterm] = a_polynomial[iterm];
  }

  // Save inside flag
  m_inside = a_inside;
}

PolynomialIF::PolynomialIF(const PolynomialIF& a_inputIF)
{
  // Copy polynomial
  int size = a_inputIF.m_polynomial.size();

  m_polynomial.resize(size);

  for (int iterm = 0; iterm < size; iterm++)
  {
    m_polynomial[iterm] = a_inputIF.m_polynomial[iterm];
  }

  // Save inside flag
  m_inside = a_inputIF.m_inside;
}

PolynomialIF::~PolynomialIF()
{
}

void PolynomialIF::GetParams(Vector<PolyTerm>& a_polynomial,
                             bool&             a_inside) const
{
  // Copy polynomial
  int size = m_polynomial.size();

  a_polynomial.resize(size);

  for (int iterm = 0; iterm < size; iterm++)
  {
    a_polynomial[iterm] = m_polynomial[iterm];
  }

  // Save inside flag
  a_inside = m_inside;
}

void PolynomialIF::SetParams(const Vector<PolyTerm>& a_polynomial,
                             const bool&             a_inside)
{
  // Copy polynomial
  int size = a_polynomial.size();

  m_polynomial.resize(size);

  for (int iterm = 0; iterm < size; iterm++)
  {
    m_polynomial[iterm] = a_polynomial[iterm];
  }

  // Save inside flag
  m_inside = a_inside;
}

Real PolynomialIF::value(const RealVect& a_point) const
{
  Real retval;

  int size = m_polynomial.size();

  // Evaluate the polynomial
  retval = 0.0;
  for (int iterm = 0; iterm < size; iterm++)
  {
    Real cur;

    cur = m_polynomial[iterm].coef;
    for (int idir = 0; idir < SpaceDim; idir++)
    {
      cur *= pow(a_point[idir],m_polynomial[iterm].powers[idir]);
    }

    retval += cur;
  }

  // Change the sign to change inside to outside
  if (!m_inside)
  {
    retval = -retval;
  }

  return retval;
}

BaseIF* PolynomialIF::newImplicitFunction() const
{
  PolynomialIF* polynomialPtr = new PolynomialIF(m_polynomial,
                                                 m_inside);

  return static_cast<BaseIF*>(polynomialPtr);
}

#include "NamespaceFooter.H"

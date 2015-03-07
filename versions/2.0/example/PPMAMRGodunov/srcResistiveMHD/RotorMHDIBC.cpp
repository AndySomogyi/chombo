#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

#include "LoHiSide.H"
#include "LoHiCenter.H"

#include "RotorMHDIBC.H"
#include "RotorMHDIBCF_F.H"

// Null constructor
RotorMHDIBC::RotorMHDIBC()
{
  m_isFortranCommonSet = false;
}

// Constructor which defines parameters used by Fortran routines
RotorMHDIBC::RotorMHDIBC(Real&       a_smallPressure,
                         const Real& a_gamma,
                         const Real& a_mu,
                         const Real& a_eta,
                         const Real& a_kappa)
{
  setFortranCommon(a_smallPressure,a_gamma,a_mu,a_eta,a_kappa);
}

RotorMHDIBC::~RotorMHDIBC()
{
}

// Sets parameters in a common block used by Fortran routines:
//   a_smallPressure - Lower limit for pressure (returned)
//   a_gamma         - Gamma for polytropic, gamma-law gas
void RotorMHDIBC::setFortranCommon(Real&       a_smallPressure,
                                   const Real& a_gamma,
				   const Real& a_mu,
				   const Real& a_eta,
				   const Real& a_kappa)
{
  CH_assert(m_isFortranCommonSet == false);

  FORT_ROTORMHDSETF(CHF_REAL(a_smallPressure),
                    CHF_CONST_REAL(a_gamma),
                    CHF_CONST_REAL(a_mu),
                    CHF_CONST_REAL(a_eta),
                    CHF_CONST_REAL(a_kappa));

  m_isFortranCommonSet = true;
}

// Factory method - this object is its own factory:
//   Return a pointer to a new PhysIBC object with m_isDefined = false (i.e.,
//   its define() must be called before it is used) and m_isFortranCommonSet
//   set to value of m_isFortranCommonset in the current (factory) object.
PhysIBC* RotorMHDIBC::new_physIBC()
{
  RotorMHDIBC* retval = new RotorMHDIBC();
  retval->m_isFortranCommonSet = m_isFortranCommonSet;

  return static_cast<PhysIBC*>(retval);
}

// Set up initial conditions
void RotorMHDIBC::initialize(LevelData<FArrayBox>& a_U)
{
  CH_assert(m_isFortranCommonSet == true);
  CH_assert(m_isDefined == true);

  DataIterator dit = a_U.boxLayout().dataIterator();

  // Iterator of all grids in this level
  for (dit.begin(); dit.ok(); ++dit)
  {
    // Storage for current grid
    FArrayBox& U = a_U[dit()];

    // Box of current grid
    Box uBox = U.box();
    //    uBox &= m_domain;

    // Set up initial condition in this grid
    FORT_ROTORMHDINITF(CHF_FRA(U),
                       CHF_CONST_REAL(m_dx),
                       CHF_BOX(uBox));
  }
}

// Set boundary fluxes
void RotorMHDIBC::primBC(FArrayBox&            a_WGdnv,
                         const FArrayBox&      a_Wextrap,
                         const FArrayBox&      a_W,
                         const int&            a_dir,
                         const Side::LoHiSide& a_side,
                         const Real&           a_time)
{
  CH_assert(m_isFortranCommonSet == true);
  CH_assert(m_isDefined == true);

  // In periodic case, this doesn't do anything
  if (!m_domain.isPeriodic(a_dir))
  {
    int lohisign;
    Box WGdnvBox = a_WGdnv.box();
    Box tmp = WGdnvBox;

    // Determine which side and thus shifting directions
    lohisign = sign(a_side);
    tmp.shiftHalf(a_dir,lohisign);

    // Is there a domain boundary next to this grid
    if (!m_domain.contains(tmp))
    {
      tmp &= m_domain;

      Box boundaryBox;

      // Find the strip of cells next to the domain boundary
      if (a_side == Side::Lo)
      {
        boundaryBox = bdryLo(tmp,a_dir);
      }
      else
      {
        boundaryBox = bdryHi(tmp,a_dir);
      }

      // Set the time varying boundary fluxes
      FORT_ROTORMHDPRIMBCF(CHF_FRA(a_WGdnv),
                           CHF_CONST_FRA(a_Wextrap),
                           CHF_CONST_REAL(a_time),
                           CHF_CONST_INT(lohisign),
                           CHF_CONST_REAL(m_dx),
                           CHF_CONST_INT(a_dir),
                           CHF_BOX(boundaryBox));
    }
  }
}

// Set boundary slopes:
//   The boundary slopes in a_dW are already set to one sided difference
//   approximations.  If this function doesn't change them they will be
//   used for the slopes at the boundaries.
void RotorMHDIBC::setBdrySlopes(FArrayBox&       a_dW,
                                const FArrayBox& a_W,
                                const int&       a_dir,
                                const Real&      a_time)
{
  CH_assert(m_isFortranCommonSet == true);
  CH_assert(m_isDefined == true);

  // In periodic case, this doesn't do anything
  if (!m_domain.isPeriodic(a_dir))
  {
    Box loBox,hiBox,centerBox,domain;
    int hasLo,hasHi;
    Box slopeBox = a_dW.box();
    slopeBox.grow(a_dir,1);

    // Generate the domain boundary boxes, loBox and hiBox, if there are
    // domain boundarys there
    loHiCenter(loBox,hasLo,hiBox,hasHi,centerBox,domain,
               slopeBox,m_domain,a_dir);

    // Set the boundary slopes if necessary
    if ((hasLo != 0) | (hasHi != 0))
    {
      FORT_ROTORMHDSLOPEBCSF(CHF_FRA(a_dW),
                             CHF_CONST_FRA(a_W),
                             CHF_CONST_REAL(m_dx),
                             CHF_CONST_INT(a_dir),
                             CHF_BOX(loBox),
                             CHF_CONST_INT(hasLo),
                             CHF_BOX(hiBox),
                             CHF_CONST_INT(hasHi));
    }
  }
}

// Adjust boundary fluxes to account for artificial viscosity
// Do nothing for this problem
void RotorMHDIBC::artViscBC(FArrayBox&       a_F,
                            const FArrayBox& a_U,
                            const FArrayBox& a_divVel,
                            const int&       a_dir,
                            const Real&      a_time)
{
}

// Set up initial conditions
void RotorMHDIBC::initializePhi(LevelData<FArrayBox>& a_U)
{
  // cout << "In initialize Phi" << endl;
  CH_assert(m_isDefined == true);
  DataIterator dit = a_U.boxLayout().dataIterator();

  // Iterator of all grids in this level
  for (dit.begin(); dit.ok(); ++dit)
  {
    // Storage for current grid
    FArrayBox& U = a_U[dit()];

    // Box of current grid
    Box uBox = U.box();
    uBox &= m_domain;

    // cout << "In initialize Phi: uBox" << uBox << endl;

    // Set up initial condition in this grid
    FORT_ROTORMHDINITPHIF(CHF_FRA(U),
                          CHF_CONST_REAL(m_dx),
                          CHF_BOX(uBox));
  }
}

// Set boundary conditions
void RotorMHDIBC::phiBC(FArrayBox&           a_U,
                        const ProblemDomain& a_domain,
                        const Real&          a_dxLevel)
{
  CH_assert(m_isDefined == true);
  // cout << "PhiBC" << endl;

  // int sign;
  Box uBox = a_U.box();
  Side::LoHiSide side;

  // Determine which side and thus shifting directions
  // if (a_side == Side::Lo)
  // {
  //   sign = -1;
  // }
  // else
  // {
  //   sign = 1;
  // }
  // cout << "phiBC: uBox=" << uBox << endl;
  // cout << "phiBC: m_domain=" << m_domain << endl;

  for (int dir1 = 0; dir1 < SpaceDim; ++dir1)
  {
    Box tmp = uBox;
    // SideIterator sit;
    // cout << "phiBC: dir" << dir1  << endl;
    // cout << "phiBC: tmp0 boundaryBox=" << tmp << endl;
    tmp.shift(dir1,1);

    // Is there a domain boundary next to this grid
    // cout << "phiBC:tmp1 boundaryBox=" << tmp << endl;
    // cout << "phiBC:tmp1 contain=" << m_domain.contains(tmp) << endl;
    if (!a_domain.contains(tmp))
    {
      // cout << "phiBC:m_domain does not contain tmp=" << tmp << endl;
      tmp &= a_domain;

      // cout << "phiBC:tmp2 boundaryBox=" << tmp << endl;
      Box boundaryBox;

      // Find the strip of cells next to the domain boundary
      boundaryBox = bdryLo(tmp,dir1);
      boundaryBox.shift(dir1,-1);
      // cout << "Lo boundaryBox=" << boundaryBox << endl;

      FORT_ROTORMHDPHIBCLOF(CHF_FRA(a_U),
                            CHF_CONST_REAL(a_dxLevel),
                            CHF_CONST_INT(dir1),
                            CHF_BOX(boundaryBox));

    }

    side=Side::Hi;
    // cout << "phiBC: dir=" << dir1 << endl;
    // cout << "phiBC: tmp0 boundaryBox=" << tmp << endl;

    tmp=uBox;
    tmp.shift(dir1,-1);

    // Is there a domain boundary next to this grid
    // cout << "phiBC:tmp1 boundaryBox=" << tmp << endl;

    if (!a_domain.contains(tmp))
    {
      tmp &= a_domain;
      // cout << "phiBC:tmp2 boundaryBox=" << tmp << endl;

      Box boundaryBox;

      // Find the strip of cells next to the domain boundary
      boundaryBox = bdryHi(tmp,dir1);
      // cout << "Hi boundaryBox=" << boundaryBox << endl;

      FORT_ROTORMHDPHIBCHIF(CHF_FRA(a_U),
                            CHF_CONST_REAL(a_dxLevel),
                            CHF_CONST_INT(dir1),
                            CHF_BOX(boundaryBox));
    }
  }
}

// Set boundary fluxes
void RotorMHDIBC::fluxBC(FArrayBox&            a_F,
                         const FArrayBox&      a_W,
                         const FArrayBox&      a_Wextrap,
                         const int&            a_dir,
                         const Side::LoHiSide& a_side,
                         const Real&           a_time)
{
  CH_assert(m_isFortranCommonSet == true);
  CH_assert(m_isDefined == true);

  int sign;
  Box FBox = a_F.box();
  Box tmp = FBox;

  // Determine which side and thus shifting directions
  if (a_side == Side::Lo)
  {
    sign = -1;
  }
  else
  {
    sign = 1;
  }

  // cout << "tmp0 boundaryBox=" << tmp << endl;
  tmp.shiftHalf(a_dir,sign);

  // Is there a domain boundary next to this grid
  // cout << "a_side=" << a_side << endl;
  // cout << "tmp1 boundaryBox=" << tmp << endl;
  if (!m_domain.contains(tmp))
  {
    tmp &= m_domain;

    // cout << "tmp2 boundaryBox=" << tmp << endl;
    Box boundaryBox;

    // Find the strip of cells next to the domain boundary
    if (a_side == Side::Lo)
    {
      boundaryBox = bdryLo(tmp,a_dir);
      // cout << "Lo boundaryBox=" << boundaryBox << endl;
    }
    else
    {
      boundaryBox = bdryHi(tmp,a_dir);
      // cout << "Hi boundaryBox=" << boundaryBox << endl;
    }

    // Shift things to all line up correctly
    boundaryBox.shiftHalf(a_dir,-sign);
    a_F.shiftHalf(a_dir,-sign);

    // Set the solid wall boundary fluxes- not for the Rotor problem
    // FORT_ROTORMHDSOLIDBCF(CHF_FRA(a_F),
    //                       CHF_CONST_FRA(a_Wextrap),
    //                       CHF_CONST_REAL(m_dx),
    //                       CHF_CONST_INT(a_dir),
    //                       CHF_BOX(boundaryBox));

    // Set the time varying boundary fluxes
    FORT_ROTORMHDBCF(CHF_FRA(a_F),
                     CHF_CONST_FRA(a_Wextrap),
    //                 CHF_CONST_FRA(a_W), RS: changed W to Wextrap
                     CHF_CONST_REAL(a_time),
                     CHF_CONST_INT(sign),
                     CHF_CONST_REAL(m_dx),
                     CHF_CONST_INT(a_dir),
                     CHF_BOX(boundaryBox));

    // FORT_ROTORMHDSOLIDBCF(CHF_FRA(a_F),
    //                       CHF_CONST_FRA(a_Wextrap),
    //                       CHF_CONST_REAL(a_time),
    //                       CHF_CONST_INT(sign),
    //                       CHF_CONST_REAL(m_dx),
    //                       CHF_CONST_INT(a_dir),
    //                       CHF_BOX(boundaryBox));

    // Shift returned fluxes to be face centered
    a_F.shiftHalf(a_dir,sign);
  }
}

// Set boundary fluxes
void RotorMHDIBC::parabolicFluxBC(FArrayBox&            a_F,
                                  FArrayBox&            a_W,
                                  const int&            a_dir,
                                  const Side::LoHiSide& a_side,
                                  const Real&           a_time)
{
  CH_assert(m_isFortranCommonSet == true);
  CH_assert(m_isDefined == true);

  int sign;
  Box FBox = a_F.box();
  Box tmp = FBox;

  // Determine which side and thus shifting directions
  if (a_side == Side::Lo)
  {
    sign = -1;
  }
  else
  {
    sign = 1;
  }

  // cout << "tmp0 boundaryBox=" << tmp << endl;
  tmp.shiftHalf(a_dir,sign);

  // Is there a domain boundary next to this grid
  // cout << "a_side=" << a_side << endl;
  // cout << "tmp1 boundaryBox=" << tmp << endl;

  if (!m_domain.contains(tmp))
  {
    tmp &= m_domain;

    // cout << "tmp2 boundaryBox=" << tmp << endl;
    Box boundaryBox;

    // Find the strip of cells next to the domain boundary
    if (a_side == Side::Lo)
    {
      boundaryBox = bdryLo(tmp,a_dir);
      // cout << "Lo boundaryBox=" << boundaryBox << endl;
    }
    else
    {
      boundaryBox = bdryHi(tmp,a_dir);
      // cout << "Hi boundaryBox=" << boundaryBox << endl;
    }

    // Shift things to all line up correctly
    boundaryBox.shiftHalf(a_dir,-sign);
    a_F.shiftHalf(a_dir,-sign);

    // cout << "Calling MHDPARABOLICFLUX" << endl;
    FORT_ROTORMHDPARABOLICFLUXBCF(CHF_FRA(a_F),
                                  CHF_FRA(a_W),
                                  CHF_CONST_REAL(a_time),
                                  CHF_CONST_INT(sign),
                                  CHF_CONST_REAL(m_dx),
                                  CHF_CONST_INT(a_dir),
                                  CHF_BOX(boundaryBox));

    // Shift returned fluxes to be face centered
    a_F.shiftHalf(a_dir,sign);
  }
}

#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

#ifndef _SCALARDIRICHLETBC_H_
#define _SCALARDIRICHLETBC_H_

#include "Box.H"
#include "FArrayBox.H"
#include "Vector.H"
#include "LoHiSide.H"
#include "GhostBC.H"


/// implements dirichlet BCs for scalars
class scalarDirichletBC  : public BoxGhostBC
{

public:
  friend class DomainGhostBC;

  /// 
  scalarDirichletBC(int            a_dir,
                    Side::LoHiSide a_side);
  
  /// 
  scalarDirichletBC(int             a_dir,
                    Side::LoHiSide  a_side, 
                    const Interval& a_comps);

  ///
  virtual ~scalarDirichletBC();

  ///
  scalarDirichletBC();
  
  /// sets Dirichlet Boundary value
  void setVal(const Real a_bcVal);


private:
  
  /// dirichlet value at edge
  Real m_BCVal;
  
  /// fills dirichletfactor with 1.0 and inhomVal with BCVal
  virtual void
  fillBCValues(FArrayBox&           a_neumfac,
               FArrayBox&           a_dirfac,
               FArrayBox&           a_inhmval,
               Real                 a_dx, 
               const ProblemDomain& a_domain) const;
  
  /// fills dirichletfactor with 1.0 and inhomVal with BCVal
  virtual void
  fillBCValues(FArrayBox& a_neumfac,
               FArrayBox& a_dirfac,
               FArrayBox& a_inhmval,
               Real       a_dx, 
               const Box& a_domain) const;

  ///
  virtual BoxGhostBC* new_boxghostbc() const;

private:
  scalarDirichletBC(const scalarDirichletBC& a_bcin) {;}

};

/// Class to set inflow BC's for scalars
class scalarInflowBC : public BoxGhostBC
{

public: 
  friend class DomainGhostBC;

  ///
  scalarInflowBC();

  /// 
  scalarInflowBC(int            a_dir,
                 Side::LoHiSide a_side);

  ///
  scalarInflowBC(int             a_dir,
                 Side::LoHiSide  a_side,
                 const Interval& a_comps);

  ///
  virtual ~scalarInflowBC();

  /// set which scalar this is
  void setScalarType(int a_scalType);
  
  ///
  int scalarType() const;

  /// 
  void isHomogeneous(bool a_isHomogeneous);
  
  /// 
  bool isHomogeneous() const;


private:

  ///
  int m_scalType;

  /// 
  bool m_isHomogeneous;

  /// fills dirichlet factor with 1.0 and inhomval with inflow profile
  virtual void 
  fillBCValues(FArrayBox&           a_neumfac,
               FArrayBox&           a_dirfac,
               FArrayBox&           a_inhmval,
               Real                 a_dx,
               const ProblemDomain& a_domain) const;


  /// fills dirichlet factor with 1.0 and inhomval with inflow profile
  virtual void 
  fillBCValues(FArrayBox& a_neumfac,
               FArrayBox& a_dirfac,
               FArrayBox& a_inhmval,
               Real       a_dx,
               const Box& a_domain) const;

  ///
  virtual BoxGhostBC* new_boxghostbc() const;
};

#endif

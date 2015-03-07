#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif
#ifndef _VCLEVELOP_H_
#define _VCLEVELOP_H_

#include <iostream>
#include <cmath>
#include <assert.h>
#include <cstdlib>
#include "ParentLevelOp.H" 
#include "FluxBox.H"

/// Variable coefficient operator class

/**
 Variable coefficient operator class in the dimensionally split
 form suitable for ADI-like line MG smoothing
 (Yuri Omelchenko, 08.20.02 )

LofPhi := sum_idir {L_idim(i)}+m_diagn(i)*phi(i)
L_idir(i) := m_lower(i)*phi(i-1)+m_centr(i)*phi(i)+m_upper(i)*phi(i+1)
 */
class VCDivOp: public ParentLevelOp
{
public:
  
  ///
  LevelOp* new_levelop() const;

  ///
  /**
     Full define function.  makes all coarse-fine
     information and sets internal variables 
  */
  void define(
              const DisjointBoxLayout& Ba,
              const DisjointBoxLayout* base_ba,
              Real  DxLevel, 
              int refratio,
              const Box& a_domain,
              bool a_homogeneousOnly=false,
              int ncomp = 1);

  ///
  /**
     Full define function.  makes all coarse-fine
     information and sets internal variables 
  */
  void define(
              const DisjointBoxLayout& Ba,
              const DisjointBoxLayout* base_ba,
              Real  DxLevel, 
              int refratio,
              const ProblemDomain& a_probdomain,
              bool a_homogeneousOnly=false,
              int ncomp = 1);
  ///
  void define(const LevelOp* opfine, int a_reftoFine);

  /// define coefficients
  void defineVC(const DisjointBoxLayout& a_grids, int a_ncomp);

  /// 
  VCDivOp();

  /// 
  ~VCDivOp();

  

  ///
  /**
     Smoother. 
     Assumes that problem has already been put in 
     residual correction form, 
     so that C/F boundary conditions are homogeneous.
  */
  virtual
  void smooth(
              LevelData<FArrayBox>& phi,
              const LevelData<FArrayBox>& rhs
              );

  ///
  void applyOpI(
                LevelData<FArrayBox>& phi,
                const LevelData<FArrayBox>* phiCoarse,
                LevelData<FArrayBox>& LOfPhi
                );

  ///
  void applyOpIcfHphys(
                       LevelData<FArrayBox>& phi,
                       const LevelData<FArrayBox>* phiCoarse,
                       LevelData<FArrayBox>& LOfPhi
                       );
  
  ///
  void applyOpH(
                LevelData<FArrayBox>& phi,
                LevelData<FArrayBox>& LOfPhi
                );

  ///
  void applyOpHcfIphys(
                       LevelData<FArrayBox>& phi,
                       LevelData<FArrayBox>& LOfPhi
                       );

  ///
  virtual
  void getFlux(
              FArrayBox& a_fineFlux,
              const FArrayBox& a_data,
	      const DataIndex& a_datInd, 
              int a_dir);
  ///
  virtual
  void levelPreconditioner(
            LevelData<FArrayBox>& a_phihat, 
            const LevelData<FArrayBox>& a_rhshat
            );

  /// compute default coefficients
  virtual
  void computeDefaultVC(Real a_dxLevel);

  /// compute coefficients through via external a_levelVC()
  virtual
  void computeVC(void (*a_levelVC)(VCDivOp &a_vcdivop));

  /// set coefficients for diagonal
  virtual
  void setDiagCoeff(const LevelData<FArrayBox>& a_diag);

  /// set coefficients for div(grad) term.
  virtual 
  void setDivCoeff(const LevelData<FluxBox>& a_coeff,
                   Real a_dxLevel);

  /// view coefficients
  virtual
  void viewOp(int idir, const Box& viewbox) const;
 
  ///
  bool isVCDefined() const;

  ///
  bool isDiagSet() const;

  ///
  bool isDivSet() const;

  ///
  void needNewLambda();

  ///
  LevelData<FArrayBox> m_lower[CH_SPACEDIM];
  ///
  LevelData<FArrayBox> m_centr[CH_SPACEDIM];
  ///
  LevelData<FArrayBox> m_upper[CH_SPACEDIM];
  ///
  LevelData<FArrayBox> m_diagn;

protected:

  ///
  virtual
  void setDefaultValues();

  ///
  virtual
  ///
  void clearMemory();

  ///
  bool m_isVCDefined;

  ///
  bool m_isDiagSet;

  ///
  bool m_isDivSet;

  ///
  bool m_isLambdaSet;

  /// block and line smoothers
  void levelGSRB(
            LevelData<FArrayBox>& a_phi, 
            const LevelData<FArrayBox>& a_rhs
            );

  ///
  void levelADG(
            LevelData<FArrayBox>& a_phi, 
            const LevelData<FArrayBox>& a_rhs,
            Real a_rho=1.0, int a_kgs=1
            );

private:
  /// scratch array (used by GSRB and ADG)
  LevelData<FArrayBox> m_lofphi;
  ///
  LevelData<FArrayBox> m_lambda;

  ///
  /** This function must be called to compute m_lambda 
      following every time new VCs are evaluated so that
      GSRB uses an updated m_lambda
  */
  void computeInverseDiag(); 

  // disallow these 
  VCDivOp(const VCDivOp&); 
  VCDivOp& operator=(const VCDivOp&); 
};

#endif

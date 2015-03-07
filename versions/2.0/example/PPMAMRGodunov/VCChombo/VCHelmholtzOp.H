#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

// VCHelmholtzOp.H
// petermc, 22 Nov 2002

#ifndef VCHELMHOLTZOP_H
#define VCHELMHOLTZOP_H

#include "VCDivOp.H"
#include "BaseHelmholtzOp.H"

/// Variable-coefficient Helmholtz operator class
/**
   This class implements the operator

   LofPhi := (alpha*I + div beta*grad)phi

   where alpha is cell-centered and beta is face-centered.
   More specifically, at point i:

   LofPhi(i) := alpha(i)*phi(i) + sum_d {L_d(i+e_d/2) - L_d(i-e_d/2)}/dx \\
   with e_d the unit vector in direction d, and the face-centered \\
   L_d(i+e_d/2) = beta[d](i+e_d/2) * (phi(i+e_d) - phi(i))/dx .
 */

class VCHelmholtzOp: public BaseHelmholtzOp
{
public:

  ///
  /**
     Full define function.  makes all coarse-fine
     information and sets internal variables 
  */
  void define(
              const DisjointBoxLayout& a_Ba,
              const DisjointBoxLayout* a_baseBa,
              Real  a_dxLevel,
              int a_refratio,
              const Box& a_domain,
              bool a_homogeneousOnly=false,
              int a_ncomp = 1);

  ///
  /**
     Full define function.  makes all coarse-fine
     information and sets internal variables 
  */
  void define(
              const DisjointBoxLayout& a_Ba,
              const DisjointBoxLayout* a_baseBa,
              Real  a_dxLevel,
              int a_refratio,
              const ProblemDomain& a_probdomain,
              bool a_homogeneousOnly=false,
              int a_ncomp = 1);

  ///
  void define(const LevelOp* opfine, int a_reftoFine);

  ///
  /**
     Smoother. 
     Assumes that problem has already been put in 
     residual correction form, 
     so that C/F boundary conditions are homogeneous.
  */
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
  /**
     Coarse / Fine interpolation operator. fills ghost 
     cells of phi.
  */
  void CFInterp(
                LevelData<FArrayBox>& phi,
                const LevelData<FArrayBox>& phiCoarse
                );

  ///
  /** homogeneous coarse/fine interpolation operator
   */
  void homogeneousCFInterp(LevelData<FArrayBox>& phi);

  ///
  bool isDefined() const;

  ///
  /**
     // This does a GSRB Pre/Conditioned BiCGStab on a level
     // for the bottom solver.  
  */
  void bottomSmoother(LevelData<FArrayBox>& phi,
                      const LevelData<FArrayBox>& rhs);

  ///
  void getFlux(
              FArrayBox& a_fineFlux,
              const FArrayBox& a_data,
	      const DataIndex& a_datInd, 
              int a_dir);
  ///
  void levelPreconditioner(
            LevelData<FArrayBox>& a_phihat, 
            const LevelData<FArrayBox>& a_rhshat
            );

  ///
  /**
     Set coefficients for diagonal, alpha, in (alpha*I + div beta*grad)phi
   */
  void setAlpha(const LevelData<FArrayBox>& a_diag);

  ///
  /**
     Set coefficients for div(grad), beta, in (alpha*I + div beta*grad)phi
   */
  void setBeta(const LevelData<FluxBox>& a_coeff,
               Real a_dxLevel);

  ///
  /**
     In BaseHelmholtzOp:
     Scale the helmholtz coeficient alpha in (alpha*I + div beta*grad)phi
     (alpha := a_scale*alpha)
  */
  void scaleAlpha(Real a_scale);

  ///
  /**
     In BaseHelmholtzOp:
     Scale the helmholtz coefficient beta in (alpha*I + div beta*grad)phi
     (beta := a_scale*beta)
  */
  void scaleBeta(Real a_scale);

  ///
  /**virtual constructor workaround
     coefficient  and bcs must be defined before called.
  */
  LevelOp* new_levelop() const;

  /// 
  /** From BaseHelmholtzOp:
      set (or reset) the physical boundary conditions
   */
  void setDomainGhostBC(const DomainGhostBC& a_dombc);

  /// set convergence metrics (passing through to bottom smoother)
  void setConvergenceMetric(Real a_metric, int a_comp);
  

  /// 
  VCHelmholtzOp() { m_vcdivop = NULL; };

  /// 
  ~VCHelmholtzOp() { delete m_vcdivop; };


protected:

  //
  VCDivOp* m_vcdivop;

private:
  
};

#endif


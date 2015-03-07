#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

#ifndef _VCAMRPOISSONOP_H_
#define _VCAMRPOISSONOP_H_

#include "AMRMultiGrid.H"
#include "REAL.H"
#include "Box.H"
#include "LevelDataOps.H"
#include "BCFunc.H"
#include "FArrayBox.H"
#include "FluxBox.H"
#include "CFIVS.H"
#include "QuadCFInterp.H"
#include "CoarseAverage.H"
#include "LevelFluxRegister.H"
#include "AMRTGA.H"

#include "NamespaceHeader.H"

///
/**
   Operator for solving variable-coefficient
   (alpha * aCoef(x) * I - beta * Div(bCoef(x) . Grad)) phi = rho
   over an AMR hierarchy.
*/
class VCAMRPoissonOp : public TGAHelmOp<LevelData<FArrayBox> >
{
public:

  ///for tga stuff
  void setAlphaAndBeta(const Real& a_alpha,
                       const Real& a_beta);

  ///for tga stuff
  void diagonalScale(LevelData<FArrayBox> & a_rhs)
  {
    for(DataIterator dit = a_rhs.dataIterator(); dit.ok(); ++dit)
      {
        a_rhs[dit()] *= (*m_aCoef)[dit()];
      }
  }
  ///also calls reset lambda
  void setCoefs(const RefCountedPtr<LevelData<FArrayBox> >& a_aCoef,
                const RefCountedPtr<LevelData<FluxBox>   >& a_bCoef,
                const Real& a_alpha,
                const Real& a_beta);

  ///should be called every time coefs are changed.
  void resetLambda();

  ///
  /**
   */
  VCAMRPoissonOp()
  {
    ;
  }

  ///
  /**
   */
  virtual ~VCAMRPoissonOp()
  {
    ;
  }

  ///
  /**
     define function for AMRLevelOp which has no finer or coarser AMR level
  */
  void define(const DisjointBoxLayout& a_grids,
              const Real&              a_dx,
              const ProblemDomain&     a_domain,
              BCHolder                   a_bc);

  ///
  /**
     define function for AMRLevelOp which has no finer AMR level
  */
  void define(const DisjointBoxLayout& a_grids,
              const DisjointBoxLayout& a_baseBAPtr,
              const Real&              a_dxLevel,
              int                      a_refRatio,
              const ProblemDomain&     a_domain,
              BCHolder                   a_bc);

  ///
  /** full define function for AMRLevelOp with both coarser and finer levels */
  void define(const DisjointBoxLayout& a_grids,
              const DisjointBoxLayout& a_gridsFiner,
              const DisjointBoxLayout& a_gridsCoarser,
              const Real&              a_dxLevel,
              int                      a_refRatio,
              int                      a_refRatioFiner,
              const ProblemDomain&     a_domain,
              BCHolder                   a_bc);

  ///
  /** full define function for AMRLevelOp with finer levels, but no coarser */
  void define(const DisjointBoxLayout& a_grids,
              const DisjointBoxLayout& a_gridsFiner,
              const Real&              a_dxLevel,
              int                      a_refRatio, //  dummy argument, send in 1
              int                      a_refRatioFiner,
              const ProblemDomain&     a_domain,
              BCHolder                   a_bc);

  /// Compute lambda once alpha, aCoef, beta, bCoef are defined
  void computeLambda();

  ///
  virtual void residual(LevelData<FArrayBox>&       a_lhs,
                        const LevelData<FArrayBox>& a_phi,
                        const LevelData<FArrayBox>& a_rhs,
                        bool                        a_homogeneous = false);

  ///
  virtual void preCond(LevelData<FArrayBox>&       a_correction,
                       const LevelData<FArrayBox>& a_residual);

  ///
  virtual void applyOp(LevelData<FArrayBox>&       a_lhs,
                       const LevelData<FArrayBox>& a_phi,
                       bool                        a_homogeneous = false);
  ///
  virtual void applyOpNoBoundary(LevelData<FArrayBox>&       a_lhs,
                                 const LevelData<FArrayBox>& a_phi);


  ///
  virtual void create(LevelData<FArrayBox>&       a_lhs,
                      const LevelData<FArrayBox>& a_rhs);

  ///
  virtual void createCoarsened(LevelData<FArrayBox>&       a_lhs,
                               const LevelData<FArrayBox>& a_rhs,
                               const int&                  a_refRat);

  ///
  virtual void assign(LevelData<FArrayBox>&       a_lhs,
                      const LevelData<FArrayBox>& a_rhs);
  ///
  virtual Real dotProduct(const LevelData<FArrayBox>& a_1,
                          const LevelData<FArrayBox>& a_2);

  ///
  virtual void incr(LevelData<FArrayBox>&       a_lhs,
                    const LevelData<FArrayBox>& a_x,
                    Real                        a_scale);

  ///
  virtual void axby(LevelData<FArrayBox>&       a_lhs,
                    const LevelData<FArrayBox>& a_x,
                    const LevelData<FArrayBox>& a_y,
                    Real                        a_a,
                    Real                        a_b);

  ///
  virtual void scale(LevelData<FArrayBox>& a_lhs,
                     const Real&           a_scale);

  ///
  virtual Real norm(const LevelData<FArrayBox>& a_x,
                    int                         a_ord);

  ///
  virtual void setToZero( LevelData<FArrayBox>& a_x);

  /*@}*/

  /**
     \name MGLevelOp functions */
  /*@{*/

  ///
  virtual void relax(LevelData<FArrayBox>&       a_e,
                     const LevelData<FArrayBox>& a_residual,
                     int                         a_iterations);

  ///
  virtual void createCoarser(LevelData<FArrayBox>&       a_coarse,
                             const LevelData<FArrayBox>& a_fine,
                             bool                        a_ghosted);

  /**
     calculate restricted residual
     a_resCoarse[2h] = I[h->2h] (rhsFine[h] - L[h](phiFine[h])
  */
  virtual void restrictResidual(LevelData<FArrayBox>&       a_resCoarse,
                                LevelData<FArrayBox>&       a_phiFine,
                                const LevelData<FArrayBox>& a_rhsFine);

  /**
     correct the fine solution based on coarse correction
     a_phiThisLevel += I[2h->h](a_correctCoarse)
  */
  virtual void prolongIncrement(LevelData<FArrayBox>&       a_phiThisLevel,
                                const LevelData<FArrayBox>& a_correctCoarse);

  /*@}*/

  /**
     \name AMRLevelOp functions */
  /*@{*/

  /** returns 1 when there are no coarser AMRLevelOp objects */
  virtual int refToCoarser()
  {
    return m_refToCoarser;
  }

  /** a_residual = a_rhs - L(a_phi, a_phiFine, a_phiCoarse) */
  virtual void AMRResidual(LevelData<FArrayBox>&              a_residual,
                           const LevelData<FArrayBox>&        a_phiFine,
                           const LevelData<FArrayBox>&        a_phi,
                           const LevelData<FArrayBox>&        a_phiCoarse,
                           const LevelData<FArrayBox>&        a_rhs,
                           bool                               a_homogeneousPhysBC,
                           AMRLevelOp<LevelData<FArrayBox> >* a_finerOp);

  /** residual assuming no more coarser AMR levels */
  virtual void AMRResidualNC(LevelData<FArrayBox>&              a_residual,
                             const LevelData<FArrayBox>&        a_phiFine,
                             const LevelData<FArrayBox>&        a_phi,
                             const LevelData<FArrayBox>&        a_rhs,
                             bool                               a_homogeneousPhysBC,
                             AMRLevelOp<LevelData<FArrayBox> >* a_finerOp);

  /** a_residual = a_rhs - L(a_phi, a_phiCoarse)  */
  virtual void AMRResidualNF(LevelData<FArrayBox>&       a_residual,
                             const LevelData<FArrayBox>& a_phi,
                             const LevelData<FArrayBox>& a_phiCoarse,
                             const LevelData<FArrayBox>& a_rhs,
                             bool                        a_homogeneousPhysBC);


  /** apply AMR operator, including coarse-fine matching */
  virtual void AMROperator(LevelData<FArrayBox>&              a_LofPhi,
                           const LevelData<FArrayBox>&        a_phiFine,
                           const LevelData<FArrayBox>&        a_phi,
                           const LevelData<FArrayBox>&        a_phiCoarse,
                           bool                               a_homogeneousPhysBC,
                           AMRLevelOp<LevelData<FArrayBox> >* a_finerOp);

  /** AMR operator assuming no more coarser AMR levels */
  virtual void AMROperatorNC(LevelData<FArrayBox>&              a_LofPhi,
                             const LevelData<FArrayBox>&        a_phiFine,
                             const LevelData<FArrayBox>&        a_phi,
                             bool                               a_homogeneousPhysBC,
                             AMRLevelOp<LevelData<FArrayBox> >* a_finerOp);

  /** AMR operator assuming no finer level */
  virtual void AMROperatorNF(LevelData<FArrayBox>&       a_LofPhi,
                             const LevelData<FArrayBox>& a_phi,
                             const LevelData<FArrayBox>& a_phiCoarse,
                             bool                        a_homogeneousPhysBC);

  /** a_resCoarse = I[h-2h]( a_residual - L(a_correction, a_coarseCorrection))
      it is assumed that a_resCoarse has already been filled in with the coarse
      version of AMRResidualNF and that this operation is free to overwrite
      in the overlap regions.
  */

  virtual void AMRRestrict(LevelData<FArrayBox>&       a_resCoarse,
                           const LevelData<FArrayBox>& a_residual,
                           const LevelData<FArrayBox>& a_correction,
                           const LevelData<FArrayBox>& a_coarseCorrection);

  /** a_correction += I[h->h](a_coarseCorrection) */
  virtual void AMRProlong(LevelData<FArrayBox>&       a_correction,
                          const LevelData<FArrayBox>& a_coarseCorrection);

  /** a_residual = a_residual - L(a_correction, a_coarseCorrection) */
  virtual void AMRUpdateResidual(LevelData<FArrayBox>&       a_residual,
                                 const LevelData<FArrayBox>& a_correction,
                                 const LevelData<FArrayBox>& a_coarseCorrection);

  ///
  /**
     compute norm over all cells on coarse not covered by finer
  */
  virtual Real AMRNorm(const LevelData<FArrayBox>& a_coarseResid,
                       const LevelData<FArrayBox>& a_fineResid,
                       const int&                  a_refRat,
                       const int&                  a_ord);

  virtual void write(const LevelData<FArrayBox>* a_data,
                     const char*                 a_filename);

  /*@}*/

  /// Identity operator constant coefficient -- if you change this call resetLambda()
  Real m_alpha;

  /// Identity operator spatially varying coefficient storage (cell-centered) --- if you change this call resetLambda()
  RefCountedPtr<LevelData<FArrayBox> > m_aCoef;

  /// Laplacian operator constant coefficient --- if you change this call resetLambda()
  Real m_beta;

  /// Laplacian operator spatially varying coefficient storage (face-centered) --- if you change this call resetLambda()
  RefCountedPtr<LevelData<FluxBox> > m_bCoef;

  /// Reciprocal of the diagonal entry of the operator matrix
  LevelData<FArrayBox> m_lambda;

  Real m_dxCrse;

  Real          m_dx;
  ProblemDomain m_domain;

protected:
  LevelDataOps<FArrayBox> m_levelOps;

  BCHolder m_bc;

  LayoutData<CFIVS> m_loCFIVS[SpaceDim];
  LayoutData<CFIVS> m_hiCFIVS[SpaceDim];

  Copier       m_exchangeCopier;
  QuadCFInterp m_interpWithCoarser;

  int m_refToCoarser;
  int m_refToFiner;

  void levelGSRB(LevelData<FArrayBox>&       a_e,
                 const LevelData<FArrayBox>& a_residual);

  void levelJacobi(LevelData<FArrayBox>&       a_phi,
                   const LevelData<FArrayBox>& a_rhs);

  void homogeneousCFInterp(LevelData<FArrayBox>& a_phif);

  void homogeneousCFInterp(LevelData<FArrayBox>& a_phif,
                           const DataIndex&      a_datInd,
                           int                   a_idir,
                           Side::LoHiSide        a_hiorlo);

  void interpOnIVSHomo(LevelData<FArrayBox>& a_phif,
                       const DataIndex&      a_datInd,
                       const int             a_idir,
                       const Side::LoHiSide  a_hiorlo,
                       const IntVectSet&     a_interpIVS);

  /// computes flux over face-centered a_facebox.
  void getFlux(FArrayBox&       a_flux,
               const FArrayBox& a_data,
               const FluxBox&   a_bCoef,
               const Box&       a_facebox,
               int              a_dir,
               int              a_ref = 1) const ;

  void reflux(const LevelData<FArrayBox>&        a_phiFine,
              const LevelData<FArrayBox>&        a_phi,
              LevelData<FArrayBox>&              a_residual,
              AMRLevelOp<LevelData<FArrayBox> >* a_finerOp);

  /// utility function which computes operator after all bc's have been set
  void computeOperatorNoBCs(LevelData<FArrayBox>&       a_lhs,
                            const LevelData<FArrayBox>& a_phi);

};

///
/**
   Factory to create VCAMRPoissonOps
*/
class VCAMRPoissonOpFactory: public AMRLevelOpFactory<LevelData<FArrayBox> >
{
public:
  VCAMRPoissonOpFactory();

  virtual ~VCAMRPoissonOpFactory()
  {
    ;
  }

  ///
  /**
     a_coarseDomain is the domain at the coarsest level.
     a_grids is the AMR  hierarchy.
     a_refRatios are the refinement ratios between levels.  The ratio lives
     with the coarser level so a_refRatios[ilev] is the ratio
     between ilev and ilev+1
     a_coarseDx is the grid spacing at the coarsest level.
     BCFunc holds the boundary conditions.
     a_alpha is the identity constant coefficient
     a_beta is the laplacian constant coefficient.
     a_aCoef is the identity spatially varying coefficient
     a_bCoef is the laplacian spatially varying coefficient.
  */
  void define(const ProblemDomain&                           a_coarseDomain,
              const Vector<DisjointBoxLayout>&               a_grids,
              const Vector<int>&                             a_refRatios,
              const Real&                                    a_coarsedx,
              BCHolder                                         a_bc,
              const Real&                                    a_alpha,
              Vector<RefCountedPtr<LevelData<FArrayBox> > >& a_aCoef,
              const Real&                                    a_beta,
              Vector<RefCountedPtr<LevelData<FluxBox> > >&   a_bCoef);

  ///
  virtual VCAMRPoissonOp* MGnewOp(const ProblemDomain& a_FineindexSpace,
                                  int                  a_depth,
                                  bool                 a_homoOnly = true);

  ///
  virtual VCAMRPoissonOp* AMRnewOp(const ProblemDomain& a_indexSpace);

  ///
  virtual int refToFiner(const ProblemDomain& a_domain) const;

  void setDefaultValues();

  int m_coefficient_average_type;

private:
  Vector<ProblemDomain>     m_domains;
  Vector<DisjointBoxLayout> m_boxes;

  Vector<Real> m_dx;
  Vector<int>  m_refRatios; // refinement to next coarser level

  BCHolder m_bc;

  Real m_alpha;
  Real m_beta;

  Vector<RefCountedPtr<LevelData<FArrayBox> > > m_aCoef;
  Vector<RefCountedPtr<LevelData<FluxBox> > >   m_bCoef;

  Vector<RefCountedPtr<LevelData<FArrayBox> > > m_lambda;
};

#include "NamespaceFooter.H"
#endif

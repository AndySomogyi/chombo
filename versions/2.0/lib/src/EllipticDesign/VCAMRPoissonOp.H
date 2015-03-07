#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

// BVS, June 18, 2003
// revised by dtg in 2006


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
#include "NamespaceHeader.H"

///
/**
   Operator for solving variable-coefficient 
   (alpha(x) I + mu*Div(Beta(x) Grad))(phi) = rho
   over an AMR hierarchy.
*/
class VCAMRPoissonOp : public AMRLevelOp<LevelData<FArrayBox> >
{
public:
  ///
  /**
   */
  VCAMRPoissonOp(){ ; }

  ///
  /**
   */
  virtual ~VCAMRPoissonOp() { ; }


  ///
  /**
     define function for AMRLevelOp which has no finer or coarser AMR level 
   */
  void define(const DisjointBoxLayout& a_grids, 
              const Real&              a_dx,
              const ProblemDomain&     a_domain, 
              BCFunc                   a_bc);     

  ///
  /**
     define function for AMRLevelOp which has no finer AMR level 
  */
  void define(const DisjointBoxLayout& a_grids,
              const DisjointBoxLayout& a_baseBAPtr,
              const Real&              a_dxLevel,
              int                      a_refRatio,
              const ProblemDomain&     a_domain,
              BCFunc                   a_bc);

  ///
  /** full define function for AMRLevelOp with both coarser and finer levels */
  void define(const DisjointBoxLayout& a_grids,
              const DisjointBoxLayout& a_gridsFiner,
              const DisjointBoxLayout& a_gridsCoarser,
              const Real&              a_dxLevel,
              int                      a_refRatio,
              int                      a_refRatioFiner,
              const ProblemDomain&     a_domain,
              BCFunc                   a_bc);

  /** full define function for AMRLevelOp with finer levels, but no coarser */
  void define(const DisjointBoxLayout& a_grids,
              const DisjointBoxLayout& a_gridsFiner,
              const Real&              a_dxLevel,
              int                      a_refRatio, //  dummy argument, send in 1
              int                      a_refRatioFiner,
              const ProblemDomain&     a_domain,
              BCFunc                   a_bc);

  ///
  virtual void residual(  LevelData<FArrayBox>& a_lhs,
                          const LevelData<FArrayBox>& a_phi,
                          const LevelData<FArrayBox>& a_rhs,
                          bool a_homogeneous = false);

  ///
  virtual void preCond(   LevelData<FArrayBox>& a_correction,
                          const LevelData<FArrayBox>& a_residual);

  ///
  virtual void applyOp(   LevelData<FArrayBox>& a_lhs,
                          const LevelData<FArrayBox>& a_phi,
                          bool a_homogeneous = false);
  ///
  virtual void create(    LevelData<FArrayBox>& a_lhs,
                          const LevelData<FArrayBox>& a_rhs);
  
  ///
  virtual void createCoarsened(    LevelData<FArrayBox>& a_lhs,
                                   const LevelData<FArrayBox>& a_rhs,
                                   const int& a_refRat);

  ///
  virtual void assign(    LevelData<FArrayBox>& a_lhs,
                          const LevelData<FArrayBox>& a_rhs) ;
  ///
  virtual Real dotProduct(const LevelData<FArrayBox>& a_1,
                          const LevelData<FArrayBox>& a_2) ;

  ///
  virtual void incr( LevelData<FArrayBox>& a_lhs,
                     const LevelData<FArrayBox>& a_x,
                     Real a_scale) ;

  ///
  virtual void axby( LevelData<FArrayBox>& a_lhs, 
                     const LevelData<FArrayBox>& a_x,
                     const LevelData<FArrayBox>& a_y,
                     Real a, Real b) ;
  
  ///
  virtual void scale(LevelData<FArrayBox>& a_lhs, const Real& a_scale) ;
  
  ///
  virtual Real norm(const LevelData<FArrayBox>& a_x, int a_ord);

  ///
  virtual void setToZero( LevelData<FArrayBox>& a_x);
  /*@}*/

  /**
     \name MGLevelOp functions */
  /*@{*/

  ///
  virtual void relax(LevelData<FArrayBox>& a_e,
                     const LevelData<FArrayBox>& a_residual,
                     int iterations);

  ///
  virtual void createCoarser(LevelData<FArrayBox>& a_coarse,
                             const LevelData<FArrayBox>& a_fine,
                             bool ghosted);
  /**
     calculate restricted residual
     a_resCoarse[2h] = I[h->2h] (rhsFine[h] - L[h](phiFine[h])
  */
  virtual void restrictResidual(LevelData<FArrayBox>& a_resCoarse,
                                LevelData<FArrayBox>& a_phiFine,
                                const LevelData<FArrayBox>& a_rhsFine);

  /**
     correct the fine solution based on coarse correction
     a_phiThisLevel += I[2h->h](a_correctCoarse)
  */
  virtual void prolongIncrement(LevelData<FArrayBox>& a_phiThisLevel,
                                const LevelData<FArrayBox>& a_correctCoarse);

  /*@}*/

  /**
     \name AMRLevelOp functions */
  /*@{*/

  /** returns 1 when there are no coarser AMRLevelOp objects */
  virtual int refToCoarser(){ return m_refToCoarser; }

  /** a_residual = a_rhs - L(a_phi, a_phiFine, a_phiCoarse) */
  virtual void AMRResidual(LevelData<FArrayBox>& a_residual,
                           const LevelData<FArrayBox>& a_phiFine,
                           const LevelData<FArrayBox>& a_phi,
                           const LevelData<FArrayBox>& a_phiCoarse,
                           const LevelData<FArrayBox>& a_rhs,
                           bool a_homogeneousPhysBC,
                           AMRLevelOp<LevelData<FArrayBox> >* a_finerOp);

  /** residual assuming no more coarser AMR levels */
  virtual void AMRResidualNC(LevelData<FArrayBox>& a_residual,
                             const LevelData<FArrayBox>& a_phiFine,
                             const LevelData<FArrayBox>& a_phi,
                             const LevelData<FArrayBox>& a_rhs,
                             bool a_homogeneousPhysBC,
                             AMRLevelOp<LevelData<FArrayBox> >* a_finerOp);

  /** a_residual = a_rhs - L(a_phi, a_phiCoarse)  */
  virtual void AMRResidualNF(LevelData<FArrayBox>& a_residual,
                             const LevelData<FArrayBox>& a_phi,
                             const LevelData<FArrayBox>& a_phiCoarse,
                             const LevelData<FArrayBox>& a_rhs,
                             bool a_homogeneousPhysBC);


  /** apply AMR operator, including coarse-fine matching */
  virtual void AMROperator(LevelData<FArrayBox>& a_LofPhi,
                           const LevelData<FArrayBox>& a_phiFine,
                           const LevelData<FArrayBox>& a_phi,
                           const LevelData<FArrayBox>& a_phiCoarse,
                           bool a_homogeneousPhysBC,
                           AMRLevelOp<LevelData<FArrayBox> >* a_finerOp);

  /** AMR operator assuming no more coarser AMR levels */
  virtual void AMROperatorNC(LevelData<FArrayBox>& a_LofPhi,
                             const LevelData<FArrayBox>& a_phiFine,
                             const LevelData<FArrayBox>& a_phi,
                             bool a_homogeneousPhysBC,
                             AMRLevelOp<LevelData<FArrayBox> >* a_finerOp);

  /** AMR operator assuming no finer level */
  virtual void AMROperatorNF(LevelData<FArrayBox>& a_LofPhi,
                             const LevelData<FArrayBox>& a_phi,
                             const LevelData<FArrayBox>& a_phiCoarse,
                             bool a_homogeneousPhysBC);

  /** a_resCoarse = I[h-2h]( a_residual - L(a_correction, a_coarseCorrection))
      it is assumed that a_resCoarse has already been filled in with the coarse
      version of AMRResidualNF and that this operation is free to overwrite
      in the overlap regions.
  */

  virtual void AMRRestrict(LevelData<FArrayBox>& a_resCoarse,
                           const LevelData<FArrayBox>& a_residual,
                           const LevelData<FArrayBox>& a_correction,
                           const LevelData<FArrayBox>& a_coarseCorrection);

  /** a_correction += I[h->h](a_coarseCorrection) */
  virtual void AMRProlong(LevelData<FArrayBox>& a_correction,
                          const LevelData<FArrayBox>& a_coarseCorrection);

  /** a_residual = a_residual - L(a_correction, a_coarseCorrection) */
  virtual void AMRUpdateResidual(LevelData<FArrayBox>& a_residual,
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

  /*@}*/

  ///public constants
  Real m_mu;

  /// coefficient storage
  RefCountedPtr<LevelData<FArrayBox> > m_alpha;
  /// face-centered coefficient storage
  RefCountedPtr<LevelData<FluxBox> > m_beta;
  Real                    m_dxCrse;

protected:
  Real                    m_dx;
  ProblemDomain           m_domain;
  LevelDataOps<FArrayBox> m_levelOps;
  BCFunc                  m_bc;
  LayoutData<CFIVS>       m_loCFIVS[SpaceDim];
  LayoutData<CFIVS>       m_hiCFIVS[SpaceDim];
  Copier                  m_exchangeCopier;
  QuadCFInterp            m_interpWithCoarser;
  int                     m_refToCoarser;
  int                     m_refToFiner;

  void levelGSRB(LevelData<FArrayBox>& a_e,  
                 const LevelData<FArrayBox>& a_residual);

  void levelJacobi(LevelData<FArrayBox>& a_phi,
                   const LevelData<FArrayBox>& a_rhs);

  void homogeneousCFInterp(LevelData<FArrayBox>& a_phif);
  void homogeneousCFInterp(LevelData<FArrayBox>& a_phif,
                           const DataIndex& a_datInd,
                           int a_idir,
                           Side::LoHiSide a_hiorlo);
  void interpOnIVSHomo(LevelData<FArrayBox>& a_phif,
                       const DataIndex& a_datInd,
                       const int a_idir,
                       const Side::LoHiSide a_hiorlo,
                       const IntVectSet& a_interpIVS);
  
  /// computes flux over face-centered a_facebox.
  void getFlux(FArrayBox&       a_flux,
               const FArrayBox& a_data,
               const FluxBox&   a_beta,
               const Box&       a_facebox,
               int              a_dir,
               int ref = 1) const ;
  
  void reflux(const LevelData<FArrayBox>& a_phiFine,
              const LevelData<FArrayBox>& a_phi,
              LevelData<FArrayBox>& residual,
              AMRLevelOp<LevelData<FArrayBox> >* a_finerOp);

  /// utility function which computes operator after all bc's have been set
  void computeOperatorNoBCs(LevelData<FArrayBox>& a_lhs,
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

  virtual ~VCAMRPoissonOpFactory(){;}

  /// 
  /**
     a_coarseDomain is the domain at the coarsest level.
     a_grids is the AMR  hierarchy.
     a_refRatios are the refinement ratios between levels.  The ratio lives with the coarser level
     so a_refRatios[ilev] is the ratio between ilev and ilev+1
     a_coarseDx is the grid spacing at the coarsest level.
     BCFunc holds the boundary conditions.
     a_alpha is the identity coefficient
     a_beta is the laplacian coefficient.
  */
  void define(const ProblemDomain& a_coarseDomain,
              const Vector<DisjointBoxLayout>& a_grids,
              const Vector<int>& a_refRatios,
              const Real&        a_coarsedx,
              BCFunc a_bc);

  // regular multigrid definition function --deprecated
  void define(const ProblemDomain& a_domain,
              const DisjointBoxLayout& a_grid,
              const Real&    a_dx,
              BCFunc a_bc,
              int maxDepth = -1);


  /// set diagonal coefficient
  void setAlpha(Vector<RefCountedPtr<LevelData<FArrayBox> > >& a_alpha);
  
  /// only sets one level
  void setAlpha(RefCountedPtr<LevelData<FArrayBox> >& a_alpha, 
                int a_level);

  ///
  void setBeta(Vector<RefCountedPtr<LevelData<FluxBox> > >& a_beta);

  ///
  void setBeta(RefCountedPtr<LevelData<FluxBox> >& a_beta,
               int a_level);
  

  /// access functions for alpha and beta
  Vector<RefCountedPtr<LevelData<FArrayBox> > >& getAlpha() {return m_alpha;}

  ///
  const Vector<RefCountedPtr<LevelData<FArrayBox> > >& getAlpha() const
  {return m_alpha;}

  /// 
  RefCountedPtr<LevelData<FArrayBox> >& getAlpha(int a_level)
  {return m_alpha[a_level];}

  ///
  const RefCountedPtr<LevelData<FArrayBox> >& getAlpha(int a_level) const
  {return m_alpha[a_level];}

  Vector<RefCountedPtr<LevelData<FluxBox> > >& getBeta()
  {return m_beta;}

  ///
  const Vector<RefCountedPtr<LevelData<FluxBox> > >& getBeta() const
  {return m_beta;}

  /// 
  RefCountedPtr<LevelData<FluxBox> >& getBeta(int a_level)
  {return m_beta[a_level];}

  ///
  const RefCountedPtr<LevelData<FluxBox> >& getBeta(int a_level) const
  {return m_beta[a_level];}


  ///
  virtual VCAMRPoissonOp*
  MGnewOp(const ProblemDomain& a_FineindexSpace,
          int depth,
          bool homoOnly = true);

  ///
  virtual  VCAMRPoissonOp* AMRnewOp(const ProblemDomain& a_indexSpace);

  ///
  virtual int refToFiner(const ProblemDomain&) const;


  int m_coefficient_average_type;


  void setDefaultValues();

private:
  Vector<ProblemDomain>     m_domains;
  Vector<DisjointBoxLayout> m_boxes;
  Vector<Real>     m_dx;
  Vector<int>       m_refRatios; // refinement to next coarser level
  BCFunc  m_bc;
  Vector<RefCountedPtr<LevelData<FArrayBox> > > m_alpha;
  Vector<RefCountedPtr<LevelData<FluxBox> > > m_beta;

};

#include "NamespaceFooter.H"
#endif

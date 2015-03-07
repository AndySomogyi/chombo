#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

#ifndef _EBVISCOUSTENSOROP_H_
#define _EBVISCOUSTENSOROP_H_

#include "REAL.H"
#include "Box.H"
#include "FArrayBox.H"
#include "Vector.H"
#include <map>
#include "RefCountedPtr.H"

#include "AMRMultiGrid.H"

#include "EBIndexSpace.H"
#include "EBCellFAB.H"
#include "EBCellFactory.H"

#include "EBLevelDataOps.H"
#include "BaseEBBC.H"
#include "AMRTGA.H"
#include "BaseDomainBC.H"
#include "CFIVS.H"
#include "EBFluxRegister.H"
#include "EBMGAverage.H"
#include "EBMGInterp.H"
#include "EBTensorCFInterp.H"
#include "EBCoarsen.H"
#include "PolyGeom.H"
#include "EBQuadCFInterp.H"
#include "EBLevelGrid.H"
#include "NamespaceHeader.H"

///
/**
   operator is alpha a I + (divF) = alpha I + beta*div(eta(grad B + grad B^T) + lambda* I div B )
   beta, lambda,  and eta are incorporated into the flux F
*/
//class EBViscousTensorOp: public AMRLevelOp<LevelData<EBCellFAB> >
class EBViscousTensorOp: public LevelTGAHelmOp<LevelData<EBCellFAB>, EBFluxFAB >
{
public:

  ///for tga to reset stuff
  virtual void setAlphaAndBeta(const Real& a_alpha,
                               const Real& a_beta)
  {
    m_alpha = a_alpha;
    m_beta = a_beta;
    defineStencils();
  }

  ///another tgaism
  virtual void diagonalScale(LevelData<EBCellFAB> & a_rhs)
  {
    EBLevelDataOps::kappaWeight(a_rhs);
    for(DataIterator dit = m_eblg.getDBL().dataIterator(); dit.ok(); ++dit)
      {
        for(int idir = 0; idir < SpaceDim; idir++)
          {
            int isrc = 0; int idst = idir; int inco = 1;
            a_rhs[dit()].mult((*m_acoef)[dit()], isrc, idst, inco);
          }
      }
  }

  ///a leveltgaism
  virtual void fillGrad(const LevelData<EBCellFAB>& a_phi);

  ///another leveltgaism
  virtual void getFlux(EBFluxFAB&                    a_flux,
                       const LevelData<EBCellFAB>&   a_data,
                       const Box&                    a_grid,
                       const DataIndex&              a_dit,
                       Real                          a_scale)
  {
    a_flux.define(m_eblg.getEBISL()[a_dit], a_grid, SpaceDim);
    a_flux.setVal(0.);
    for(int idir = 0; idir < SpaceDim; idir++)
      {
        FArrayBox&       regFlux =       (FArrayBox&) a_flux[idir ].getSingleValuedFAB();
        const FArrayBox& regData = (const FArrayBox&) a_data[a_dit].getSingleValuedFAB();
        const FArrayBox& regGrad = (const FArrayBox&) m_grad[a_dit].getSingleValuedFAB();
        Box faceBox = surroundingNodes(a_grid, idir);
        getFlux(regFlux, regData, regGrad, faceBox, idir, a_dit);
      }
  }

  /// apply operator without any boundary or coarse-fine boundary conditions and no finer level
  virtual void applyOpNoBoundary(LevelData<EBCellFAB>&       a_opPhi,
                                 const LevelData<EBCellFAB>& a_phi)
  {
    s_turnOffBCs = true;
    applyOp(a_opPhi, a_phi, true);
    s_turnOffBCs = false;
  }

//functions used by the wider world

  ///
  /**
   */
  EBViscousTensorOp(const EBLevelGrid &                                a_eblgFine,
                    const EBLevelGrid &                                a_eblg,
                    const EBLevelGrid &                                a_eblgCoar,
                    const EBLevelGrid &                                a_eblgCoarMG,
                    const Real&                                        a_alpha,
                    const Real&                                        a_beta,
                    const RefCountedPtr<LevelData<EBCellFAB> >&        a_acoef,
                    const RefCountedPtr<LevelData<EBFluxFAB> >&        a_eta,
                    const RefCountedPtr<LevelData<EBFluxFAB> >&        a_lambda,
                    const RefCountedPtr<LevelData<BaseIVFAB<Real> > >& a_etaIrreg,
                    const RefCountedPtr<LevelData<BaseIVFAB<Real> > >& a_lambdaIrreg,
                    const Real&                                        a_dx,
                    const Real&                                        a_dxCoar,
                    const int&                                         a_refToFine,
                    const int&                                         a_refToCoar,
                    const RefCountedPtr<ViscousBaseDomainBC>&          a_domainBC,
                    const RefCountedPtr<ViscousBaseEBBC>&              a_ebBC,
                    const bool&                                        a_hasMGObjects,
                    const IntVect&                                     a_ghostCellsPhi,
                    const IntVect&                                     a_ghostCellsRHS);

  ///
  virtual ~EBViscousTensorOp();

  ///
  /** a_residual = a_rhs - L(a_phiFine, a_phi)   no coaser AMR level*/
  void AMRResidualNC(LevelData<EBCellFAB>&       a_residual,
                     const LevelData<EBCellFAB>& a_phiFine,
                     const LevelData<EBCellFAB>& a_phi,
                     const LevelData<EBCellFAB>& a_rhs,
                     bool a_homogeneousBC,
                     AMRLevelOp<LevelData<EBCellFAB> >* a_finerOp);


  ///
  /** apply AMR operator   no coaser AMR level*/
  void AMROperatorNC(LevelData<EBCellFAB>&       a_LofPhi,
                     const LevelData<EBCellFAB>& a_phiFine,
                     const LevelData<EBCellFAB>& a_phi,
                     bool a_homogeneousBC,
                     AMRLevelOp<LevelData<EBCellFAB> >* a_finerOp);

  ///
  /**
   */
  virtual void residual(LevelData<EBCellFAB>&       a_residual,
                        const LevelData<EBCellFAB>& a_phi,
                        const LevelData<EBCellFAB>& a_rhs,
                        bool                        a_homogeneousPhysBC=false);

  ///
  /**
   */
  virtual void preCond(LevelData<EBCellFAB>&       a_opPhi,
                       const LevelData<EBCellFAB>& a_phi);


  ///
  /**
     this is the linearop function.  CFBC is set to homogeneous.  phic is null
  */
  virtual void applyOp(LevelData<EBCellFAB>&             a_opPhi,
                       const LevelData<EBCellFAB>&       a_phi,
                       bool                              a_homogeneousPhysBC);

  ///
  /**
   */
  virtual void create(LevelData<EBCellFAB>&       a_lhs,
                      const LevelData<EBCellFAB>& a_rhs);

  ///
  virtual void createCoarsened(LevelData<EBCellFAB>&       a_lhs,
                               const LevelData<EBCellFAB>& a_rhs,
                               const int&                  a_refRat);

  ///
  Real
  AMRNorm(const LevelData<EBCellFAB>& a_coarResid,
          const LevelData<EBCellFAB>& a_fineResid,
          const int& a_refRat,
          const int& a_ord);


  ///
  /**
   */
  virtual void assign(LevelData<EBCellFAB>&       a_lhs,
                      const LevelData<EBCellFAB>& a_rhs);

  ///
  /**
   */
  virtual Real dotProduct(const LevelData<EBCellFAB>& a_1,
                          const LevelData<EBCellFAB>& a_2);

  ///
  /**
   */
  virtual void incr(LevelData<EBCellFAB>&       a_lhs,
                    const LevelData<EBCellFAB>& a_x,
                    Real                        a_scale);

  ///
  /**
   */
  virtual void axby(LevelData<EBCellFAB>&       a_lhs,
                    const LevelData<EBCellFAB>& a_x,
                    const LevelData<EBCellFAB>& a_y,
                    Real                        a_a,
                    Real                        a_b);

  ///
  /**
   */
  virtual void scale(LevelData<EBCellFAB>& a_lhs,
                     const Real&           a_scale);

  ///
  /**
   */
  virtual Real norm(const LevelData<EBCellFAB>& a_rhs,
                    int                         a_ord);

  ///
  /**
   */
  virtual void setToZero(LevelData<EBCellFAB>& a_lhs);

  ///
  /**
   */
  virtual void setVal(LevelData<EBCellFAB>& a_lhs, const Real& a_value);

  ///
  /**
   */
  virtual void createCoarser(LevelData<EBCellFAB>&       a_coarse,
                             const LevelData<EBCellFAB>& a_fine,
                             bool                        a_ghosted);

  ///
  /**
   */
  virtual void relax(LevelData<EBCellFAB>&       a_e,
                     const LevelData<EBCellFAB>& a_residual,
                     int                         a_iterations);

  ///
  /**
     Calculate restricted residual:
     a_resCoarse[2h] = I[h->2h] (a_rhsFine[h] - L[h](a_phiFine[h]))
  */
  virtual void restrictResidual(LevelData<EBCellFAB>&       a_resCoarse,
                                LevelData<EBCellFAB>&       a_phiFine,
                                const LevelData<EBCellFAB>& a_rhsFine);

  ///
  /**
     Correct the fine solution based on coarse correction:
     a_phiThisLevel += I[2h->h] (a_correctCoarse)
  */
  virtual void prolongIncrement(LevelData<EBCellFAB>&       a_phiThisLevel,
                                const LevelData<EBCellFAB>& a_correctCoarse);

  ///
  /** Refinement ratio between this level and coarser level.
      Returns 1 when there are no coarser AMRLevelOp objects */
  virtual int refToCoarser();

  ///
  /** Refinement ratio between this level and coarser level.
      Returns 1 when there are no coarser AMRLevelOp objects */
  virtual int refToFiner();

  ///
  /** a_residual = a_rhs - L(a_phi, a_phiFine, a_phiCoarse) */
  virtual void AMRResidual(LevelData<EBCellFAB>& a_residual,
                           const LevelData<EBCellFAB>& a_phiFine,
                           const LevelData<EBCellFAB>& a_phi,
                           const LevelData<EBCellFAB>& a_phiCoarse,
                           const LevelData<EBCellFAB>& a_rhs,
                           bool a_homogeneousBC,
                           AMRLevelOp<LevelData<EBCellFAB> >* a_finerOp);

  ///
  /** a_residual = a_rhs - L(a_phi, a_phiCoarse)  */
  virtual void AMRResidualNF(LevelData<EBCellFAB>& a_residual,
                             const LevelData<EBCellFAB>& a_phi,
                             const LevelData<EBCellFAB>& a_phiCoarse,
                             const LevelData<EBCellFAB>& a_rhs,
                             bool a_homogeneousBC);


  ///
  /** a_residual = a_rhs - L(a_phi, a_phiFine, a_phiCoarse) */
  virtual void AMROperator(LevelData<EBCellFAB>& a_LofPhi,
                           const LevelData<EBCellFAB>& a_phiFine,
                           const LevelData<EBCellFAB>& a_phi,
                           const LevelData<EBCellFAB>& a_phiCoarse,
                           bool a_homogeneousBC,
                           AMRLevelOp<LevelData<EBCellFAB> >* a_finerOp);

  ///
  /** a_residual = a_rhs - L(a_phi, a_phiCoarse)  */
  virtual void AMROperatorNF(LevelData<EBCellFAB>& a_LofPhi,
                             const LevelData<EBCellFAB>& a_phi,
                             const LevelData<EBCellFAB>& a_phiCoarse,
                             bool a_homogeneousBC);

  ///
  /** a_resCoarse = I[h-2h] (a_residual - L(a_correction, a_coarseCorrection)) */
  virtual void AMRRestrict(LevelData<EBCellFAB>& a_resCoarse,
                           const LevelData<EBCellFAB>& a_residual,
                           const LevelData<EBCellFAB>& a_correction,
                           const LevelData<EBCellFAB>& a_coarseCorrection);

  ///
  /** a_correction += I[2h->h](a_coarseCorrection) */
  virtual void AMRProlong(LevelData<EBCellFAB>&       a_correction,
                          const LevelData<EBCellFAB>& a_coarseCorrection);

  ///
  /** a_residual = a_residual - L(a_correction, a_coarseCorrection) */
  virtual void AMRUpdateResidual(LevelData<EBCellFAB>&       a_residual,
                                 const LevelData<EBCellFAB>& a_correction,
                                 const LevelData<EBCellFAB>& a_coarseCorrection);

  ///
  /** a_residual = a_residual - L(a_correction, a_coarseCorrection) */
  static void
  getDivergenceStencil(VoFStencil&      a_divStencil,
                       const FaceIndex& a_face,
                       const DataIndex& a_dit,
                       const Real     & a_dx,
                       const EBLevelGrid& a_eblg);

  ///
  /** a_residual = a_residual - L(a_correction, a_coarseCorrection) */
  static void
  getGradientStencil(VoFStencil&  a_gradStencil,
                     int a_ivar,
                     int a_diffDir,
                     const FaceIndex& a_face,
                     const DataIndex& a_dit,
                     const Real     & a_dx,
                     const EBLevelGrid& a_eblg);
protected:

  //internal functions not usable without significant context
  //  void fillGrad(const LevelData<EBCellFAB>&       a_phi);

  void incrOpRegularDir(EBCellFAB&             a_lhs,
                        const EBCellFAB&       a_phi,
                        const bool&            a_homogeneous,
                        const int&             a_dir,
                        const DataIndex&       a_datInd);

  void applyOpIrregular(EBCellFAB&             a_lhs,
                        const EBCellFAB&       a_phi,
                        const bool&            a_homogeneous,
                        const DataIndex&       a_datInd);

  void getFlux(FArrayBox&                    a_flux,
               const FArrayBox&              a_phi,
               const FArrayBox&              a_gradPhi,
               const Box&                    a_faceBox,
               const int&                    a_idir,
               const DataIndex&              a_datInd);
  void
  homogeneousCFInterp(EBCellFAB&             a_phi,
                      const DataIndex&       a_datInd,
                      const int&             a_idir,
                      const Side::LoHiSide&  a_hiorlo);
  void
  cfinterp(const LevelData<EBCellFAB>&       a_phi,
           const LevelData<EBCellFAB>&       a_phiCoarse);

  void reflux(const LevelData<EBCellFAB>& a_phiFine,
              const LevelData<EBCellFAB>& a_phi,
              LevelData<EBCellFAB>& residual,
              AMRLevelOp<LevelData<EBCellFAB> >* a_finerOp);

  void homogeneousCFInterp(LevelData<EBCellFAB>&   a_phi);

  void  defineStencils();

  void getVoFStencil(VoFStencil&      a_vofStencil,
                     const VolIndex&  a_vof,
                     const DataIndex& a_dit,
                     int              a_ivar);

  void getFluxStencil(VoFStencil&      a_fluxStencil,
                      const FaceIndex& a_face,
                      const DataIndex& a_dit,
                      int a_ivar);

  void getDivergenceStencil(VoFStencil&  a_gradStencil,
                            const FaceIndex& a_face,
                            const DataIndex& a_dit);

  void getGradientStencil(VoFStencil&  a_gradStencil,
                          int a_ivar,
                          int a_diffDir,
                          const FaceIndex& a_face,
                          const DataIndex& a_dit);


  void gsrbColor(LevelData<EBCellFAB>&       a_phi,
                 const LevelData<EBCellFAB>& a_lph,
                 const LevelData<EBCellFAB>& a_rhs,
                 const IntVect&              a_color);

  void getFaceCenteredFluxStencil(VoFStencil&      a_fluxStencil,
                                  const FaceIndex& a_face,
                                  const DataIndex& a_dit,
                                  int a_ivar);

  void cellGrad(EBCellFAB&             a_gradPhi,
                const  EBCellFAB&      a_phi,
                const Box&             a_grid,
                const DataIndex&       a_datInd);

  void
  incrementFRCoar(const LevelData<EBCellFAB>& a_phiFine,
                  const LevelData<EBCellFAB>& a_phi);


  void
  incrementFRFine(const LevelData<EBCellFAB>& a_phiFine,
                  const LevelData<EBCellFAB>& a_phi,
                  AMRLevelOp<LevelData<EBCellFAB> >* a_finerOp);

  //member data, bloated though it may be
  static bool s_turnOffBCs;
  //locations of coarse fine interfaces
  LayoutData<CFIVS>                            m_loCFIVS[CH_SPACEDIM];
  LayoutData<CFIVS>                            m_hiCFIVS[CH_SPACEDIM];

  //grid information
  EBLevelGrid                                  m_eblgFine;
  EBLevelGrid                                  m_eblg;
  EBLevelGrid                                  m_eblgCoar;
  EBLevelGrid                                  m_eblgCoarMG;

  //operator coefficients
  Real                                         m_alpha;
  Real                                         m_beta;
  RefCountedPtr<LevelData<EBCellFAB> >         m_acoef;
  RefCountedPtr<LevelData<EBFluxFAB> >         m_eta;
  RefCountedPtr<LevelData<EBFluxFAB> >         m_lambda;
  RefCountedPtr<LevelData<BaseIVFAB<Real> > >  m_etaIrreg;
  RefCountedPtr<LevelData<BaseIVFAB<Real> > >  m_lambdaIrreg;


  //flux register with finer level
  EBFluxRegister                               m_fluxReg;

  //grid spacing
  //dxcoar not implied by ref ratio--we could be in some MG level where dx > dxcoar
  //needed for homogeneous cf interp
  Real                                         m_dx;
  Real                                         m_dxCoar;

  //refinement ratios and whether this object has multigrid objects
  bool                                         m_hasFine;
  bool                                         m_hasCoar;
  int                                          m_refToFine;
  int                                          m_refToCoar;
  bool                                         m_hasMGObjects;

  //ghost cell information needed for ebstencil optimizations
  IntVect                                      m_ghostCellsPhi;
  IntVect                                      m_ghostCellsRHS;

  //stencil to apply irregular operator
  LayoutData<RefCountedPtr<EBStencil> >        m_stencil[CH_SPACEDIM];

  //relaxation coefficent
  LevelData<EBCellFAB>                         m_relCoef;

  //gradient of solution at cell centers
  LevelData<EBCellFAB>                         m_grad;

  //for irregular cell iteration (includes buffer around multivalued cells)
  LayoutData<VoFIterator >                     m_vofIter;
  LayoutData<VoFIterator >                     m_vofMult;
  //for domain boundary conditions at ir regular cells
  LayoutData<VoFIterator >                     m_vofIterDomLo[CH_SPACEDIM];
  LayoutData<VoFIterator >                     m_vofIterDomHi[CH_SPACEDIM];

  //for inhomogeneous cf interpolation
  RefCountedPtr<EBTensorCFInterp>              m_interpWithCoarser;

  //averaging operators
  EBMGAverage                                  m_ebAverage;
  EBMGAverage                                  m_ebAverageMG;

  //interpolation operators
  EBMGInterp                                   m_ebInterp;
  EBMGInterp                                   m_ebInterpMG;

  //boundary condition operators
  RefCountedPtr<ViscousBaseDomainBC>           m_domainBC;
  RefCountedPtr<ViscousBaseEBBC>               m_ebBC;
  Vector<IntVect>                              m_colors;
private:
  ///weak construction is bad
  EBViscousTensorOp()
  {
    MayDay::Error("invalid operator");
  }

  //copy constructor and operator= disallowed for all the usual reasons
  EBViscousTensorOp(const EBViscousTensorOp& a_opin)
  {
    MayDay::Error("invalid operator");
  }

  void operator=(const EBViscousTensorOp& a_opin)
  {
    MayDay::Error("invalid operator");
  }
};


#include "NamespaceFooter.H"
#endif

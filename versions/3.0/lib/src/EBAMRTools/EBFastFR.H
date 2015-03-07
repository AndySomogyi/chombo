#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

#ifndef _EBFASTFR_H_
#define _EBFASTFR_H_

#include "REAL.H"
#include "Vector.H"
#include "EBCellFAB.H"
#include "EBFaceFAB.H"
#include "EBISLayout.H"
#include "EBISBox.H"
#include "IntVectSet.H"
#include "CFStencil.H"
#include "LoHiSide.H"
#include "LevelData.H"
#include "EBLevelGrid.H"
#include "LayoutData.H"
#include "LevelFluxRegister.H"
#include "NamespaceHeader.H"
class EBIndexSpace;

///EBFastFR-A class to encapsulate a levels worth of flux registers.
/**
   A EBFastFR  handles all the data choreography
   necessary to create a levels worth of flux registers.
   This only does the elliptic flux register stuff for now (no interacting with
   redistribution, no RZ, no regular and regular separate).
*/
class EBFastFR
{
public:

  ///
  /**
     Default constructor.  Leaves object undefined.
  */
  EBFastFR();

  ///
  /**
     Full constructor.  Calls the define function which creates
     a levels worth of flux registers.
  */
  EBFastFR(const EBLevelGrid&       a_eblgFine,
           const EBLevelGrid&       a_eblgCoar,
           const int&               a_nref,
           const int&               a_nvar);

  ///
  ~EBFastFR();

  ///
  /**
     Full define function.  Creates a levels worth of flux registers.
  */
  void
  define(const EBLevelGrid&       a_eblgFine,
         const EBLevelGrid&       a_eblgCoar,
         const int&               a_nref,
         const int&               a_nvar);

  ///
  /**
     Initialize values of registers  to zero.
  */
  void
  setToZero();


  ///
  void
  incrementCoarseBoth(const EBFaceFAB&      a_coarseFlux,
                      const Real&           a_scale,
                      const DataIndex&      a_coarseDataIndex,
                      const Interval&       a_variables,
                      const int&            a_dir,
                      const Side::LoHiSide& a_sd);


  ///
  void
  incrementFineBoth(const EBFaceFAB&      a_fineFlux,
                    const Real&           a_scale,
                    const DataIndex&      a_fineDataIndex,
                    const Interval&       a_variables,
                    const int&            a_dir,
                    const Side::LoHiSide& a_sd);



  ///
  void reflux(LevelData<EBCellFAB>& a_uCoarse,
              const Interval&       a_variables,
              const Real&           a_scale);


  ///
  /**
     Increments mass array  with left-over mass
     from reflux divergence.   this is to test this pig.
     Ordinarily this mass would go into redistribution.   this
     is used in the test suite to  check whether a constant
     flux refluxed and then unrefluxed (this function) ends up
     with a constant solution.
     Correction at each cell = (1-kappa)refluxCorrection.
  */
  void incrementDensityArray(LevelData<EBCellFAB>& a_coarMass,
                             const Interval&       a_variables,
                             const Real&           a_scale);

  ///
  bool isDefined() const;

protected:

  static int index(int a_dir, Side::LoHiSide a_side);

  void incrementFineIrreg(const EBFaceFAB&      a_fineFlux,
                          const Real&           a_newScale,
                          const DataIndex&      a_fineDatInd,
                          const Interval&       a_variables,
                          const int&            a_dir,
                          const Side::LoHiSide& a_sd);

  void incrementCoarIrreg(const EBFaceFAB&       a_coarFlux,
                          const Real&            a_scale,
                          const DataIndex&       a_coarDatInd,
                          const Interval&        a_variables,
                          const int&             a_dir,
                          const Side::LoHiSide&  a_sd);

  void  irregSetToZero();

  void cacheOldSolution(const LevelData<EBCellFAB>& a_uCoar,
                        const Interval&             a_variables);

  void restoreOldSolution(LevelData<EBCellFAB>&       a_uCoar,
                          const Interval&             a_variables);

  void irregReflux(LevelData<EBCellFAB>& a_uCoar,
                   const Interval&       a_variables,
                   const Real&           a_scale);

  LevelFluxRegister m_levelFluxReg;
  bool m_isDefined;

  EBLevelGrid       m_eblgFine;
  EBLevelGrid       m_eblgCoar;
  int               m_refRat;
  int               m_nComp;

  //ebisl of coarse buffers
  EBISLayout m_ebislBufCoar[2*SpaceDim];
  //ebisl of fine buffers
  EBISLayout m_ebislBufFine[2*SpaceDim];

  //map from input coarse grid to coarse buffers
  LayoutData<Vector<DataIndex>  > m_coarIndexMap[SpaceDim*2];

  LayoutData<bool>  m_noEBCFine[SpaceDim*2];
  LayoutData<bool>  m_noEBCCoar[SpaceDim*2];

  //grids for fine buffers
  DisjointBoxLayout m_bufGridsFine[SpaceDim*2];
  //grids for coarse buffers
  DisjointBoxLayout m_bufGridsCoar[SpaceDim*2];

  //fine buffers
  LevelData<BaseIVFAB<Real> > m_regsCoar[SpaceDim*2];
  //coarse buffers
  LevelData<BaseIVFAB<Real> > m_regsFine[SpaceDim*2];

  //sets for buffers
  LayoutData<IntVectSet>  m_cfivsFine[SpaceDim*2];
  LayoutData<IntVectSet>  m_cfivsCoar[SpaceDim*2];

  LayoutData<VoFIterator> m_vofitBufCoar[SpaceDim*2];
  LayoutData<VoFIterator> m_vofitBufFine[SpaceDim*2];
  LayoutData<VoFIterator> m_vofitCoar;
  //sets for coarse buffers
  //scratch space for refluxing
  LevelData<BaseIVFAB<Real> > m_scratchc[SpaceDim*2];
  LayoutData<BaseIVFAB<Real> > m_saveoldc;
  //nrefdmo is the number of fine faces per coar face
  //this is intrinsically dimension-dependent
  Real m_nrefdmo;
  bool m_hasEBCF;
private:
  //for all the usual reasons,
  //there is no copy constructor for this class
  //and there is no operator= for this class
  //(there is no rule six).
  void operator= (const EBFastFR& out)
  {
    MayDay::Error("invalid operator");
  }
  EBFastFR(const EBFastFR& out)
  {
    MayDay::Error("invalid operator");
  }
};
#include "NamespaceFooter.H"
#endif

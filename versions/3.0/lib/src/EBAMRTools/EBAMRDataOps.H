#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

#ifndef _EBAMRDATAOPS_H_
#define _EBAMRDATAOPS_H_

#include "EBLevelDataOps.H"
#include "FaceIterator.H"
#include "EBFluxFAB.H"
#include "EBCoarseAverage.H"
#include "EBCoarsen.H"
#include "EBQuadCFInterp.H"
#include "EBPWLFillPatch.H"
#include "EBPWQuadFillPatch.H"
#include "EBIndexSpace.H"
#include "EBAMRIO.H"
#include "EBLevelGrid.H"
#include "NamespaceHeader.H"

#define EBAMRDATAOPS_INTERIORREGVOFS   0x1
#define EBAMRDATAOPS_BOUNDARYREGVOFS   0x2
#define EBAMRDATAOPS_INTERIORIRREGVOFS 0x4
#define EBAMRDATAOPS_BOUNDARYIRREGVOFS 0x8

#define EBAMRDATAOPS_ALLVOFS         (EBAMRDATAOPS_INTERIORREGVOFS   | \
                                        EBAMRDATAOPS_BOUNDARYREGVOFS   | \
                                        EBAMRDATAOPS_INTERIORIRREGVOFS | \
                                        EBAMRDATAOPS_BOUNDARYIRREGVOFS )

#define EBAMRDATAOPS_REGULARVOFS     (EBAMRDATAOPS_INTERIORREGVOFS | \
                                        EBAMRDATAOPS_BOUNDARYREGVOFS )

#define EBAMRDATAOPS_IRREGULARVOFS   (EBAMRDATAOPS_INTERIORIRREGVOFS | \
                                        EBAMRDATAOPS_BOUNDARYIRREGVOFS )

#define EBAMRDATAOPS_BOUNDARYVOFS    (EBAMRDATAOPS_BOUNDARYREGVOFS | \
                                        EBAMRDATAOPS_IRREGULARVOFS   )
///
/**
   Repository of functions for EB Hierarchies.
 */
class EBAMRDataOps
{
public:
  EBAMRDataOps()
  {
  }

  static void
  getErrorFromCoarseAndFine(Vector< LevelData<EBCellFAB>* >&           a_errorCoar,
                            const Vector< LevelData<EBCellFAB>* >&     a_solnCoar,
                            const Vector< DisjointBoxLayout >&         a_gridsCoar,
                            const Vector< EBISLayout >&                a_ebislCoar,
                            const ProblemDomain&                       a_level0DomainCoar,
                            const Vector< LevelData<EBCellFAB>* >&     a_solnFine,
                            const Vector< DisjointBoxLayout >&         a_gridsFine,
                            const Vector< EBISLayout >&                a_ebislFine,
                            const ProblemDomain&                       a_level0DomainFine,
                            const Vector<int>&                         a_refRat);

  ~EBAMRDataOps(){;}

  ///
  /**
     count vofs, boxes, regular cells, etc.
   */
  static int countVoF(const Vector<DisjointBoxLayout>& a_dbl,
                      const Vector<EBISLayout>&        a_ebisl,
                      const Vector<ProblemDomain>&     a_domain);

  ///
  /**
     average a cell centered velocity vector (with SpaceDim components) to face centers
   */
  static void averageCellToFacesMAC(Vector<LevelData<EBFluxFAB>* >&        a_dataFlux,
                                    const Vector<LevelData<EBCellFAB>* >&  a_dataCell,
                                    const Vector<EBLevelGrid>&             a_eblg,
                                    const Vector<int>&                     a_refRatio);


  ///
  /**
     average cell centered data to all face centers for component a_comp
   */
  static void averageCellToFaces(Vector<LevelData<EBFluxFAB>* >&        a_dataFlux,
                                 const Vector<LevelData<EBCellFAB>* >&  a_dataCell,
                                 const Vector<EBLevelGrid>&             a_eblg,
                                 const Vector<int>&                     a_refRatio,
                                 const int&                             a_comp);

  ///
  /**
     average cell centered data to all face centers for component a_comp
   */
  static void averageCellToFaces(Vector<RefCountedPtr<LevelData<EBFluxFAB> > >& a_dataFlux,
                                 const Vector<LevelData<EBCellFAB>* >&          a_dataCell,
                                 const Vector<EBLevelGrid>&                     a_eblg,
                                 const Vector<int>&                             a_refRatio,
                                 const int&                                     a_comp);


  ///
  /**
   */
  static void quadCFInterpAll(Vector<LevelData<EBCellFAB>* >&  a_data,
                              const Vector<DisjointBoxLayout>& a_dbl,
                              const Vector<EBISLayout>&        a_ebisl,
                              const Vector<ProblemDomain>&     a_domain,
                              const Vector<int>&               a_refRatio);

  ///
  /**
   */
  static void quadCFInterpAll(Vector<LevelData<EBCellFAB>* >&  a_data,
                              const Vector<EBLevelGrid>&       a_eblg,
                              const Vector<int>&               a_refRatio);


  ///
  /**
   */
  static void quadCFInterpOne(Vector<LevelData<EBCellFAB>* >&  a_data,
                              const Vector<EBLevelGrid>&       a_eblg,
                              const Vector<int>&               a_refRatio,
                              const int&                       a_fineLevel);

  ///
  /**
   */
  static void quadCFInterpOne(Vector<LevelData<EBCellFAB>* >&  a_data,
                              const Vector<DisjointBoxLayout>& a_dbl,
                              const Vector<EBISLayout>&        a_ebisl,
                              const Vector<ProblemDomain>&     a_domain,
                              const Vector<int>&               a_refRatio,
                              const int&                       a_fineLevel);


  ///
  /**
   */
  static void quadCFInterpOne(LevelData<EBCellFAB>&       a_dataFine,
                              const LevelData<EBCellFAB>& a_dataCoar,
                              const DisjointBoxLayout&    a_dblFine,
                              const DisjointBoxLayout&    a_dblCoar,
                              const EBISLayout&           a_ebislFine,
                              const EBISLayout&           a_ebislCoar,
                              const ProblemDomain&        a_domainCoar,
                              const int&                  a_refRatio);


  ///
  /**
   */
  static void pwlFillPatchAll(Vector<LevelData<EBCellFAB>* >&  a_data,
                              const Vector<EBLevelGrid>&       a_eblg,
                              const Vector<int>&               a_refRatio);
  ///
  /**
   */
  static void pwlFillPatchAll(Vector<LevelData<EBCellFAB>* >&  a_data,
                              const Vector<DisjointBoxLayout>& a_dbl,
                              const Vector<EBISLayout>&        a_ebisl,
                              const Vector<ProblemDomain>&     a_domain,
                              const Vector<int>&               a_refRatio);


  ///
  /**
   */
  static void pwlFillPatchOne(LevelData<EBCellFAB>&       a_dataFine,
                              LevelData<EBCellFAB>&       a_dataCoar,
                              const DisjointBoxLayout&    a_dblFine,
                              const DisjointBoxLayout&    a_dblCoar,
                              const EBISLayout&           a_ebislFine,
                              const EBISLayout&           a_ebislCoar,
                              const ProblemDomain&        a_domainCoar,
                              const int&                  a_refRatioCoar);


  ///
  /**
   */
  static void pwqFillPatchAll(Vector<LevelData<EBCellFAB>* >&  a_data,
                              const Vector<EBLevelGrid>&       a_eblg,
                              const Vector<int>&               a_refRatio);


  ///
  /**
   */
  static void pwqFillPatchOne(LevelData<EBCellFAB>& a_dataFine,
                              LevelData<EBCellFAB>& a_dataCoar,
                              const EBLevelGrid&    a_eblgFine,
                              const EBLevelGrid&    a_eblgCoar,
                              const int&                  a_refRatioCoar);


  ///
  /**
   */
  static void exchangeAll(Vector<LevelData<EBCellFAB>* >& a_phi);

  ///
  /**
   */
  static void exchangeCorners(Vector<LevelData<EBCellFAB>* >& a_data,
                              const ProblemDomain&            a_domain);

  ///
  /**
   */
  static void exchangeComp(Vector<LevelData<EBCellFAB>* >& a_data,
                           const int&                      a_comp);


  ///
  /**
   */
  static void exchangeAll(Vector<LevelData<EBFluxFAB>* >& a_phi);


  ///
  /**
   */
  static void exchangeComp(Vector<LevelData<EBFluxFAB>* >& a_data,
                           const int&                      a_comp);


  ///
  /**
   */
  static void exchangeAll(Vector<RefCountedPtr<LevelData<EBCellFAB> > >& a_phi);


  ///
  /**
   */
  static void exchangeComp(Vector<RefCountedPtr<LevelData<EBCellFAB> > >& a_data,
                           const int&                      a_comp);


  ///
  /**
   */
  static void exchangeAll(Vector<RefCountedPtr<LevelData<EBFluxFAB> > >& a_phi);


  ///
  /**
   */
  static void exchangeComp(Vector<RefCountedPtr<LevelData<EBFluxFAB> > >& a_data,
                           const int&                      a_comp);


  ///
  /**
   */
  static void coarsenDown(Vector<LevelData<EBCellFAB>* >&  a_data,
                          const Vector<EBLevelGrid>&       a_eblg,
                          const Vector<int>&               a_refRatio);


  ///
  /**
   */
  static void averageDown(Vector<LevelData<EBCellFAB>* >&  a_data,
                          const Vector<EBISLayout>&        a_ebisl,
                          const Vector<DisjointBoxLayout>& a_dbl,
                          const Vector<ProblemDomain>&     a_domain,
                          const Vector<int>&               a_refRatio);

  ///
  /**
   */
  static void averageDown(Vector<LevelData<EBCellFAB>* >&  a_data,
                          const Vector<EBLevelGrid>&       a_eblg,
                          const Vector<int>&               a_refRatio);


  ///
  /**
   */
  static void averageDown(LevelData<EBCellFAB>&       a_dataCoar,
                          const LevelData<EBCellFAB>& a_dataFine,
                          const EBISLayout&           a_ebislCoar,
                          const EBISLayout&           a_ebislFine,
                          const DisjointBoxLayout&    a_dblCoar,
                          const DisjointBoxLayout&    a_dblFine,
                          const ProblemDomain&        a_domainCoar,
                          const int&                  a_refRatio);


  static Real subtractOffMean(Vector<LevelData<EBCellFAB>* >&  a_data,
                              const Vector<EBLevelGrid>&       a_eblg,
                              const Vector<int> &              a_refRat);

  static void averageDown(Vector<LevelData<EBFluxFAB> * >& a_data,
                          const Vector<EBLevelGrid>&       a_eblg,
                          const Vector<int>&               a_refRatio);
  ///
  /**
   */
  static void averageDown(Vector<LevelData<EBFluxFAB>* >&  a_data,
                          const Vector<EBISLayout>&        a_ebisl,
                          const Vector<DisjointBoxLayout>& a_dbl,
                          const Vector<ProblemDomain>&     a_domain,
                          const Vector<int>&               a_refRatio);


  ///
  /**
   */
  static void averageDown(Vector<RefCountedPtr<LevelData<EBCellFAB> > >& a_data,
                          const Vector<EBISLayout>&                      a_ebisl,
                          const Vector<DisjointBoxLayout>&               a_dbl,
                          const Vector<ProblemDomain>&                   a_domain,
                          const Vector<int>&                             a_refRatio);


  ///
  /**
   */
  static void averageDown(Vector<RefCountedPtr<LevelData<EBFluxFAB> > >& a_data,
                          const Vector<EBISLayout>&                      a_ebisl,
                          const Vector<DisjointBoxLayout>&               a_dbl,
                          const Vector<ProblemDomain>&                   a_domain,
                          const Vector<int>&                             a_refRatio);

  ///
  /**
   */
  static void averageDown(Vector<RefCountedPtr<LevelData<EBFluxFAB> > >& a_data,
                          const Vector<EBLevelGrid>&                     a_eblg,
                          const Vector<int>&                             a_refRatio);


  ///
  /**
   */
  static void setCoveredAMRVal(Vector<LevelData<EBCellFAB>* >& a_data,
                               const Vector<EBISLayout>&       a_ebisl,
                               const Vector<int>&              a_refRat,
                               const Real&                     a_value);

  static void setCoveredAMRVal(Vector<LevelData<EBCellFAB>* >& a_data,
                               const Vector<EBLevelGrid>&      a_eblg,
                               const Vector<int>&              a_refRat,
                               const Real&                     a_value);


  ///
  /**
   */
  static void setCoveredVal(Vector<LevelData<EBCellFAB>* >&a_data,
                            const Real& a_value);


  ///
  /**
   */
  static void deleteAMRData(Vector<LevelData<EBCellFAB>* >&  a_amrData);


  ///
  /**
   */
  static void defineAMRData(Vector<LevelData<EBCellFAB>* >&  a_amrData,
                            const Vector<EBLevelGrid>&       a_eblg,
                            const IntVect&                   a_ghosts,
                            const int&                       a_nComp,
                            const int&                       a_numLevels);
  ///
  /**
   */
  static void defineAMRData(Vector<LevelData<EBFluxFAB>* >&  a_amrData,
                            const Vector<EBLevelGrid>&       a_eblg,
                            const IntVect&                   a_ghosts,
                            const int&                       a_nComp,
                            const int&                       a_numLevels);


  ///
  /**
   */
  static void defineAMRData(Vector<RefCountedPtr<LevelData<EBCellFAB> > >& a_amrData,
                            const Vector<EBLevelGrid>&                     a_eblg,
                            const IntVect&                                 a_ghosts,
                            const int&                                     a_nComp,
                            const int&                                     a_numLevels);


  ///
  /**
   */
  static void defineAMRData(Vector<RefCountedPtr<LevelData<EBFluxFAB> > >& a_amrData,
                            const Vector<EBLevelGrid>&                     a_eblg,
                            const IntVect&                                 a_ghosts,
                            const int&                                     a_nComp,
                            const int&                                     a_numLevels);


  ///
  /**
   */
  static void scale(Vector<LevelData<EBFluxFAB>* >& a_lhs,
                    const Real&                     a_scale) ;


  ///
  /**
   */
  static void scale(Vector<LevelData<EBCellFAB>* >& a_lhs,
                    const Real&                     a_scale) ;


  ///
  /**
   */
  static void scale(Vector<LevelData<EBCellFAB>* >& a_lhs,
                    const Real&                     a_scale,
                    const int&                      a_comp);


  ///
  /**
   */
  static void setToZero(Vector<LevelData<EBCellFAB>* >& a_result);


  ///
  /**
   */
  static void setVal(Vector<LevelData<EBCellFAB>* >& a_result,
                     const Real&                     a_value);


  ///
  /**
   */
  static void setVal(Vector<LevelData<EBCellFAB>* >& a_result,
                     const Real&                     a_value,
                     const int&                      a_comp);


  ///
  /**
   */
  static void setToZero(Vector<LevelData<EBFluxFAB>* >& a_result);


  ///
  /**
   */
  static void setVal(Vector<LevelData<EBFluxFAB>* >& a_result,
                     const Real&                     a_value);


  ///
  /**
   */
  static void assign(Vector<LevelData<EBCellFAB>* >&       a_lhs,
                     const Vector<LevelData<EBCellFAB>* >& a_rhs);


  ///
  /**
   */
  static void assign(Vector<RefCountedPtr<LevelData<EBCellFAB> > >&       a_lhs,
                     const Vector<LevelData<EBCellFAB>* >&                a_rhs);


  ///
  /**
   */
  static void assign(Vector<LevelData<EBFluxFAB>* >&       a_lhs,
                     const Vector<LevelData<EBFluxFAB>* >& a_rhs);


  ///
  /**
   */
  static void assign(Vector<RefCountedPtr<LevelData<EBFluxFAB> > >& a_lhs,
                     const Vector<LevelData<EBFluxFAB>* >&          a_rhs);


  ///
  /**
   */
  static void assign(Vector<LevelData<EBCellFAB>* >&       a_to,
                     const Vector<LevelData<EBCellFAB>* >& a_from,
                     const Interval&                       a_toInterval,
                     const Interval&                       a_fromInterval);


  ///
  /**
   */
  static void assign(Vector<RefCountedPtr<LevelData<EBCellFAB> > >& a_to,
                     const Vector<LevelData<EBCellFAB>* >&          a_from,
                     const Interval&                                a_toInterval,
                     const Interval&                                a_fromInterval);


  ///
  /**
   */
  static void incr(Vector<LevelData<EBCellFAB>* >& a_lhs,
                   const Vector<LevelData<EBCellFAB>* >&   a_rhs,
                   const Real& a_scale);



  ///
  /**
   */
  static void incr(Vector<LevelData<EBCellFAB>* >& a_lhs,
                   const Real& a_scale);


  ///
  /**
   */
  static void axby(Vector<LevelData<EBCellFAB>* >&       a_lhs,
                   const Vector<LevelData<EBCellFAB>* >& a_x,
                   const Vector<LevelData<EBCellFAB>* >& a_y,
                   const Real& a,
                   const Real& b);


  ///
  /**
   */
  static void axby(Vector<LevelData<EBCellFAB>* >&       a_lhs,
                   const Vector<LevelData<EBCellFAB>* >& a_x,
                   const Vector<LevelData<EBCellFAB>* >& a_y,
                   const Real& a_a,
                   const Real& a_b,
                   const int&  a_lhsComp,
                   const int&  a_xComp,
                   const int&  a_yComp);


  ///
  /**
   */
  static  void sum(Vector<LevelData<EBCellFAB>* >&       a_result,
                   const Vector<LevelData<EBCellFAB>* >& a_in1,
                   const Vector<LevelData<EBCellFAB>* >& a_in2);


  ///
  /**
   */
  static  void addConstant(Vector<LevelData<EBCellFAB>* >& a_data,
                           const Real&                     a_constant);


  ///
  /**
   */
  static  void product(Vector<LevelData<EBCellFAB>* >&       a_result,
                       const Vector<LevelData<EBCellFAB>* >& a_in1,
                       const Vector<LevelData<EBCellFAB>* >& a_in2);


  ///
  /**
   */
  static  void product(Vector<LevelData<EBCellFAB>* >&       a_result,
                       const Vector<LevelData<EBCellFAB>* >& a_in1,
                       const Vector<LevelData<EBCellFAB>* >& a_in2,
                       const int&                            a_rComp,
                       const int&                            a_1Comp,
                       const int&                            a_2Comp);


  ///
  /**
   */
  static void divideVectorByScalar(Vector<LevelData<EBCellFAB>* >&       a_vectorOut,
                                   const Vector<LevelData<EBCellFAB>* >& a_vectorIn,
                                   const Vector<LevelData<EBCellFAB>* >& a_scalar);

  ///
  /**
   */
  static  void divide(Vector<LevelData<EBCellFAB>* >&       a_result,
                      const Vector<LevelData<EBCellFAB>* >& a_in1,
                      const Vector<LevelData<EBCellFAB>* >& a_in2);


  ///
  /**
   */
  static  void divide(Vector<LevelData<EBCellFAB>* >&       a_result,
                      const Vector<LevelData<EBCellFAB>* >& a_in1,
                      const Vector<LevelData<EBCellFAB>* >& a_in2,
                      const int&                            a_rComp,
                      const int&                            a_1Comp,
                      const int&                            a_2Comp);


  ///
  /**
   */
  static  void product(Vector<LevelData<EBFluxFAB>* >&       a_result,
                       const Vector<LevelData<EBFluxFAB>* >& a_in1,
                       const Vector<LevelData<EBFluxFAB>* >& a_in2);


  ///
  /**
   */
  static  void kappaWeight(Vector<LevelData<EBCellFAB>* >& a_data);


  ///
  /**
   */
  static  void kappaScale(Vector<LevelData<EBCellFAB>* >& a_data,
                          const Real&                     a_scale);


  ///
  /**
   */
  static  Real kappaNorm(Real&                                 a_volume,
                         const Vector<LevelData<EBCellFAB>* >& a_data,
                         int                                   a_which,
                         const Vector<ProblemDomain>&          a_domain,
                         int                                   a_p=2);


  ///
  /**
   */
  static  Real noKappaNorm(Real&                                 a_volume,
                           const Vector<LevelData<EBCellFAB>* >& a_data,
                           int                                   a_which,
                           const Vector<ProblemDomain>&          a_domain,
                           int                                   a_p=2);


  ///
  /**
   */
  static  Real kappaDotProduct(Real&                                 a_volume,
                               const Vector<LevelData<EBCellFAB>* >& a_data1,
                               const Vector<LevelData<EBCellFAB>* >& a_data2,
                               int                                   a_which,
                               const Vector<ProblemDomain>&          a_domain);


  ///
  /**
   */
  static  Real noKappaDotProduct(Real&                                 a_volume,
                                 const Vector<LevelData<EBCellFAB>* >& a_data1,
                                 const Vector<LevelData<EBCellFAB>* >& a_data2,
                                 int                                   a_which,
                                 const Vector<ProblemDomain>&          a_domain);



  ///
  /**
   */
  static void setMaxMin(Vector<LevelData<EBCellFAB>* >& a_data,
                        const Real&                     a_maxVal,
                        const Real&                     a_minVal,
                        const int&                      a_comp);

  ///
  /**
   */
  static void getMaxMin(Real&                                 a_maxVal,
                        Real&                                 a_minVal,
                        const Vector<LevelData<EBCellFAB>* >& a_data,
                        const int&                            a_comp);


  ///
  /**
   */
  static Real subtractOffMean(Vector<LevelData<EBCellFAB>* >&  a_data,
                              const Vector<DisjointBoxLayout>& a_grids,
                              const Vector<EBISLayout>&        a_ebisl,
                              const Vector<int> &              a_refRat);


  ///
  /**
   */
  static void checkThisData(const Vector<LevelData<EBCellFAB>* >& a_data,
                            const string&                         a_name,
                            const IntVect&                        a_iv1 = IntVect::Zero,
                            const IntVect&                        a_iv2 = IntVect::Zero,
                            const Real&                           a_shift = 0.0);


  ///
  /**
   */
  static bool checkNANINF(const Vector<LevelData<EBCellFAB>* >& a_data,
                          const IntVect&                        a_iv1 = IntVect::Zero,
                          const IntVect&                        a_iv2 = IntVect::Zero,
                          const Real&                           a_shift = 0.0);


  ///
  /**
   */
  static void checkThisData(const Vector<LevelData<EBFluxFAB>* >& a_data,
                            const string&                         a_name);


  ///
  /**
   */
  static bool checkNANINF(const Vector<LevelData<EBFluxFAB>* >& a_data);



  ///
  /**
   */
  static  Real sum(const Vector<LevelData<EBCellFAB>* >& a_data,
                   const Vector<DisjointBoxLayout>&      a_grids,
                   const Vector<EBISLayout>&             a_ebisl,
                   const Vector<int> &                   a_refRat,
                   int   a_comp,
                   bool  a_mutiplyByKappa);

protected:

};
#include "NamespaceFooter.H"
#endif

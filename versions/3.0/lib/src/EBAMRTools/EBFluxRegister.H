#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

#ifndef _EBFLUXREGISTER_H_
#define _EBFLUXREGISTER_H_

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
#include "LayoutData.H"
#include "NamespaceHeader.H"
class EBIndexSpace;

class EBCoarToFineRedist;
class EBCoarToCoarRedist;

///EBFluxRegister-A class to encapsulate a levels worth of flux registers.
/**
   A EBFluxRegister  handles all the data choreography
   necessary to create a levels worth of flux registers.
 */
class EBFluxRegister
{
public:

  ///
  /**
     Default constructor.  Leaves object undefined.
  */
  EBFluxRegister();

  ///
  /**
     Full constructor.  Calls the define function which creates
     a levels worth of flux registers.
  */
  EBFluxRegister(const DisjointBoxLayout& a_dblFine,
                 const DisjointBoxLayout& a_dblCoar,
                 const EBISLayout&        ebislFine,
                 const EBISLayout&        a_ebislCoar,
                 const Box&               a_domainCoar,
                 const int&               a_nref,
                 const int&               a_nvar,
                 const EBIndexSpace*      ebisPtr);

  ///
  ~EBFluxRegister();

  ///
  /**
     Full define function.  Creates a levels worth of flux registers.
  */
  void
  define(const DisjointBoxLayout& a_dblFine,
         const DisjointBoxLayout& a_dblCoar,
         const EBISLayout&        a_ebislFine,
         const EBISLayout&        a_ebislCoar,
         const ProblemDomain&     a_domainCoar,
         const int&               a_nref,
         const int&               a_nvar,
         const EBIndexSpace*      ebisPtr);

  ///
  /**
      Initialize values of registers  to zero.
  */
  void
  setToZero();

  ///
  /**
     increments the register with data from coarseFlux, multiplied by scale.
     Increment coarse buffer with the area-weighted sum
     of the fluxes on the faces.
     buf += -scale*sum(areaFrac*flux)
     CoarseFlux must contain the coarse fluxes in the dir direction
     for the grid   m_coarseLayout[coarseDataIndex].
  */
  void
  incrementCoarseRegular(const EBFaceFAB&   a_coarseFlux,
                         const Real&        a_scale,
                         const DataIndex&   a_coarseDataIndex,
                         const Interval&    a_variables,
                         const int&         a_dir);

  ///
  void
  incrementCoarseBoth(const EBFaceFAB&   a_coarseFlux,
                      const Real&        a_scale,
                      const DataIndex&   a_coarseDataIndex,
                      const Interval&    a_variables,
                      const int&         a_dir);

  void incrementCoarse(const EBFaceFAB&   a_coarseFlux,
                       const Real&        a_scale,
                       const DataIndex&   a_coarseDataIndex,
                       const Interval&    a_varSrc,
                       const Interval&    a_varDst,
                       const int&         a_dir)
  {
    CH_assert(a_varSrc==a_varDst);
    incrementCoarseBoth(a_coarseFlux, a_scale, a_coarseDataIndex, a_varSrc, a_dir);
  }

  ///
  /**
     Here is the story on flux registers for RZ.
     (1)  You do need to multiply the radius into the flux before you send
     them to the flux register increment functions.  See LevelGodunov calls
     EBPatchGodunov::assembleFluxReg which  does this:

     newflux = radius*godunovflux(chf_ix[i;j;k], iv)
     if(  ((iv .eq. CMOMX).and.(idir.eq.0)).or.
     &        ((iv .eq. CMOMY).and.(idir.eq.1))) then
     newflux = newflux +
     &           radius*godunovflux(chf_ix[i;j;k], CPRES)
     endif

     fluxregflux(chf_ix[i;j;k], iv) = newflux

     You will have to do something similar.

     (2)  The math of what happens in the RZ functions of flux register is
     in section 5.3 of the ChomboDoc/EBAMRGodunov document.   I will
     here (and in  EBFluxRegister.H) explain what each function is doing.
     In all of this kvol = (1/2*cellvol)*\int(r^2 nr dl), the integral
     volume fraction of the cell.

     (2.1) EBFluxRegister::incrementCoarse{reg,irreg}lRZ,
     coarsebuffer = dt*area*fluxregflux*volFrac/(kvol*rad)
     The time step is the scale sent in.
     coarReg += dt_c*volFrac*areaFrac*fluxregflux/(kvol*rad_cell_center);

     (2.2) EBFluxRegister::incrementFine{reg,irreg}RZ.  f subscripts are
     from the fine level, c from the coarse.
     fineBuffer += dt_f*volFrac_coar*areaFrac_f*flux/(kvol_c*rad_c_cell*nref^{D-1});

  */
  void
  incrementCoarseRegulRZ(const EBFaceFAB&  a_coarseFlux,
                         const Real&       a_scale,
                         const DataIndex&  a_coarseDataIndex,
                         const Interval&   a_variables,
                         const int&        a_dir,
                         const Real&       a_dx);

  ///
  /**
     increments the register with data from coarseFlux, multiplied by scale.
     Increment coarse buffer with the area-weighted sum
     of the fluxes on the faces.
     buf += -scale*sum(areaFrac*flux)
     CoarseFlux must contain the coarse fluxes in the dir direction
     for the grid   m_coarseLayout[coarseDataIndex].
  */
  void
  incrementCoarseIrregular(const BaseIFFAB<Real>& a_coarseFlux,
                           const Real&            a_scale,
                           const DataIndex&       a_coarseDataIndex,
                           const Interval&        a_variables,
                           const int&             a_dir);

  ///
  /**
     For when you have an EBFaceFAB lying around instead of a BaseIFFAB
   */
  void incrementCoarseIrregular(const EBFaceFAB&       a_coarFlux,
                                const Real&            a_scale,
                                const DataIndex&       a_coarDatInd,
                                const Interval&        a_variables,
                                const int&             a_dir);
  ///
  /**
     Here is the story on flux registers for RZ.
     (1)  You do need to multiply the radius into the flux before you send
     them to the flux register increment functions.  See LevelGodunov calls
     EBPatchGodunov::assembleFluxReg which  does this:

     newflux = radius*godunovflux(chf_ix[i;j;k], iv)
     if(  ((iv .eq. CMOMX).and.(idir.eq.0)).or.
     &        ((iv .eq. CMOMY).and.(idir.eq.1))) then
     newflux = newflux +
     &           radius*godunovflux(chf_ix[i;j;k], CPRES)
     endif

     fluxregflux(chf_ix[i;j;k], iv) = newflux

     You will have to do something similar.

     (2)  The math of what happens in the RZ functions of flux register is
     in section 5.3 of the ChomboDoc/EBAMRGodunov document.   I will
     here (and in  EBFluxRegister.H) explain what each function is doing.
     In all of this kvol = (1/2*cellvol)*\int(r^2 nr dl), the integral
     volume fraction of the cell.

     (2.1) EBFluxRegister::incrementCoarse{reg,irreg}lRZ,
     coarsebuffer = dt*area*fluxregflux*volFrac/(kvol*rad)
     The time step is the scale sent in.
     coarReg += dt_c*volFrac*areaFrac*fluxregflux/(kvol*rad_cell_center);

     (2.2) EBFluxRegister::incrementFine{reg,irreg}RZ.  f subscripts are
     from the fine level, c from the coarse.
     fineBuffer += dt_f*volFrac_coar*areaFrac_f*flux/(kvol_c*rad_c_cell*nref^{D-1});

  */
  void
  incrementCoarseIrregulRZ(const BaseIFFAB<Real>& a_coarseFlux,
                           const Real&            a_scale,
                           const DataIndex&       a_coarseDataIndex,
                           const Interval&        a_variables,
                           const int&             a_dir,
                           const Real&            a_dx);

  ///
  /**
     increments the register with data from fineFlux, multiplied by scale.
     fineFlux must contain the coarse fluxes in the dir direction
     for the grid m_fineLayout[fineDataIndex].
     This operation is local.
  */
  void
  incrementFineRegular(const EBFaceFAB&      a_fineFlux,
                       const Real&           a_scale,
                       const DataIndex&      a_fineDataIndex,
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

  void
  incrementFine(const EBFaceFAB&      a_fineFlux,
                const Real&           a_scale,
                const DataIndex&      a_fineDataIndex,
                const Interval&       a_varSrc,
                const Interval&       a_varDst,
                const int&            a_dir)
  {
    CH_assert(a_varSrc == a_varDst);
    for(SideIterator sit; sit.ok(); ++sit)
      {
        incrementFineBoth(a_fineFlux, a_scale, a_fineDataIndex, a_varSrc, a_dir, sit());
      }
  }




  ///
  /**
     Here is the story on flux registers for RZ.
     (1)  You do need to multiply the radius into the flux before you send
     them to the flux register increment functions.  See LevelGodunov calls
     EBPatchGodunov::assembleFluxReg which  does this:

     newflux = radius*godunovflux(chf_ix[i;j;k], iv)
     if(  ((iv .eq. CMOMX).and.(idir.eq.0)).or.
     &        ((iv .eq. CMOMY).and.(idir.eq.1))) then
     newflux = newflux +
     &           radius*godunovflux(chf_ix[i;j;k], CPRES)
     endif

     fluxregflux(chf_ix[i;j;k], iv) = newflux

     You will have to do something similar.

     (2)  The math of what happens in the RZ functions of flux register is
     in section 5.3 of the ChomboDoc/EBAMRGodunov document.   I will
     here (and in  EBFluxRegister.H) explain what each function is doing.
     In all of this kvol = (1/2*cellvol)*\int(r^2 nr dl), the integral
     volume fraction of the cell.

     (2.1) EBFluxRegister::incrementCoarse{reg,irreg}lRZ,
     coarsebuffer = dt*area*fluxregflux*volFrac/(kvol*rad)
     The time step is the scale sent in.
     coarReg += dt_c*volFrac*areaFrac*fluxregflux/(kvol*rad_cell_center);

     (2.2) EBFluxRegister::incrementFine{reg,irreg}RZ.  f subscripts are
     from the fine level, c from the coarse.
     fineBuffer += dt_f*volFrac_coar*areaFrac_f*flux/(kvol_c*rad_c_cell*nref^{D-1});

  */
  void
  incrementFineRegulRZ(const EBFaceFAB&      a_fineFlux,
                       const Real&           a_scale,
                       const DataIndex&      a_fineDataIndex,
                       const Interval&       a_variables,
                       const int&            a_dir,
                       const Side::LoHiSide& a_sd,
                       const Real&           a_dx);

  ///
  /**
     increments the register with data from fineFlux, multiplied by scale.
     fineFlux must contain the coarse fluxes in the dir direction
     for the grid m_fineLayout[fineDataIndex].
     This operation is local.
  */
  void
  incrementFineIrregular(const BaseIFFAB<Real>& a_fineFlux,
                         const Real&            a_scale,
                         const DataIndex&       a_fineDataIndex,
                         const Interval&        a_variables,
                         const int&             a_dir,
                         const Side::LoHiSide&  a_sd);

  ///
  /**
     For when you have a EBFaceFAB lying around instead of a BaseIFFAB
  */
  void
  incrementFineIrregular(const EBFaceFAB&       a_fineFlux,
                         const Real&            a_scale,
                         const DataIndex&       a_fineDataIndex,
                         const Interval&        a_variables,
                         const int&             a_dir,
                         const Side::LoHiSide&  a_sd);

  ///
  void
  incrementFineBothRZ(const EBFaceFAB&      a_fineFlux,
                      const Real&           a_scale,
                      const DataIndex&      a_fineDataIndex,
                      const Interval&       a_variables,
                      const int&            a_dir,
                      const Side::LoHiSide& a_sd,
                      const Real&       a_dx);
  ///
  void
  incrementCoarseBothRZ(const EBFaceFAB&   a_coarseFlux,
                        const Real&        a_scale,
                        const DataIndex&   a_coarseDataIndex,
                        const Interval&    a_variables,
                        const int&         a_dir,
                        const Real&        a_dx);
  ///
  /**
     Here is the story on flux registers for RZ.
     (1)  You do need to multiply the radius into the flux before you send
     them to the flux register increment functions.  See LevelGodunov calls
     EBPatchGodunov::assembleFluxReg which  does this:

     newflux = radius*godunovflux(chf_ix[i;j;k], iv)
     if(  ((iv .eq. CMOMX).and.(idir.eq.0)).or.
     &        ((iv .eq. CMOMY).and.(idir.eq.1))) then
     newflux = newflux +
     &           radius*godunovflux(chf_ix[i;j;k], CPRES)
     endif

     fluxregflux(chf_ix[i;j;k], iv) = newflux

     You will have to do something similar.

     (2)  The math of what happens in the RZ functions of flux register is
     in section 5.3 of the ChomboDoc/EBAMRGodunov document.   I will
     here (and in  EBFluxRegister.H) explain what each function is doing.
     In all of this kvol = (1/2*cellvol)*\int(r^2 nr dl), the integral
     volume fraction of the cell.

     (2.1) EBFluxRegister::incrementCoarse{reg,irreg}lRZ,
     coarsebuffer = dt*area*fluxregflux*volFrac/(kvol*rad)
     The time step is the scale sent in.
     coarReg += dt_c*volFrac*areaFrac*fluxregflux/(kvol*rad_cell_center);

     (2.2) EBFluxRegister::incrementFine{reg,irreg}RZ.  f subscripts are
     from the fine level, c from the coarse.
     fineBuffer += dt_f*volFrac_coar*areaFrac_f*flux/(kvol_c*rad_c_cell*nref^{D-1});

  */
  void
  incrementFineIrregulRZ(const BaseIFFAB<Real>& a_fineFlux,
                         const Real&            a_scale,
                         const DataIndex&       a_fineDataIndex,
                         const Interval&        a_variables,
                         const int&             a_dir,
                         const Side::LoHiSide&  a_sd,
                         const Real&            a_dx);

  ///
  /**
     increments uCoarse with the reflux divergence of the
     contents of the flux register.
     This operation is global and blocking.
  */
  void reflux(LevelData<EBCellFAB>& a_uCoarse,
              const Interval&       a_variables,
              const Real&           a_scale);

  ///
  /**
     increments uCoarse with the reflux divergence of the
     contents of the flux register.
     beta is a pointwise multiplier to scale
     This operation is global and blocking.
  */
  void reflux(LevelData<EBCellFAB>& a_uCoarse,
              const Interval&       a_variables,
              const Real&           a_scale,
              const LevelData<EBCellFAB>& a_beta);

  ///
  void reflux(LevelData<EBCellFAB>& a_uCoarse,
              const Interval&       a_variables,
              const Real&           a_scale,
              const int&            a_refluxDir);

  ///
  void reflux(LevelData<EBCellFAB>& a_uCoarse,
              const Interval&       a_variables,
              const Real&           a_scale,
              const int&            a_refluxDir,
              const LevelData<EBCellFAB>& a_beta);

  ///
  /**
     Increments redistribution register with left-over mass
     from reflux divergence.
  */
  void incrementRedistRegister(EBCoarToFineRedist& a_register,
                               const Interval&     a_variables,
                               const Real&         a_scale);

  ///
  /**
     Increments redistribution register with left-over mass
     from reflux divergence.
  */
  void incrementRedistRegister(EBCoarToCoarRedist& a_register,
                               const Interval&     a_variables,
                               const Real&         a_scale);

  ///
  void incrementRedistRegistRZ(EBCoarToCoarRedist& a_register,
                               const Interval&     a_variables,
                               const Real&         a_dx);

  ///
  /**
     Increments mass array  with left-over mass
     from reflux divergence.
  */
  void incrementDensityArray(LevelData<EBCellFAB>& a_coarMass,
                             const Interval&       a_variables,
                             const Real&           a_scale);

  ///
  bool isDefined() const;

  ///
  void dumpCoar(const int& a_idir,
                const Side::LoHiSide& a_sd);
  ///
  void dumpFine(const int& a_idir,
                const Side::LoHiSide& a_sd);

  ///
  void dumpBIVF(const BaseIVFAB<Real>& reg,
                const EBISBox& a_ebisBox);
protected:

  static int index(int a_dir, Side::LoHiSide a_side);

  //internal use only
  void setDefaultValues();

  // both fabCoarse and fabFine have fabs that are
  // defined in the coarse index space.  fabCoarse
  // can ONLY be indexed into by the DataIndex stored in
  // coarToCoarMap.  fabFine is indexed by the same indexing
  // used in the fine DisjointBoxLayout.

  //input fine grid
  DisjointBoxLayout    m_gridsFine;
  //input coarse grid
  DisjointBoxLayout    m_gridsCoar;

  //grids for fine buffers
  DisjointBoxLayout m_bufGridsFine[SpaceDim*2];
  //grids for coarse buffers
  DisjointBoxLayout m_bufGridsCoar[SpaceDim*2];

  //fine buffers
  LevelData<BaseIVFAB<Real> > m_regsCoar[SpaceDim*2];
  //coarse buffers
  LevelData<BaseIVFAB<Real> > m_regsFine[SpaceDim*2];

  //scratch space for refluxing
  LevelData<BaseIVFAB<Real> > m_scratchc[SpaceDim*2];

  //sets for fine buffers
  LayoutData<IntVectSet> m_cfivsFine[SpaceDim*2];
  //sets for coarse buffers
  LayoutData<IntVectSet> m_cfivsCoar[SpaceDim*2];

  //ebisl of input fine grid
  EBISLayout m_ebislFine;
  //ebisl of input coar grid
  EBISLayout m_ebislCoar;

  //ebisl of coarse buffers
  EBISLayout m_ebislBufCoar[2*SpaceDim];
  //ebisl of fine buffers
  EBISLayout m_ebislBufFine[2*SpaceDim];

  //map from input coarse grid to coarse buffers
  LayoutData<Vector<DataIndex>  > m_coarIndexMap[SpaceDim*2];

  bool m_isDefined;
  int m_nComp;
  int m_refRat;
  ProblemDomain m_domainCoar;

private:
  //for all the usual reasons,
  //there is no copy constructor for this class
  //and there is no operator= for this class
  //(there is no rule six).
  void operator= (const EBFluxRegister&)
  {
    MayDay::Error("invalid operator");
  }
  EBFluxRegister(const EBFluxRegister&)
  {
    MayDay::Error("invalid operator");
  }
};
#include "NamespaceFooter.H"
#endif

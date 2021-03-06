#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

#ifndef _EBPatchPolytropic_H_
#define _EBPatchPolytropic_H_
#include "EBPatchGodunov.H"
#include "Stencils.H"
#include "CoveredFaceFAB.H"
#include "BaseIFFAB.H"

///
/**
   An EBPatchGodunov for gamma-law gas dynamics.
*/
class EBPatchPolytropic : public EBPatchGodunov
{
public:
  ///
  /**
   */
  EBPatchPolytropic();

  ///
  /**
   */
  virtual ~EBPatchPolytropic();

  ///
  /**
   */
  void
  normalPred(EBCellFAB&       a_primLo,
             EBCellFAB&       a_primHi,
             const EBCellFAB& a_primState,
             const EBCellFAB& a_slopePrim,
             const Real&      a_scale,
             const int&       a_dir,
             const Box&       a_box);

  ///
  /**
   */
  void
  computeEBIrregFlux(BaseIVFAB<Real>&  a_ebIrregFlux,
                     const EBCellFAB&  a_primState,
                     const EBCellFAB   a_slopePrim[SpaceDim],
                     const IntVectSet& a_irregIVS,
                     const EBCellFAB&  a_flux);

  ///
  void
  computeBoundaryPress(BaseIVFAB<Real>&  a_pstar,
                       BaseIVFAB<Real>&  a_vnorm,
                       BaseIVFAB<Real>&  a_normDotModiano,
                       BaseIVFAB<Real>&  a_normError,
                       BaseIVFAB<Real>&  a_pstarNoVNorm,
                       BaseIVFAB<Real>&  a_vDotTrueNorm,
                       const RealVect&   a_modianoAxis,
                       const RealVect&   a_modianoCorner,
                       const EBCellFAB&  a_primState,
                       const EBCellFAB&  a_nonConsDivF,
                       const EBCellFAB   a_slopePrim[SpaceDim],
                       const IntVectSet& a_irregIVS,
                       const EBCellFAB&  a_flux);

  ///
  void
  computeBoundaryWDivU(BaseIVFAB<Real>&  a_pstar,
                       BaseIVFAB<Real>&  a_vnorm,
                       BaseIVFAB<Real>&  a_normDotModiano,
                       const RealVect&   a_modianoAxis,
                       const EBCellFAB&  a_primState,
                       const EBCellFAB&  a_nonConsDivF,
                       const EBCellFAB   a_slopePrim[SpaceDim],
                       const IntVectSet& a_irregIVS,
                       const EBCellFAB&  a_flux);

  void 
  primitivesAndDivergences(EBCellFAB&          a_nonConsDivF,
                           EBCellFAB&          a_consState,
                           EBFluxFAB&          a_facePrim,
                           BaseIVFAB<Real>     a_coveredPrimMinu[SpaceDim],
                           BaseIVFAB<Real>     a_coveredPrimPlus[SpaceDim],
                           Vector<VolIndex>    a_coveredFaceMinu[SpaceDim],
                           Vector<VolIndex>    a_coveredFacePlus[SpaceDim],
                           EBFluxFAB&          a_flux,
                           BaseIVFAB<Real>&    a_ebIrregFlux,
                           BaseIVFAB<Real>&    a_nonConservativeDivergence,
                           const EBCellFAB&    a_flattening,
                           const EBCellFAB&    a_source,
                           const Box&          a_box,
                           const IntVectSet&   a_ivsSmall,
                           const DataIndex&    a_dit,
                           bool                a_verbose);


  void 
  getRHSTemperature(EBCellFAB&          a_rhsTemp,
                    EBFluxFAB&          a_facePrim,
                    BaseIVFAB<Real>     a_coveredPrimMinu[SpaceDim],
                    BaseIVFAB<Real>     a_coveredPrimPlus[SpaceDim],
                    const Real&         a_specHeat,
                    const Box&          a_box,
                    const IntVectSet&   a_ivsSmall);



  void 
  getRHSVelocity(EBCellFAB&          a_rhsTemp,
                 EBFluxFAB&          a_facePrim,
                 BaseIVFAB<Real>     a_coveredPrimMinu[SpaceDim],
                 BaseIVFAB<Real>     a_coveredPrimPlus[SpaceDim],
                 const Box&          a_box,
                 const IntVectSet&   a_ivsSmall);


  ///
  /**
   */
  Real
  getMaxWaveSpeed(const EBCellFAB& a_consState,
                  const Box& a_box);

  ///
  /**
     Compute the primitive state given the conserved state.
     W = W(U).
  */
  void
  consToPrim(EBCellFAB&       a_primState,
             const EBCellFAB& a_conState,
             const Box&       a_box,
             int              a_logflag,
             bool             a_verbose=false);

  void
  consToPrim(BaseIVFAB<Real>&       a_primState,
             const BaseIVFAB<Real>& a_conState,
             const IntVectSet&      a_ivs);

#ifdef CH_USE_HDF5
  virtual void expressions(HDF5HeaderData& a_expressions);
#endif

  ///
  void  getCoveredValuesPrim(Vector<Real>& a_covValues);

  ///
  void  getCoveredValuesCons(Vector<Real>& a_covValues);

  ///
  /**
   */
  void
  consToPrim(BaseIVFAB<Real>&  a_primState,
             const EBCellFAB&  a_conState,
             const IntVectSet& a_ivs);

  void
  setRZSource(EBCellFAB&       a_source,
              const EBCellFAB& a_consState,
              const Box&       a_box);

  void
  setSource(EBCellFAB&       a_source,
            const EBCellFAB& a_consState,
            const Box&       a_box);

  void
  assembleFluxReg(EBFaceFAB&       a_fluxRegFlux,
                  const EBFaceFAB& a_godunovFlux,
                  const int&       a_idir,
                  const Box&       a_cellBox);

  void
  assembleFluxIrr(BaseIFFAB<Real>&       a_fluxRegFlux,
                  const BaseIFFAB<Real>& a_godunovFlux,
                  const int&             a_idir,
                  const Box&             a_cellBox,
                  const IntVectSet&      a_set);

  ///
  /**
   */
  void
  primToCons(EBCellFAB&       a_primState,
             const EBCellFAB& a_conState,
             const Box&       a_box);

  ///
  /**
   */
  void
  primToCons(BaseIVFAB<Real>&       a_primState,
             const BaseIVFAB<Real>& a_conState,
             const IntVectSet&      a_ivs);

  ///

  /**
     Given input left and right states, compute a suitably-upwinded
     flux (e.g. by solving a Riemann problem), as in
  */
  void
  riemann(EBFaceFAB&       a_flux,
          const EBCellFAB& a_primLeft,
          const EBCellFAB& a_primRight,
          const int&       a_dir,
          const Box&       a_box);

  ///
  /**
     Given input left and right states, compute a suitably-upwinded
     flux (e.g. by solving a Riemann problem).
  */
  void
  riemann(BaseIVFAB<Real>&        a_coveredFlux,
          const BaseIVFAB<Real>&  a_extendedState,
          const EBCellFAB&        a_primState,
          const Vector<VolIndex>& a_region,
          const int&              a_dir,
          const Side::LoHiSide&   a_sd,
          const Box&       a_box);

  /**
     Given input left and right states, compute a suitably-upwinded
     flux (e.g. by solving a Riemann problem), as in
  */
  void
  getFlux(EBFaceFAB&       a_flux,
          const EBFaceFAB& a_prim,
          const int&       a_dir,
          const Box&       a_box);

  ///
  /**
     Given the value of the flux, update the conserved quantities as in
     and modify  in place the flux for the purpose of passing it to a
     LevelFluxRegister.
     \\
     consstate_i += a_scale*(flux_i-1/2 - flux_i+1/2)
  */
  void
  updateCons(EBCellFAB&              a_consState,
             const EBFaceFAB&        a_flux,
             const BaseIVFAB<Real>&  a_coveredFluxMinu,
             const BaseIVFAB<Real>&  a_coveredFluxPlus,
             const Vector<VolIndex>& a_coveredFaceMinu,
             const Vector<VolIndex>& a_coveredFacePlus,
             const int&              a_dir,
             const Box&              a_box,
             const Real&             a_scale);

  void
  updateConsRZ(EBCellFAB&              a_consState,
             const EBFaceFAB&        a_flux,
             const BaseIVFAB<Real>&  a_coveredFluxMinu,
             const BaseIVFAB<Real>&  a_coveredFluxPlus,
             const Vector<VolIndex>& a_coveredFaceMinu,
             const Vector<VolIndex>& a_coveredFacePlus,
             const int&              a_dir,
             const Box&              a_box,
             const Real&             a_scale);


  ///
  /**
   */
  void
  setCoveredConsVals(EBCellFAB& a_consState);

  ///
  /**
     Return the names of the variables.  A default
     implementation is provided that puts in generic names.
  */
  Vector<string> primNames();
  Vector<string> primNamesNoLog();

  ///
  /**
     Return number of components for primitive variables.
  */
  int numPrimitives() const;

  ///
  /**
     Returns number of components for flux variables.
  */
  int numFluxes() const;

  ///
  /**
     Returns number of components for conserved variables.
  */
  int numConserved() const;

  ///
  /**
   */
  int numSlopes() const;

  ///
  /**
     Return true if the application is using flattening.
  */
  bool usesFlattening() const;

  ///
  /**
     Return true if the application is using artificial viscosity.
  */
  bool usesArtificialViscosity() const;

  ///
  /**
     Return true if you are using fourth-order slopes.
     Return false if you are using second-order slopes.
  */
  bool usesFourthOrderSlopes() const;

  ///
  /**
     Returns the interval of component indices in the primitive variable
     EBCellFAB for the velocities.
  */
  Interval velocityInterval() const;

  ///
  /**
     Returns the component index for the density.
  */
  virtual int densityIndex() const;

  ///
  /**
     Returns the component index for the pressure. Called only if flattening is used.
  */
  int pressureIndex() const;

  ///
  /**
     Return the names of the variables.  A default
     implementation is provided that puts in generic names.
  */
  Vector<string> stateNames();

  ///
  /**
     Returns the component index for the bulk modulus, used as a
     normalization to measure shock strength in flattening.
     Called only if flattening is used.
  */
  int bulkModulusIndex() const;

  ///
  /**
     Returns value of artificial viscosity. Called only if artificial
     viscosity is being used.
  */
  Real artificialViscosityCoefficient() const;


  void
  floorPrimitives(EBCellFAB&  a_primState,
                  const Box&  a_box);

  void
  floorConserved(EBCellFAB&  a_consState,
                 const Box&  a_box);


  void
  floorConserved(BaseIVFAB<Real>&  a_consState,
                 const IntVectSet&  a_set);

  void
  floorPrimitives(BaseIVFAB<Real>&  a_primState,
                  const IntVectSet&  a_set);


  /**
   */
  bool isDefined() const;

  ///
  /**
   */
  void
  setGamma(const Real& a_gamma);

  ///
  /**
   */
  Real getGamma() const;

  ///
  /**
   */
  void
  nonconservativeDivergence(EBCellFAB&             a_divF,
                            const EBFluxFAB&       a_flux,
                            const BaseIVFAB<Real>  a_coveredFluxMinu[SpaceDim],
                            const BaseIVFAB<Real>  a_coveredFluxPlus[SpaceDim],
                            const Vector<VolIndex> a_coveredFaceMinu[SpaceDim],
                            const Vector<VolIndex> a_coveredFacePlus[SpaceDim],
                            const Box&             a_box);


  ///
  /**
   */
  void
  consUndividedDivergence(BaseIVFAB<Real>&       a_divF,
                          const BaseIFFAB<Real>  a_centroidFlux[SpaceDim],
                          const BaseIVFAB<Real>& a_ebIrregFlux,
                          const IntVectSet&      a_ivs);

  void doRZCoords(bool a_doRZCoords);

  void
  nonconservativeDivergenRZ(EBCellFAB&             a_divF,
                            const EBFluxFAB&       a_flux,
                            const BaseIVFAB<Real>  a_coveredFluxMinu[SpaceDim],
                            const BaseIVFAB<Real>  a_coveredFluxPlus[SpaceDim],
                            const Vector<VolIndex> a_coveredFaceMinu[SpaceDim],
                            const Vector<VolIndex> a_coveredFacePlus[SpaceDim],
                            const Box&             a_box);

  void
  consUndividedDivergenRZ(BaseIVFAB<Real>&       a_divF,
                          const BaseIFFAB<Real>  a_centroidFlux[SpaceDim],
                          const BaseIVFAB<Real>& a_ebIrregFlux,
                          const IntVectSet&      a_ivs);

private:

  bool m_isGammaSet;
  Real m_gamma;
  bool m_doRZCoords;

private:
  //disallowed for all the usual reasons
  void operator=(const EBPatchPolytropic& a_input)
  {
    MayDay::Error("invalid operator");
  }
  EBPatchPolytropic(const EBPatchPolytropic& a_input)
  {
    MayDay::Error("invalid operator");
  }

};
#endif

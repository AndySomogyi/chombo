#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

#ifndef _EBAMRCNS_H_
#define _EBAMRCNS_H_

#include "AMRLevel.H"
#include "EBCellFAB.H"
#include "EBLevelTGA.H"
#include "BaseIVFAB.H"
#include "LevelData.H"
#include "EBLevelRedist.H"
#include "EBCoarToCoarRedist.H"
#include "EBPatchGodunov.H"
#include "EBPatchGodunovFactory.H"
#include "EBCoarToFineRedist.H"
#include "EBFineToCoarRedist.H"
#include "EBCoarseAverage.H"
#include "EBPWLFineInterp.H"
#include "EBFluxRegister.H"
#include "EBLevelCNS.H"
#include "Box.H"
#include "BiCGStabSolver.H"
#include "EBAMRCNSParams.H"
#include "IntVectSet.H"
#include "AMRMultiGrid.H"
#include "EBPatchPolytropicFactory.H"
#include "AMRTGA.H"
#include "BaseLevelTGA.H"
#include "Vector.H"
#include "DisjointBoxLayout.H"
#include "EBQuadCFInterp.H"
#include "NamespaceHeader.H"


///
/**
 */
class EBAMRCNS : public AMRLevel
{
public:

  EBAMRCNS(const EBAMRCNSParams&                          a_params,
           const RefCountedPtr<EBPatchPolytropicFactory>&    a_godFactory)
  {
    m_params                = a_params;
    m_ebPatchGodunovFactory = a_godFactory;
  }

  virtual ~EBAMRCNS()
  {
  }

  virtual void define(AMRLevel*            a_coarser_level_ptr,
                      const ProblemDomain& a_problem_domain,
                      int                  a_level,
                      int                  a_ref_ratio);



  virtual Real advance();

  virtual void postTimeStep();

  virtual void tagCells(IntVectSet& a_tags);

  virtual void tagCellsInit(IntVectSet& a_tags);

  virtual void regrid(const Vector<Box>& a_new_grids);

  virtual void initialGrid(const Vector<Box>& a_new_grids);

  virtual void initialData();

  virtual void postInitialize();

  virtual Real computeDt();

  virtual Real computeInitialDt();
  LevelData<EBCellFAB>& getStateNew();
  LevelData<EBCellFAB>& getStateOld();
  Real getDt() const;
  EBISLayout getEBISLayout() const;
  void fillConsAndPrim(LevelData<EBCellFAB>& a_data) ;
  static int s_NewPlotFile;
  void getVelocity(LevelData<EBCellFAB>&       a_velo,
                   const LevelData<EBCellFAB>& a_state) const;

  void getTemperature(LevelData<EBCellFAB>&       a_temp,
                      const LevelData<EBCellFAB>& a_state) const;

  void evalOperator(LevelData<EBCellFAB>&              a_source,
                    const LevelData<EBCellFAB>&        a_primitive,
                    const LevelData<EBCellFAB>&        a_primiCoar,
                    TGAHelmOp<LevelData<EBCellFAB> >*  a_op,
                    Real a_alpha, Real a_beta);

#ifdef CH_USE_HDF5
  virtual void writePlotHeaderOld    (HDF5Handle& a_handle) const;
  virtual void writePlotLevelOld     (HDF5Handle& a_handle) const;
  virtual void writePlotHeader       (HDF5Handle& a_handle) const;
  virtual void writePlotLevel        (HDF5Handle& a_handle) const;
  virtual void writeCheckpointHeader (HDF5Handle& a_handle) const;
  virtual void writeCheckpointLevel  (HDF5Handle& a_handle) const;
  virtual void readCheckpointHeader  (HDF5Handle& a_handle);
  virtual void readCheckpointLevel   (HDF5Handle& a_handle);
#endif

  void sumConserved(Real& a_sumcons,
                    const int& a_ivar) const;


  virtual void syncWithFineLevel();

protected:
  void addULaplUToRHSTemperature();

  void getHalfState(LevelData<EBCellFAB>& a_stateInt);
  EBAMRCNSParams m_params;

  EBLevelGrid    m_eblg;

  LevelData<EBCellFAB>         m_stateOld;
  LevelData<EBCellFAB>         m_stateNew;
  LevelData<EBCellFAB>         m_rhsVelo;
  LevelData<EBCellFAB>         m_rhsTemp;
  LevelData<BaseIVFAB<Real> >  m_massDiff;
  LayoutData<IntVectSet>       m_sets;

  EBCoarToFineRedist m_ebCoarToFineRedist;
  EBCoarToCoarRedist m_ebCoarToCoarRedist;
  EBFineToCoarRedist m_ebFineToCoarRedist;

  int m_nComp;
  int m_nPrim;
  int m_nGhost;
  Vector<string> m_stateNames;
  Vector<string> m_primNames;
  RealVect m_dx;
  Real m_dtNew;
  bool m_hasCoarser, m_hasFiner;
  EBCoarseAverage                                               m_ebCoarseAverage;
  EBPWLFineInterp                                               m_ebFineInterp;
  EBFluxRegister                                                m_divFFluxRegister;
  EBFluxRegister                                                m_veloFluxRegister;
  EBFluxRegister                                                m_tempFluxRegister;
  RefCountedPtr<EBPatchGodunov>                                 m_ebPatchGodunov;
  RefCountedPtr<EBPatchPolytropicFactory>                       m_ebPatchGodunovFactory;
  EBLevelCNS                                                    m_ebLevelGodunov;
  EBLevelRedist                                                 m_ebLevelRedist;

  RefCountedPtr<LevelData<EBCellFAB> >                          m_acoTemp;
  RefCountedPtr<LevelData<EBCellFAB> >                          m_acoVelo;
  RefCountedPtr<LevelData<EBFluxFAB> >                          m_eta;
  RefCountedPtr<LevelData<EBFluxFAB> >                          m_bcoTemp;
  RefCountedPtr<LevelData<EBFluxFAB> >                          m_lambda;
  RefCountedPtr<LevelData<BaseIVFAB<Real> > >                   m_etaIrreg;
  RefCountedPtr<LevelData<BaseIVFAB<Real> > >                   m_bcoTempIrreg;
  RefCountedPtr<LevelData<BaseIVFAB<Real> > >                   m_lambdaIrreg;
  RefCountedPtr<EBQuadCFInterp>                                 m_quadCFI;

  static RefCountedPtr<AMRLevelOpFactory<LevelData<EBCellFAB> > >  s_tempFactory;
  static RefCountedPtr<AMRLevelOpFactory<LevelData<EBCellFAB> > >  s_veloFactory;
  static RefCountedPtr<AMRMultiGrid<     LevelData<EBCellFAB> > >  s_tempSolver;
  static RefCountedPtr<AMRMultiGrid<     LevelData<EBCellFAB> > >  s_veloSolver;
  static BiCGStabSolver<LevelData<EBCellFAB> >                     s_botSolver;
  static RefCountedPtr< EBLevelTGA >                               s_veloLevTGA;
  static RefCountedPtr< EBLevelTGA >                               s_tempLevTGA;
  void fillCoefficients();
  void    defineSolvers();

  void coarseFineIncrement();
  void postTimeStepRefluxRedistDance();
private:

  Real calculateMass() const;

  EBAMRCNS* getCoarserLevel() const;
  EBAMRCNS* getFinerLevel() const;
  void getTGADiffusion(LevelData<EBCellFAB>& a_diffusion,  
                       LevelData<EBCellFAB>& a_phiold,
                       LevelData<EBCellFAB>& a_rhs, 
                       RefCountedPtr< EBLevelTGA> & a_levTGA,
                       bool a_isVarVel);

  //this one is explicit and also does reweighting of the 
  //redistribution object in the case of mass weighting
  void redistributeDensity();

  //adds redistribution energy to m_rhsTemp
  void redistributeEnergy();

  //adds redistribution momentum to m_rhsVelo
  void redistributeMomentum();

  void fluxDivergence(          LevelData<EBCellFAB>& a_divergeF,
                                LevelData<EBCellFAB>& a_halfStCo);
  void advanceDensity(    const LevelData<EBCellFAB>& a_divergeF);
  void advanceVelocity(   const LevelData<EBCellFAB>& a_divergeF);
  void advanceTemperature(const LevelData<EBCellFAB>& a_divergeF,
                          const LevelData<EBCellFAB>& a_halfStCo);
  void hyperbolicSource(LevelData<EBCellFAB>& a_source,
                        const LevelData<EBCellFAB>& a_halfSt,
                        const LevelData<EBCellFAB>& a_halfCo);

  void levelSetup();
  //weak construction is an
  //evil that must be rooted out.
  //verily, it is a character flaw.
  EBAMRCNS()
  {
    MayDay::Error("invalid operator");
  }
  //disallowed for all the usual reasons
  void operator=(const EBAMRCNS& a_input)
  {
    MayDay::Error("invalid operator");
  }
  EBAMRCNS(const EBAMRCNS& a_input)
  {
    MayDay::Error("invalid operator");
  }

};

#include "NamespaceFooter.H"
#endif

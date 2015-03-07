#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

#ifndef _EBAMRNoSubcycleOld_H_
#define _EBAMRNoSubcycleOld_H_

#include "LevelData.H"
#include "EBCellFAB.H"
#include "EBISLayout.H"
#include "EBIBCFactory.H"
#include "AMRTGA.H"
#include "EBCompositeMACProjector.H"
#include "EBCompositeCCProjector.H"
#include "EBAMRPoissonOp.H"
#include "EBLevelAdvect.H"
#include "EBCoarseAverage.H"
#include "EBSimpleSolver.H"
#include "EBAMRNoSubcycle.H"

#include "UsingNamespace.H"

/// Class to manage non-subcycled advance for EBAMRINS.
/**
*/
class EBAMRNoSubcycleOld
{

public:

  ///
  /**
       Note that the EBAMRNoSubcycleOld class defines its own grid hierarchy
       Initial and boundary conditions come through the EBIBCFactory.
       If fixedGrids != NULL, regridding gets turned off and the input
       grids are used
   */
  EBAMRNoSubcycleOld(const AMRParameters&             a_params,
                  const EBIBCFactory&              a_IBC,
                  const ProblemDomain&             a_coarsestDomain,
                  Real                             a_viscosity);

  ///
  void
  setSolverParams(int  a_numSmooth,
                  int  a_itermax,
                  int  a_mgcycle,
                  Real a_hang,
                  Real a_tolerance);

  void
  computeAdvectionVelocities(Vector<LevelData<EBFluxFAB> *>&                       a_advVel,
                             Vector<LayoutData< Vector< BaseIVFAB<Real> * > >* >&  a_coveredAdvVelLo,
                             Vector<LayoutData< Vector< BaseIVFAB<Real> * > >* >&  a_coveredAdvVelHi,
                             const Vector<LevelData<EBCellFAB> *>&                 a_veloOld);

  void
  computeAdvectiveDerivative(Vector<LevelData<EBCellFAB>* >                     &  a_uDotDelS,
                             Vector<LevelData<EBFluxFAB>* >                     &  a_macAdvVel,
                             Vector<LevelData<EBFluxFAB>* >                     &  a_macScalar,
                             Vector<LayoutData< Vector< BaseIVFAB<Real> * > >* >&  a_coveredAdvVelLo,
                             Vector<LayoutData< Vector< BaseIVFAB<Real> * > >* >&  a_coveredAdvVelHi,
                             Vector<LayoutData< Vector< BaseIVFAB<Real> * > >* >&  a_coveredScalarLo,
                             Vector<LayoutData< Vector< BaseIVFAB<Real> * > >* >&  a_coveredScalarHi,
                             bool a_nonConsOnly = false,
                             bool a_consOnly    = false);
  //fixes timestep to input
  void useFixedDt(Real a_dt);

  /// destructor.
  ~EBAMRNoSubcycleOld();

  /// set up grids, etc for fixed-hiearchy run
  void setupForFixedHierarchyRun(const Vector<Vector<Box> >& a_vectBoxes);

  /// set up grids, etc for AMR run
  void setupForAMRRun();

  /// set up for restart
  void setupForRestart(const string& a_restartFile);

  /// advance solution.  returns time of final solution
  Real run(Real a_maxTime, int a_maxStep);

  ///for convergence tests
  Vector<DisjointBoxLayout>& getGrids()
  {
    return m_grids;
  }
  ///for convergence tests
  Vector<EBISLayout>& getEBISLayouts()
  {
    return m_ebisl;
  }

  ///for convergence tests
  Vector<LevelData<EBCellFAB>* >& getVeloNew()
  {
    return m_veloNew;
  }
  ///for convergence tests
  Vector<LevelData<EBCellFAB>* >& getGphiNew()
  {
    return m_gphiNew;
  }

  ///for convergence tests
  Vector<LevelData<EBCellFAB>* >& getVeloOld()
  {
    return m_veloOld;
  }
  ///for convergence tests
  Vector<LevelData<EBCellFAB>* >& getGphiOld()
  {
    return m_gphiOld;
  }
  ///for convergence tests
  Vector<LevelData<EBCellFAB>* >& getPresOld()
  {
    return m_presOld;
  }
  ///for convergence tests
  Vector<LevelData<EBCellFAB>* >& getPresNew()
  {
    return m_presNew;
  }
  ///for convergence tests
  Real getDt()
  {
    return m_dt;
  }


  Vector< LayoutData< Vector< Vector<VolIndex> > >* >&  getCoveredFaceLo() {return m_coveredFaceLitLo;}
  Vector< LayoutData< Vector< Vector<VolIndex> > >* >&  getCoveredFaceHi() {return m_coveredFaceLitHi;}
  Vector< LayoutData< Vector< IntVectSet > >* >      &  getCoveredSetsLo() {return m_coveredSetsLitLo;}
  Vector< LayoutData< Vector< IntVectSet > >* >      &  getCoveredSetsHi() {return m_coveredSetsLitHi;}

  Vector<Real>& getDx()     {return m_dx;}
  Vector<ProblemDomain> & getDomain() {return m_domain;}
  EBCompositeMACProjector&  getMACProj() {return *m_macProjector;}
  EBCompositeCCProjector&    getCCProj() {return *m_ccProjector;}
protected:


  /// write final checkpoint and plot files
  void conclude();

  /// accumulate number of points updated
  void pointsUpdated();


#ifdef CH_USE_HDF5
  ///
  void writeCheckpointFile();

  ///
  void readCheckpointFile(const string& a_restartFile);

  ///
  void writePlotFile();

#endif

  void filter( Vector<LevelData<EBCellFAB>* >&   a_veloc);
  void advance();
  void advanceVelocity();
  void computeAdvectiveDerivative(Vector<LevelData<EBCellFAB>* >& a_uGradS,
                                  Vector<LevelData<EBCellFAB>* >& a_scalOld,
                                  bool                            a_reallyVelocity);


  /// things to do after a timestep
  void postTimeStep();

  /// create tags
  void tagCells(Vector<IntVectSet>& a_tags);

  /// create tags
  void tagCellsLevel(IntVectSet& a_tags, int a_level);

  /// regrid
  void regrid();

  ///
  /**
     Perform any post-regridding ops.   This includes the final smoothing and
     projecting the velocity and reinitializing the pressure.
   */
  void postRegrid();

  ///
  /**
     Perform any pre-regrdding ops.   This is vel = (I-mu*lapl)vel.
   */
  void preRegrid();

  ///
  void initialGrid(const Vector<DisjointBoxLayout>& a_new_grids);

  /// initialize data
  void initialData();

  /// things do do after initialization
  void postInitialize();

  void defineGrids(const Vector<Vector<Box> >& a_vectBoxes);
  void defineEBISLs();
  void defineNew();
  void defineOld();
  void defineProjections();
  void allocateDataHolders();
  /// compute dt
  Real computeDt();

  /// compute dt with initial data
  Real computeInitialDt();

  // compute vorticity on level a_level
  void computeVorticity(LevelData<EBCellFAB>& a_vorticity, int a_level) ;


  /// set up the AMR hierarchy, etc for a_vectGrids.
  void setup();

  /// regrid  once new boxes have been generated
  void regrid(const Vector<Vector<Box> >& a_new_grids);

  /// move new-time state into old-time state
  void swapOldAndNewStates();

  void initialGrid(const Vector<Vector<Box> >& a_vectBoxes);

public:
  void
  extrapolateScalarCol(Vector<LevelData<EBFluxFAB>* >                     &  a_macScalar,
                       Vector<LayoutData< Vector< BaseIVFAB<Real> * > >* >&  a_coveredMacLo,
                       Vector<LayoutData< Vector< BaseIVFAB<Real> * > >* >&  a_coveredMacHi,
                       const RefCountedPtr<EBPhysIBCFactory>              &  a_advectBC,
                       const Vector<LevelData<EBFluxFAB>* >               &  a_advectiveVel,
                       const Vector<LevelData<EBCellFAB>*>*                  a_sourceTerm,
                       const Vector<LevelData<EBCellFAB>* >               &  a_cellScalar,
                       const Vector<LevelData<EBCellFAB>* >               &  a_cellVelocity);

  //fills temporary stuff.  need to do this for some convergence tests.  do not do this unless you are sure
  //you know the context
  void allocateTemporaries();
  //delete temporary stuff.  need to do this for some convergence tests.  do not do this unless you are sure
  //you know the context
  void deleteTemporaries();
protected:
  void copyNewToOld(int a_ilev);
  void copyOldToNew(int a_ilev);

  void applyEBAMROp(Vector<LevelData<EBCellFAB>* >&       a_lap,
                    Vector<LevelData<EBCellFAB>* >&       a_phi,
                    Vector<LevelData<EBCellFAB>* >&       a_zero,
                    int                                   a_velComp);

  void
  viscousSourceForAdvect(Vector<LevelData<EBCellFAB>* >&       a_source,
                         Vector<LevelData<EBCellFAB>* >&       a_vel,
                         Vector<LevelData<EBCellFAB>* >&       a_zero,
                         int                                   a_velComp);
  void
  extrapolateToCoveredFaces(Vector<LayoutData< Vector< BaseIVFAB<Real> * > >* >&  a_coveredMacLo,
                            Vector<LayoutData< Vector< BaseIVFAB<Real> * > >* >&  a_coveredMacHi,
                            Vector<LevelData<EBFluxFAB>* >&                       a_macOpen,
                            Vector<LevelData<EBCellFAB>* >&                       a_cellOpen,
                            int                                                   a_idir);

  void averageDown(Vector<LevelData<EBCellFAB>* >&  a_data);
  void averageDown(Vector<LevelData<EBFluxFAB>* >&  a_data);

  //parmaters
  AMRParameters          m_params;
  EBIBC*                 m_ibc;

  //grid description
  Vector<ProblemDomain>     m_domain;
  Vector<Real>              m_dx;
  Vector<DisjointBoxLayout> m_grids;
  Vector<EBISLayout>        m_ebisl;
  Vector<EBLevelGrid>       m_eblg;

  //the actual data
  Vector<LevelData<EBCellFAB>* > m_veloOld;
  Vector<LevelData<EBCellFAB>* > m_veloNew;
  Vector<LevelData<EBCellFAB>* > m_presOld;
  Vector<LevelData<EBCellFAB>* > m_presNew;
  Vector<LevelData<EBCellFAB>* > m_gphiOld;
  Vector<LevelData<EBCellFAB>* > m_gphiNew;

  EBCompositeCCProjector*        m_ccProjector;
  EBCompositeMACProjector*       m_macProjector;
  RefCountedPtr<AMRMultiGrid<LevelData<EBCellFAB> > >  m_solver[SpaceDim];
  RefCountedPtr<AMRTGA<LevelData<EBCellFAB> > >  m_tgaSolver[SpaceDim];
  BiCGStabSolver<LevelData<EBCellFAB> >          m_bottomSolver;
  EBSimpleSolver                                 m_bottomSolverSimp;
  Vector<RefCountedPtr<EBQuadCFInterp> >    m_quadCFI;
  Vector<RefCountedPtr<EBCoarseAverage> >   m_aveOper;
  Vector<RefCountedPtr<EBCoarseAverage> >   m_aveSpac;

  //  finest extant level
  int m_finestLevel;
  int m_curStep;

  //  time stepping information
  Real m_time;
  Real m_dt;
  Real m_prescribedDt;
  Real m_viscosity;

  long long m_pointsUpdated;

  //switches that i would rather lose but cannot
  bool m_useFixedDt;

  bool m_isSetup;
  //if viscosity == 0 this turns off viscous solvers, etc
  bool m_viscousCalc;

  bool m_doRestart;

  bool m_advanceGphiOnly;


  //The following is the ugly part of this stuff.
  //Covered data (advection velocities, etc) are defined with this stuff.
  //These are used to define covered data and to iterate over it.
  //The actual data is defined on the fly to avoid confusion (no really)
  //and to save on memory since all these things are transitory.
  Vector< LayoutData< Vector< Vector<VolIndex> > >* >  m_coveredFaceLitLo;
  Vector< LayoutData< Vector< Vector<VolIndex> > >* >  m_coveredFaceLitHi;
  Vector< LayoutData< Vector< IntVectSet > >* >        m_coveredSetsLitLo;
  Vector< LayoutData< Vector< IntVectSet > >* >        m_coveredSetsLitHi;

  //defines the above stuff
  void defineIrregularData();

  //temporaries.  only live during advance.  it makes memory management easier.  really.
  // compute advection velocities
  Vector<LevelData<EBCellFAB>* >                       m_uDotDelU;
  Vector<LevelData<EBCellFAB>* >                       m_uDotDelS;
  Vector<LevelData<EBCellFAB>* >                       m_cellScratch;
  Vector<LevelData<EBCellFAB>* >                       m_cellScratc1;
  Vector<LevelData<EBCellFAB>* >                       m_cellScratc2;
  Vector<LevelData<EBCellFAB>* >                       m_cellScratc3;

  Vector<LevelData<EBFluxFAB>* >                       m_advVel;
  Vector<LevelData<EBFluxFAB>* >                       m_macGradient;
  Vector<LevelData<EBFluxFAB>* >                       m_macScratch1;
  Vector<LevelData<EBFluxFAB>* >                       m_macScratch2;

  Vector<LayoutData< Vector< BaseIVFAB<Real> * > >* >  m_coveredAdvVelLo;
  Vector<LayoutData< Vector< BaseIVFAB<Real> * > >* >  m_coveredAdvVelHi;
  Vector<LayoutData< Vector< BaseIVFAB<Real> * > >* >  m_coveredScratchLo;
  Vector<LayoutData< Vector< BaseIVFAB<Real> * > >* >  m_coveredScratchHi;

  Vector<RefCountedPtr<EBLevelAdvect > > m_ebLevAd;
private:

  // weak construction is bad
  EBAMRNoSubcycleOld()
  {
    MayDay::Error("invalid constructor");
  }

  //we disallow copy construction and assignment of objects
  //which contain pointered data
  EBAMRNoSubcycleOld(const EBAMRNoSubcycleOld& a_input)
  {
    MayDay::Error("invalid constructor");
  }

  //we disallow copy construction and assignment of objects
  //which contain pointered data
  void operator=(const EBAMRNoSubcycleOld& a_input)
  {
    MayDay::Error("invalid constructor");
  }

};

#endif

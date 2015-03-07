#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include <cmath>

#include "parstream.H"
#include "ParmParse.H"
#include "PolyGeom.H"

#include "DebugOut.H"
#include "Box.H"
#include "Vector.H"
#include "IntVectSet.H"
#include "EBCellFAB.H"
#include "DisjointBoxLayout.H"
#include "LayoutIterator.H"
#include "CH_HDF5.H"
#include "SPMD.H"
#include "EBLoadBalance.H"
#include "ProblemDomain.H"
#include "BoxIterator.H"
#include "AMRMultiGrid.H"
#include "EBAMRIO.H"
#include "BaseIVFactory.H"
#include "EBViscousTensorOpFactory.H"
#include "EBConductivityOpFactory.H"

#include "AMRLevel.H"
#include "EBAMRCNS.H"
#include "EBCellFactory.H"
#include "BaseIVFactory.H"
#include "VoFIterator.H"
#include "EBPatchGodunovF_F.H"
#include "EBPatchPolytropicF_F.H"
#include "EBIndexSpace.H"
#include "EBArith.H"
#include "EBLGIntegrator.H"
#include "EBLevelDataOps.H"
#include "NamespaceHeader.H"

int                                                      EBAMRCNS::s_NewPlotFile  = 0;
RefCountedPtr<AMRLevelOpFactory<LevelData<EBCellFAB> > > EBAMRCNS::s_tempFactory  =  RefCountedPtr<AMRLevelOpFactory<LevelData<EBCellFAB> > >();
RefCountedPtr<AMRLevelOpFactory<LevelData<EBCellFAB> > > EBAMRCNS::s_veloFactory  =  RefCountedPtr<AMRLevelOpFactory<LevelData<EBCellFAB> > >();
RefCountedPtr<EBLevelTGA >                               EBAMRCNS::s_veloLevTGA   =  RefCountedPtr<EBLevelTGA >();
RefCountedPtr<EBLevelTGA >                               EBAMRCNS::s_tempLevTGA   =  RefCountedPtr<EBLevelTGA >();
RefCountedPtr<AMRMultiGrid<     LevelData<EBCellFAB> > > EBAMRCNS::s_tempSolver   =  RefCountedPtr<AMRMultiGrid<     LevelData<EBCellFAB> > >();
RefCountedPtr<AMRMultiGrid<     LevelData<EBCellFAB> > > EBAMRCNS::s_veloSolver   =  RefCountedPtr<AMRMultiGrid<     LevelData<EBCellFAB> > >();
BiCGStabSolver<LevelData<EBCellFAB> >                    EBAMRCNS::s_botSolver    = BiCGStabSolver<LevelData<EBCellFAB> >();
void
EBAMRCNS::
fillCoefficients()
{
  EBCellFactory fact(m_eblg.getEBISL());
  LevelData<EBCellFAB> halfSt(m_eblg.getDBL(),m_nComp, 4*IntVect::Unit, fact);
  getHalfState(halfSt);
  for (DataIterator dit = m_eblg.getDBL().dataIterator(); dit.ok(); ++dit)
    {
      (*m_eta     )    [dit()].setVal(m_params.m_viscosityMu);
      (*m_bcoTemp )    [dit()].setVal(m_params.m_thermalCond);
      (*m_etaIrreg)    [dit()].setVal(m_params.m_viscosityMu);
      (*m_lambda     ) [dit()].setVal(m_params.m_viscosityLa);
      (*m_lambdaIrreg) [dit()].setVal(m_params.m_viscosityLa);
      (*m_bcoTempIrreg)[dit()].setVal(m_params.m_thermalCond);
      Interval srcInt(CRHO, CRHO);
      Interval dstInt(0, 0);
      const Box& region = m_eblg.getDBL().get(dit());
      (*m_acoVelo)[dit()].copy(region, dstInt, region, halfSt[dit()], srcInt);
      (*m_acoTemp)[dit()].copy(region, dstInt, region, halfSt[dit()], srcInt);
      (*m_acoTemp)[dit()] *= m_params.m_specHeatCv;
    }
}
void
EBAMRCNS::
defineSolvers()
{
  if(m_params.m_doDiffusion)
    {
      Vector<AMRLevel*> hierarchy = AMRLevel::getAMRLevelHierarchy();
      int nlevels = hierarchy.size();
      Vector<int>                                           refRat(nlevels);
      Vector<EBLevelGrid>                                   eblgs(nlevels);
      Vector<DisjointBoxLayout>                             grids(nlevels);
      Vector<RefCountedPtr<LevelData<EBCellFAB> >        >  aVelo(nlevels);
      Vector<RefCountedPtr<LevelData<EBFluxFAB> >        >  eta(nlevels);
      Vector<RefCountedPtr<LevelData<EBCellFAB> >        >  acoTemp(nlevels);
      Vector<RefCountedPtr<LevelData<EBFluxFAB> >        >  bcoTemp(nlevels);
      Vector<RefCountedPtr<LevelData<EBFluxFAB> >        >  lambda(nlevels);
      Vector<RefCountedPtr<LevelData<BaseIVFAB<Real> > > >  etaIrreg(nlevels);
      Vector<RefCountedPtr<LevelData<BaseIVFAB<Real> > > >  bcoTempIrreg(nlevels);
      Vector<RefCountedPtr<LevelData<BaseIVFAB<Real> > > >  lambdaIrreg(nlevels);
      Vector<RefCountedPtr<EBQuadCFInterp> >                quadCFI(nlevels);

      EBAMRCNS* coarsestLevel = (EBAMRCNS*)(hierarchy[0]);
      Real           lev0Dx      = (coarsestLevel->m_dx[0]);
      ProblemDomain lev0Dom      =  coarsestLevel->m_eblg.getDomain();


      for(int ilev = 0; ilev < nlevels; ilev++)
        {
          EBAMRCNS* cnsLevel = (EBAMRCNS*)(hierarchy[ilev]);
          cnsLevel->fillCoefficients();

          eblgs       [ilev] = cnsLevel->m_eblg;
          grids       [ilev] = cnsLevel->m_eblg.getDBL();
          refRat      [ilev] = cnsLevel->m_ref_ratio;
          aVelo       [ilev] = cnsLevel->m_acoVelo;
          acoTemp     [ilev] = cnsLevel->m_acoTemp;
          eta         [ilev] = cnsLevel->m_eta;
          etaIrreg    [ilev] = cnsLevel->m_etaIrreg;
          lambda      [ilev] = cnsLevel->m_lambda;
          lambdaIrreg [ilev] = cnsLevel->m_lambdaIrreg;
          quadCFI     [ilev] = cnsLevel->m_quadCFI;
          bcoTemp     [ilev] = cnsLevel->m_bcoTemp;
          bcoTempIrreg[ilev] = cnsLevel->m_bcoTempIrreg;
        }

      //alpha, beta get replaced in tga solves
      Real alpha = 1;
      Real beta = 1;
      IntVect  giv    = 4*IntVect::Unit;
      RealVect origin =  RealVect::Zero;
      int numSmooth, numMG, maxIter, mgverb;
      Real tolerance, hang, normThresh;
      ParmParse pp("amrmultigrid");
      pp.get("num_smooth", numSmooth);
      pp.get("num_mg",     numMG);
      pp.get("hang_eps",   hang);
      pp.get("norm_thresh",normThresh);
      pp.get("tolerance",  tolerance);
      pp.get("max_iter",   maxIter);
      pp.get("verbosity",  mgverb);

      s_veloFactory =
        RefCountedPtr<AMRLevelOpFactory<LevelData<EBCellFAB> > >
        ((AMRLevelOpFactory<LevelData<EBCellFAB> >*)
         (new EBViscousTensorOpFactory(eblgs, alpha, beta, aVelo, eta,
                                       lambda, etaIrreg, lambdaIrreg,lev0Dx, refRat,
                                       m_params.m_doBCVelo, m_params.m_ebBCVelo,giv, giv)));
      s_tempFactory =
        RefCountedPtr<AMRLevelOpFactory<LevelData<EBCellFAB> > >
        ((AMRLevelOpFactory<LevelData<EBCellFAB> >*)
         (new EBConductivityOpFactory(eblgs, quadCFI,
                                      alpha, beta, acoTemp, bcoTemp, bcoTempIrreg,
                                      lev0Dx, refRat, m_params.m_doBCTemp, m_params.m_ebBCTemp,
                                      giv, giv)));


      s_veloSolver = RefCountedPtr<AMRMultiGrid< LevelData<EBCellFAB> > >( new AMRMultiGrid< LevelData<EBCellFAB> > ());
      s_tempSolver = RefCountedPtr<AMRMultiGrid< LevelData<EBCellFAB> > >( new AMRMultiGrid< LevelData<EBCellFAB> > ());
  
      s_veloSolver->define(lev0Dom, *s_veloFactory, &s_botSolver, nlevels);
      s_tempSolver->define(lev0Dom, *s_tempFactory, &s_botSolver, nlevels);


      s_veloSolver->setSolverParameters(numSmooth, numSmooth, numSmooth, numMG, maxIter, tolerance, hang, normThresh);
      s_tempSolver->setSolverParameters(numSmooth, numSmooth, numSmooth, numMG, maxIter, tolerance, hang, normThresh);
      s_veloSolver->m_verbosity = mgverb;
      s_tempSolver->m_verbosity = mgverb;
      s_botSolver.m_verbosity   = mgverb-3;
      //  s_botSolver.m_numRestarts = 0;

      s_tempLevTGA = RefCountedPtr<EBLevelTGA>(new EBLevelTGA(grids, refRat, lev0Dom, s_tempFactory, s_tempSolver));
      s_veloLevTGA = RefCountedPtr<EBLevelTGA>(new EBLevelTGA(grids, refRat, lev0Dom, s_veloFactory, s_veloSolver));
      s_tempLevTGA->setEBLG(eblgs);
      s_veloLevTGA->setEBLG(eblgs);
    }
}

/***************************/
void
EBAMRCNS::
define(AMRLevel*            a_coarser_level_ptr,
       const ProblemDomain& a_problem_domain,
       int                  a_level,
       int                  a_ref_ratio)
{
  if(m_params.m_verbosity >= 3)
    {
      pout() << "EBAMRCNS::define, level=" << a_level << endl;
    }

  EBPatchGodunov::useConservativeSource(true);
  m_isDefined = true;
  AMRLevel::define(a_coarser_level_ptr,
                   a_problem_domain,
                   a_level,
                   a_ref_ratio);

  if (a_coarser_level_ptr != NULL)
    {
      EBAMRCNS* amrg_ptr =
        dynamic_cast<EBAMRCNS*>(a_coarser_level_ptr);
      if (amrg_ptr == NULL)
        {
          pout() << "EBAMRG::define:cannot cast  to derived class"
                 << endl;
          MayDay::Error();
        }
      m_params = amrg_ptr->m_params;
    }
  Real dxSize  = m_params.m_domainLength/m_problem_domain.domainBox().size(0);
  m_dx = dxSize*RealVect::Unit;
  m_nGhost = 4;

  if(m_ebPatchGodunov != NULL)
    delete m_ebPatchGodunov;

  m_ebPatchGodunov = RefCountedPtr<EBPatchGodunov>(m_ebPatchGodunovFactory->create());
  m_ebPatchGodunov->define(m_problem_domain, m_dx);

  m_nComp      = m_ebPatchGodunov->numConserved();
  m_nPrim      = m_ebPatchGodunov->numPrimitives();
  m_stateNames = m_ebPatchGodunov->stateNames();
  m_primNames  = m_ebPatchGodunov->primNames();
  m_ref_ratio  = a_ref_ratio;
}
/***************************/
Real
EBAMRCNS::
advance()
{
  EBPatchGodunov::s_whichLev = m_level;
  if(m_params.m_verbosity >= 3)
    {
      pout() << " in EBAMRCNS advance for level " << m_level << endl;
    }
  m_stateNew.copyTo(m_stateNew.interval(),
                    m_stateOld,
                    m_stateOld.interval());


  EBCellFactory fact(m_eblg.getEBISL());
  LevelData<EBCellFAB> halfStCo;
  LevelData<EBCellFAB>  divergeF(m_eblg.getDBL(),m_nComp, 4*IntVect::Unit, fact);
  fluxDivergence(       divergeF, halfStCo);
  advanceDensity(       divergeF);

  if(m_params.m_doDiffusion)//otherwise, everything gets done in advance density
    {
      //now can define the solvers with half state density
      defineSolvers();
      advanceVelocity(      divergeF);
      advanceTemperature(   divergeF, halfStCo);
    }

  //deal with time and time step
  Real maxWaveSpeed = m_ebLevelGodunov.getMaxWaveSpeed(m_stateOld);
  Real new_dt = m_params.m_cfl*m_dx[0]/maxWaveSpeed;
  m_time += m_dt;

  //save stable timestep to save computational effort
  m_dtNew = new_dt;

  return new_dt;
}
/***************************/
void
EBAMRCNS::
getVelocity(LevelData<EBCellFAB>&       a_velocity,
            const LevelData<EBCellFAB>& a_state) const
{
  Interval velInterv(0, SpaceDim-1);
  Interval momInterv(CMOMX, CMOMX+SpaceDim-1);
  for(DataIterator dit = m_eblg.getDBL().dataIterator(); dit.ok(); ++dit)
    {
      const Box& region = m_eblg.getDBL().get(dit());
      a_velocity[dit()].copy(region, velInterv, region, a_state[dit()], momInterv);
      for(int idir = 0; idir < SpaceDim; idir++)
        {
          a_velocity[dit()].divide(a_state[dit()], CRHO, idir, 1);
        }
    }
}
/***************************/
void
EBAMRCNS::
getTemperature(LevelData<EBCellFAB>&       a_temperature,
               const LevelData<EBCellFAB>& a_state) const
{
  for(DataIterator dit = m_eblg.getDBL().dataIterator(); dit.ok(); ++dit)
    {
      const Box& region = m_eblg.getDBL().get(dit());
      BaseFab<Real>& regTemp = a_temperature[dit()].getSingleValuedFAB();
      const BaseFab<Real>& regStat = a_state[dit()].getSingleValuedFAB();
      FORT_GETTEMPERATURE(CHF_FRA1(regTemp,0),
                          CHF_CONST_FRA(regStat),
                          CHF_BOX(region),
                          CHF_CONST_REAL(m_params.m_specHeatCv));
      IntVectSet ivsMulti = m_eblg.getEBISL()[dit()].getMultiCells(region);
      for(VoFIterator vofit(ivsMulti, m_eblg.getEBISL()[dit()].getEBGraph()); vofit.ok(); ++vofit)
        {
          Vector<Real> cvec(CNUM);
          for(int ivar = 0; ivar< CNUM; ivar++)
            {
              cvec[ivar] = a_state[dit()](vofit(), ivar);
            }
          Real temper;
          FORT_POINTGETTEMPERATURE(CHF_REAL(temper),
                                   CHF_VR(cvec),
                                   CHF_CONST_REAL(m_params.m_specHeatCv));

          a_temperature[dit()](vofit(), 0)  = temper;
        }
    }
}
/***************************/
void
EBAMRCNS::
evalOperator(LevelData<EBCellFAB>&              a_source,
             const LevelData<EBCellFAB>&        a_primitive,
             const LevelData<EBCellFAB>&        a_primiCoar,
             TGAHelmOp<LevelData<EBCellFAB> >*  a_op,
             Real a_alpha, Real a_beta)
{
  LevelData<EBCellFAB> zero;
  a_op->create(zero, a_source);
  a_op->setToZero(zero);
  a_op->setAlphaAndBeta(a_alpha, a_beta);
  if(!m_hasCoarser)
    {
      a_op->residual(     a_source, a_primitive,              zero, false);
    }
  else
    {
      a_op->AMRResidualNF(a_source, a_primitive, a_primiCoar, zero, false);
    }
  a_op->scale(a_source, -1.0);
}
/***************************/
void
EBAMRCNS::
getHalfState(LevelData<EBCellFAB>& a_stateInt)
{
  //interpolate state to n+1/2
  Real told = 0; Real tnew = 1; Real time = 0.5;
  EBArith::timeInterpolate(a_stateInt, m_stateOld, m_stateNew,
                           m_eblg.getDBL(), time, told, tnew);
}
/***************************/
void
EBAMRCNS::
hyperbolicSource(LevelData<EBCellFAB>&       a_source,
                 const LevelData<EBCellFAB>& a_halfSt,
                 const LevelData<EBCellFAB>& a_halfCo)
{
  EBPatchGodunov::useConservativeSource(true);
  EBLevelDataOps::setVal(a_source, 0.0);
  if(m_params.m_doDiffusion)
    {
      EBCellFactory fact(m_eblg.getEBISL());
      LevelData<EBCellFAB>  momeSour(m_eblg.getDBL(), SpaceDim,4*IntVect::Unit, fact);
      LevelData<EBCellFAB>  enerSour(m_eblg.getDBL(),        1,4*IntVect::Unit, fact);
      LevelData<EBCellFAB>  velocity(m_eblg.getDBL(), SpaceDim,4*IntVect::Unit, fact);
      LevelData<EBCellFAB>  temperat(m_eblg.getDBL(),        1,4*IntVect::Unit, fact);
      LevelData<EBCellFAB>  veloCoar, tempCoar;
      getVelocity(   velocity,  a_halfSt);
      getTemperature(temperat,  a_halfSt);
      if(m_hasCoarser)
        {
          EBAMRCNS* coarPtr = getCoarserLevel();
          EBCellFactory factCoar(coarPtr->m_eblg.getEBISL());
          veloCoar.define(coarPtr->m_eblg.getDBL(), SpaceDim,4*IntVect::Unit, factCoar);
          tempCoar.define(coarPtr->m_eblg.getDBL(),        1,4*IntVect::Unit, factCoar);
          getVelocity(   veloCoar, a_halfCo);
          getTemperature(tempCoar, a_halfCo);
        }

      TGAHelmOp<LevelData<EBCellFAB> >* ebvtop = (TGAHelmOp<LevelData<EBCellFAB> >*)s_veloFactory->AMRnewOp(m_eblg.getDomain());
      TGAHelmOp<LevelData<EBCellFAB> >* amrpop = (TGAHelmOp<LevelData<EBCellFAB> >*)s_tempFactory->AMRnewOp(m_eblg.getDomain());

      evalOperator(momeSour, velocity, veloCoar, ebvtop, 0.0, 1.0);
      evalOperator(enerSour, temperat, tempCoar, amrpop, 0.0, 1.0);

      //AMRNewOp is a new so we need to clean up.
      delete ebvtop;
      delete amrpop;

      Interval momeInterv(CMOMX, CMOMX+SpaceDim-1);
      Interval enerInterv(CENG,  CENG);
      Interval vectInterv(0, SpaceDim-1);
      Interval scalInterv(0, 0);
      for(DataIterator dit = a_source.dataIterator(); dit.ok(); ++dit)
        {
          const Box& region = m_eblg.getDBL().get(dit());
          a_source[dit()].copy(region, momeInterv, region, momeSour[dit()], vectInterv);
          a_source[dit()].copy(region, enerInterv, region, enerSour[dit()], scalInterv);
        }
    }
}
/***************************/
void
EBAMRCNS::
fluxDivergence( LevelData<EBCellFAB>& a_divergeF,
                LevelData<EBCellFAB>& a_halfStCo)
{
  IntVect ivGhost = 4*IntVect::Unit;
  EBCellFactory fact(m_eblg.getEBISL());
  LevelData<EBCellFAB> source(m_eblg.getDBL(),m_nComp, ivGhost, fact);
  LevelData<EBCellFAB> halfSt(m_eblg.getDBL(),m_nComp, ivGhost, fact);
  //set up arguments to step
  //undefined lfr in case we need it
  EBFluxRegister lfr;
  //undefined leveldata in case we need it
  const LevelData<EBCellFAB> ld;
  //set arguments to dummy arguments and
  //then fix if they are available
  EBFluxRegister* coarFR = &lfr;
  EBFluxRegister* fineFR = &lfr;
  const LevelData<EBCellFAB>* coarDataOld = &ld;
  const LevelData<EBCellFAB>* coarDataNew = &ld;

  Real told = 0.0;
  Real tnew = 0.0;
  if(m_hasCoarser)
    {
      EBAMRCNS* coarPtr = getCoarserLevel();
      EBCellFactory factCoar(coarPtr->m_eblg.getEBISL());
      a_halfStCo.define(coarPtr->m_eblg.getDBL(), m_nComp, ivGhost, factCoar);
      //recall that my flux register goes between my
      //level and the next finer level
      coarFR = &coarPtr->m_divFFluxRegister;
      coarDataOld = &coarPtr->m_stateOld;
      coarDataNew = &coarPtr->m_stateNew;
      tnew = coarPtr->m_time;
      told = tnew - coarPtr->m_dt;
      //time should never be greater than the newest coarse
      //time.  time might be very slightly smaller than
      //told because of the above subtraction.
      Real eps = 1.0e-10;
      if((m_time > tnew) || (m_time < (told - eps)))
        {
          MayDay::Error("out of bounds time input to AMRCNS");
        }
      //correct for said floating-point nastiness
      m_time = Max(m_time, told);
      EBArith::timeInterpolate(a_halfStCo, *coarDataOld, *coarDataNew,
                               coarPtr->m_eblg.getDBL(), m_time, told, tnew);
    }
  if(m_hasFiner)
    {
      //recall that my flux register goes between my
      //level and the next finer level
      fineFR = &m_divFFluxRegister;
    }

#ifndef NDEBUG
  if(!m_hasCoarser && (m_params.m_verbosity > 1))
    {
      Real summass;
      int densityIndex = m_ebPatchGodunov->densityIndex();
      sumConserved(summass,  densityIndex);

      pout() << "sum mass = " << summass << endl;
    }
#endif

  hyperbolicSource(source,  m_stateOld, a_halfStCo);
  // flux register manipulation happens in levelgodunov
  //the rhs for the solves for temperature and velocity 
  //also come out of here since we need the half time
  //values  of the primitive state.    Not the most intuitive
  //interface, I know.
  m_ebLevelGodunov.divergeF(a_divergeF,
                            m_rhsVelo,
                            m_rhsTemp,
                            m_massDiff,
                            *fineFR,
                            *coarFR,
                            m_stateOld,
                            source,
                            *coarDataOld,
                            *coarDataNew,
                            m_time,
                            told,
                            tnew,
                            m_dt);

  coarseFineIncrement();
}
void
EBAMRCNS::
coarseFineIncrement()
{
  Interval interv(0, m_nComp-1);
  // increment redistribution register between this level
  // and next coarser level
  {
    CH_TIME("coarse-fine guano");
    if(m_hasCoarser)
      {
        for(DataIterator dit = m_eblg.getDBL().dataIterator(); dit.ok(); ++dit)
          {
            m_ebFineToCoarRedist.increment(m_massDiff[dit()], dit(), interv);
          }
      }

    //initialize redistribution register between this level
    // and next finer level
    //this includes re-redistribution registers
    if(m_hasFiner)
      {

        m_ebCoarToFineRedist.setToZero();
        m_ebCoarToCoarRedist.setToZero();
        for(DataIterator dit = m_eblg.getDBL().dataIterator(); dit.ok(); ++dit)
          {
            BaseIVFAB<Real>& massDiffFAB = m_massDiff[dit()];
            m_ebCoarToFineRedist.increment(massDiffFAB, dit(), interv);
            m_ebCoarToCoarRedist.increment(massDiffFAB, dit(), interv);
          }
      }
  }
}
/***************************/
void
EBAMRCNS::
redistributeMomentum()
{
  //first redistribute into the rhs for the solve
  Interval srcInterv(CMOMX, CMOMX+SpaceDim-1);
  Interval dstInterv(0, SpaceDim-1);
  m_ebLevelRedist.redistribute(m_rhsVelo, srcInterv, dstInterv);
}
/***************************/
void
EBAMRCNS::
redistributeEnergy()
{
  //first redistribute into the rhs for the solve
  Interval srcInterv(CENG, CENG);
  Interval dstInterv(0, 0);
  m_ebLevelRedist.redistribute(m_rhsTemp, srcInterv, dstInterv);
}
/***************************/
void
EBAMRCNS::
redistributeDensity()
{
  Interval consInterv(0, m_nComp-1);
  m_ebLevelRedist.setToZero();
  if(m_params.m_useMassRedist)
    {
      //if use mass weighting, need to
      //fix weights of redistribution object
      int densityIndex = m_ebPatchGodunov->densityIndex();
      m_stateNew.exchange(Interval(0, m_nComp-1));
      m_ebLevelRedist.resetWeights(m_stateNew, densityIndex);
    }
  for (DataIterator dit = m_eblg.getDBL().dataIterator(); dit.ok(); ++dit)
    {
      m_ebLevelRedist.increment(m_massDiff[dit()], dit(), consInterv);
    }

  if(m_params.m_doDiffusion)
    {
      //first redistribute density explictly
      Interval densInterv(CRHO, CRHO);
      m_ebLevelRedist.redistribute(m_stateNew, densInterv);
    }
  else
    {
      //no diffusion means that we can do everything explicitly
      m_ebLevelRedist.redistribute(m_stateNew, consInterv);
    }
}
/***************************/
void
EBAMRCNS::
advanceDensity(const LevelData<EBCellFAB>& a_divergeF)
{
  if(m_params.m_doDiffusion)
    {
      int bigStart = CRHO;
      int litStart = 0;
      int numComp  = 1;
      int srcComp, dstComp;
      //thankfully, we have something that can be advanced explicitly
      for (DataIterator dit = m_eblg.getDBL().dataIterator(); dit.ok(); ++dit)
        {

          //begin debug
          //if((m_level == 1)  && m_eblg.getDBL().get(dit()).contains(EBPatchGodunov::s_debugIV))
          //   {
          //     pout() << "here" << endl;
          //   }
          //end debug
          EBCellFAB dtDivergeF(m_eblg.getEBISL()[dit()], m_eblg.getDBL().get(dit()), 1);
          dtDivergeF.setVal(0.);
          srcComp = bigStart; dstComp = litStart;
          dtDivergeF.plus(a_divergeF[dit()], srcComp, dstComp, numComp);
          dtDivergeF *= m_dt;
          srcComp = litStart; dstComp = bigStart;
          m_stateNew[dit()].minus(dtDivergeF, srcComp, dstComp, numComp);
        }
    }
  else
    {
      //advance everything explicitly
      for (DataIterator dit = m_eblg.getDBL().dataIterator(); dit.ok(); ++dit)
        {
          EBCellFAB dtDivergeF(m_eblg.getEBISL()[dit()], m_eblg.getDBL().get(dit()), a_divergeF.nComp());
          //begin debug
          //if((m_level == 1)  && m_eblg.getDBL().get(dit()).contains(EBPatchGodunov::s_debugIV))
          //   {
          //     pout() << "here" << endl;
          //   }
          //end debug
          dtDivergeF.setVal(0.);
          dtDivergeF += a_divergeF[dit()];
          dtDivergeF *= m_dt;
          m_stateNew[dit()] -= dtDivergeF;
        }
    }
  //this also resets the weights of the redistribution object in the case of 
  //mass weighting.  redistributes everything in the case of no diffusion
  redistributeDensity(); 
}

void
EBAMRCNS::
getTGADiffusion(LevelData<EBCellFAB>& a_diffusion,  
                LevelData<EBCellFAB>& a_phiold,
                LevelData<EBCellFAB>& a_rhs, 
                RefCountedPtr< EBLevelTGA> & a_levTGA,
                bool a_isVarVel)
{
  //use TGA to advance velocity
  EBFluxRegister*       coarPhiFRPtr = NULL;
  EBFluxRegister*       finePhiFRPtr = NULL;
  LevelData<EBCellFAB>* phiCoaOldPtr = NULL;
  LevelData<EBCellFAB>* phiCoaNewPtr = NULL;
  Real tCoarOld = 0;
  Real tCoarNew = 0;
  int ncomp = 1;
  if(a_isVarVel) ncomp = SpaceDim;
  if(m_hasCoarser)
    {
      EBAMRCNS* coarCNS = getCoarserLevel();
      if(a_isVarVel)
        coarPhiFRPtr = &coarCNS->m_veloFluxRegister;
      else
        coarPhiFRPtr = &coarCNS->m_tempFluxRegister;

      const  EBLevelGrid& ceblg = coarCNS->m_eblg;
      EBCellFactory cfact(ceblg.getEBISL());
      phiCoaNewPtr = new LevelData<EBCellFAB>(ceblg.getDBL(), ncomp, 4*IntVect::Unit, cfact);
      phiCoaOldPtr = new LevelData<EBCellFAB>(ceblg.getDBL(), ncomp, 4*IntVect::Unit, cfact);

      if(a_isVarVel)
        {
          coarCNS->getVelocity(*phiCoaOldPtr, coarCNS->m_stateOld);
          coarCNS->getVelocity(*phiCoaNewPtr, coarCNS->m_stateNew);
        }
      else
        {
          coarCNS->getTemperature(*phiCoaOldPtr, coarCNS->m_stateOld);
          coarCNS->getTemperature(*phiCoaNewPtr, coarCNS->m_stateNew);
        }
      tCoarNew = coarCNS->m_time;
      tCoarOld = tCoarNew - coarCNS->m_dt;
    }
  if(m_hasFiner)
    {
      if(a_isVarVel)
        finePhiFRPtr = &m_veloFluxRegister;
      else
        finePhiFRPtr = &m_tempFluxRegister;
    }
  //this should be compute diffusion
  a_levTGA->computeDiffusion(a_diffusion, a_phiold, a_rhs,
                             finePhiFRPtr, coarPhiFRPtr,
                             phiCoaOldPtr, phiCoaNewPtr,
                             m_time, tCoarOld, tCoarNew,
                             m_dt, m_level);

  if(m_hasCoarser)
    {
      delete phiCoaOldPtr;
      delete phiCoaNewPtr;
    }
}
/***************************/
void
EBAMRCNS::
advanceVelocity(const LevelData<EBCellFAB>& a_divergeF)
{
  if(m_params.m_doDiffusion) //otherwise, everything is done in advanceDensity
    {
      EBCellFactory fact(m_eblg.getEBISL());

      LevelData<EBCellFAB>  velOld(m_eblg.getDBL(), SpaceDim, 4*IntVect::Unit, fact);
      LevelData<EBCellFAB>  lapVel(m_eblg.getDBL(), SpaceDim, 4*IntVect::Unit, fact);
      getVelocity(   velOld, m_stateOld);

      //this adds redistribution mass to m_rhsVelo
      redistributeMomentum();
      getTGADiffusion(lapVel, velOld, m_rhsVelo, s_veloLevTGA, true);    // 

      int litStart = 0;
      int bigStart = CMOMX;
      int numComp = SpaceDim;
      int srcComp, dstComp;
      for (DataIterator dit = m_eblg.getDBL().dataIterator(); dit.ok(); ++dit)
        {
          //advance hyperbolic bit of solution
          EBCellFAB dtDivergeF(m_eblg.getEBISL()[dit()], m_eblg.getDBL().get(dit()), SpaceDim);
          dtDivergeF.setVal(0.);
          srcComp = bigStart; dstComp = litStart;
          dtDivergeF.plus(a_divergeF[dit()], srcComp, dstComp, numComp);
          dtDivergeF *= m_dt;

          srcComp = litStart; dstComp = bigStart;
          m_stateNew[dit()].minus(dtDivergeF, srcComp, dstComp, numComp);
          //now add in diffusion
          //first need to multiply by dt
          lapVel[dit()] *= m_dt;
          //lapvel now holds dt*lapl(vel).  add it in.
          m_stateNew[dit()].plus(lapVel[dit()], srcComp, dstComp, numComp);
        }
    }
}
/***/
void
EBAMRCNS::
addULaplUToRHSTemperature()
{
  EBCellFactory factoryNew(m_eblg.getEBISL());
  IntVect ivGhost = m_nGhost*IntVect::Unit;
  LevelData<EBCellFAB> stateInt(m_eblg.getDBL(),m_nComp , ivGhost, factoryNew);
  LevelData<EBCellFAB> velocity(m_eblg.getDBL(),SpaceDim, ivGhost, factoryNew);
  LevelData<EBCellFAB> laplVelo(m_eblg.getDBL(),SpaceDim, ivGhost, factoryNew);
  getHalfState(stateInt);
  getVelocity(velocity,  stateInt);


  LevelData<EBCellFAB>* veloCoaOldPtr = NULL;
  LevelData<EBCellFAB>* veloCoaNewPtr = NULL;
  LevelData<EBCellFAB>* veloCoar;
  Real tCoarOld, tCoarNew;
  //get coarse velocity at right time
  if(m_hasCoarser)
    {
      EBAMRCNS* coarCNS = getCoarserLevel();
      const  EBLevelGrid& ceblg = coarCNS->m_eblg;
      EBCellFactory cfact(ceblg.getEBISL());
      veloCoaNewPtr = new LevelData<EBCellFAB>(ceblg.getDBL(), SpaceDim, 4*IntVect::Unit, cfact);
      veloCoaOldPtr = new LevelData<EBCellFAB>(ceblg.getDBL(), SpaceDim, 4*IntVect::Unit, cfact);
      veloCoar      = new LevelData<EBCellFAB>(ceblg.getDBL(), SpaceDim, 4*IntVect::Unit, cfact);

      coarCNS->getVelocity(*veloCoaOldPtr, coarCNS->m_stateOld);
      coarCNS->getVelocity(*veloCoaNewPtr, coarCNS->m_stateNew);
      tCoarNew = coarCNS->m_time;
      tCoarOld = tCoarNew - coarCNS->m_dt;

      EBArith::timeInterpolate(*veloCoar, *veloCoaOldPtr, *veloCoaNewPtr,
                               ceblg.getDBL(), m_time, tCoarOld, tCoarNew);
    }

  //evaluate lapl vel
  TGAHelmOp<LevelData<EBCellFAB> >* ebvtop = (TGAHelmOp<LevelData<EBCellFAB> >*)s_veloFactory->AMRnewOp(m_eblg.getDomain());
  evalOperator(laplVelo, velocity, *veloCoar, ebvtop, 0.0, 1.0);

  //now need to take the dot product of velocity and lapl velocity and add it to the rhs
  LevelData<EBCellFAB> uLaplU(m_eblg.getDBL(),1, ivGhost, factoryNew);
  for(DataIterator dit = m_eblg.getDBL().dataIterator(); dit.ok(); ++dit)
    {
      uLaplU[dit()].setVal(0.);
      for(int idir = 0; idir < SpaceDim; idir++)
        {
          int isrc = idir; int idst = 0; int inc = 1;
          uLaplU[dit()].plus(velocity[dit()], isrc, idst, inc);
          uLaplU[dit()].mult(laplVelo[dit()], isrc, idst, inc);
        }
      m_rhsTemp[dit()] += uLaplU[dit()];
    }
  
  //clean up
  delete ebvtop; //AMRNewOp is a new so we need to clean it up.   I suppose I could use reclaim to further obfuscate.
  if(m_hasCoarser)
    {
      delete veloCoar;
      delete veloCoaNewPtr;
      delete veloCoaOldPtr;
    }
}
/**/
void
EBAMRCNS::
advanceTemperature(const LevelData<EBCellFAB>& a_divergeF,
                   const LevelData<EBCellFAB>& a_halfCo)
{
  if(m_params.m_doDiffusion) //otherwise, everything is done in advanceDensity
    {
      EBCellFactory fact(m_eblg.getEBISL());

      LevelData<EBCellFAB>  tempOld(m_eblg.getDBL(), SpaceDim, 4*IntVect::Unit, fact);
      LevelData<EBCellFAB>  lapTemp(m_eblg.getDBL(), SpaceDim, 4*IntVect::Unit, fact);
      getTemperature(   tempOld, m_stateOld);
      getTemperature(   lapTemp, m_stateNew);

      //the hyperbolics could not supply the u lapl u term to the rhs 
      //of temperature because it did not have the velocity at the right time
      //or centering.   so we have to add it now.
      addULaplUToRHSTemperature();

      //this adds redistribution mass to m_rhsTemp
      redistributeEnergy();

      getTGADiffusion(lapTemp, tempOld, m_rhsTemp, s_tempLevTGA, false);

      int litStart = 0;
      int bigStart = CENG;
      int numComp = 1;
      int srcComp, dstComp;
      for(DataIterator dit = m_eblg.getDBL().dataIterator(); dit.ok(); ++dit)
        {
          //advance hyperbolic bit of solution
          EBCellFAB dtDivergeF(m_eblg.getEBISL()[dit()], m_eblg.getDBL().get(dit()), 1);
          dtDivergeF.setVal(0.);
          srcComp = bigStart; dstComp = litStart;
          dtDivergeF.plus(a_divergeF[dit()], srcComp, dstComp, numComp);
          dtDivergeF *= m_dt;

          srcComp = litStart; dstComp = bigStart;
          m_stateNew[dit()].minus(dtDivergeF, srcComp, dstComp, numComp);
          //now add in diffusion
          //first need to multiply by dt
          lapTemp[dit()] *= m_dt;
          //lapTemp now holds dt*lapl(Temp).  add it in.
          m_stateNew[dit()].plus(lapTemp[dit()], srcComp, dstComp, numComp);
        }
    }
}

/***************************/
void
EBAMRCNS::
postTimeStep()
{
  if(m_params.m_verbosity >= 3)
    {
      pout() << " in EBAMRCNS postTimeStep for level " << m_level << endl;
    }
  if (m_hasFiner)
    {
      Interval interv(0, m_nComp-1);
      EBAMRCNS* finePtr = getFinerLevel();
      finePtr->m_ebCoarseAverage.average(m_stateNew,
                                         finePtr->m_stateNew,
                                         interv);
    }
  //this does the refluxing and redistribution evil dance
  postTimeStepRefluxRedistDance();
}
void
EBAMRCNS::
postTimeStepRefluxRedistDance()
{
  //this does the refluxing and redistribution evil dance
  Interval interv(0, m_nComp-1);
  if(m_params.m_doDiffusion)
    {
      MayDay::Error("this stuff needs to be implicit");
    }
  if(m_hasCoarser)
    {
      if(m_params.m_doSmushing)
        {
          //redistibute to coarser level
          EBAMRCNS* coarPtr = getCoarserLevel();
          //if use mass weighting, need to
          //fix weights of redistribution object
          if(m_params.m_useMassRedist)
            {
              int densevar = m_ebPatchGodunov->densityIndex();
              m_stateNew.exchange(Interval(0, m_nComp-1));
              m_ebFineToCoarRedist.resetWeights(coarPtr->m_stateNew, densevar);
            }
          m_ebFineToCoarRedist.redistribute(coarPtr->m_stateNew, interv);
        }
      m_ebFineToCoarRedist.setToZero();
    }
  if (m_hasFiner)
    {
      EBAMRCNS* finePtr = getFinerLevel();
      // reflux from finer level solution
      CH_assert(Abs(m_dx[0]-m_dx[1])<1.e-9);
      Real scale = -1.0/m_dx[0];
      m_divFFluxRegister.reflux(m_stateNew, interv, scale);
      //the flux register must modify the redistribution
      //registers
      m_divFFluxRegister.incrementRedistRegister(m_ebCoarToFineRedist,
                                                 interv, scale);

      m_divFFluxRegister.incrementRedistRegister(m_ebCoarToCoarRedist,
                                                 interv, scale);

      if(m_params.m_doSmushing)
        {
          //if use mass weighting, need to
          //fix weights of redistribution object
          if(m_params.m_useMassRedist)
            {
              int densevar = m_ebPatchGodunov->densityIndex();
              m_stateNew.exchange(Interval(0, m_nComp-1));
              m_ebCoarToFineRedist.resetWeights(m_stateNew, densevar);
              m_ebCoarToCoarRedist.resetWeights(m_stateNew, densevar);
            }
          //redistibute to finer level
          m_ebCoarToFineRedist.redistribute(finePtr->m_stateNew, interv);
          //do the re redistirubtion
          m_ebCoarToCoarRedist.redistribute(m_stateNew, interv);
        }

      m_ebCoarToFineRedist.setToZero();
      m_ebCoarToCoarRedist.setToZero();
    }
}
/****************************/
void
EBAMRCNS::
tagCells(IntVectSet& a_tags)
{
  if (m_params.m_verbosity >= 3)
    {
      pout() << "EBAMRCNS::tagCells for level " << m_level << endl;
    }

  // Create tags based on undivided gradient of density
  IntVectSet localTags;

  vector<Real> problo(SpaceDim);
  // If there is a coarser level interpolate undefined ghost cells
  //only interpolate the density
  int densityIndex = m_ebPatchGodunov->densityIndex();
  Interval intervDensity(densityIndex, densityIndex);
  EBCellFactory factory(m_eblg.getEBISL());
  int nCons = m_ebPatchGodunov->numConserved();
  LevelData<EBCellFAB> consTemp(m_eblg.getDBL(), nCons, IntVect::Unit, factory);
  Interval consInterv(0, nCons-1);
  m_stateNew.copyTo(consInterv, consTemp, consInterv);
  if (m_hasCoarser)
    {
      const EBAMRCNS* coarCNS = getCoarserLevel();
      int refRatCrse = coarCNS->refRatio();
      int nghost = 1;
      EBPWLFillPatch patcher(m_eblg.getDBL(),
                             coarCNS->m_eblg.getDBL(),
                             m_eblg.getEBISL(),
                             coarCNS->m_eblg.getEBISL(),
                             coarCNS->m_eblg.getDomain().domainBox(),
                             refRatCrse, m_nComp, nghost);

      Real coarTimeOld = 0.0;
      Real coarTimeNew = 1.0;
      Real fineTime    = 0.0;
      patcher.interpolate(consTemp,
                          coarCNS->m_stateOld,
                          coarCNS->m_stateNew,
                          coarTimeOld,
                          coarTimeNew,
                          fineTime,
                          intervDensity);
    }
  consTemp.exchange(intervDensity);

  // Compute undivided gradient
  for (DataIterator dit = m_eblg.getDBL().dataIterator(); dit.ok(); ++dit)
    {
      const Box& b = m_eblg.getDBL().get(dit());
      const EBISBox& ebisBox = m_eblg.getEBISL()[dit()];
      EBCellFAB gradFab(ebisBox, b, SpaceDim);
      const EBCellFAB& stateFab = consTemp[dit()];
      BaseFab<Real>& regGradFab = gradFab.getSingleValuedFAB();
      const BaseFab<Real>& regStateFab = stateFab.getSingleValuedFAB();

      for (int idir = 0; idir < SpaceDim; ++idir)
        {
          const Box bCenter = b & grow(m_eblg.getDomain(),-BASISV(idir));
          const Box bLo     = b & adjCellLo(bCenter,idir);
          const int hasLo = ! bLo.isEmpty();
          const Box bHi     = b & adjCellHi(bCenter,idir);
          const int hasHi = ! bHi.isEmpty();

          FORT_GETRELATIVEGRAD(CHF_FRA1(regGradFab,idir),
                               CHF_CONST_FRA1(regStateFab,0),
                               CHF_CONST_INT(idir),
                               CHF_BOX(bLo),
                               CHF_CONST_INT(hasLo),
                               CHF_BOX(bHi),
                               CHF_CONST_INT(hasHi),
                               CHF_BOX(bCenter));

          //do one-sided diffs where necessary at irregular cells.
          IntVectSet ivsIrreg = ebisBox.getIrregIVS(b);
          for(VoFIterator vofit(ivsIrreg, ebisBox.getEBGraph());
              vofit.ok(); ++vofit)
            {
              const VolIndex& vof = vofit();
              const IntVect&  iv = vof.gridIndex();
              //one-sided diffs on domain bndry
              bool onLeftDomain = iv[idir] == m_eblg.getDomain().domainBox().smallEnd(idir);
              bool onRighDomain = iv[idir] == m_eblg.getDomain().domainBox().bigEnd(idir);
              bool hasFacesLeft = (ebisBox.numFaces(vof, idir, Side::Lo) > 0) && !onLeftDomain;
              bool hasFacesRigh = (ebisBox.numFaces(vof, idir, Side::Hi) > 0) && !onRighDomain;

              Real valCent = stateFab(vof, densityIndex);
              Real dpl = 0.;
              Real dpr = 0.;
              Real dpc = 0.;

              //compute one-sided diffs where you have them
              if(hasFacesLeft)
                {
                  Vector<FaceIndex> facesLeft =
                    ebisBox.getFaces(vof, idir, Side::Lo);
                  Real valLeft = 0.0;
                  for(int iface = 0; iface <facesLeft.size(); iface++)
                    {
                      VolIndex vofLeft = facesLeft[iface].getVoF(Side::Lo);
                      valLeft += stateFab(vofLeft, densityIndex);
                    }
                  valLeft /= Real(facesLeft.size());
                  dpl = valCent - valLeft;
                }
              if(hasFacesRigh)
                {
                  Vector<FaceIndex> facesRigh =
                    ebisBox.getFaces(vof, idir, Side::Hi);
                  Real valRigh = 0.0;
                  for(int iface = 0; iface <facesRigh.size(); iface++)
                    {
                      VolIndex vofRigh = facesRigh[iface].getVoF(Side::Hi);
                      valRigh += stateFab(vofRigh, densityIndex);
                    }
                  valRigh /= Real(facesRigh.size());
                  dpr = valRigh - valCent;
                }
              if(hasFacesLeft && hasFacesRigh)
                {
                  dpc = 0.5*(dpl+dpr);
                }
              else if(!hasFacesLeft && !hasFacesRigh)
                {
                  dpc = 0.0;
                }
              else if(hasFacesLeft && !hasFacesRigh)
                {
                  dpc = dpl;
                }
              else if(hasFacesRigh && !hasFacesLeft)
                {
                  dpc = dpr;
                }

              gradFab(vof, idir) = dpc/valCent;
            }
        }

      EBCellFAB gradMagFab(ebisBox, b, 1);
      BaseFab<Real>& regGradMagFab = gradMagFab.getSingleValuedFAB();
      FORT_MAGNITUDE(CHF_FRA1(regGradMagFab,0),
                     CHF_CONST_FRA(regGradFab),
                     CHF_BOX(b));

      //pointwise op so just have to iterate over multivalued cells
      IntVectSet ivsMulti = ebisBox.getMultiCells(b);
      for(VoFIterator vofit(ivsMulti, ebisBox.getEBGraph());
          vofit.ok(); ++vofit)
        {
          const VolIndex& vof = vofit();
          Real mag = 0.0;
          for(int idir = 0; idir < SpaceDim; idir++)
            {
              Real graddir = gradFab(vof, idir);
              mag += graddir*graddir;
            }
          mag = sqrt(mag);
          gradMagFab(vof, 0) = mag;
        }

      // Tag where gradient exceeds threshold

      IntVectSet ivsTot(b);
      for(VoFIterator vofit(ivsTot, ebisBox.getEBGraph());
          vofit.ok(); ++vofit)
        {
          const VolIndex& vof = vofit();
          const IntVect& iv = vof.gridIndex();
          if (gradMagFab(vof, 0) >= m_params.m_refineThresh)
            {
              localTags |= iv;
            }
        }

      bool tagAllIrregular = false;
      ParmParse pp;
      if(pp.contains("tag_all_irregular"))
        {
          pp.get("tag_all_irregular", tagAllIrregular);
        }
      if(tagAllIrregular)
        {
          //refine all irregular cells
          //this is probably not ideal.
          IntVectSet irregIVS = ebisBox.getIrregIVS(b);
          localTags |= irregIVS;
        }
    }

  localTags.grow(m_params.m_tagBufferSize);

  // Need to do this in two steps unless a IntVectSet::operator &=
  // (ProblemDomain) operator is defined
  Box localTagsBox = localTags.minBox();
  localTagsBox &= m_eblg.getDomain();
  localTags &= localTagsBox;
  a_tags = localTags;
}
/***************************/
void
EBAMRCNS::
tagCellsInit(IntVectSet& a_tags)
{
  if (m_params.m_verbosity >= 3)
    {
      pout() << "EBAMRCNS::tagCellsInit for level " << m_level << endl;
    }

  // Since tags are calculated using only current time step data, use
  // the same tagging function for initialization and for regridding.
  tagCells(a_tags);
}
/***************************/
void
EBAMRCNS::
regrid(const Vector<Box>& a_new_grids)
{
  //first save old data
  // save data for later copy
  //not using m_eblg.getEBISL() because it gets wiped later
  const EBIndexSpace* const ebisPtr = Chombo_EBIS::instance();
  EBISLayout ebislOld;
  ebisPtr->fillEBISLayout(ebislOld, m_eblg.getDBL(), m_eblg.getDomain().domainBox(), m_nGhost);
  EBCellFactory factoryOld(ebislOld);
  IntVect ivGhost = m_nGhost*IntVect::Unit;
  LevelData<EBCellFAB> stateSaved(m_eblg.getDBL(), m_nComp, ivGhost, factoryOld);
  Interval interv(0,m_nComp-1);
  stateSaved.define(m_eblg.getDBL(), m_nComp, ivGhost, factoryOld);
  m_stateNew.copyTo(interv, stateSaved, interv);

  //create grids and ebis layouts
  //the ebisl has to know about the fact that we really have
  //four ghost cells.
  m_level_grids = a_new_grids;
  Vector<int> proc_map;

  EBLoadBalance(proc_map,a_new_grids, m_eblg.getDomain().domainBox());

  DisjointBoxLayout grids(a_new_grids, proc_map);

  m_eblg.define(grids, m_problem_domain, m_nGhost, ebisPtr);

  if(m_params.m_verbosity >= 3)
    {
      pout() << " in EBAMRCNS regrid for level " << m_level << endl;
    }



  // set up data structures
  levelSetup();

  // interpolate to coarser level
  if (m_hasCoarser)
    {
      EBAMRCNS* coarPtr = getCoarserLevel();
      m_ebFineInterp.interpolate(m_stateNew,
                                 coarPtr->m_stateNew,
                                 interv);
    }

  // copy from old state
  stateSaved.copyTo(interv,m_stateNew, interv);
  defineSolvers();
}
/***************************/
void
EBAMRCNS::
initialGrid(const Vector<Box>& a_new_grids)
{
  if(m_params.m_verbosity >= 3)
    {
      pout() << " in EBAMRCNS initialGrid for level " << m_level << endl;
    }

  const EBIndexSpace* const ebisPtr = Chombo_EBIS::instance();
  m_level_grids = a_new_grids;

  // load balance and create boxlayout
  Vector<int> proc_map;
  EBLoadBalance(proc_map,a_new_grids, m_problem_domain.domainBox());
  if(m_params.m_verbosity >= 3)
    {
      pout() << " just loadbalanced " << m_level << endl;
    }

  DisjointBoxLayout grids(a_new_grids,proc_map);
  m_eblg.define(grids, m_problem_domain, m_nGhost, ebisPtr);
  if(m_params.m_verbosity >= 3)
    {
      pout() << "EBAMRCNS::initialgrid grids " << endl;
      DisjointBoxLayout dbl = m_eblg.getDBL();
      dumpDBL(&dbl);
    }

  // set up data structures
  levelSetup();
  defineSolvers();
}
/***************************/
void
EBAMRCNS::
initialData()
{

  if(m_params.m_verbosity >= 3)
    {
      pout() << " in EBAMRCNS initialData for level " << m_level << endl;
    }

  const EBPhysIBC* const ebphysIBCPtr =
    m_ebPatchGodunov->getEBPhysIBC();

  //initialize both new and old states to
  //be the same thing
  ebphysIBCPtr->initialize(m_stateNew, m_eblg.getEBISL());
  ebphysIBCPtr->initialize(m_stateOld, m_eblg.getEBISL());
}
/***************************/
void
EBAMRCNS::
postInitialize()
{
}
void
EBAMRCNS::
syncWithFineLevel()
{
  if(m_params.m_verbosity >= 3)
    {
      pout() << " in EBAMRCNS syncWithFineLevel for level " << m_level << endl;
    }
  //stuff that needs to be setup from the finer
  //level.  A bunch of objects depend upon the layouts
  //from both levels and the finer level changes more
  //often from regrid so this needs to be called from the finer
  //level
  CH_assert(m_hasFiner);
  if(m_hasFiner)
    {
      EBAMRCNS* finePtr = getFinerLevel();
      int nRefFine = refRatio();
      const EBLevelGrid& finer_eblg = finePtr->m_eblg;
      const DisjointBoxLayout& finer_dbl = finer_eblg.getDBL();
      const EBISLayout& finer_ebisl      = finer_eblg.getEBISL();
      // maintain flux registers
      m_divFFluxRegister.define(finer_dbl,
                                m_eblg.getDBL(),
                                finer_ebisl,
                                m_eblg.getEBISL(),
                                m_eblg.getDomain().domainBox(),
                                nRefFine,
                                m_nComp, Chombo_EBIS::instance());

      m_veloFluxRegister.define(finer_dbl,
                                m_eblg.getDBL(),
                                finer_ebisl,
                                m_eblg.getEBISL(),
                                m_eblg.getDomain().domainBox(),
                                nRefFine,
                                SpaceDim, Chombo_EBIS::instance());

      m_tempFluxRegister.define(finer_dbl,
                                m_eblg.getDBL(),
                                finer_ebisl,
                                m_eblg.getEBISL(),
                                m_eblg.getDomain().domainBox(),
                                nRefFine,
                                1, Chombo_EBIS::instance());

      //define fine to coarse redistribution object
      //for now set to volume weighting
      m_ebCoarToFineRedist.define(finer_eblg, m_eblg, nRefFine , m_nComp, 1);
      //define coarse to coarse redistribution object
      m_ebCoarToCoarRedist.define(finer_eblg, m_eblg, nRefFine , m_nComp, 1);

      //set all the registers to zero
      m_ebCoarToFineRedist.setToZero();
      m_ebCoarToCoarRedist.setToZero();
      m_divFFluxRegister.setToZero();
      m_veloFluxRegister.setToZero();
      m_tempFluxRegister.setToZero();
    }

}
/***************************/
Real
EBAMRCNS::
computeDt()
{
  return m_dtNew;
}
/***************************/
Real
EBAMRCNS::
computeInitialDt()
{
  Real maxwavespeed =  m_ebLevelGodunov.getMaxWaveSpeed(m_stateNew);
  CH_assert(Abs(m_dx[0]-m_dx[1])<1.e-9);
  Real newDT = m_initial_dt_multiplier * m_dx[0] /maxwavespeed;

  return newDT;
}
/***************************/
void
EBAMRCNS::
sumConserved(Real& a_sumcons,
             const int& a_ivar) const
{
  Real sumgrid = 0;
  for(DataIterator dit= m_eblg.getDBL().dataIterator(); dit.ok(); ++dit)
    {
      const EBCellFAB& state = m_stateNew[dit()];
      Box thisBox = m_eblg.getDBL().get(dit());
      IntVectSet uberIVS(thisBox);
      const EBISBox& ebisBox = m_eblg.getEBISL()[dit()];
      for(VoFIterator vofit(uberIVS, ebisBox.getEBGraph());
          vofit.ok(); ++vofit)
        {
          const VolIndex& vof = vofit();
          Real consVal = state(vof, a_ivar);
          Real volFrac = ebisBox.volFrac(vof);
          Real volume = volFrac;
          sumgrid += consVal*volume;
        }
    }

  Vector<Real> all_sum;
  gather(all_sum,sumgrid,uniqueProc(SerialTask::compute));
  Real sumallgrid = 0.;
  if (procID() == uniqueProc(SerialTask::compute))
    {
      for (int i = 0; i < all_sum.size(); ++i)
        {
          sumallgrid += all_sum[i];
        }
    }
  broadcast(sumallgrid,uniqueProc(SerialTask::compute));
  a_sumcons = sumallgrid;
}
/***************************/
EBAMRCNS*
EBAMRCNS::
getCoarserLevel() const
{
  EBAMRCNS* retval = NULL;
  if(m_coarser_level_ptr != NULL)
    {
      retval = dynamic_cast <EBAMRCNS*> (m_coarser_level_ptr);

      if (retval == NULL)
        {
          pout() << "EBAMRG::getCoarserLevel: dynamic cast failed"
                 << endl;
          MayDay::Error();
        }
    }
  return retval;
}

/***************************/
EBAMRCNS*
EBAMRCNS::
getFinerLevel() const
{
  EBAMRCNS* retval = NULL;
  if(m_finer_level_ptr != NULL)
    {
      retval = dynamic_cast <EBAMRCNS*> (m_finer_level_ptr);

      if (retval == NULL)
        {
          pout() << "EBAMRG::getFinerLevel: dynamic cast failed"
                 << endl;
          MayDay::Error();
        }
    }
  return retval;
}
/***************************/
void
EBAMRCNS::
levelSetup()
{
  if(m_params.m_verbosity >= 3)
    {
      pout() << " in EBAMRCNS levelSetup for level " << m_level << endl;
    }

  EBCellFactory factoryNew(m_eblg.getEBISL());
  IntVect ivGhost = m_nGhost*IntVect::Unit;
  m_stateNew.define(m_eblg.getDBL(),m_nComp, ivGhost, factoryNew);
  m_stateOld.define(m_eblg.getDBL(),m_nComp, ivGhost, factoryNew);
  m_rhsVelo.define( m_eblg.getDBL(),SpaceDim,ivGhost, factoryNew);
  m_rhsTemp.define( m_eblg.getDBL(),       1,ivGhost, factoryNew);
  m_sets.define(m_eblg.getDBL());
  for(DataIterator dit = m_eblg.getDBL().dataIterator(); dit.ok(); ++dit)
    {
      m_sets[dit()] = m_eblg.getEBISL()[dit()].getIrregIVS(m_eblg.getDBL().get(dit()));
    }
  EBCellFactory       cellFact(m_eblg.getEBISL());
  EBFluxFactory       fluxFact(m_eblg.getEBISL());
  BaseIVFactory<Real> bivfFact(m_eblg.getEBISL(), m_sets);
  m_acoVelo     = RefCountedPtr< LevelData<EBCellFAB> >       (new LevelData<EBCellFAB>       (m_eblg.getDBL(), 1, IntVect::Unit, cellFact));
  m_acoTemp     = RefCountedPtr< LevelData<EBCellFAB> >       (new LevelData<EBCellFAB>       (m_eblg.getDBL(), 1, IntVect::Unit, cellFact));
  m_eta         = RefCountedPtr< LevelData<EBFluxFAB> >       (new LevelData<EBFluxFAB>       (m_eblg.getDBL(), 1, IntVect::Unit, fluxFact));
  m_lambda      = RefCountedPtr< LevelData<EBFluxFAB> >       (new LevelData<EBFluxFAB>       (m_eblg.getDBL(), 1, IntVect::Unit, fluxFact));
  m_bcoTemp     = RefCountedPtr< LevelData<EBFluxFAB> >       (new LevelData<EBFluxFAB>       (m_eblg.getDBL(), 1, IntVect::Unit, fluxFact));
  m_etaIrreg    = RefCountedPtr< LevelData<BaseIVFAB<Real> > >(new LevelData<BaseIVFAB<Real> >(m_eblg.getDBL(), 1, IntVect::Unit, bivfFact));
  m_bcoTempIrreg= RefCountedPtr< LevelData<BaseIVFAB<Real> > >(new LevelData<BaseIVFAB<Real> >(m_eblg.getDBL(), 1, IntVect::Unit, bivfFact));
  m_lambdaIrreg = RefCountedPtr< LevelData<BaseIVFAB<Real> > >(new LevelData<BaseIVFAB<Real> >(m_eblg.getDBL(), 1, IntVect::Unit, bivfFact));


  EBAMRCNS* coarPtr = getCoarserLevel();
  EBAMRCNS* finePtr = getFinerLevel();

  m_hasCoarser = (coarPtr != NULL);
  m_hasFiner   = (finePtr != NULL);

  //define redistribution object for this level
  //for now set to volume weighting
  m_ebLevelRedist.define(m_eblg.getDBL(),
                         m_eblg.getEBISL(),
                         m_eblg.getDomain(),
                         m_nComp);

  if (m_hasCoarser)
    {
      int nRefCrse = m_coarser_level_ptr->refRatio();
      const EBLevelGrid& coEBLG = coarPtr->m_eblg;

      m_ebCoarseAverage.define(m_eblg.getDBL(),
                               coEBLG.getDBL(),
                               m_eblg.getEBISL(),
                               coEBLG.getEBISL(),
                               coEBLG.getDomain().domainBox(),
                               nRefCrse,
                               m_nComp, Chombo_EBIS::instance());
      m_ebFineInterp.define(m_eblg.getDBL(),
                            coEBLG.getDBL(),
                            m_eblg.getEBISL(),
                            coEBLG.getEBISL(),
                            coEBLG.getDomain().domainBox(),
                            nRefCrse,
                            m_nComp);

      // maintain levelgodunov
      m_ebLevelGodunov.define(m_eblg.getDBL(),
                              coEBLG.getDBL(),
                              m_eblg.getEBISL(),
                              coEBLG.getEBISL(),
                              m_eblg.getDomain(),
                              nRefCrse,
                              m_dx,
                              m_params,
                              m_ebPatchGodunov,
                              m_hasCoarser,
                              m_hasFiner);

      //define fine to coarse redistribution object
      //for now set to volume weighting
      {
        CH_TIME("fineToCoar_defs");
        m_ebFineToCoarRedist.define(m_eblg, coEBLG,
                                    nRefCrse, m_nComp, 
                                    m_params.m_redistRad);
                                   
        m_ebFineToCoarRedist.setToZero();
      }

      int nvarQuad = 1; //temperature
      //no EBCF crossing for now because viscous tensor op does not have it
      bool ebCF = false;
      m_quadCFI = RefCountedPtr<EBQuadCFInterp>
        (new EBQuadCFInterp(m_eblg.getDBL(),
                            coEBLG.getDBL(),
                            m_eblg.getEBISL(),
                            coEBLG.getEBISL(),
                            coEBLG.getDomain(),
                            nRefCrse, nvarQuad,
                            (*m_eblg.getCFIVS()),
                            Chombo_EBIS::instance(),
                            ebCF));

      coarPtr->syncWithFineLevel();
    }
  else
    {
      m_quadCFI = RefCountedPtr<EBQuadCFInterp>(new EBQuadCFInterp());
      m_ebLevelGodunov.define(m_eblg.getDBL(),
                              DisjointBoxLayout(),
                              m_eblg.getEBISL(),
                              EBISLayout(),
                              m_eblg.getDomain(),
                              m_ref_ratio,
                              m_dx,
                              m_params,
                              m_ebPatchGodunov,
                              m_hasCoarser,
                              m_hasFiner);
    }

  m_sets.define(m_eblg.getDBL());
  for(DataIterator dit = m_eblg.getDBL().dataIterator();
      dit.ok(); ++dit)
    {
      Box thisBox = m_eblg.getDBL().get(dit());
      m_sets[dit()] = m_eblg.getEBISL()[dit()].getIrregIVS(thisBox);
    }
  BaseIVFactory<Real> factory(m_eblg.getEBISL(), m_sets);
  //the ghost cells on the mass redistribution array
  //are tied to the redistribution radius
  m_massDiff.define(m_eblg.getDBL(), m_nComp, ivGhost, factory);
  for(DataIterator dit = m_eblg.getDBL().dataIterator(); dit.ok(); ++dit)
    {
      m_massDiff[dit()].setVal(0.0);
    }
}
/***************************/
LevelData<EBCellFAB>&
EBAMRCNS::
getStateNew()
{
  return m_stateNew;
}
/***************************/
LevelData<EBCellFAB>&
EBAMRCNS::
getStateOld()
{
  return m_stateOld;
}
/***************************/
Real
EBAMRCNS::
getDt() const
{
  return m_dt;
}
/***************************/
EBISLayout
EBAMRCNS::
getEBISLayout() const
{
  EBISLayout ebisl = m_eblg.getEBISL();
  return ebisl;
}
/***************************/
void
EBAMRCNS::
fillConsAndPrim(LevelData<EBCellFAB>& a_data)
{
  if (m_params.m_verbosity >= 3)
    {
      pout() << "EBAMRCNS::fillConsAndPrim" << endl;
    }
  //fill output data including a bunch of geometric crud
  int nCons = m_ebPatchGodunov->numConserved();
  int nPrim = m_ebPatchGodunov->numPrimitives() ;
  int consAndPrim = nCons + nPrim;

  EBCellFactory ebcellfact(m_eblg.getEBISL());
  a_data.define(m_eblg.getDBL(), consAndPrim, IntVect::Zero, ebcellfact);

  for (DataIterator dit = m_eblg.getDBL().dataIterator(); dit.ok(); ++dit)
    {
      const EBISBox& ebisbox = m_eblg.getEBISL()[dit()];
      const Box& grid = m_eblg.getDBL().get(dit());
      const EBCellFAB& consfab = m_stateNew[dit()];
      EBCellFAB primfab(ebisbox, grid, nPrim);
      //cfivs, time and timestep fake and not used here
      Real faket = 1.0;
      int logflag = 0;
      IntVectSet emptyivs;
      m_ebPatchGodunov->setValidBox(grid, ebisbox, emptyivs, faket, faket);
      m_ebPatchGodunov->consToPrim(primfab, consfab, grid,logflag);

      EBCellFAB& outputfab = a_data[dit()];

      Interval consIntervSrc(0, nCons-1);
      Interval consIntervDst(0, nCons-1);
      Interval primIntervSrc(0, nPrim-1);
      Interval primIntervDst(nCons, consAndPrim-1);

      // copy regular data
      outputfab.copy(grid, consIntervDst,  grid, consfab, consIntervSrc);
      outputfab.copy(grid, primIntervDst,  grid, primfab, primIntervSrc);
      Real coveredVal = -10.0;
      for(int ivar = 0; ivar < consAndPrim; ivar++)
        {
          outputfab.setInvalidData(coveredVal, ivar);
        }

    }//end loop over grids
}
/***************************/
#ifdef CH_USE_HDF5
/***************************/
void
EBAMRCNS::
writeCheckpointHeader(HDF5Handle& a_handle) const
{
  if (m_params.m_verbosity >= 3)
    {
      pout() << "EBAMRCNS::writeCheckpointHeader" << endl;
    }

  //stuff in non-eb checkpoint header is already
  //set in the define function, such as
  // Setup the number of components
  // Setup the component names
  //so i will eliminate the middleman
}
/***************************/
void
EBAMRCNS::
readCheckpointHeader(HDF5Handle& a_handle)
{
  if (m_params.m_verbosity >= 3)
    {
      pout() << "EBAMRCNS::readCheckpointHeader" << endl;
    }

  //stuff in non-eb checkpoint header is already
  //set in the define function, such as
  // Setup the number of components
  // Setup the component names
  // So i will eliminate the middleman.
}
/***************************/
void
EBAMRCNS::
writeCheckpointLevel(HDF5Handle& a_handle) const
{
  if (m_params.m_verbosity >= 3)
    {
      pout() << "EBAMRCNS::writeCheckpointLevel" << endl;
    }

  //reasons for deviations from non-eb stuff
  // the tag buffer size is set by the factory.
  // dx is set in the define function
  // the problem domain is set in the define function

  // Setup the level string
  char levelStr[20];
  sprintf(levelStr,"%d",m_level);
  const std::string label = std::string("level_") + levelStr;

  a_handle.setGroup(label);

  // Setup the level header information
  HDF5HeaderData header;

  header.m_int ["ref_ratio"]       = m_ref_ratio;
  header.m_real["dt"]              = m_dt;

  if (m_params.m_verbosity >= 3)
    {
      pout() << "EBAMRCNS::writeCheckpointLevel 2" << endl;
    }

  // Write the header for this level
  header.writeToFile(a_handle);

  // Write the data for this level
  if (m_params.m_verbosity >= 3)
    {
      pout() << "EBAMRCNS::writeCheckpointLevel 3" << endl;
    }
  write(a_handle,m_eblg.getDBL());
  if (m_params.m_verbosity >= 3)
    {
      pout() << "EBAMRCNS::writeCheckpointLevel 4" << endl;
    }
  write(a_handle,m_stateOld,"dataOld");
  if (m_params.m_verbosity >= 3)
    {
      pout() << "EBAMRCNS::writeCheckpointLevel 5" << endl;
    }
  write(a_handle,m_stateNew,"dataNew");
  if (m_params.m_verbosity >= 3)
    {
      pout() << "EBAMRCNS::writeCheckpointLevel 6" << endl;
    }
}
/***************************/
/***************************/
void
EBAMRCNS::
readCheckpointLevel(HDF5Handle& a_handle)
{
  CH_assert(m_isDefined);
  if (m_params.m_verbosity >= 3)
    {
      pout() << "EBAMRCNS::readCheckpointLevel" << endl;
    }

  // Setup the level string
  char levelStr[20];
  sprintf(levelStr,"%d",m_level);
  const std::string label = std::string("level_") + levelStr;

  // Read the header for this level
  a_handle.setGroup(label);

  HDF5HeaderData header;
  header.readFromFile(a_handle);

  // Get the refinement ratio
  if (header.m_int.find("ref_ratio") == header.m_int.end())
    {
      MayDay::Error("::readCheckpointLevel: file does not contain ref_ratio");
    }
  m_ref_ratio = header.m_int["ref_ratio"];

  //reasons for deviations from non-eb stuff
  // the tag buffer size is set by the factory.
  // dx is set in the define function
  // the problem domain is set in the define function

  // Get dt
  if (header.m_real.find("dt") == header.m_real.end())
    {
      MayDay::Error("readCheckpointLevel: file does not contain dt");
    }
  m_dt = header.m_real["dt"];

  // Get the grids
  Vector<Box> vboxGrids;
  const int gridStatus = read(a_handle, vboxGrids);
  if (gridStatus != 0)
    {
      MayDay::Error("readCheckpointLevel: file has no grids");
    }

  Vector<int> proc_map;
  EBLoadBalance(proc_map,vboxGrids, m_eblg.getDomain().domainBox());

  const EBIndexSpace* const ebisPtr = Chombo_EBIS::instance();
  DisjointBoxLayout grids(vboxGrids,proc_map);
  m_eblg.define(grids, m_problem_domain, m_nGhost, ebisPtr);

  //this keeps the order of the AMRLevel m_level_grids
  //consistent with m_eblg.getDBL()
  LayoutIterator lit = m_eblg.getDBL().layoutIterator();
  for (lit.begin(); lit.ok(); ++lit)
    {
      Box b = m_eblg.getDBL().get(lit());
      m_level_grids.push_back(b);
    }

  EBCellFactory factoryNew(m_eblg.getEBISL());
  //m_nghost is set in define function
  IntVect ivGhost = m_nGhost*IntVect::Unit;
  m_stateNew.define(m_eblg.getDBL(),m_nComp, ivGhost, factoryNew);
  m_stateOld.define(m_eblg.getDBL(),m_nComp, ivGhost, factoryNew);
  m_rhsVelo.define( m_eblg.getDBL(),SpaceDim,ivGhost, factoryNew);
  m_rhsTemp.define( m_eblg.getDBL(),       1,ivGhost, factoryNew);

  //  Interval vars(0, m_nComp-1);
  //the false says to not redefine the data
  int dataStatusNew = read<EBCellFAB>(a_handle,
                                      m_stateNew,
                                      "dataNew",
                                      m_eblg.getDBL(),
                                      Interval(),
                                      false);

  int dataStatusOld = read<EBCellFAB>(a_handle,
                                      m_stateOld,
                                      "dataOld",
                                      m_eblg.getDBL(),
                                      Interval(),
                                      false);

  if ((dataStatusNew != 0) || (dataStatusOld != 0))
    {
      MayDay::Error("file does not contain state data");
    }
  // Set up data structures
  levelSetup();
}

/***************************/
void
EBAMRCNS::
writePlotHeader(HDF5Handle& a_handle) const
{
  if(s_NewPlotFile == 0)
    {
      writePlotHeaderOld(a_handle);
      return;
    }
  Vector<int> refRatios;
  const EBAMRCNS* current = this;
  int nlevs = 0;
  while(current != NULL){
    refRatios.push_back(current->refRatio());
    nlevs++;
    current = (const EBAMRCNS*)(current-> m_finer_level_ptr);
  }

  headerEBFile(a_handle, nlevs, refRatios,
               m_eblg.getDomain().domainBox(), m_dx, m_stateNew.ghostVect());

  writeCellCenteredNames(a_handle, m_stateNames);
  a_handle.setGroup("/Expressions");
  HDF5HeaderData expressions;
  EBPatchGodunov& patchgod = (EBPatchGodunov&)(*m_ebPatchGodunov);
  patchgod.expressions(expressions);
  expressions.writeToFile(a_handle);

}
/***************************/
void
EBAMRCNS::
writePlotHeaderOld(HDF5Handle& a_handle) const
{
  if (m_params.m_verbosity >= 3)
    {
      pout() << "EBAMRCNS::writePlotHeader" << endl;
    }

  HDF5HeaderData header;
  // Setup the number of components
  //have to add in a lot of geometric crap.
  // 3 norms + 6 area fracs + 1 distance + 1 volFrac
  // =  11 extra components
  //forces 3d
  int nCons = m_ebPatchGodunov->numConserved();
  int nPrim = m_ebPatchGodunov->numPrimitives() ;
  int consAndPrim = nCons + nPrim;

  int indexVolFrac = consAndPrim;
  int indexAreaFrac = indexVolFrac+1;
  int indexNormal = indexAreaFrac+ 2*SpaceDim;
  int indexDist = indexNormal+SpaceDim;
  int nCompTotal = indexDist+1;
  CH_assert(nCompTotal == consAndPrim + 3*SpaceDim+2);

  Vector<string> names(nCompTotal);

  for (int i = 0; i < nCons; i++)
    {
      names[i] = m_stateNames[i];
    }
  for (int i = 0; i < nPrim; i++)
    {
      names[nCons + i] = m_primNames[i];
    }

  string volFracName("fraction-0");
  Vector<string> normName(3);
  Vector<string> areaName(6);
  string distName("distance-0");

  normName[0] = "xnormal-0";
  normName[1] = "ynormal-0";
  normName[2] = "znormal-0";

  areaName[0] = "xAreafractionLo-0";
  areaName[1] = "xAreafractionHi-0";
  areaName[2] = "yAreafractionLo-0";
  areaName[3] = "yAreafractionHi-0";
  areaName[4] = "zAreafractionLo-0";
  areaName[5] = "zAreafractionHi-0";

  names[indexVolFrac] = volFracName;

  for (int i=0; i < 2*SpaceDim; i++)
    {
      names[indexAreaFrac+i] = areaName[i];
    }

  for (int i=0; i < SpaceDim; i++)
    {
      names[indexNormal+i] = normName[i];
    }

  names[indexDist] = distName;

  //now output this into the hdf5 handle
  header.m_int["num_components"] = nCompTotal;
  // Setup the component names
  char compStr[30];
  for (int comp = 0; comp < nCompTotal; ++comp)
    {
      sprintf(compStr,"component_%d",comp);
      header.m_string[compStr] = names[comp];
    }

  // Write the header to the file
  header.writeToFile(a_handle);

  if (m_params.m_verbosity >= 4)
    {
      pout() << header << endl;
    }
}

/***************************/
void EBAMRCNS::writePlotLevel(HDF5Handle& a_handle) const
{
  if(s_NewPlotFile == 0)
    {
      writePlotLevelOld(a_handle);
      return;
    }
  if (m_params.m_verbosity >= 3)
    {
      pout() << "EBAMRCNS::writePlotLevel " << m_level<< endl;
    }
  writeCellCentered(a_handle, m_level, &m_stateNew);
}

/***************************/
void EBAMRCNS::writePlotLevelOld(HDF5Handle& a_handle) const
{

  if (m_params.m_verbosity >= 3)
    {
      pout() << "EBAMRCNS::writePlotLevel" << endl;
    }
  //fill output data including a bunch of geometric crud
  int nCons = m_ebPatchGodunov->numConserved();
  int nPrim = m_ebPatchGodunov->numPrimitives() ;
  int consAndPrim = nCons + nPrim;
  int indexVolFrac = consAndPrim;
  int indexAreaFrac = indexVolFrac+1;
  int indexNormal = indexAreaFrac+ 2*SpaceDim;
  int indexDist = indexNormal+SpaceDim;
  int nCompTotal = indexDist+1;
  CH_assert(nCompTotal == consAndPrim + 3*SpaceDim+2);

  Vector<Real> coveredValuesCons(nCons, -10.0);
  Vector<Real> coveredValuesPrim(nPrim, -10.0);

  Vector<Real> coveredValues;
  coveredValues.append(coveredValuesCons);
  coveredValues.append(coveredValuesPrim);

  LevelData<FArrayBox> fabData(m_eblg.getDBL(), nCompTotal, IntVect::Zero);

#ifdef CH_MPI
  MPI_Barrier(Chombo_MPI::comm);
#endif
  for (DataIterator dit = m_eblg.getDBL().dataIterator(); dit.ok(); ++dit)
    {
      const EBISBox& ebisbox = m_eblg.getEBISL()[dit()];
      const Box& grid = m_eblg.getDBL().get(dit());
      const EBCellFAB& consfab = m_stateNew[dit()];
      EBCellFAB primfab(ebisbox, grid, nPrim);
      //cfivs, time and timestep fake and not used here
      Real faket = 1.0;
      IntVectSet emptyivs;
      ParmParse pp;
      int logflag = 1;
      if(pp.contains("logflag"))
        {
          pp.get("logflag", logflag);
        }
      EBPatchGodunov& patchgod = (EBPatchGodunov&)(*m_ebPatchGodunov);
      patchgod.setValidBox(grid, ebisbox, emptyivs, faket, faket);
      patchgod.consToPrim(primfab, consfab, grid, logflag);

      FArrayBox& currentFab = fabData[dit()];

      // copy regular data
      currentFab.copy(consfab.getSingleValuedFAB(),0,0,nCons);
      currentFab.copy(primfab.getSingleValuedFAB(),0,nCons,nPrim);

      // set default volume fraction
      currentFab.setVal(1.0,indexVolFrac);

      // set default area fractions
      for (int i=0; i < 2*SpaceDim; i++)
        {
          currentFab.setVal(1.0,indexAreaFrac+i);
        }

      // set default normal
      for (int i=0; i < SpaceDim; i++)
        {
          currentFab.setVal(0.0,indexNormal+i);
        }

      // set default distance of EB from corner
      currentFab.setVal(0.0,indexDist);

      // set special values
      // iterate through the current grid
      // NOTE:  this is probably an inefficient way to do this
      for (BoxIterator bit(grid); bit.ok(); ++bit)
        {
          const IntVect& iv = bit();
          // set special values for covered cells
          if (ebisbox.isCovered(iv))
            {
              for (int icomp = 0; icomp < consAndPrim; icomp++)
                {
                  Real cval = coveredValues[icomp];

                  currentFab(iv,icomp) = cval;
                }
              // volume fraction is zero
              currentFab(iv,indexVolFrac) = 0.0;

              // area fractions are zero
              for (int i=0; i < 2*SpaceDim; i++)
                {
                  currentFab(iv,indexAreaFrac+i) = 0.0;
                }
            }

          // set special values for irregular cells
          if (ebisbox.isIrregular(iv))
            {
              Vector<VolIndex> vofs = ebisbox.getVoFs(iv);
              Real volFrac = ebisbox.volFrac(vofs[0]);
              RealVect normal = ebisbox.normal(vofs[0]);

              // set volume fraction
              currentFab(iv,indexVolFrac) = volFrac;

              // set area fractions--use only the first face you find
              for (int i=0; i < SpaceDim; i++)
                {
                  Vector<FaceIndex> faces;

                  faces = ebisbox.getFaces(vofs[0],i,Side::Lo);
                  if (faces.size() == 0)
                    {
                      currentFab(iv,indexAreaFrac+2*i) = 0.0;
                    }
                  else
                    {
                      currentFab(iv,indexAreaFrac+2*i) =
                        ebisbox.areaFrac(faces[0]);
                    }

                  faces = ebisbox.getFaces(vofs[0],i,Side::Hi);
                  if (faces.size() == 0)
                    {
                      currentFab(iv,indexAreaFrac+2*i+1) = 0.0;
                    }
                  else
                    {
                      currentFab(iv,indexAreaFrac+2*i+1) =
                        ebisbox.areaFrac(faces[0]);
                    }
                }

              // set normal
              for (int i=0; i < SpaceDim; i++)
                {
                  currentFab(iv,indexNormal+i) = normal[i];
                }

              // set distance unless the length of the normal is zero
              Real length = PolyGeom::dot(normal,normal);

              if (length > 0)
                {
                  Real dist = PolyGeom::computeAlpha(volFrac,normal)*m_dx[0];
                  currentFab(iv,indexDist) = -dist;
                }
            } //end if(isIrregular)
        }//end loop over cells
    }//end loop over grids

#ifdef CH_MPI
  MPI_Barrier(Chombo_MPI::comm);
#endif
  // Setup the level string
  char levelStr[20];
  sprintf(levelStr,"%d",m_level);
  const std::string label = std::string("level_") + levelStr;

  a_handle.setGroup(label);

  // Setup the level header information
  HDF5HeaderData header;

  header.m_int ["ref_ratio"]   = m_ref_ratio;
  header.m_real["dx"]          = m_dx[0];
  header.m_real["dt"]          = m_dt;
  header.m_real["time"]        = m_time;
  header.m_box ["prob_domain"] = m_eblg.getDomain().domainBox();

  // Write the header for this level
  header.writeToFile(a_handle);

  if (m_params.m_verbosity >= 4)
    {
      pout() << header << endl;
    }

  // Write the data for this level
  write(a_handle,fabData.boxLayout());
  write(a_handle,fabData,"data");
}

#endif
#include "NamespaceFooter.H"

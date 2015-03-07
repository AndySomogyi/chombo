#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

#ifndef _BASELEVELTGA_H__
#define _BASELEVELTGA_H__

#include <iostream>
#include <math.h>
#include "SPACE.H"
#include <stdlib.h>
#include <REAL.H>
#include <Box.H>
#include <DisjointBoxLayout.H>
#include <LevelData.H>
#include <ProblemDomain.H>
#include "AMRTGA.H"
#include "NamespaceHeader.H"


/// Implements TGA algorithm to update solution to heat equation
/** The ResistivitySolver class implements the Runge-Kutta-based approach
    to solving the heat equation due to Twizell, Gumel, and Arigu.
*/
template <class T, class TFlux, class TFR>
class BaseLevelTGA
{

public:

  LevelTGAHelmOp<T, TFlux>* newOp(const ProblemDomain& a_indexSpace,
                                  RefCountedPtr<AMRLevelOpFactory<T > >&   a_opFact)
  {
    LevelTGAHelmOp<T, TFlux>* retval = (LevelTGAHelmOp<T,TFlux>*)
      (a_opFact->AMRnewOp(a_indexSpace));

    return retval;
  }

  /// full constructor
  /**
   */
  BaseLevelTGA(const Vector<DisjointBoxLayout>&            a_grids,
               const Vector<int>&                          a_refRat,
               const ProblemDomain&                        a_level0Domain,
               RefCountedPtr<AMRLevelOpFactory< T > >&     a_opFact,
               const RefCountedPtr<AMRMultiGrid<T > >&     a_solver);

  /// destructor
  virtual ~BaseLevelTGA()
  {
  }


  /// do solve
  /** update ALL components of phi
   */
  void updateSoln(T&       a_phiNew,
                  T&       a_phiOld,
                  T&       a_src,
                  LevelData<TFlux>&       a_flux,
                  TFR*                    a_FineFluxRegPtr,
                  TFR*                    a_CrseFluxRegPtr,
                  const T* a_crsePhiOldPtr,
                  const T* a_crsePhiNewPtr,
                  Real                    oldTime,
                  Real                    crseOldTime,
                  Real                    crseNewTime,
                  Real                    dt,
                  int                     a_level,
                  bool                    a_zeroPhi = true,
                  int                     a_fluxStartComponent = 0);


  /// computes time-centered diffusive term for explicit updating
  /** compute time-centered L(phi) for use in subsequent update
      operations.  In this case, we do solve for phiNew, then
      subtract source and old phi back out to get L(phi).
  */
  void computeDiffusion(T&       a_DiffusiveTerm,
                        T&       a_phiOld,
                        T&       a_src,
                        LevelData<TFlux>&       a_flux,
                        TFR*                    a_FineFluxRegPtr,
                        TFR*                    a_crseFluxRegPtr,
                        const T* a_crsePhiOldPtr,
                        const T* a_crsePhiNewPtr,
                        Real                    a_oldTime,
                        Real                    a_crseOldTime,
                        Real                    a_crseNewTime,
                        Real                    a_dt,
                        int                     a_level);


protected:

  void incrementFlux(LevelData<TFlux>& a_diffusiveFlux,
                     T& a_phi,
                     int               a_level,
                     Real              a_mu,
                     Real              a_dt,
                     Real              a_sign,
                     bool              a_setToZero);

  ///
  /**
     solves helm(a_phi) = a_rhs, with a_phiC on coarser level;
     (I - a_mu*a_dt*laplacian)(a_phi) = a_rhs
   */
  void solveHelm(T& a_phi,
                 T& a_phiC,
                 T& a_rhs,
                 int               a_level,
                 Real              a_mu,
                 Real              a_dt,
                 bool              a_zeroPhi = true);

  ///
  /**
     returns a_ans = helm(a_phi) with a_phiC pointer to coarser level
                   = (I + a_mu*a_dt*laplacian)(a_phi)

     the a_applyBC argument specifies both coarse-fine inhomogeneous BCs and domain
     boundary inhomogeneous BCs.  Do them all, inhomogeneously, or do none of them.
   */
  void applyHelm(T&          a_ans,
                 const T&    a_phi,
                 const T*    a_phiC,
                 int                        a_level,
                 Real                       a_mu,
                 Real                       a_dt,
                 bool                       a_applyBC);

  ///
  /**
     sets alpha and beta in helmholtz operator
   */
  void resetSolverAlphaAndBeta(const Real& a_alpha,
                               const Real& a_beta);

  ///
  /**
     returns a_data at time a_time, by interpolating from
     a_oldData at time a_oldTime, and
     a_newData at time a_newTime.
   */
  void timeInterp(T&       a_data,
                  const T& a_oldData,
                  const T& a_newData,
                  Real                    a_time,
                  Real                    a_oldTime,
                  Real                    a_newTime,
                  int                     a_Lev);

  //app dependent
  virtual void setSourceGhostCells(T&        a_src,
                                   const DisjointBoxLayout& a_dbl) = 0;



  Vector<DisjointBoxLayout>                             m_grids ;
  Vector<int>                                           m_refRat;
  ProblemDomain                                         m_level0Domain;
  Vector<RefCountedPtr<LevelTGAHelmOp<T,TFlux> > >  m_ops;
  RefCountedPtr<AMRMultiGrid<T > >       m_solver;
  Real m_mu1, m_mu2, m_mu3, m_mu4, m_r1;

private:

  void operator=(const BaseLevelTGA& a_opin)
  {
    MayDay::Error("invalid operator");
  }
  BaseLevelTGA(BaseLevelTGA& a_opin)
  {
    MayDay::Error("invalid operator");
  }
  /// default constructor
  BaseLevelTGA()
  {
    MayDay::Error("weak construction is bad");
  }

};

template <class T, class TFlux, class TFR>
BaseLevelTGA<T, TFlux, TFR>::
BaseLevelTGA(const Vector<DisjointBoxLayout>&            a_grids,
             const Vector<int>&                          a_refRat,
             const ProblemDomain&                        a_level0Domain,
             RefCountedPtr<AMRLevelOpFactory< T > >&     a_opFact,
             const RefCountedPtr<AMRMultiGrid<T > >&     a_solver)
{
  m_grids  =       a_grids;
  m_refRat =       a_refRat;
  m_level0Domain = a_level0Domain;
  m_solver =       a_solver;

  Real tgaEpsilon = 1.e-12;
  Real a = 2.0 - sqrt(2.0) - tgaEpsilon;
  m_mu1 = (a - sqrt(pow(a,2) - 4.0*a + 2.0))/2.0 ;
  m_mu2 = (a + sqrt(pow(a,2) - 4.0*a + 2.0))/2.0 ;
  m_mu3 = (1.0 - a);
  m_mu4 = 0.5 - a;

  Real discr = sqrt(a*a - 4.0*a + 2.0);
  m_r1 = (2.0*a - 1.0)/(a + discr);

  m_ops.resize(a_grids.size());
  ProblemDomain domLev = a_level0Domain;
  for (int ilev = 0; ilev < m_ops.size(); ilev++)
    {
      m_ops[ilev] = RefCountedPtr<LevelTGAHelmOp<T,TFlux> >
        (newOp(domLev, a_opFact));
      if (ilev < m_ops.size() - 1)
        {
          domLev.refine(a_refRat[ilev]);
        }
    }
}

template <class T, class TFlux, class TFR>
void
BaseLevelTGA<T, TFlux, TFR>::
resetSolverAlphaAndBeta(const Real& a_alpha,
                        const Real& a_beta)
{
  Vector<MGLevelOp<T >* > ops = m_solver->getAllOperators();
  for(int iop = 0; iop < ops.size(); iop++)
    {
      LevelTGAHelmOp<T,TFlux>* helmop = (LevelTGAHelmOp<T,TFlux>*) ops[iop];
      helmop->setAlphaAndBeta(a_alpha, a_beta);
    }
  for(int iop = 0; iop < m_ops.size(); iop++)
    {
      m_ops[iop]->setAlphaAndBeta(a_alpha, a_beta);
    }
}
template <class T, class TFlux, class TFR>
void
BaseLevelTGA<T, TFlux, TFR>::
applyHelm(T&         a_ans,
          const T&   a_phi,
          const T*   a_phiC,
          int                           a_level,
          Real                          a_mu,
          Real                          a_dt,
          bool                          a_applyBC)
{
  Real factor  = a_mu*a_dt;
  m_ops[a_level]->setAlphaAndBeta(1.0, factor);
  T zero;
  m_ops[a_level]->create(zero, a_ans);
  m_ops[a_level]->setToZero(zero);

  // set a_ans = helm(a_phi)
  //           = (I + factor*laplacian)(a_phi)
  if(a_applyBC)
    {
      if( (a_phiC == NULL) || (a_level==0))
        {
          m_ops[a_level]->applyOp(a_ans, a_phi, false);
        }
      else
        {
          m_ops[a_level]->AMROperatorNF(a_ans, a_phi, *a_phiC, false);
        }
    }
  else
    {
      m_ops[a_level]->applyOpNoBoundary(a_ans, a_phi);
    }
}

/***/
template <class T, class TFlux, class TFR>
void
BaseLevelTGA<T, TFlux, TFR>::
timeInterp(T&       a_data,
           const T& a_oldData,
           const T& a_newData,
           Real                        a_time,
           Real                        a_oldTime,
           Real                        a_newTime,
           int                         a_lev)
{
  Real eps = 1.0e-10;
  CH_assert(a_newTime >= a_oldTime);
  Real diff = (a_newTime - a_oldTime);
  m_ops[a_lev]->setToZero(a_data);
  if(diff < eps)
    {
      //no real time advance and don't want to divide by zero
      m_ops[a_lev]->incr(a_data, a_oldData, 1.0);
    }
  else
    {
      Real factor = (a_time-a_oldTime)/(a_newTime - a_oldTime);
      m_ops[a_lev]->incr(a_data, a_oldData, 1.0-factor);
      m_ops[a_lev]->incr(a_data, a_newData, factor);
    }
}

/***/
template <class T, class TFlux, class TFR>
void
BaseLevelTGA<T, TFlux, TFR>::
solveHelm(T&       a_phi,
          T&       a_phiC,
          T&       a_rhs,
          int                     a_level,
          Real                    a_mu,
          Real                    a_dt,
          bool                    a_zeroPhi)
{
  if (a_zeroPhi)
    {
      m_ops[a_level]->setToZero(a_phi);
    }
  Vector<T* > phi(m_grids.size(), NULL);
  Vector<T* > rhs(m_grids.size(), NULL);
  phi[a_level] = &a_phi;
  rhs[a_level] = &a_rhs;
  if(a_level > 0)
    {
      phi[a_level-1] = &a_phiC;
    }

  Real factor  = -a_dt*a_mu;
  resetSolverAlphaAndBeta(1.0, factor);

  m_solver->solve(phi, rhs, a_level, a_level, a_zeroPhi);
  int solverExitStatus = m_solver->m_exitStatus;
  if (solverExitStatus==2 || solverExitStatus==4 || solverExitStatus==6)
    {
      // These status codes correspond to the cases in which
      // norm is neither small enough nor reduced enough.
      // Either we've reached the maximum number of iterations,
      // or we've hung.
      pout() << "BaseLevelTGA:: WARNING: solver exitStatus == "
             << solverExitStatus << std::endl;
    }
}

template <class T, class TFlux, class TFR>
void
BaseLevelTGA<T, TFlux, TFR>::
incrementFlux(LevelData<TFlux>&      a_diffusiveFlux,
              T&      a_phi,
              int                    a_level,
              Real                   a_mu,
              Real                   a_dt,
              Real                   a_sign,
              bool                   a_setToZero)
{
  Real factor  = a_sign*a_dt*a_mu;
  m_ops[a_level]->setAlphaAndBeta(1.0, factor);

  // increment flux
  TFlux tempFlux;
  m_ops[a_level]->fillGrad(a_phi);
  for (DataIterator  dit = a_phi.dataIterator(); dit.ok(); ++dit)
    {
      TFlux& thisFlux = a_diffusiveFlux[dit];
      if(a_setToZero)
        {
          thisFlux.setVal(0.0);
        }

      m_ops[a_level]->getFlux(tempFlux, a_phi, m_grids[a_level][dit], dit(), 1.0);
      thisFlux += tempFlux;
    }
}

template <class T, class TFlux, class TFR>
void
BaseLevelTGA<T, TFlux, TFR>::
updateSoln(T&       a_phiNew,
           T&       a_phiOld,
           T&       a_src,
           LevelData<TFlux>&       a_flux,
           TFR*                    a_fineFluxRegPtr,
           TFR*                    a_crseFluxRegPtr,
           const T* a_crsePhiOldPtr,
           const T* a_crsePhiNewPtr,
           Real a_oldTime,
           Real a_crseOldTime,
           Real a_crseNewTime,
           Real a_dt,
           int a_level,
           bool a_zeroPhi,
           int a_fluxStartComponent)
{
  int ncomp = a_phiNew.nComp();
  Interval intervBase(0, ncomp-1);
  Interval intervFlux(a_fluxStartComponent, a_fluxStartComponent + ncomp-1);

  CH_assert(a_level >= 0);
  CH_assert(a_level <  m_grids.size());
  CH_assert((a_level == 0) || (a_crsePhiOldPtr != NULL));
  CH_assert((a_level == 0) || (a_crsePhiNewPtr != NULL));
  CH_assert(a_crseNewTime >= a_crseOldTime);
  CH_assert(a_dt >= 0.);

  T rhst, srct;
  m_ops[a_level]->create(rhst, a_src);
  m_ops[a_level]->create(srct, a_phiNew);
  m_ops[a_level]->setToZero(srct);
  m_ops[a_level]->setToZero(rhst);
  m_ops[a_level]->incr(srct, a_src, 1.0);
  T coarseData;

  if ((a_crsePhiOldPtr != NULL) && (a_level > 0))
    {
      m_ops[a_level-1]->create(coarseData, *a_crsePhiOldPtr);
      setSourceGhostCells(srct, m_grids[a_level]);
    }

  applyHelm(rhst, srct, NULL, a_level, m_mu4, a_dt, false); //'false', ie, use the extrap data

  m_ops[a_level]->scale(rhst, a_dt);
  incrementFlux(a_flux, srct, a_level, a_dt*m_mu4, a_dt, -1.0, true);

  // first need to compute coarse-level BC at this level's old time
  if (a_level > 0)
    {
      timeInterp(coarseData,
                 *a_crsePhiOldPtr,
                 *a_crsePhiNewPtr,
                 a_oldTime,
                 a_crseOldTime,
                 a_crseNewTime,
                 a_level-1);

    }


  applyHelm(a_phiNew, a_phiOld, &coarseData, a_level, m_mu3, a_dt, true); //'true' apply CF and domain BC

  incrementFlux(a_flux, a_phiOld, a_level, m_mu3, a_dt, -1., false);

  m_ops[a_level]->incr(rhst, a_phiNew, 1.0);


  // now construct coarse-level BC for intermediate solve
  // intermediate solution will be at time = oldTime + (1-r1)dt
  if(a_level > 0)
    {
      Real intermediateTime = a_oldTime + (1-m_r1)*a_dt;

      timeInterp(coarseData,
                 *a_crsePhiOldPtr,
                 *a_crsePhiNewPtr,
                 intermediateTime,
                 a_crseOldTime,
                 a_crseNewTime,
                 a_level-1);

    }

  solveHelm(a_phiNew, coarseData, rhst, a_level, m_mu2, a_dt, a_zeroPhi);

  incrementFlux(a_flux, a_phiNew, a_level, m_mu2, a_dt, -1.0, false);

  m_ops[a_level]->assign(rhst, a_phiNew);


  // now construct CF-BC for final solve
  if (a_level > 0)
    {
      Real newTime = a_oldTime + a_dt;
      timeInterp(coarseData,
                 *a_crsePhiOldPtr,
                 *a_crsePhiNewPtr,
                 newTime,
                 a_crseOldTime,
                 a_crseNewTime,
                 a_level-1);
    }



  solveHelm(a_phiNew, coarseData, rhst, a_level, m_mu1, a_dt, a_zeroPhi);

  incrementFlux(a_flux, a_phiNew, a_level, m_mu1, a_dt, -1.0, false);
  // now increment flux registers -- note that because of the way
  // we defined the fluxes, the dt multiplier is already in the
  // flux
  if ((a_fineFluxRegPtr != NULL) && (a_level > 0))
    {
      Real fluxMult = 1.0;
      for (DataIterator dit = m_grids[a_level].dataIterator(); dit.ok(); ++dit)
        {
          TFlux& thisFlux = a_flux[dit];
          for(int dir=0; dir<SpaceDim; ++dir)
            {
              a_fineFluxRegPtr->incrementCoarse(thisFlux[dir],
                                                fluxMult, dit(),
                                                intervBase, // source
                                                intervFlux, // dest
                                                dir);
            }
        }
    } // end if there is a finer-level

  if ((a_crseFluxRegPtr != NULL) && (a_level > 0))
    {
      Real fluxMult = 1.0;

      for (DataIterator dit = m_grids[a_level].dataIterator(); dit.ok(); ++dit)
        {
          TFlux& thisFlux = a_flux[dit];
          for(int dir=0; dir<SpaceDim; ++dir)
            {
              a_crseFluxRegPtr->incrementFine(thisFlux[dir],
                                              fluxMult, dit(),
                                              intervBase, // source
                                              intervFlux, // dest
                                              dir);
            }
        }
    } // end if there is a coarser level
}

// in this function, we call the updateSoln function, compute
// the update, then subtract off the extra pieces to return the
// diffusive part of the update
template <class T, class TFlux, class TFR>
void
BaseLevelTGA<T, TFlux, TFR>::
computeDiffusion(T&       a_diffusiveTerm,
                 T&       a_phiOld,
                 T&       a_src,
                 LevelData<TFlux>&       a_flux,
                 TFR*                    a_fineFluxRegPtr,
                 TFR*                    a_crseFluxRegPtr,
                 const T* a_crsePhiOldPtr,
                 const T* a_crsePhiNewPtr,
                 Real a_oldTime,
                 Real a_crseOldTime,
                 Real a_crseNewTime,
                 Real a_dt,
                 int a_level)
{
  // first compute updated solution
  T tempSoln;

  m_ops[a_level]->create(tempSoln, a_phiOld);
  m_ops[a_level]->setToZero(tempSoln);

  updateSoln(tempSoln, a_phiOld, a_src, a_flux,
             a_fineFluxRegPtr,
             a_crseFluxRegPtr,
             a_crsePhiOldPtr,
             a_crsePhiNewPtr,
             a_oldTime,
             a_crseOldTime,
             a_crseNewTime, a_dt, a_level);

  // now subtract everything off to leave us with diffusive term
  m_ops[a_level]->incr(tempSoln, a_phiOld, -1.0);
  m_ops[a_level]->scale(tempSoln, 1.0/a_dt);

  //now multiply by a if there is an a
  m_ops[a_level]->diagonalScale(tempSoln);

  // and finally, subtract off a_src
  m_ops[a_level]->incr(tempSoln, a_src, -1.0);

  // what's left should be the time-centered diffusive part of the update
  m_ops[a_level]->assign(a_diffusiveTerm, tempSoln);
}

#include "NamespaceFooter.H"
#endif

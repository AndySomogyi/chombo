#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

// NodeBiCGStabSmoother.cpp
// adapted from NodeBiCGStabSmoother by DFMartin, Sun, May 5, 2002
// petermc, 5 June 2002

#include "NodeBiCGStabSmoother.H"
#include "LayoutIterator.H"
#include "Norms.H"
#include "NodeDotProduct.H"
// using std::cout;
// using std::cerr;
using std::endl;

// ---------------------------------------------------------
NodeBiCGStabSmoother::NodeBiCGStabSmoother()
{
  m_maxIter = 23;
  m_solverTol = 1.0e-7;
  m_small = 1.0e-60;
  m_converge_small = 1.0e-2;
  // m_verbose = false;
  m_verbose = true;
}


// ---------------------------------------------------------
NodeBiCGStabSmoother::~NodeBiCGStabSmoother()
{
}


// ---------------------------------------------------------
void
NodeBiCGStabSmoother::setMaxIter(int a_max_iter)
{
  m_maxIter = a_max_iter;
}


// ---------------------------------------------------------
void
NodeBiCGStabSmoother::setSolverTol(Real a_solverTol)
{
  m_solverTol = a_solverTol;
}


// ---------------------------------------------------------
NodeBaseBottomSmoother*
NodeBiCGStabSmoother::new_bottomSmoother() const
{
  NodeBiCGStabSmoother* newsmoother = new NodeBiCGStabSmoother();
  newsmoother->setVerbose(m_verbose);

  if (newsmoother == NULL)
    {
      MayDay::Error("Out of Memory in NodeBiCGStabSmoother::new_bottomSmoother");
    }
  return static_cast<NodeBaseBottomSmoother*> (newsmoother);
}



// ---------------------------------------------------------
/***********************/
// This does a GSRB Pre/Conditioned BiCGStab on a level
// for the bottom solver.
/***********************/
void
NodeBiCGStabSmoother::doBottomSmooth(LevelData<NodeFArrayBox>& a_phi,
                                     const LevelData<NodeFArrayBox>& a_rhs,
                                     NodeLevelOp* a_levelopPtr)
{
  // Preconditioned BiConjugate Gradient Stabilized Method
  // Barrett et al., "Templates for the Solution of Linear Systems:
  // Building Blocks for Iterative Methods", Fig. 2.10 (p. 27).

  // For solve M * yhat = y, we set yhat = y and run gsrb(yhat, resid) twice.

  bool alreadyRestarted = false;
  Real small2 = 1.0e-8;
  int ncomp = a_rhs.nComp();

 CH_assert (ncomp == a_phi.nComp());

  const DisjointBoxLayout grids = a_rhs.getBoxes();
 CH_assert (grids == a_phi.getBoxes());

  LevelData<NodeFArrayBox> corremf(grids, ncomp, IntVect::Unit);
  LevelData<NodeFArrayBox> residmf(grids, ncomp, IntVect::Zero);
  LevelData<NodeFArrayBox> rtwidmf(grids, ncomp, IntVect::Zero);
  LevelData<NodeFArrayBox> shatmf(grids, ncomp, IntVect::Unit);
  LevelData<NodeFArrayBox> phatmf(grids, ncomp, IntVect::Unit);
  LevelData<NodeFArrayBox> pmf(grids, ncomp, IntVect::Zero);
  LevelData<NodeFArrayBox> vmf(grids, ncomp, IntVect::Zero);
  LevelData<NodeFArrayBox> tmf(grids, ncomp, IntVect::Zero);

  DataIterator dit = residmf.dataIterator();

  // Initialize corremf = 0, residmf = 0, vmf = 0.
  for (dit.begin(); dit.ok(); ++dit)
    {
      corremf[dit()].getFab().setVal(0.);
      residmf[dit()].getFab().setVal(0.);
      vmf[dit()].getFab().setVal(0.);

      // these are necessary to prevent uninitialized memory reads
      // petermc, 20 May 2002
      rtwidmf[dit()].getFab().setVal(0.);
      shatmf[dit()].getFab().setVal(0.);
      phatmf[dit()].getFab().setVal(0.);
      pmf[dit()].getFab().setVal(0.);
      tmf[dit()].getFab().setVal(0.);
    }

  // Set residmf = a_rhs - L(a_phi).
  a_levelopPtr->residualH(residmf, a_phi, a_rhs);

  // Added by petermc, 28 March 2003:  need these for norms, so that
  // you're taking norms over the valid nodes only.
  const ProblemDomain domain(a_levelopPtr->m_domain);
  const LayoutData< Vector<IntVectSet> >& IVSVext(a_levelopPtr->m_IVSVext);
  const Interval intvl = residmf.interval();

  // Real initnorm = maxnorm(residmf, residmf.interval(), m_verbose);
  Real initnorm = maxnorm(residmf, domain, IVSVext, intvl, m_verbose);
  Real oldnorm = initnorm;
  Real rnorm = initnorm;

  // Set rtwidmf = residmf.
  for (dit.begin(); dit.ok(); ++dit)
    rtwidmf[dit()].copy(residmf[dit()]);

  Vector<Real> rho;
  Vector<Real> alpha;
  // Vector<Real> beta;
  Vector<Real> omega;
  //dummy omega[0], alpha[0] (not used)
  omega.push_back(0.0);
  alpha.push_back(0.0);
  bool finished = false;
  for (int iter = 1; (iter <= m_maxIter && !finished); iter++)
    {
      // DotProductNodes() is in NodeFArrayBox.cpp
      // Real rhodot = DotProductNodes(residmf, rtwidmf, grids);
      Real rhodot = DotProductNodes(residmf, rtwidmf, domain, IVSVext, intvl);
      if (Abs(rhodot) < m_small)
        {
          if (m_verbose)
            pout() << "Returning from NodeBiCGStabSmoother::bottomSmoother because dot product too small; iteration " << iter << endl;
          for (dit.begin(); dit.ok(); ++dit)
            a_phi[dit()].getFab() += corremf[dit()].getFab();

          return;
        }
      rho.push_back(rhodot);

      if (iter == 1 || alreadyRestarted)
        {
          // Set pmf = residmf.
          for (dit.begin(); dit.ok(); ++dit)
            pmf[dit()].copy(residmf[dit()]);
        }
      else // iter >= 2
        {
          Real betaim1=(rho[iter-1]/rho[iter-2])*(alpha[iter-1]/omega[iter-1]);
          // beta.push_back(betaim1);
          Real omegaim1 = omega[iter-1];
          // Set pmf = residmf + betaim1 * (pmf - (omegaim1 * vmf)).
          //p(i) = r(i-1) + beta(i-1)*(p(i-1) -omega(i-1)*v(i-1))
          for (dit.begin(); dit.ok(); ++dit)
            {
              FArrayBox& pfab = pmf[dit()].getFab();
              FArrayBox& vfab = vmf[dit()].getFab();
              FArrayBox& rfab = residmf[dit()].getFab();

              Box fabbox(surroundingNodes(grids.get(dit())));
              FArrayBox temp(fabbox, ncomp);
              temp.copy(vfab);
              temp *= -omegaim1;
              temp += pfab;
              temp *= betaim1;
              temp += rfab;
              pfab.copy(temp);
            }
        }

      //solve M(phat) = p(i)
      a_levelopPtr->levelPreconditioner(phatmf, pmf);

      //v(i) = A(phat)
      // vmf = L(phatmf)
      a_levelopPtr->applyOpH(vmf, phatmf);
      //alpha(i) = rho(i-1)/(rtwid^T v(i))
      // Real denom = DotProductNodes(vmf, rtwidmf, grids);
      Real denom = DotProductNodes(vmf, rtwidmf, domain, IVSVext, intvl);
      if (Abs(denom) > small2 * rho[iter-1])
        {
          Real alphai = rho[iter-1] / denom;
          alpha.push_back(alphai);

          // increment residual and solution
          for (dit.begin(); dit.ok(); ++dit)
            {
              FArrayBox& rfab = residmf[dit()].getFab();
              FArrayBox& vfab = vmf[dit()].getFab();
              FArrayBox& phatfab = phatmf[dit()].getFab();
              FArrayBox& corrfab = corremf[dit()].getFab();

              rfab.plus(vfab, -alphai);
              corrfab.plus(phatfab, alphai);
            }
        }
      else
        {
          alpha.push_back(0.0);

          // if correction was bad, don't increment anything,
          // but set residual to zero (forces recomputation of resid)
          for (dit.begin(); dit.ok(); ++dit)
            residmf[dit()].getFab().setVal(0.0);
        }

      // compute norm of residual
      // Real norms = maxnorm(residmf, residmf.interval(), m_verbose);
      Real norms = maxnorm(residmf, domain, IVSVext, intvl, m_verbose);
      if (m_verbose)
        pout() << "NodeBiCGStabSmoother::doBottomSmooth() norms = "
               << norms << endl;

      // If we skip second update step, still need to store something
      // in omega.  try using 0
      Real omegai = 0.;

      //check norm of s.
      //if small enough, step to end (force recomputation of residual, etc.)
      // if not "done" here, do second correction step
      if (norms > m_solverTol * initnorm)
        {
          // solve Mshat = s
          a_levelopPtr->levelPreconditioner(shatmf, residmf);

          //t = A(shat)
          a_levelopPtr->applyOpH(tmf, shatmf);

          //w(i) = tTresid/tTt
          // Real numerw = DotProductNodes(tmf, residmf, grids);
          Real numerw = DotProductNodes(tmf, residmf, domain, IVSVext, intvl);
          // Real denomw = DotProductNodes(tmf, tmf, grids);
          Real denomw = DotProductNodes(tmf, tmf, domain, IVSVext, intvl);

          if (Abs(denomw) > m_small)
            omegai = numerw / denomw;

          //corr(i) = corr(i-1) + omega(i)*shat
          //resid(i) = resid(i) - omega(i)*t

          for (dit.begin(); dit.ok(); ++dit)
            {
              FArrayBox& corrfab = corremf[dit()].getFab();
              FArrayBox& shatfab = shatmf[dit()].getFab();
              FArrayBox& tfab = tmf[dit()].getFab();
              FArrayBox& residfab = residmf[dit()].getFab();

              // increment correction: corr = corr+omega*shat
              corrfab.plus(shatfab, omegai);
              residfab.plus(tfab, -omegai);
            } // end loop over grids

          // reset this just in case (if we made it this
          // far, then we can assume that we're not caught in a
          // loop
          // petermc, 29 May 2002, actually don't reset this;
          // allow only one restart
          // alreadyRestarted = false;
        } // end second update loop

      omega.push_back(omegai);

      // now check residual.
      // rnorm = maxnorm(residmf, residmf.interval(), m_verbose);
      rnorm = maxnorm(residmf, domain, IVSVext, intvl, m_verbose);
      if (m_verbose && (rnorm > 0))
        {
          pout() << "NodeBiCGStabSmoother::doBottomSmooth():  after "
                 << iter << " iterations, residnorm = "
                 << rnorm << endl;
        }

      // if residnorm is small enough (or if correction is
      // small enough), then recompute residual from scratch
      // and check again (to avoid roundoff-drift issues).  if
      // we _still_ think we're done, then exit
      if ((rnorm <= m_solverTol*initnorm) ||
          (Abs(omegai) <= m_solverTol) ||
          (rnorm > oldnorm * (1.0 - m_converge_small)))
        {
          // recompute residual from scratch -- first increment
          // a_phi with correction
          for (dit.begin(); dit.ok(); ++dit)
            {
              FArrayBox& phifab = a_phi[dit()].getFab();
              FArrayBox& corrfab = corremf[dit()].getFab();

              phifab += corrfab;
              corrfab.setVal(0.0);
            }

          // Set residmf = a_rhs - L(a_phi).
          a_levelopPtr->residualH(residmf, a_phi, a_rhs);

          // now recompute residual norm and check for convergence
          oldnorm = rnorm;
          // rnorm = maxnorm(residmf, residmf.interval(), m_verbose);
          rnorm = maxnorm(residmf, domain, IVSVext, intvl, m_verbose);

          // if we've already restarted and haven't gotten
          // anywhere, then we've stalled and should probably
          // exit.  also, if newly-computed residual still indicates
          // convergence, then we're done.
          if (alreadyRestarted || (rnorm <= m_solverTol*initnorm))
            {
              finished = true;
            }
          else
            {
              // otherwise, go back and try again with new residual
              if (m_verbose)
                pout() << "NodeBiCGStabSmoother::doBottomSmooth() -- restarting iterations" << endl;
              // also reset rtwid
              for (dit.begin(); dit.ok(); ++dit)
                {
                  rtwidmf[dit()].copy(residmf[dit()]);
                }
              // try doing this recursively?

              finished = false;
              alreadyRestarted = true;
            }

          if (finished && m_verbose)
            pout() << "NodeBiCGStabSmoother::doBottomSmooth() final Norm = "
                   << rnorm << endl;
        } // end if we're recomputing the residual
    } // end main loop over BiCGStab iterations
}

#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

// NodeCGSmoother.cpp
// adapted from NodeCGSmoother by DFMartin, Sun, May 5, 2002
// petermc, 5 June 2002

#include "NodeCGSmoother.H"
#include "LayoutIterator.H"
#include "Norms.H"
#include "NodeDotProduct.H"
// using std::cout;
// using std::cerr;
using std::endl;

// ---------------------------------------------------------
NodeCGSmoother::NodeCGSmoother()
{
  m_maxIter = 23;
  m_solverTol = 1.0e-7;
  m_small = 1.0e-60;
  m_converge_small = 1.0e-2;
  // m_verbose = false;
  m_verbose = true;
}


// ---------------------------------------------------------
NodeCGSmoother::~NodeCGSmoother()
{
}


// ---------------------------------------------------------
void
NodeCGSmoother::setMaxIter(int a_max_iter)
{
  m_maxIter = a_max_iter;
}


// ---------------------------------------------------------
void
NodeCGSmoother::setSolverTol(Real a_solverTol)
{
  m_solverTol = a_solverTol;
}


// ---------------------------------------------------------
NodeBaseBottomSmoother*
NodeCGSmoother::new_bottomSmoother() const
{
  NodeCGSmoother* newsmoother = new NodeCGSmoother();
  newsmoother->setVerbose(m_verbose);

  if (newsmoother == NULL)
    {
      MayDay::Error("Out of Memory in NodeCGSmoother::new_bottomSmoother");
    }
  return static_cast<NodeBaseBottomSmoother*> (newsmoother);
}



// ---------------------------------------------------------
/***********************/
// This does a GSRB Pre/Conditioned CG on a level
// for the bottom solver.
/***********************/
void
NodeCGSmoother::doBottomSmooth(LevelData<NodeFArrayBox>& a_phi,
                               const LevelData<NodeFArrayBox>& a_rhs,
                               NodeLevelOp* a_levelopPtr)
{
  // Preconditioned BiConjugate Gradient Stabilized Method
  // Barrett et al., "Templates for the Solution of Linear Systems:
  // Building Blocks for Iterative Methods", Fig. 2.10 (p. 15).

  // For solve M * yhat = y, we set yhat = y and run gsrb(yhat, resid) twice.

  // We assume the operator and preconditioner are both NEGATIVE definite.

  bool alreadyRestarted = false;
  Real small2 = 1.0e-8;
  int ncomp = a_rhs.nComp();

 CH_assert (ncomp == a_phi.nComp());

  const DisjointBoxLayout grids = a_rhs.getBoxes();
 CH_assert (grids == a_phi.getBoxes());

  LevelData<NodeFArrayBox> residmf(grids, ncomp, IntVect::Zero);
  LevelData<NodeFArrayBox> zmf(grids, ncomp, IntVect::Unit);
  LevelData<NodeFArrayBox> pmf(grids, ncomp, IntVect::Unit);
  LevelData<NodeFArrayBox> qmf(grids, ncomp, IntVect::Zero);

  DataIterator dit = residmf.dataIterator();

  // Initialize everything to zero to prevent uninitialized memory reads.
  for (dit.begin(); dit.ok(); ++dit)
    {
      residmf[dit()].getFab().setVal(0.);
      zmf[dit()].getFab().setVal(0.);
      pmf[dit()].getFab().setVal(0.);
      qmf[dit()].getFab().setVal(0.);
    }

  // Set residmf = a_rhs - L(a_phi).
  a_levelopPtr->residualH(residmf, a_phi, a_rhs);
  // added 18 Jun 2002
  // for (dit.begin(); dit.ok(); ++dit)
  // residmf[dit()].getFab().negate();

  // Added by petermc, 28 March 2003:  need these for norms, so that
  // you're taking norms over the valid nodes only.
  const ProblemDomain domain(a_levelopPtr->m_domain);
  const LayoutData< Vector<IntVectSet> >& IVSVext(a_levelopPtr->m_IVSVext);
  const Interval intvl = residmf.interval();

  // Real initnorm = maxnorm(residmf, residmf.interval(), m_verbose);
  Real initnorm = maxnorm(residmf, domain, IVSVext, intvl, m_verbose);
  Real rnorm = initnorm;
  Real oldnorm = initnorm;

  Vector<Real> rho;
  bool finished = false;
  for (int iter = 1; (iter <= m_maxIter && !finished); iter++)
    {
      //solve M(z(i-1)) = r(i-1)
      a_levelopPtr->levelPreconditioner(zmf, residmf);
      // added 18 Jun 2002
      // for (dit.begin(); dit.ok(); ++dit)
      // zmf[dit()].getFab().negate();

      // DotProductNodes() is in NodeFArrayBox.cpp
      // Real rhodot = DotProductNodes(residmf, zmf, grids);
      Real rhodot = DotProductNodes(residmf, zmf, domain, IVSVext, intvl);
      // added 18 Jun 2002
      if (rhodot > 0.)
        pout() << "Warning: in NodeCGSmoother::doBottomSmooth(): levelPreconditioner is not negative definite "
               << "at iteration " << iter << " with rhodot=" << rhodot << endl;

      Real abs_rhodot = Abs( rhodot );
      if (abs_rhodot == 0.0) return ; //[NOTE: added by <dbs> 6/03]
      if (abs_rhodot < m_small)
        {
          pout() << "Warning: in NodeCGSmoother::doBottomSmooth(): CG failed at iteration "
                 << iter << ", rhodot=" << rhodot << " is less than threshold " << m_small << endl;
          // for (dit.begin(); dit.ok(); ++dit)
          // a_phi[dit()].getFab() += corremf[dit()].getFab();
          return;
        }
      rho.push_back(rhodot);

      if (iter == 1 || alreadyRestarted)
        {
          // Set pmf = zmf.
          for (dit.begin(); dit.ok(); ++dit)
            pmf[dit()].copy(zmf[dit()]);
        }
      else // iter >= 2
        {
          Real betaim1 = rho[iter-1] / rho[iter-2];
          // Set pmf = zmf + betaim1 * pmf
          // p(i) = z(i-1) + beta(i-1)*p(i-1)
          for (dit.begin(); dit.ok(); ++dit)
            {
              FArrayBox& pfab = pmf[dit()].getFab();
              FArrayBox& zfab = zmf[dit()].getFab();

              pfab *= betaim1;
              pfab += zfab;
            }
        }

      // q(i) = A(p(i))
      // qmf = L(pmf)
      a_levelopPtr->applyOpH(qmf, pmf);
      // added 18 Jun 2002
      // for (dit.begin(); dit.ok(); ++dit)
      // qmf[dit()].getFab().negate();

      //alpha(i) = rho(i-1)/(p(i)^T q(i))
      // Real denom = DotProductNodes(pmf, qmf, grids);
      Real denom = DotProductNodes(pmf, qmf, domain, IVSVext, intvl);
      // added 18 Jun 2002
      if (denom > 0.)
        pout() << "Warning: in NodeCGSmoother::doBottomSmooth(): operator is not negative definite "
               << "at iteration " << iter << " with denom=" << denom << endl;

      // this condition is equivalent to |alpha(i)| < 1/small2
      if (Abs(denom) > small2 * Abs(rho[iter-1]))
        {
          Real alphai = rho[iter-1] / denom;

          // increment residual and solution
          for (dit.begin(); dit.ok(); ++dit)
            {
              FArrayBox& phifab = a_phi[dit()].getFab();
              FArrayBox& pfab = pmf[dit()].getFab();
              FArrayBox& qfab = qmf[dit()].getFab();
              FArrayBox& rfab = residmf[dit()].getFab();

              phifab.plus(pfab, alphai);
              rfab.plus(qfab, -alphai);
            }
        }
      else
        {
          // if correction was bad, don't increment anything,
          // but set residual to zero (forces recomputation of resid)
          for (dit.begin(); dit.ok(); ++dit)
            residmf[dit()].getFab().setVal(0.0);
        }

      // now check residual.
      oldnorm = rnorm;
      // rnorm = maxnorm(residmf, residmf.interval(), m_verbose);
      rnorm = maxnorm(residmf, domain, IVSVext, intvl, m_verbose);
      if (m_verbose && (rnorm > 0))
        {
          pout() << "Info: in NodeCGSmoother::doBottomSmooth(): after "
               << iter << " iterations, residnorm = "
               << rnorm << endl;
        }

      // if residnorm is small enough (or if correction is
      // small enough), then recompute residual from scratch
      // and check again (to avoid roundoff-drift issues).  if
      // we _still_ think we're done, then exit
      if (rnorm <= m_solverTol*initnorm ||
          rnorm > oldnorm * (1.0 - m_converge_small) )
        {
          // recompute residual from scratch
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
                pout() << "Info: in NodeCGSmoother::doBottomSmooth(): restarting iterations" << endl;
              // try doing this recursively?

              finished = false;
              alreadyRestarted = true;
            }
        }

      if (m_verbose)
        {
          if (finished)
            pout() << "Info: in NodeCGSmoother::doBottomSmooth(): final Norm = "
                   << rnorm << endl;
          else if (iter == m_maxIter)
            pout() << "Warning: in NodeCGSmoother::doBottomSmooth(): not converged after " << iter
                   << " iterations, resid/initRes = " << rnorm/initnorm
                   << endl;
        }
    } // end main loop over CG iterations
}

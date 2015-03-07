#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

#include "MACProjectorF_F.H"
#include "LevelMACProjector.H"
#include "computeSum.H"
#include "BoxIterator.H"

/*****/
void
LevelMACProjector::
project(LevelData<FluxBox>&         a_velocity,
        LevelData<FluxBox>&         a_gradPhi,
        LevelData<FArrayBox>&       a_phi,
        const LevelData<FArrayBox>* a_phiCoarse)
{
  CH_assert(a_velocity.nComp() == 1);
  CH_assert(a_gradPhi.nComp() == 1);

  macEnforceVelocityBC( a_velocity,  m_grids,  m_domain, m_dx, m_bcVel);

  macDivergence(m_rhs  , a_velocity, m_grids,  m_domain, m_dx);

#ifndef NDEBUG
  Real sumRHS = computeSum(m_rhs, NULL, 2, m_dx);
  pout() << "LevelMACProjector::project sum of rhs = " << sumRHS << endl;
#endif

  //solve kappa lapl phi = kappa div(u).
  solve(a_phi, m_rhs, a_phiCoarse);

  //compute gradient and subtract off of velocity.   at domain boundaries gradient is set to zero
  //so boundary values are not changed.  exchange on phi done internally
  macGradient(          a_gradPhi,    a_phi, m_grids,  m_domain, m_dx);
  macEnforceGradientBC( a_gradPhi,    a_phi, m_grids,  m_domain, m_dx, m_bcGrad);

  for(DataIterator dit = m_grids.dataIterator(); dit.ok(); ++dit)
    {
      for(int idir = 0; idir < SpaceDim; idir++)
        {
          a_velocity[dit()][idir] -= a_gradPhi[dit()][idir];
        }
    }

  //enforce boundary conditons again
  macEnforceVelocityBC( a_velocity,  m_grids,  m_domain, m_dx, m_bcVel);
}
void
LevelMACProjector::
solve(LevelData<FArrayBox>&            a_phi,
      const LevelData<FArrayBox>&      a_rhs,
      const LevelData<FArrayBox>*      a_phiCoarse)
{
  Vector<LevelData<FArrayBox> *> rhsHier(m_numLevels, NULL), phiHier(m_numLevels, NULL);
  rhsHier[m_lbase] = &((LevelData<FArrayBox>&)(a_rhs));
  phiHier[m_lbase] = &a_phi;
  if(m_lbase > 0)
    {
      CH_assert(a_phiCoarse != NULL);
      phiHier[m_lbase-1] = (LevelData<FArrayBox>*)(a_phiCoarse);
    }
  for(DataIterator dit = a_phi.dataIterator(); dit.ok(); ++dit)
    {
      a_phi[dit()].setVal(0.);
    }
  m_solver->solve(phiHier, rhsHier, m_lbase, m_lbase);
}

/*****/
void
macGradient(LevelData<FluxBox>&          a_gradPhi,
            const LevelData<FArrayBox>&  a_phi,
            const DisjointBoxLayout&     a_grids,
            const ProblemDomain&         a_domain,
            const Real&                  a_dx)
{
  Interval interv(0,0);
  //cast to non-const so we can do exchange.
  //none of the valid data is changed so const-nitude
  //is preserved at least in spirit.
  LevelData<FArrayBox>& phi = (LevelData<FArrayBox>&) a_phi;
  phi.exchange(interv);
  for(DataIterator dit = a_grids.dataIterator(); dit.ok(); ++dit)
    {

      macGradient(a_gradPhi[dit()], a_phi[dit()],
                  a_grids.get(dit()),a_domain, a_dx);

    }
}
/****/
void
macGradient(FluxBox&             a_gradPhi,
            const FArrayBox&     a_phi,
            const Box&           a_box,
            const ProblemDomain& a_domain,
            const Real&          a_dx)
{
  a_gradPhi.setVal(0.);
  for(int idir = 0; idir < SpaceDim; idir++)
    {
      //do interior faces.  domain boundary faces are boundary conditions
      Box interiorBox = a_box;
      interiorBox.grow(idir, 1);
      interiorBox &= a_domain;
      interiorBox.grow(idir, -1);

      Box faceBox = surroundingNodes(interiorBox, idir);

      FORT_MACGRADPHI(CHF_FRA1(a_gradPhi[idir], 0),
                      CHF_CONST_FRA1(a_phi, 0),
                      CHF_CONST_INT(idir),
                      CHF_CONST_REAL(a_dx),
                      CHF_BOX(faceBox));

    }
}

/*****/
void
macEnforceGradientBC(LevelData<FluxBox>&          a_gradPhi,
                     const LevelData<FArrayBox>&  a_phi,
                     const DisjointBoxLayout&     a_grids,
                     const ProblemDomain&         a_domain,
                     const Real&                  a_dx,
                     GradientBCFunction           a_bcFunc)
{
  Interval interv(0,0);
  //cast to non-const so we can do exchange.
  //none of the valid data is changed so const-nitude
  //is preserved at least in spirit.
  LevelData<FArrayBox>& phi = (LevelData<FArrayBox>&) a_phi;
  phi.exchange(interv);
  for(DataIterator dit = a_grids.dataIterator(); dit.ok(); ++dit)
    {
      macEnforceGradientBC(a_gradPhi[dit()],
                           a_phi[dit()],
                           a_grids.get(dit()),
                           a_domain, a_dx, a_bcFunc);
    }
}
/****/
void
macEnforceGradientBC(FluxBox&                    a_gradPhi,
                     const FArrayBox&            a_phi,
                     const Box   &               a_box,
                     const ProblemDomain&        a_domain,
                     const Real&                 a_dx,
                     GradientBCFunction          a_bcFunc)
{
  for(int idir = 0; idir < SpaceDim; idir++)
    {
      if(!a_domain.isPeriodic(idir))
        {
          for(SideIterator sit; sit.ok(); ++sit)
            {
              //check if domain and grid have the same bounds on sit() side of idir
              IntVect domBounds = a_domain.domainBox().sideEnd(sit());
              IntVect boxBounds =                a_box.sideEnd(sit());

              if(domBounds[idir] == boxBounds[idir])
                {
                  //get cells on the other side of the domain
                  Box bFace = adjCellBox(a_box, idir, sit(), 1);
                  //shift to the domain faces
                  bFace.shiftHalf(idir, -sign(sit()));
                  for(BoxIterator bit(bFace); bit.ok(); ++bit)
                    {
                      a_gradPhi[idir](bit(), 0) = a_bcFunc(bit(), a_dx, a_phi, idir, sit());
                    }
                } // end if(we are on domain)
            } //end loop over sides
        }//end check if !periodic
    } //end loop over directions
}
/*****/
void
macDivergence(LevelData<FArrayBox>&        a_divVel,
              const LevelData<FluxBox>&    a_velocity,
              const DisjointBoxLayout&     a_grids,
              const ProblemDomain&         a_domain,
              const Real&                  a_dx)
{
  Interval interv(0,0);
  LevelData<FluxBox>& vel = (LevelData<FluxBox>&) a_velocity;
  vel.exchange(interv);
  for(DataIterator dit = a_grids.dataIterator(); dit.ok(); ++dit)
    {
      macDivergence(a_divVel[dit()],
                    a_velocity[dit()],
                    a_grids.get(dit()),
                    a_domain, a_dx);

    }
}
/*****/
void
macDivergence(FArrayBox&             a_divVel,
              const FluxBox&         a_velo,
              const Box&             a_box,
              const ProblemDomain&   a_domain,
              const Real&            a_dx)
{
  //set the divergence initially to zero
  //then loop through directions and increment the divergence
  //with each directions flux difference.
  a_divVel.setVal(0.0);
  int ncons = a_divVel.nComp();
  for(int idir = 0; idir < SpaceDim; idir++)
    {
      FORT_MACDIVERGEF( CHF_BOX(a_box),
                        CHF_FRA(a_divVel),
                        CHF_CONST_FRA(a_velo[idir]),
                        CHF_CONST_INT(idir),
                        CHF_CONST_INT(ncons),
                        CHF_CONST_REAL(a_dx));
    }
}
/*****/
void
macEnforceVelocityBC(LevelData<FluxBox>&                a_velocity,
                     const DisjointBoxLayout&           a_grids,
                     const ProblemDomain&               a_domain,
                     const Real&                        a_dx,
                     VelocityBCFunction                 a_bcfunc)
{
  for(DataIterator dit = a_grids.dataIterator(); dit.ok(); ++dit)
    {
      macEnforceVelocityBC(a_velocity[dit()],
                           a_grids.get(dit()),
                           a_domain, a_dx, a_bcfunc);
    }
}
/*****/
void
macEnforceVelocityBC(FluxBox&                a_velocity,
                     const Box&              a_box,
                     const ProblemDomain&    a_domain,
                     const Real&             a_dx,
                     VelocityBCFunction      a_bcfunc)
{
  for(int idir = 0; idir < SpaceDim; idir++)
    {
      if(!a_domain.isPeriodic(idir))
        {
          for(SideIterator sit; sit.ok(); ++sit)
            {
              //check if domain and grid have the same bounds on sit() side of idir
              IntVect domBounds = a_domain.domainBox().sideEnd(sit());
              IntVect boxBounds =                a_box.sideEnd(sit());

              if(domBounds[idir] == boxBounds[idir])
                {
                  //get cells on the other side of the domain
                  Box bFace = adjCellBox(a_box, idir, sit(), 1);
                  //shift to the domain faces
                  bFace.shiftHalf(idir, -sign(sit()));
                  for(BoxIterator bit(bFace); bit.ok(); ++bit)
                    {
                      a_velocity[idir](bit(), 0) = a_bcfunc(bit(), a_dx, a_velocity, idir, sit());
                    }
                } // end if(we are on domain)
            } //end loop over sides
        }//end check if !periodic
    } //end loop over directions
}

/*****/

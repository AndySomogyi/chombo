#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "LevelCCProjector.H"
#include "QuadCFInterp.H"
#include "CCProjectorF_F.H"

/***/
void
LevelCCProjector::
divergence(LevelData<FArrayBox>&       a_divergence,
           LevelData<FArrayBox>&       a_velocity,
           LevelData<FArrayBox>*       a_velCoarse)
{
  if(m_lbase > 0)
    {
      CH_assert(a_velCoarse != NULL);
      DisjointBoxLayout coarGrids = a_velCoarse->disjointBoxLayout();
      QuadCFInterp interpolator(m_grids, &coarGrids, m_dx, m_refToCoar, SpaceDim, m_domain);
      
      interpolator.coarseFineInterp(a_velocity, *a_velCoarse);
    }

  a_velocity.exchange(Interval(0, SpaceDim-1));

  //average velocity to faces
  averageCellsToFaces(m_macVel, a_velocity, m_grids, m_domain);

  macDivergence(a_divergence, m_macVel,     m_grids, m_domain, m_dx);
}

void
LevelCCProjector::
project(LevelData<FArrayBox>&       a_velocity,
        LevelData<FArrayBox>&       a_gradient,
        LevelData<FArrayBox>&       a_phi,
        const LevelData<FArrayBox>* a_velCoarse,
        const LevelData<FArrayBox>* a_phiCoarse)
{
  if(m_lbase > 0)
    {
      CH_assert(a_velCoarse != NULL);
      DisjointBoxLayout coarGrids = a_phiCoarse->disjointBoxLayout();
      QuadCFInterp interpolator(m_grids, &coarGrids, m_dx, m_refToCoar, SpaceDim, m_domain);
      
      interpolator.coarseFineInterp(a_velocity, *a_velCoarse);
    }

  a_velocity.exchange(Interval(0, SpaceDim-1));
  //average velocity to faces
  averageCellsToFaces(m_macVel, a_velocity, m_grids, m_domain);

  //use mac projection to make divergence-free
  m_macProj->project(m_macVel, m_macGph, a_phi, a_phiCoarse);

  // extrapolate gradient to domain faces so that
  // we get reasonable answers near domain boundary
  extrapolateToDomainBoundaries(  m_macGph, m_grids, m_domain);
  
  //put the gradient into cell centers
  averageFacesToCells(a_gradient, m_macGph, m_grids, m_domain);

  //subtract cell-centered gradient off input vel
  for(DataIterator dit  = m_grids.dataIterator(); dit.ok(); ++dit)
    {
      a_velocity[dit()] -= a_gradient[dit()];
    }
}
/***/
void
averageCellsToFaces(LevelData<FluxBox>&           a_macVeloc,
                    const LevelData<FArrayBox>&   a_cellVeloc,
                    const DisjointBoxLayout&      a_grids,
                    const ProblemDomain&          a_domain)
{
  for(DataIterator dit = a_grids.dataIterator(); dit.ok(); ++dit)
    {
      averageCellsToFaces(a_macVeloc [dit()],
                          a_cellVeloc[dit()],
                          a_grids.get(dit()),
                          a_domain);
    }
}
/***/
void
averageCellsToFaces(FluxBox&               a_macVel,
                    const FArrayBox&       a_cellVel,
                    const Box&             a_grid,
                    const ProblemDomain&   a_domain)
{
  for(int faceDir = 0; faceDir < SpaceDim; faceDir++)
    {
      //only want non-domain boundary faces
      Box faceBox = a_grid;
      faceBox.grow(faceDir, 1);
      faceBox  &= a_domain;
      faceBox.grow(faceDir, -1);
      faceBox.surroundingNodes(faceDir);
      FORT_CCPAVECELLTOFACE(CHF_FRA1(a_macVel[faceDir], 0),
                            CHF_CONST_FRA1(a_cellVel, faceDir),
                            CHF_CONST_INT(faceDir),
                            CHF_BOX(faceBox));
      
    }
}
/***/
void
averageFacesToCells(LevelData<FArrayBox>&         a_cellData,
                    const LevelData<FluxBox>&     a_macData,
                    const DisjointBoxLayout&      a_grids,
                    const ProblemDomain&          a_domain)
{
  for(DataIterator dit = a_grids.dataIterator(); dit.ok(); ++dit)
    {
      averageFacesToCells(a_cellData [dit()],
                          a_macData  [dit()],
                          a_grids.get(dit()),
                          a_domain);
    }
}
/***/
void
averageFacesToCells(FArrayBox&             a_cellData,
                    const FluxBox&         a_fluxData,
                    const Box&             a_grid,
                    const ProblemDomain&   a_domain)
{
  for(int faceDir = 0; faceDir < SpaceDim; faceDir++)
    {
      FORT_CCPAVEFACETOCELL( CHF_FRA1(      a_cellData,faceDir),
                             CHF_CONST_FRA1(a_fluxData[faceDir], 0),
                             CHF_CONST_INT(faceDir),
                             CHF_BOX(a_grid));
    }
}
/***/
void
extrapolateToDomainBoundaries(LevelData<FluxBox>&      a_macData,
                              const DisjointBoxLayout& a_grids,
                              const ProblemDomain&     a_domain)
{
  for(DataIterator dit = a_grids.dataIterator(); dit.ok(); ++dit)
    {
      for(int faceDir = 0; faceDir < SpaceDim; faceDir++)
        {
          for(SideIterator sit; sit.ok(); ++sit)
            {
              extrapolateToDomainBoundaries(a_macData  [dit()],
                                            a_grids.get(dit()),
                                            a_domain, faceDir, sit());
            }
        }
    }
}
/***/
void
extrapolateToDomainBoundaries(FluxBox&               a_faceData,
                              const Box&             a_grid,
                              const ProblemDomain&   a_domain,
                              int                    a_faceDir,
                              Side::LoHiSide         a_side)
{
  IntVect ivSideGrid =   a_grid.sideEnd(a_side);
  IntVect ivSideDom  = a_domain.domainBox().sideEnd(a_side);
  if(ivSideGrid[a_faceDir] == ivSideDom[a_faceDir])
    {
      //create face-centered box along domain boundary
      //and interior to domain
      Box sideBox = adjCellBox(a_grid, a_faceDir, a_side, 1);
      int ishift = -sign(a_side);
      sideBox.shiftHalf(a_faceDir, ishift);
      CH_assert(a_grid.size(a_faceDir) > 3);
      FORT_CCPEXTRAPTODOMFACE( CHF_FRA1(a_faceData[a_faceDir], 0),
                               CHF_CONST_INT(a_faceDir),
                               CHF_CONST_INT(ishift),
                               CHF_BOX(sideBox));
    }
}
/***/




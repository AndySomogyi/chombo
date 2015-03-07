#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

#include "LayoutIterator.H"
#include "ComputeEnergy.H"
#include "ComputeEnergyF_F.H"


////////////////////////////////////////////////////////////////

Real computeEnergy(const Vector<LevelData<FArrayBox>* >& a_phi,
                   const Vector<LevelData<FArrayBox>* >& a_pi,
                   const Vector<int>&                    a_nRefFine,
                   Real                                  a_dxCrse,
                   Real                                  a_lambda)
{
  int numLevels = a_phi.size();
  Real sum, sumLevel;

  sum = 0.0;

  // loop over levels
  for (int lev = 0; lev < numLevels; lev++) {
    //in case there are extra levels which are not defined
    if (a_phi[lev] != NULL)
    {
      CH_assert(a_phi[lev]->isDefined());
      CH_assert(a_pi [lev]->isDefined());

      const DisjointBoxLayout* finerGridsPtr;

      if (lev < numLevels-1)
      {
        finerGridsPtr = &(a_phi[lev+1]->getBoxes());
      }
      else
      {
        finerGridsPtr = NULL;
      }
      // compute on current level with zeros in covered cells
      sumLevel = computeEnergy(*(a_phi[lev]),*(a_pi[lev]),finerGridsPtr,a_nRefFine[lev],
                               a_dxCrse, a_lambda);
      sum += sumLevel;

      a_dxCrse = a_dxCrse/a_nRefFine[lev];
    }
  }

  return sum;
}

////////////////////////////////////////////////////////////////

Real computeEnergy(const LevelData<FArrayBox>& a_phi,
                   const LevelData<FArrayBox>& a_pi,
                   const DisjointBoxLayout*    a_finerGridsPtr,
                   int                         a_nRefFine,
                   Real                        a_dx,
                   Real                        a_lambda)
{
  Real sum, sumLevel;
  sum = 0.0;

  const DisjointBoxLayout& levelGrids = a_phi.getBoxes();

  // compute the energy in box, overwrite the cells
  // that are covered by finer levels with zero, then sum
  for (DataIterator dit = a_phi.dataIterator(); dit.ok(); ++dit)
  {
    FArrayBox energy(levelGrids.get(dit()),1);

    FORT_COMPUTEENERGYINBOX(CHF_FRA1(energy,0),
                            CHF_CONST_FRA(a_phi[dit()]),
                            CHF_CONST_FRA(a_pi[dit()]),
                            CHF_BOX(levelGrids.get(dit())),
                            CHF_CONST_REAL(a_lambda), CHF_CONST_REAL(a_dx));

    if (a_finerGridsPtr != NULL)
    {
      // now loop over fine boxes and set covered regions to 0
      for (LayoutIterator litFine = a_finerGridsPtr->layoutIterator(); litFine.ok(); ++litFine)
      {
        Box coveredBox(a_finerGridsPtr->get(litFine()));
        coveredBox.coarsen(a_nRefFine);
        coveredBox &= levelGrids.get(dit());
        if (!coveredBox.isEmpty())
        {
          energy.setVal(0.0,coveredBox,0);
        }
      }  // end loop over fine-grid boxes
    } // end if there is a finer level

    sumLevel = energy.sum(0); //[NOTE: redundant var to simplify debugging]
    sum += sumLevel;
  } // end loop over this level's grids


  // do broadcast/gather thing here
#ifdef CH_MPI
  Real recv;
  int result = MPI_Allreduce(&sum, &recv, 1, MPI_CH_REAL,
                             MPI_SUM, Chombo_MPI::comm);
  if (result != MPI_SUCCESS)
    {
      MayDay::Error("Sorry, but I had a communication error in computeSum");
    }

  sum = recv;
#endif
  return sum;
}

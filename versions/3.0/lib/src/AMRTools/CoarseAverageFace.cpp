#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

#include "CoarseAverageFace.H"
#include "AverageFaceF_F.H"
#include "NamespaceHeader.H"

// ----------------------------------------------------------
CoarseAverageFace::CoarseAverageFace() :m_isDefined(false)
{
}


// ----------------------------------------------------------
CoarseAverageFace::~CoarseAverageFace()
{
}

// ----------------------------------------------------------
CoarseAverageFace::CoarseAverageFace(const DisjointBoxLayout& a_fineGrids,
                                     int a_nComp, int a_nRef)
  : m_isDefined(false)
{
  define(a_fineGrids, a_nComp, a_nRef);
}

// ----------------------------------------------------------
void
CoarseAverageFace::define(const DisjointBoxLayout& a_fineGrids,
                          int a_nComp, int a_nRef)
{
  m_nRef = a_nRef;

  DisjointBoxLayout coarsened_fine_domain;
  coarsen(coarsened_fine_domain, a_fineGrids, m_nRef);
  m_coarsenedFineData.define(coarsened_fine_domain, a_nComp);

  m_isDefined = true;
}

// ----------------------------------------------------------
bool
CoarseAverageFace::isDefined() const
{
  return m_isDefined;
}

// ----------------------------------------------------------
// this function is shamelessly based on the ANAG CoarseAverage
// (cell-centered) version
void
CoarseAverageFace::averageToCoarse(LevelData<FluxBox>& a_coarseData,
                                   const LevelData<FluxBox>& a_fineData)
{
  computeAverages(a_coarseData, a_fineData, arithmetic);
}


// ----------------------------------------------------------
// this function is shamelessly based on the ANAG CoarseAverage
// (cell-centered) version
void
CoarseAverageFace::averageToCoarseHarmonic(LevelData<FluxBox>& a_coarseData,
                                           const LevelData<FluxBox>& a_fineData)
{
  computeAverages(a_coarseData, a_fineData, harmonic);
}


void
CoarseAverageFace::computeAverages(LevelData<FluxBox>& a_coarseData,
                                   const LevelData<FluxBox>& a_fineData,
                                   int a_averageType)
{

  CH_assert(isDefined());

  DataIterator dit = a_fineData.dataIterator();
  for (dit.reset(); dit.ok(); ++dit)
    {
      FluxBox& coarsenedFine = m_coarsenedFineData[dit()];
      const FluxBox& fine = a_fineData[dit()];

      // coarsen from the entire fine grid onto the entire coarse grid
      averageGridData(coarsenedFine, fine, a_averageType);
    }

  m_coarsenedFineData.copyTo(m_coarsenedFineData.interval(),
                             a_coarseData, a_coarseData.interval());

}


// ----------------------------------------------------------
void
CoarseAverageFace::averageGridData(FluxBox& a_coarsenedFine,
                                   const FluxBox& a_fine,
                                   int a_averageType) const
{

  for (int dir=0; dir<SpaceDim; dir++)
    {
      FArrayBox& coarseFab = a_coarsenedFine[dir];
      const FArrayBox& fineFab = a_fine[dir];

      const Box& coarseBox = coarseFab.box();

      // set up refinement box
      int boxHi = m_nRef-1;
      IntVect hiVect(D_DECL6(boxHi,boxHi,boxHi,
                             boxHi,boxHi,boxHi));
      // don't want to index at all in dir direction --
      // instead, want to just march along face.
      hiVect.setVal(dir,0);
      IntVect loVect(D_DECL6(0,0,0,0,0,0));
      Box refBox(loVect, hiVect);

      if (a_averageType == arithmetic)
        {
          FORT_AVERAGEFACE( CHF_FRA(coarseFab),
                            CHF_CONST_FRA(fineFab),
                            CHF_BOX(coarseBox),
                            CHF_CONST_INT(dir),
                            CHF_CONST_INT(m_nRef),
                            CHF_BOX(refBox));
        }
      else if (a_averageType == harmonic)
        {
          FORT_AVERAGEFACEHARMONIC( CHF_FRA(coarseFab),
                                    CHF_CONST_FRA(fineFab),
                                    CHF_BOX(coarseBox),
                                    CHF_CONST_INT(dir),
                                    CHF_CONST_INT(m_nRef),
                                    CHF_BOX(refBox));
        }
      else
        {
          MayDay::Error("CoarseAverageFace::averageGridData -- bad averageType");
        }
    }

}

#include "NamespaceFooter.H"

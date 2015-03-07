#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif
// Dan Martin, Fri, Jan 14, 2000

#ifndef _CoarseAverageFace_H_
#define _CoarseAverageFace_H_

#include "REAL.H"
#include "FluxBox.H"
#include "LevelData.H"
#include "NamespaceHeader.H"


/// replaces edge-centered coarse-level data w/ averaged fine-level data
/** This class replaces edge-centered data on a coarse level of refinement
    with the average of the finer-level data which overlays the edge.
    This class is similar to CoarseAverage */
class CoarseAverageFace
{
public:
  /// Default constructor
  CoarseAverageFace();

  /// defining constructor
  CoarseAverageFace(const DisjointBoxLayout& a_fineGrids, 
		    int a_nComp, int a_nRef);

  /// destructor
  ~CoarseAverageFace();

  /// defines the object
  void define(const DisjointBoxLayout& a_fineGrids, 
	      int a_nComp, int a_nRef);
  
  ///
  bool isDefined() const;

  /// averages fine-level data to coarse level
  void averageToCoarse(LevelData<FluxBox>& a_coarse_data, 
		       const LevelData<FluxBox>& a_fine_data);

  /// averages fine-level data to coarse level using harmonic averaging
  void averageToCoarseHarmonic(LevelData<FluxBox>& a_coarse_data, 
                               const LevelData<FluxBox>& a_fine_data);

protected:

  ///
  enum averageType {arithmetic = 0,
                    harmonic, 
                    NUM_AVERAGE_TYPES};

  bool m_isDefined;

  int m_nRef;

  // work array for coarsening of fine data, same "shape" as fine data
  LevelData<FluxBox> m_coarsenedFineData;

  /** utility function called by both averageToCoarse and 
      averageCoarseHarmonic (to avoid code duplication)
  */
  void computeAverages(LevelData<FluxBox>& a_coarse_data, 
		       const LevelData<FluxBox>& a_fine_data,
                       int a_averageType);


  /// averages entire single grid data from fine->crse
  void averageGridData(FluxBox& a_coarsenedFine, 
                       const FluxBox& fine,
                       int a_averageType) const;

  
};

#include "NamespaceFooter.H"
#endif
  




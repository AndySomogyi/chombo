#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

#ifndef _EBPatch_Polytropic_Factory_H_
#define _EBPatch_Polytropic_Factory_H_

#include "EBPatchPolytropic.H"
#include "EBPatchGodunovFactory.H"

///
/**
 */
class EBPatchPolytropicFactory: public EBPatchGodunovFactory
{
public:
  ///
  /**
   */
  virtual ~EBPatchPolytropicFactory();

  ///
  /**
   */
  EBPatchGodunov* create() const;


  ///
  /**
   */
  EBPatchPolytropicFactory(const EBPhysIBCFactory* const a_bc,
                           const Real&                   a_gamma,
                           const bool&                   a_useFourthOrderSlopes,
                           const bool&                   a_useFlattening,
                           const bool&                   a_useArtificialVisc,
                           const bool&                   a_useLimiting,
                           const bool&                   a_doRZCoords);

protected:
  const EBPhysIBCFactory* m_bcFactoryPtr;
  Real              m_gamma;
  bool              m_useFourthOrderSlopes;
  bool              m_useFlattening;
  bool              m_useArtificialVisc;
  bool              m_doRZCoords;
  bool              m_useLimiting;

private:

  //disallowed for all the usual reasons
  EBPatchPolytropicFactory()
  {
    MayDay::Error("invalid operator");
  }
  void operator=(const EBPatchPolytropicFactory& a_input)
  {
    MayDay::Error("invalid operator");
  }
  EBPatchPolytropicFactory(const EBPatchPolytropicFactory& a_input)
  {
    MayDay::Error("invalid operator");
  }
};
#endif

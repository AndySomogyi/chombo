#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

#ifndef _EBAMRCNSFACTORY_H_
#define _EBAMRCNSFACTORY_H_

#include "EBAMRCNS.H"
#include "AMRLevelFactory.H"
#include "EBPatchGodunovFactory.H"
#include "EBPatchPolytropicFactory.H"
#include "NamespaceHeader.H"

///
/**
 */
class EBAMRCNSFactory : public AMRLevelFactory
{
public:
  ///
  AMRLevel* new_amrlevel() const;

  ///
  EBAMRCNSFactory(const EBAMRCNSParams&                       a_params,
                  const RefCountedPtr<EBPatchPolytropicFactory> a_factory)
  {
    m_params  = a_params;
    m_factory = a_factory;
  }

  ///
  virtual ~EBAMRCNSFactory()
  {;}

protected:
  EBAMRCNSParams                       m_params;
  RefCountedPtr<EBPatchPolytropicFactory> m_factory;

private:
  //weak construction is bad
  EBAMRCNSFactory()
  {
    MayDay::Error("invalid operator");
  }
  //not disallowing operator= or copy construction because they make sense here.
};

#include "NamespaceFooter.H"
#endif

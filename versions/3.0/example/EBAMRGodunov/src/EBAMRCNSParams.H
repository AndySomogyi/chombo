#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

#ifndef _EBAMRCNSParams_H_
#define _EBAMRCNSParams_H_

#include "Vector.H"
#include "RealVect.H"
#include "REAL.H"
#include "BaseDomainBC.H"
#include "BaseEBBC.H"
#include "NamespaceHeader.H"

///
class EBAMRCNSParams
{
public:
  bool                                m_doDiffusion;
  bool                                m_slipBoundaries;
  bool                                m_useMassRedist;
  bool                                m_doSmushing;
  int                                 m_redistRad;
  int                                 m_tagBufferSize;
  int                                 m_verbosity;
  Real                                m_domainLength;
  Real                                m_initialDtMultiplier;
  Real                                m_cfl;
  Real                                m_refineThresh;
  Real                                m_specHeatCv;
  Real                                m_thermalCond;
  Real                                m_viscosityMu;
  Real                                m_viscosityLa;
  RefCountedPtr<BaseDomainBCFactory>  m_doBCTemp;
  RefCountedPtr<BaseDomainBCFactory>  m_doBCVelo;
  RefCountedPtr<BaseEBBCFactory>      m_ebBCVelo;
  RefCountedPtr<BaseEBBCFactory>      m_ebBCTemp;
  EBAMRCNSParams()
  {
    m_initialDtMultiplier = 0.1;
    m_cfl                 = 0.5;
    m_redistRad           = 1;
    m_domainLength        = 1.0;
    m_refineThresh        = -123456789;
    m_tagBufferSize       = 1;
    m_verbosity           = 4;
    m_useMassRedist       = true;
    m_doSmushing          = true;
    m_viscosityMu         = 1.0e-3;
    m_viscosityLa         = 1.0e-3;
    m_specHeatCv          = 1.0e-3;
    m_thermalCond         = 1.0e-3;
  }
  ~EBAMRCNSParams()
  {;}
};


#include "NamespaceFooter.H"
#endif

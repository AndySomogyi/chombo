#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif


#include "EBConstantCFInterp.H"
#include "EBConstantCFInterpF_F.H"
#include "NamespaceHeader.H"


/////////////
EBConstantCFInterp::~EBConstantCFInterp()
{
}
//////////////
EBConstantCFInterp::
EBConstantCFInterp(const DisjointBoxLayout&  a_dbl,
                   const EBISLayout&         a_ebisl,
                   const ProblemDomain&      a_domain,
                   const IntVect&            a_nGhost)
{

  m_domain         = a_domain;
  m_dbl            = a_dbl;
  m_ebisl          = a_ebisl;
  m_nGhost         = a_nGhost;
  m_cornerCopier.define(a_dbl, a_dbl, m_domain, a_nGhost, true);
}


////
void
EBConstantCFInterp::
interpolate(LevelData<EBCellFAB>&   a_phi)
{
  CH_assert(m_nGhost == a_phi.ghostVect());
  CH_assert(a_phi.nComp() == 1);
  CH_assert( a_phi.ghostVect() >= IntVect::Unit);

  a_phi.exchange(Interval(0,0));

  for (DataIterator dit = m_dbl.dataIterator(); dit.ok(); ++dit)
    {
      for (int idir = 0; idir < SpaceDim; idir++)
        {
          for (SideIterator sit; sit.ok(); sit.next())
            {
              interpolate(a_phi[dit()],dit(),idir,sit());
            }
        }
    }
  Interval interv(0,0);
  a_phi.exchange(interv);
  //now exchange into corner cells covered by other ghost cells.
  a_phi.exchange(interv, m_cornerCopier);
}
/////////
void
EBConstantCFInterp::
interpolate(EBCellFAB&            a_phi,
            const DataIndex&      a_datInd,
            int                   a_idir,
            Side::LoHiSide        a_hiorlo)
{
  CH_assert((a_idir >= 0) && (a_idir  < SpaceDim));
  CH_assert((a_hiorlo == Side::Lo )||(a_hiorlo == Side::Hi ));
  CH_assert(a_phi.nComp() == 1);

  Box phiBox = a_phi.getRegion();

  int isign = sign(a_hiorlo);
  int nghost = m_nGhost[a_idir];
  CH_assert(nghost >= 1);
  int ishift = -(nghost+1)*isign;
  //want first row of valid cells of phi box
  Box validBox = adjCellBox(phiBox, a_idir, a_hiorlo, 1);
  validBox.shift(a_idir, ishift);
  BaseFab<Real>& regPhi = a_phi.getSingleValuedFAB();
  //this will break when the  EB crosses the coarse-fine interface.
  FORT_REGCONSTANTINTERP(CHF_FRA1(regPhi,0),
                         CHF_BOX(validBox),
                         CHF_CONST_INT(a_idir),
                         CHF_CONST_INT(nghost),
                         CHF_CONST_INT(isign));

}
////
#include "NamespaceFooter.H"

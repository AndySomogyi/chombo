#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

#include "EBViscousTensorOp.H"
#include "EBAMRPoissonOpF_F.H"
#include "EBAMRPoissonOp.H"
#include "InterpF_F.H"
#include "EBArith.H"
#include "ViscousTensorOp.H"
#include "ViscousTensorOpF_F.H"
#include "EBViscousTensorOpF_F.H"
#include "AMRPoissonOpF_F.H"
#include "NamespaceHeader.H"

bool EBViscousTensorOp::s_turnOffBCs = false; //needs to default to false

/*****/
EBViscousTensorOp::
EBViscousTensorOp(const EBLevelGrid &                                a_eblgFine,
                  const EBLevelGrid &                                a_eblg,
                  const EBLevelGrid &                                a_eblgCoar,
                  const EBLevelGrid &                                a_eblgCoarMG,
                  const Real&                                        a_alpha,
                  const Real&                                        a_beta,
                  const RefCountedPtr<LevelData<EBCellFAB> >&        a_acoef,
                  const RefCountedPtr<LevelData<EBFluxFAB> >&        a_eta,
                  const RefCountedPtr<LevelData<EBFluxFAB> >&        a_lambda,
                  const RefCountedPtr<LevelData<BaseIVFAB<Real> > >& a_etaIrreg,
                  const RefCountedPtr<LevelData<BaseIVFAB<Real> > >& a_lambdaIrreg,
                  const Real&                                        a_dx,
                  const Real&                                        a_dxCoar,
                  const int&                                         a_refToFine,
                  const int&                                         a_refToCoar,
                  const RefCountedPtr<ViscousBaseDomainBC>&          a_domainBC,
                  const RefCountedPtr<ViscousBaseEBBC>&              a_ebBC,
                  const bool&                                        a_hasMGObjects,
                  const IntVect&                                     a_ghostCellsPhi,
                  const IntVect&                                     a_ghostCellsRHS)
{
  EBArith::getMultiColors(m_colors);
  m_eblgFine      =  a_eblgFine;
  m_eblg          =  a_eblg;
  m_eblgCoar      =  a_eblgCoar;
  m_eblgCoarMG    =  a_eblgCoarMG;
  m_alpha         =  a_alpha;
  m_beta          =  a_beta;
  m_acoef         =  a_acoef;
  m_eta           =  a_eta;
  m_etaIrreg      =  a_etaIrreg;
  m_lambda        =  a_lambda;
  m_lambdaIrreg   =  a_lambdaIrreg;
  m_dx            =  a_dx;
  m_dxCoar        =  a_dxCoar;
  m_refToFine     =  a_refToFine;
  m_refToCoar     =  a_refToCoar;
  m_domainBC      =  a_domainBC;
  m_ebBC          =  a_ebBC;
  m_hasMGObjects  =  a_hasMGObjects;
  m_ghostCellsPhi =  a_ghostCellsPhi;
  m_ghostCellsRHS =  a_ghostCellsRHS;

  m_hasCoar = a_eblgCoar.isDefined();
  m_hasFine = a_eblgFine.isDefined();
  if(m_hasCoar)
    {
      ProblemDomain domainCoar = coarsen(m_eblg.getDomain(), m_refToCoar);
      m_interpWithCoarser = RefCountedPtr<EBTensorCFInterp>(new EBTensorCFInterp(m_eblg.getDBL(),  m_eblgCoar.getDBL(),
                                                                                 m_eblg.getEBISL(),m_eblgCoar.getEBISL(),
                                                                                 domainCoar, m_refToCoar, SpaceDim, m_dx));

      for (int idir = 0; idir < SpaceDim; idir++)
        {
          m_loCFIVS[idir].define(m_eblg.getDBL());
          m_hiCFIVS[idir].define(m_eblg.getDBL());

          for (DataIterator dit = m_eblg.getDBL().dataIterator();  dit.ok(); ++dit)
            {
              m_loCFIVS[idir][dit()].define(m_eblg.getDomain(), m_eblg.getDBL().get(dit()),
                                            m_eblg.getDBL(), idir,Side::Lo);
              m_hiCFIVS[idir][dit()].define(m_eblg.getDomain(), m_eblg.getDBL().get(dit()),
                                            m_eblg.getDBL(), idir,Side::Hi);
            }
        }

      //if this fails, then the AMR grids violate proper nesting.
      ProblemDomain domainCoarsenedFine;
      DisjointBoxLayout dblCoarsenedFine;

      int maxBoxSize = 32;
      bool dumbool;
      bool hasCoarser = EBAMRPoissonOp::getCoarserLayouts(dblCoarsenedFine,
                                                          domainCoarsenedFine,
                                                          m_eblg.getDBL(),
                                                          m_eblg.getEBISL(),
                                                          m_eblg.getDomain(),
                                                          m_refToCoar,
                                                          m_eblg.getEBIS(),
                                                          maxBoxSize, dumbool);

      //should follow from coarsenable
      if(hasCoarser)
        {
          EBLevelGrid eblgCoarsenedFine(dblCoarsenedFine, domainCoarsenedFine, 4, Chombo_EBIS::instance());
          m_ebInterp.define( m_eblg.getDBL(),     m_eblgCoar.getDBL(),
                             m_eblg.getEBISL(), m_eblgCoar.getEBISL(),
                             domainCoarsenedFine, m_refToCoar, SpaceDim,
                             m_eblg.getEBIS(),     a_ghostCellsPhi);
          m_ebAverage.define(m_eblg.getDBL(),   eblgCoarsenedFine.getDBL(),
                             m_eblg.getEBISL(), eblgCoarsenedFine.getEBISL(),
                             domainCoarsenedFine, m_refToCoar, SpaceDim,
                             m_eblg.getEBIS(), a_ghostCellsRHS);

        }
    }
  if(m_hasMGObjects)
    {
      int mgRef = 2;
      m_eblgCoarMG = a_eblgCoarMG;

      m_ebInterpMG.define( m_eblg.getDBL(),     m_eblgCoarMG.getDBL(),
                           m_eblg.getEBISL(), m_eblgCoarMG.getEBISL(),
                          m_eblgCoarMG.getDomain(), mgRef, SpaceDim,
                           m_eblg.getEBIS(),   a_ghostCellsPhi);
      m_ebAverageMG.define(m_eblg.getDBL(),     m_eblgCoarMG.getDBL(),
                           m_eblg.getEBISL(), m_eblgCoarMG.getEBISL(),
                           m_eblgCoarMG.getDomain() , mgRef, SpaceDim,
                           m_eblg.getEBIS(),   a_ghostCellsRHS);

    }

  defineStencils();
}
/*****/
void
EBViscousTensorOp::
defineStencils()
{
  int ncomp = SpaceDim;
  Real safety = 0.49;
  m_vofIter.define(m_eblg.getDBL());
  m_vofMult.define(m_eblg.getDBL());
  LayoutData<BaseIVFAB<VoFStencil> > slowStencil(m_eblg.getDBL());

  EBCellFactory ebcellfact(m_eblg.getEBISL());
  m_relCoef.define(m_eblg.getDBL(), SpaceDim, IntVect::Zero, ebcellfact);
  m_grad.define(m_eblg.getDBL(), SpaceDim*SpaceDim, IntVect::Unit, ebcellfact);
  //define regular relaxation coefficent
  for(DataIterator dit = m_eblg.getDBL().dataIterator(); dit.ok(); ++dit)
    {
      const Box& grid = m_eblg.getDBL().get(dit());
      const EBCellFAB& acofab = (*m_acoef)[dit()];
      //initialize lambda = alpha*acoef
      m_relCoef[dit()].setVal(0.);
      for(int idir = 0; idir < SpaceDim; idir++)
        {
          m_relCoef[dit()].plus(acofab, 0, idir, 1);
        }
      m_relCoef[dit()]*= m_alpha;

      BaseFab<Real>& regRel =   m_relCoef[dit()].getSingleValuedFAB();

      for(int idir = 0; idir < SpaceDim; idir++)
        {
          BaseFab<Real>& regEta = (*m_eta)   [dit()][idir].getSingleValuedFAB();
          BaseFab<Real>& regLam = (*m_lambda)[dit()][idir].getSingleValuedFAB();
          FORT_DECRINVRELCOEFVTOP(CHF_FRA(regRel),
                                  CHF_FRA(regEta),
                                  CHF_FRA(regLam),
                                  CHF_CONST_REAL(m_beta),
                                  CHF_BOX(grid),
                                  CHF_REAL(m_dx),
                                  CHF_INT(idir),
                                  CHF_INT(ncomp));
        }

      //now invert so lambda = stable lambda for variable coef lapl
      //(according to phil, this is the correct lambda)
      FORT_INVERTLAMBDAVTOP(CHF_FRA(regRel),
                            CHF_REAL(safety),
                            CHF_BOX(grid),
                            CHF_INT(ncomp));
    }
  EBLevelDataOps::setToZero(m_grad);

  Box sideBoxLo[SpaceDim];
  Box sideBoxHi[SpaceDim];
  for(int ivar = 0; ivar < SpaceDim; ivar++)
    {
      int idir = ivar;
      Box domainBox = m_eblg.getDomain().domainBox();
      sideBoxLo[idir] = adjCellLo(domainBox, idir, 1);
      sideBoxLo[idir].shift(idir,  1);
      sideBoxHi[idir] = adjCellHi(domainBox, idir, 1);
      sideBoxHi[idir].shift(idir, -1);
      m_stencil[ivar].define(m_eblg.getDBL());
      m_vofIterDomLo[ivar].define( m_eblg.getDBL()); // vofiterator cache for domain lo
      m_vofIterDomHi[ivar].define( m_eblg.getDBL()); // vofiterator cache for domain hi
    }
  //first define the iterators
  for(DataIterator dit = m_eblg.getDBL().dataIterator(); dit.ok(); ++dit)
    {
      const EBISBox& ebisBox = m_eblg.getEBISL()[dit()];
      const Box&     grid = m_eblg.getDBL().get(dit());
      //need to grow the irregular set by one near multivalued cells
      IntVectSet ivsIrreg = ebisBox.getIrregIVS(grid);
      IntVectSet ivsMulti = ebisBox.getMultiCells(grid);
      ivsMulti.grow(1);
      ivsMulti &= grid;

      IntVectSet ivsComplement(grid);
      ivsComplement -= ivsMulti;
      ivsComplement -= ivsIrreg;
      //ivscomplement now contains the complement of the cells we need;

      IntVectSet ivsSten(grid);
      ivsSten -= ivsComplement;

      m_vofIter[dit()].define(ivsSten, ebisBox.getEBGraph()   );
      m_vofMult[dit()].define(ivsMulti, ebisBox.getEBGraph()   );
      for(int idir = 0; idir < SpaceDim; idir++)
        {
          IntVectSet loIrreg = ivsSten;
          IntVectSet hiIrreg = ivsSten;
          loIrreg &= sideBoxLo[idir];
          hiIrreg &= sideBoxHi[idir];
          m_vofIterDomLo[idir][dit()].define(loIrreg,ebisBox.getEBGraph());
          m_vofIterDomHi[idir][dit()].define(hiIrreg,ebisBox.getEBGraph());
        }
      slowStencil[dit()].define(ivsSten, ebisBox.getEBGraph(), 1);

    }

  m_domainBC->setCoef(m_eblg,   m_beta,      m_eta,      m_lambda);
  m_ebBC->setCoef(    m_eblg,   m_beta,      m_etaIrreg, m_lambdaIrreg);
  m_ebBC->define(    *m_eblg.getCFIVS(), 1.0);

  //now define the stencils for each variable
  for(int ivar = 0; ivar < SpaceDim; ivar++)
    {
      LayoutData<BaseIVFAB<VoFStencil> >* fluxStencil = m_ebBC->getFluxStencil(ivar);
      for(DataIterator dit = m_eblg.getDBL().dataIterator(); dit.ok(); ++dit)
        {
          const EBISBox& ebisBox = m_eblg.getEBISL()[dit()];
          for(m_vofIter[dit()].reset(); m_vofIter[dit()].ok(); ++m_vofIter[dit()])
            {
              const VolIndex& vof = m_vofIter[dit()]();
              VoFStencil& vofsten = slowStencil[dit()](vof, 0);
              getVoFStencil(vofsten, vof, dit(), ivar);
              if(fluxStencil != NULL)
                {
                  BaseIVFAB<VoFStencil>& stenFAB = (*fluxStencil)[dit()];
                  //only real irregular cells are going to have EB fluxes
                  //but we are working on a bigger set
                  if(stenFAB.getIVS().contains(vof.gridIndex()))
                    {
                      VoFStencil fluxStencilPt = stenFAB(vof, 0);
                      Real   betaPt =   m_beta;
                      Real bndryArea = ebisBox.bndryArea(vof);
                      //this might need to get multiplied by the -normal[idir]
                      Real factor = betaPt*bndryArea/m_dx;
                      fluxStencilPt *= factor;
                      slowStencil[dit()](vof, 0) += fluxStencilPt;
                    }

                }

              Real diagWeight = EBArith::getDiagWeight(slowStencil[dit()](vof, 0), vof);
              m_relCoef[dit()](vof, ivar) = safety/(diagWeight);
            }

          m_stencil[ivar][dit()] = RefCountedPtr<EBStencil>(new EBStencil(m_vofIter[dit()].getVector(),  slowStencil[dit()],
                                                                          m_eblg.getDBL().get(dit()), m_eblg.getEBISL()[dit()],
                                                                          m_ghostCellsPhi, m_ghostCellsRHS, ivar));
        }
    }
}
/*****/
/* generate vof stencil as a divergence of flux stencils */
/***/
void
EBViscousTensorOp::
getVoFStencil(VoFStencil&      a_vofStencil,
              const VolIndex&  a_vof,
              const DataIndex& a_dit,
              int             a_ivar)
{
  const EBISBox& ebisBox = m_eblg.getEBISL()[a_dit];
  a_vofStencil.clear();
  for(int idir = 0; idir < SpaceDim; idir++)
    {
      for(SideIterator sit; sit.ok(); ++sit)
        {
          int isign = sign(sit());
          Vector<FaceIndex> faces = ebisBox.getFaces(a_vof, idir, sit());
          for(int iface = 0; iface < faces.size(); iface++)
            {
              VoFStencil fluxStencil;
              getFluxStencil(fluxStencil, faces[iface], a_dit, a_ivar);
              Real areaFrac = ebisBox.areaFrac(faces[iface]);
              fluxStencil *= Real(isign)*areaFrac/m_dx;
              a_vofStencil += fluxStencil;
            }
        }
    }

  //recall that eta and lambda are part of the flux.
  //mulitply divergence of flux by beta
  Real betaPt = m_beta;
  a_vofStencil *= betaPt;
  //debug
  //a_vofStencil.clear();
  //end debug
  //so now add in kappa*alpha*I
  Real kappa = ebisBox.volFrac(a_vof);
  Real acoefPt = (*m_acoef)[a_dit](a_vof,0);
  a_vofStencil.add(a_vof, kappa*m_alpha*acoefPt, a_ivar);
}
/***/
/// stencil for flux computation.   the truly ugly part of this computation
/// beta and eta are multiplied in here
/****/
void
EBViscousTensorOp::
getFluxStencil(VoFStencil&      a_fluxStencil,
               const FaceIndex& a_face,
               const DataIndex& a_dit,
               int a_ivar)
{
  //need to do this by interpolating to centroids
  //so get the stencil at each face center and add with
  //interpolation weights
  FaceStencil interpSten = EBArith::getInterpStencil(a_face,
                                                     (*m_eblg.getCFIVS())[a_dit],
                                                     m_eblg.getEBISL()[a_dit],
                                                     m_eblg.getDomain());

  a_fluxStencil.clear();
  for(int isten = 0; isten < interpSten.size(); isten++)
    {
      const FaceIndex& face = interpSten.face(isten);
      const Real&    weight = interpSten.weight(isten);
      VoFStencil faceCentSten;
      getFaceCenteredFluxStencil(faceCentSten, face, a_dit, a_ivar);
      faceCentSten *= weight;
      a_fluxStencil += faceCentSten;
    }
  //debug
  //a_fluxStencil.clear();
  //if(a_face.direction() == 1)  getFaceCenteredFluxStencil(a_fluxStencil, a_face, a_dit, a_ivar);
  //end debug

}
void
EBViscousTensorOp::
getFaceCenteredFluxStencil(VoFStencil&      a_fluxStencil,
                           const FaceIndex& a_face,
                           const DataIndex& a_dit,
                           int a_ivar)
{
  int faceDir= a_face.direction();
  VoFStencil  gBTranSten, gBNormSten;
  //the flux is lambda I diverge(B) + eta(gradB + gradB tran)
  //so the ivar flux for the faceDir direction is
  // lambda*delta(ivar, faceDir)*divergence(B) + eta*(partial(B_ivar, faceDir) + partial(B_faceDir, ivar))

  getGradientStencil(gBNormSten, a_ivar , faceDir, a_face, a_dit);
  getGradientStencil(gBTranSten, faceDir, a_ivar , a_face, a_dit);
  Real etaFace = (*m_eta)[a_dit][faceDir](a_face, 0);
  gBNormSten *=    etaFace;
  gBTranSten *=    etaFace;

  a_fluxStencil += gBNormSten;
  a_fluxStencil += gBTranSten;
  if(a_ivar == faceDir)
    {
      Real lambdaFace = (*m_lambda)[a_dit][faceDir](a_face, 0);
      VoFStencil divergeSten;
      getDivergenceStencil(divergeSten, a_face, a_dit);
      divergeSten   *= lambdaFace;
      a_fluxStencil += divergeSten;
    }
}
/***/
void
EBViscousTensorOp::
getGradientStencil(VoFStencil&  a_gradStencil,
                   int a_ivar,
                   int a_diffDir,
                   const FaceIndex& a_face,
                   const DataIndex& a_dit)
{
  getGradientStencil(a_gradStencil, a_ivar, a_diffDir, a_face, a_dit, m_dx, m_eblg);
}

void
EBViscousTensorOp::
getGradientStencil(VoFStencil&  a_gradStencil,
                   int a_ivar,
                   int a_diffDir,
                   const FaceIndex& a_face,
                   const DataIndex& a_dit,
                   const Real     & a_dx,
                   const EBLevelGrid& a_eblg)
{
  a_gradStencil.clear();
  if((a_face.direction() == a_diffDir) && (!a_face.isBoundary()))
    {
      a_gradStencil.add(a_face.getVoF(Side::Hi),  1.0/a_dx, a_ivar);
      a_gradStencil.add(a_face.getVoF(Side::Lo), -1.0/a_dx, a_ivar);
    }
  else
    {
      int numGrad = 0;
      for(SideIterator sit; sit.ok(); ++sit)
        {
          const VolIndex& vofSide = a_face.getVoF(sit());
          if(a_eblg.getDomain().contains(vofSide.gridIndex()))
            {
              VoFStencil stenSide;
              IntVectSet& cfivs = (*a_eblg.getCFIVS())[a_dit];
              EBArith::getFirstDerivStencil(stenSide, vofSide,
                                            a_eblg.getEBISL()[a_dit],
                                            a_diffDir, a_dx, &cfivs, a_ivar);
              a_gradStencil += stenSide;
              numGrad++;
            }
        }
      if(numGrad > 1)
        {
          a_gradStencil *= 1.0/Real(numGrad);
        }
    }
}

/***/
void
EBViscousTensorOp::
getDivergenceStencil(VoFStencil&      a_divStencil,
                     const FaceIndex& a_face,
                     const DataIndex& a_dit)
{
  getDivergenceStencil(a_divStencil, a_face, a_dit, m_dx, m_eblg);
}
void
EBViscousTensorOp::
getDivergenceStencil(VoFStencil&      a_divStencil,
                     const FaceIndex& a_face,
                     const DataIndex& a_dit,
                     const Real     & a_dx,
                     const EBLevelGrid& a_eblg)
{
  a_divStencil.clear();
  for(int diffDir = 0; diffDir < SpaceDim; diffDir++)
    {
      VoFStencil diffStencilDir;
      //difference direction and variable are the same
      //because that is the definition of the divergence.
      getGradientStencil(diffStencilDir, diffDir, diffDir, a_face, a_dit, a_dx, a_eblg);
      a_divStencil += diffStencilDir;
    }
}
/***/
EBViscousTensorOp::
~EBViscousTensorOp()
{
}
/***/
void
EBViscousTensorOp::
AMRResidualNC(LevelData<EBCellFAB>&       a_residual,
              const LevelData<EBCellFAB>& a_phiFine,
              const LevelData<EBCellFAB>& a_phi,
              const LevelData<EBCellFAB>& a_rhs,
              bool a_homogeneousBC,
              AMRLevelOp<LevelData<EBCellFAB> >* a_finerOp)
{
  //dummy. there is no coarse when this is called
  CH_assert(a_residual.ghostVect() == m_ghostCellsRHS);
  CH_assert(a_rhs.ghostVect() == m_ghostCellsRHS);
  LevelData<EBCellFAB> phiC;
  AMRResidual(a_residual, a_phiFine, a_phi, phiC, a_rhs, a_homogeneousBC, a_finerOp);
}
/***/
void
EBViscousTensorOp::
AMROperatorNC(LevelData<EBCellFAB>&       a_LofPhi,
              const LevelData<EBCellFAB>& a_phiFine,
              const LevelData<EBCellFAB>& a_phi,
              bool a_homogeneousBC,
              AMRLevelOp<LevelData<EBCellFAB> >* a_finerOp)
{
  //dummy. there is no coarse when this is called
  CH_assert(a_LofPhi.ghostVect() == m_ghostCellsRHS);
  LevelData<EBCellFAB> phiC;
  AMROperator(a_LofPhi, a_phiFine, a_phi, phiC,
              a_homogeneousBC, a_finerOp);
}
/*****/
void
EBViscousTensorOp::
residual(LevelData<EBCellFAB>&       a_residual,
         const LevelData<EBCellFAB>& a_phi,
         const LevelData<EBCellFAB>& a_rhs,
         bool                        a_homogeneousBC)
{
  //this is a multigrid operator so only homogeneous CF BC
  //and null coar level
  CH_assert(a_residual.ghostVect() == m_ghostCellsRHS);
  CH_assert(a_phi.ghostVect() == m_ghostCellsPhi);
  applyOp(a_residual,a_phi, a_homogeneousBC);
  incr(a_residual, a_rhs, -1.0);
  scale(a_residual, -1.0);
}

/*****/
void
EBViscousTensorOp::
preCond(LevelData<EBCellFAB>&       a_phi,
        const LevelData<EBCellFAB>& a_rhs)
{
  relax(a_phi, a_rhs, 40);
}



/*****/
void
EBViscousTensorOp::
create(LevelData<EBCellFAB>&       a_lhs,
       const LevelData<EBCellFAB>& a_rhs)
{
  int ncomp = a_rhs.nComp();
  EBCellFactory ebcellfact(m_eblg.getEBISL());
  a_lhs.define(m_eblg.getDBL(), ncomp, a_rhs.ghostVect(), ebcellfact);
}

/*****/
void
EBViscousTensorOp::
createCoarsened(LevelData<EBCellFAB>&       a_lhs,
                const LevelData<EBCellFAB>& a_rhs,
                const int&                  a_refRat)
{
  int ncomp = a_rhs.nComp();
  IntVect ghostVect = a_rhs.ghostVect();

  CH_assert(m_eblg.getDBL().coarsenable(a_refRat));

  //fill ebislayout
  DisjointBoxLayout dblCoarsenedFine;
  coarsen(dblCoarsenedFine, m_eblg.getDBL(), a_refRat);

  EBISLayout ebislCoarsenedFine;
  IntVect ghostVec = a_rhs.ghostVect();
  //const EBIndexSpace* const ebisPtr = Chombo_EBIS::instance();
  ProblemDomain coarDom = coarsen(m_eblg.getDomain(), a_refRat);
  m_eblg.getEBIS()->fillEBISLayout(ebislCoarsenedFine, dblCoarsenedFine, coarDom , ghostVec[0]);
  if(m_refToCoar > 1)
    {
      ebislCoarsenedFine.setMaxRefinementRatio(m_refToCoar, m_eblg.getEBIS());
    }

  //create coarsened data
  EBCellFactory ebcellfactCoarsenedFine(ebislCoarsenedFine);
  a_lhs.define(dblCoarsenedFine, ncomp,ghostVec, ebcellfactCoarsenedFine);
}

/*****/
Real
EBViscousTensorOp::
AMRNorm(const LevelData<EBCellFAB>& a_coarResid,
        const LevelData<EBCellFAB>& a_fineResid,
        const int& a_refRat,
        const int& a_ord)
{
  int ncomp = a_coarResid.nComp();
  Interval interv(0, ncomp-1);
  const DisjointBoxLayout& coarGrids = a_coarResid.disjointBoxLayout();
  const DisjointBoxLayout& fineGrids = a_fineResid.disjointBoxLayout();
  CH_assert(coarGrids == m_eblg.getDBL());

  //create temp and zero out under finer grids
  EBCellFactory ebcellfact(m_eblg.getEBISL());
  LevelData<EBCellFAB> coarTemp(coarGrids, ncomp, IntVect::Zero, ebcellfact);
  a_coarResid.copyTo(interv, coarTemp, interv);

  for(DataIterator dit = coarGrids.dataIterator(); dit.ok(); ++dit)
    {
      EBCellFAB& coarTempFAB = coarTemp[dit()];
      const EBISBox& ebisBox = m_eblg.getEBISL()[dit()];

      LayoutIterator litFine = fineGrids.layoutIterator();
      for (litFine.reset(); litFine.ok(); ++litFine)
        {
          Box overlayBox = coarTempFAB.box();
          Box coarsenedGrid = coarsen(fineGrids[litFine()], a_refRat);

          overlayBox &= coarsenedGrid;
          if (!overlayBox.isEmpty())
            {
              BaseFab<Real>& regToZeroFAB = coarTempFAB.getSingleValuedFAB();
              FORT_AMRPZEROSUB(CHF_FRA(regToZeroFAB),
                               CHF_BOX(overlayBox),
                               CHF_INT(ncomp));

              IntVectSet ivsZero = ebisBox.getMultiCells(overlayBox);

              for(VoFIterator vofit(ivsZero, ebisBox.getEBGraph()); vofit.ok(); ++vofit)
                {
                  for(int ivar =0; ivar < ncomp; ivar++)
                    {
                      coarTempFAB(vofit(), ivar) = 0.0;
                    }
                }
            }
        }
    }
  //return norm of temp
  return norm(coarTemp, a_ord);
}

/*****/
void
EBViscousTensorOp::
assign(LevelData<EBCellFAB>&       a_lhs,
       const LevelData<EBCellFAB>& a_rhs)
{
  EBLevelDataOps::assign(a_lhs,a_rhs);
}

/*****/
Real
EBViscousTensorOp::
dotProduct(const LevelData<EBCellFAB>& a_1,
           const LevelData<EBCellFAB>& a_2)
{
  ProblemDomain domain;
  Real volume;

  return EBLevelDataOps::kappaDotProduct(volume,a_1,a_2,EBLEVELDATAOPS_ALLVOFS,domain);
}

/*****/
void
EBViscousTensorOp::
incr(LevelData<EBCellFAB>&       a_lhs,
     const LevelData<EBCellFAB>& a_x,
     Real                        a_scale)
{
  EBLevelDataOps::incr(a_lhs,a_x,a_scale);
}

/*****/
void
EBViscousTensorOp::
axby(LevelData<EBCellFAB>&       a_lhs,
     const LevelData<EBCellFAB>& a_x,
     const LevelData<EBCellFAB>& a_y,
     Real                        a_a,
     Real                        a_b)
{
  EBLevelDataOps::axby(a_lhs,a_x,a_y,a_a,a_b);
}

/*****/
void
EBViscousTensorOp::
scale(LevelData<EBCellFAB>& a_lhs,
      const Real&           a_scale)
{
  EBLevelDataOps::scale(a_lhs,a_scale);
}

/*****/
Real
EBViscousTensorOp::
norm(const LevelData<EBCellFAB>& a_rhs,
     int                         a_ord)
{
  CH_TIME("EBAMRPoissonOp::norm");
  ProblemDomain domain;
  Real volume, sum;
  //multilevel linearop just wants the vol weighted sum
  IntVectSet ivsExclude;
  int comp = 0;
  RealVect vectDx = m_dx*RealVect::Unit;
  EBArith::volWeightedSum(sum, volume, a_rhs, m_eblg.getDBL(), m_eblg.getEBISL(), vectDx, ivsExclude, comp, a_ord);
  return sum;
}

/*****/
void
EBViscousTensorOp::
setToZero(LevelData<EBCellFAB>& a_lhs)
{
  EBLevelDataOps::setToZero(a_lhs);
}

/*****/
void
EBViscousTensorOp::
setVal(LevelData<EBCellFAB>& a_lhs, const Real& a_value)
{
  EBLevelDataOps::setVal(a_lhs, a_value);
}

void
EBViscousTensorOp::
createCoarser(LevelData<EBCellFAB>&       a_coar,
              const LevelData<EBCellFAB>& a_fine,
              bool                        a_ghosted)
{
  const DisjointBoxLayout& dbl = m_eblgCoarMG.getDBL();
  ProblemDomain coarDom = coarsen(m_eblg.getDomain(), 2);

  int nghost = a_fine.ghostVect()[0];
  EBISLayout coarEBISL;

  const EBIndexSpace* const ebisPtr = Chombo_EBIS::instance();
  ebisPtr->fillEBISLayout(coarEBISL,
                          dbl, coarDom, nghost);

  EBCellFactory ebcellfact(coarEBISL);
  a_coar.define(dbl, SpaceDim,a_fine.ghostVect(),ebcellfact);
}

/*****/
void
EBViscousTensorOp::
relax(LevelData<EBCellFAB>&       a_phi,
      const LevelData<EBCellFAB>& a_rhs,
      int                         a_iterations)
{
  CH_assert(a_phi.isDefined());
  CH_assert(a_rhs.isDefined());
  CH_assert(a_phi.ghostVect() >= IntVect::Unit);
  CH_assert(a_phi.nComp() == a_rhs.nComp());
  LevelData<EBCellFAB> lphi;
  create(lphi, a_rhs);
  // do first red, then black passes
  for (int whichIter =0; whichIter < a_iterations; whichIter++)
    {
      for(int icolor = 0; icolor < m_colors.size(); icolor++)
        {
          homogeneousCFInterp(a_phi);

          //after this lphi = L(phi)
          //this call contains bcs and exchange
          applyOp(  lphi,  a_phi, true);
          gsrbColor(a_phi, lphi, a_rhs, m_colors[icolor]);
        }
    }
}

void
EBViscousTensorOp::
homogeneousCFInterp(LevelData<EBCellFAB>&   a_phi)
{
  if(m_hasCoar)
    {
      EBCellFactory factCoarse(m_eblgCoar.getEBISL());
      LevelData<EBCellFAB> phiCoarse(m_eblgCoar.getDBL(), a_phi.nComp(), a_phi.ghostVect(), factCoarse);
      EBLevelDataOps::setToZero(phiCoarse);
      cfinterp(a_phi, phiCoarse);
    }
}
void
EBViscousTensorOp::
homogeneousCFInterp(EBCellFAB&            a_phi,
                    const DataIndex&      a_datInd,
                    const int&            a_idir,
                    const Side::LoHiSide& a_hiorlo)
{
  if(m_hasCoar)
    {
      const CFIVS* cfivsPtr = NULL;

      if (a_hiorlo == Side::Lo)
        {
          cfivsPtr = &m_loCFIVS[a_idir][a_datInd];
        }
      else
        {
          cfivsPtr = &m_hiCFIVS[a_idir][a_datInd];
        }

      const IntVectSet& interpIVS = cfivsPtr->getFineIVS();
      if (cfivsPtr->isPacked() )
        {
          const int ihiorlo = sign(a_hiorlo);
          FORT_INTERPHOMO(CHF_FRA(a_phi.getSingleValuedFAB()),
                          CHF_BOX(cfivsPtr->packedBox()),
                          CHF_CONST_REAL(m_dx),
                          CHF_CONST_REAL(m_dxCoar),
                          CHF_CONST_INT(a_idir),
                          CHF_CONST_INT(ihiorlo));

        }
      else
        {
          if (!interpIVS.isEmpty())
            {
              const EBISBox&  ebisBox = m_eblg.getEBISL()[a_datInd];
              Real dxcoar = m_dxCoar;
              Real dxfine = m_dx;
              EBArith::interpolateCFH(a_phi, a_idir, a_hiorlo,
                                      ebisBox,  dxfine, dxcoar,
                                      interpIVS);

            }
        }
    }
}
/*****/
void
EBViscousTensorOp::
gsrbColor(LevelData<EBCellFAB>&       a_phi,
          const LevelData<EBCellFAB>& a_lph,
          const LevelData<EBCellFAB>& a_rhs,
          const IntVect&              a_color)
{

  const DisjointBoxLayout& dbl = a_phi.disjointBoxLayout();
  for(DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit)
    {
      Box dblBox  = dbl.get(dit());
      BaseFab<Real>&       regPhi =     a_phi[dit()].getSingleValuedFAB();
      const BaseFab<Real>& regLph =     a_lph[dit()].getSingleValuedFAB();
      const BaseFab<Real>& regRhs =     a_rhs[dit()].getSingleValuedFAB();
      const BaseFab<Real>& regRel = m_relCoef[dit()].getSingleValuedFAB();
      IntVect loIV = dblBox.smallEnd();
      IntVect hiIV = dblBox.bigEnd();

      for (int idir = 0; idir < SpaceDim; idir++)
        {
          if (loIV[idir] % 2 != a_color[idir])
            {
              loIV[idir]++;
            }
        }

      if (loIV <= hiIV)
        {
          int ncomp = SpaceDim;
          Box coloredBox(loIV, hiIV);
          FORT_GSRBVTOP(CHF_FRA(regPhi),
                        CHF_CONST_FRA(regLph),
                        CHF_CONST_FRA(regRhs),
                        CHF_CONST_FRA(regRel),
                        CHF_BOX(coloredBox),
                        CHF_CONST_INT(ncomp));
        }

      for(m_vofMult[dit()].reset(); m_vofMult[dit()].ok(); ++m_vofMult[dit()])
        {
          const VolIndex& vof = m_vofIter[dit()]();
          const IntVect& iv = vof.gridIndex();

          bool doThisVoF = true;
          for (int idir = 0; idir < SpaceDim; idir++)
            {
              if (iv[idir] % 2 != a_color[idir])
                {
                  doThisVoF = false;
                  break;
                }
            }

          if(doThisVoF)
            {
              for(int ivar = 0; ivar < SpaceDim; ivar++)
                {
                  Real resid = a_rhs[dit()](vof, ivar) - a_lph[dit()](vof, ivar);
                  a_phi[dit()](vof, ivar) += m_relCoef[dit()](vof, ivar)*resid;
                }
            }
        }
    }
}
/*****/
void
EBViscousTensorOp::
restrictResidual(LevelData<EBCellFAB>&       a_resCoar,
                 LevelData<EBCellFAB>&       a_phi,
                 const LevelData<EBCellFAB>& a_rhs)
{
  LevelData<EBCellFAB> res;
  bool homogeneous = true;

  create(res, a_rhs);

  // Get the residual on the fine grid
  residual(res,a_phi,a_rhs,homogeneous);

  // now use our nifty averaging operator
  Interval variables(0, SpaceDim-1);
  if(m_hasMGObjects)
    {
      m_ebAverageMG.average(a_resCoar, res, variables);
    }
  else
    {
      m_ebAverage.average(a_resCoar, res, variables);
    }
}

/*****/
void
EBViscousTensorOp::
prolongIncrement(LevelData<EBCellFAB>&       a_phi,
                 const LevelData<EBCellFAB>& a_cor)
{
  Interval vars(0, SpaceDim-1);
  if(m_hasMGObjects)
    {
      m_ebInterpMG.pwcInterp(a_phi, a_cor, vars);
    }
  else
    {
      m_ebInterp.pwcInterp(a_phi, a_cor, vars);
    }
}

/*****/
int
EBViscousTensorOp::
refToCoarser()
{
  return m_refToCoar;
}

/*****/
int
EBViscousTensorOp::
refToFiner()
{
  return m_refToFine;
}

/*****/
void
EBViscousTensorOp::
AMRResidual(LevelData<EBCellFAB>& a_residual,
            const LevelData<EBCellFAB>& a_phiFine,
            const LevelData<EBCellFAB>& a_phi,
            const LevelData<EBCellFAB>& a_phiCoarse,
            const LevelData<EBCellFAB>& a_rhs,
            bool a_homogeneousBC,
            AMRLevelOp<LevelData<EBCellFAB> >* a_finerOp)
{
  //fillgrad is called in applyop
  this->cfinterp(a_phi, a_phiCoarse);

  applyOp(a_residual, a_phi, a_homogeneousBC);

  if (a_finerOp != NULL)
    {
      reflux(a_phiFine, a_phi,  a_residual, a_finerOp);
    }
  incr(a_residual, a_rhs, -1.0);
  scale(a_residual, -1.0);
}

/*****/
void
EBViscousTensorOp::
AMRResidualNF(LevelData<EBCellFAB>& a_residual,
              const LevelData<EBCellFAB>& a_phi,
              const LevelData<EBCellFAB>& a_phiCoarse,
              const LevelData<EBCellFAB>& a_rhs,
              bool a_homogeneousBC)
{
  this->cfinterp(a_phi, a_phiCoarse);
  this->residual(a_residual, a_phi, a_rhs, a_homogeneousBC ); //apply boundary conditions
}

/*****/
void
EBViscousTensorOp::
reflux(const LevelData<EBCellFAB>&        a_phiFine,
       const LevelData<EBCellFAB>&        a_phi,
       LevelData<EBCellFAB>&              a_residual,
       AMRLevelOp<LevelData<EBCellFAB> >* a_finerOp)
{
  Interval interv(0,SpaceDim-1);
  CH_assert(a_phiFine.nComp() == SpaceDim);
  int ncomp = SpaceDim;

  if(!m_fluxReg.isDefined())
    {
      m_fluxReg.define(m_eblgFine.getDBL(),     m_eblg.getDBL(),
                       m_eblgFine.getEBISL(), m_eblg.getEBISL(),
                       m_eblg.getDomain(),
                       m_refToFine, ncomp, m_eblg.getEBIS());
    }
  m_fluxReg.setToZero();

  incrementFRCoar(a_phiFine, a_phi);

  incrementFRFine(a_phiFine, a_phi, a_finerOp);

  for(int idir = 0; idir < SpaceDim; idir++)
    {
      Real scale = m_beta/m_dx;
      m_fluxReg.reflux(a_residual, interv, scale, idir);
    }
}

/****/
void
EBViscousTensorOp::
incrementFRCoar(const LevelData<EBCellFAB>& a_phiFine,
                const LevelData<EBCellFAB>& a_phi)
{
  int ncomp = SpaceDim;
  Interval interv(0,SpaceDim-1);

  for(DataIterator dit = m_eblg.getDBL().dataIterator(); dit.ok(); ++dit)
    {
      const EBISBox& ebisBox = m_eblg.getEBISL()[dit()];
      for (int idir = 0; idir < SpaceDim; idir++)
        {
          //only want interior faces
          Box cellBox = m_eblg.getDBL().get(dit());
          cellBox.grow(idir, 1);
          cellBox &= m_eblg.getDomain();
          cellBox.grow(idir, -1);
#ifndef NDEBUG
          IntVectSet ivsIrreg = ebisBox.getIrregIVS(cellBox);
          if(!ivsIrreg.isEmpty())
            {
              pout() << "EBViscousTensorOp not wired for EBCF yet " << endl;
              MayDay::Error();
            }
#endif          
          Box faceBox = surroundingNodes(cellBox, idir);

          EBFaceFAB coarflux(ebisBox, cellBox, idir, ncomp);
          FArrayBox& regFlux       =       (FArrayBox&)     coarflux.getSingleValuedFAB();
          const FArrayBox& regPhi  = (const FArrayBox&) a_phi[dit()].getSingleValuedFAB();
          const FArrayBox& regGrad = (const FArrayBox&)m_grad[dit()].getSingleValuedFAB();
          getFlux(regFlux, regPhi,  regGrad, faceBox, idir, dit());
          //not wired for EBCF yet.   Not sure what to do with the humongous stencil.
          Real scale = m_beta;
          m_fluxReg.incrementCoarseRegular(coarflux, scale, dit(), interv, idir);
        }
    }
}
/***/
void
EBViscousTensorOp::
incrementFRFine(const LevelData<EBCellFAB>& a_phiFine,
                const LevelData<EBCellFAB>& a_phi,
                AMRLevelOp<LevelData<EBCellFAB> >* a_finerOp)
{
  int ncomp = SpaceDim;
  Interval interv(0,SpaceDim-1);
  EBViscousTensorOp& finerEBAMROp = (EBViscousTensorOp& )(*a_finerOp);

  //ghost cells of phiFine need to be filled
  LevelData<EBCellFAB>& phiFine = (LevelData<EBCellFAB>&) a_phiFine;
  finerEBAMROp.cfinterp(phiFine, a_phi);
  finerEBAMROp.fillGrad(phiFine);
  phiFine.exchange(interv);

  DataIterator ditf = a_phiFine.dataIterator();
  for (ditf.reset(); ditf.ok(); ++ditf)
    {
      const EBISBox& ebisBoxFine = m_eblgFine.getEBISL()[ditf()];
      const EBCellFAB& phiFine = a_phiFine[ditf()];
      const EBCellFAB& gradFine = finerEBAMROp.m_grad[ditf()];

      for (int idir = 0; idir < SpaceDim; idir++)
        {
          SideIterator sit;
          for (sit.begin(); sit.ok(); sit.next())
            {
              //only want interior faces
              Box cellBox = m_eblgFine.getDBL().get(ditf());
              cellBox.grow(idir, 1);
              cellBox &= m_eblgFine.getDomain();
              cellBox.grow(idir, -1);
              Box faceBox = surroundingNodes(cellBox, idir);
#ifndef NDEBUG
              IntVectSet ivsIrreg = ebisBoxFine.getIrregIVS(cellBox);
              if(!ivsIrreg.isEmpty())
                {
                  pout() << "EBViscousTensorOp not wired for EBCF yet " << endl;
                  MayDay::Error();
                }
#endif          

              EBFaceFAB fluxFine(ebisBoxFine, cellBox, idir, ncomp);

              FArrayBox&       regFlux = (FArrayBox&)      fluxFine.getSingleValuedFAB();
              const FArrayBox& regPhi  = (const FArrayBox&) phiFine.getSingleValuedFAB();
              const FArrayBox& regGrad = (const FArrayBox&)gradFine.getSingleValuedFAB();
              finerEBAMROp.getFlux(regFlux, regPhi,  regGrad, faceBox, idir, ditf());
              Real scale = m_beta;
              //not wired for EBCF yet.   Not sure what to do with the
              //humongous stencil.
              m_fluxReg.incrementFineRegular(fluxFine, scale, ditf(), interv,  idir,    sit());
            }
        }
    }
}
/*****/
void
EBViscousTensorOp::
AMROperator(LevelData<EBCellFAB>& a_LofPhi,
            const LevelData<EBCellFAB>& a_phiFine,
            const LevelData<EBCellFAB>& a_phi,
            const LevelData<EBCellFAB>& a_phiCoarse,
            bool a_homogeneousBC,
            AMRLevelOp<LevelData<EBCellFAB> >* a_finerOp)
{
  //fillgrad is called in applyop
  cfinterp(a_phi, a_phiCoarse);

  applyOp(a_LofPhi, a_phi, a_homogeneousBC);
  if (a_finerOp != NULL)
    {
      reflux(a_phiFine, a_phi,  a_LofPhi, a_finerOp);
    }
}

/*****/
void
EBViscousTensorOp::
AMROperatorNF(LevelData<EBCellFAB>& a_LofPhi,
              const LevelData<EBCellFAB>& a_phi,
              const LevelData<EBCellFAB>& a_phiCoarse,
              bool a_homogeneousBC)
{
  //fillgrad is called in applyop
  cfinterp(a_phi, a_phiCoarse);

  //apply boundary conditions in applyOp
  this->applyOp(a_LofPhi, a_phi, a_homogeneousBC );
}
/*****/
void
EBViscousTensorOp::
AMRRestrict(LevelData<EBCellFAB>&       a_resCoar,
            const LevelData<EBCellFAB>& a_residual,
            const LevelData<EBCellFAB>& a_correction,
            const LevelData<EBCellFAB>& a_coarseCorr)
{
  LevelData<EBCellFAB> res;

  EBCellFactory ebcellfactTL(m_eblg.getEBISL());
  IntVect ghostVec = a_residual.ghostVect();

  res.define(m_eblg.getDBL(), SpaceDim, ghostVec, ebcellfactTL);
  EBLevelDataOps::setVal(res, 0.0);

  cfinterp(a_correction, a_coarseCorr);
  //API says that we must average(a_residual - L(correction, coarCorrection))
  applyOp(res, a_correction,  true);
  incr(res, a_residual, -1.0);
  scale(res,-1.0);

  //use our nifty averaging operator
  Interval variables(0, SpaceDim-1);
  m_ebAverage.average(a_resCoar, res, variables);
}
void
EBViscousTensorOp::
AMRProlong(LevelData<EBCellFAB>&       a_correction,
           const LevelData<EBCellFAB>& a_coarseCorr)
{
  Interval variables(0, SpaceDim-1);
  m_ebInterp.pwcInterp(a_correction, a_coarseCorr, variables);
}

/*****/
void
EBViscousTensorOp::
AMRUpdateResidual(LevelData<EBCellFAB>&       a_residual,
                  const LevelData<EBCellFAB>& a_correction,
                  const LevelData<EBCellFAB>& a_coarseCorrection)
{
  LevelData<EBCellFAB> r;
  this->create(r, a_residual);
  this->AMRResidualNF(r, a_correction, a_coarseCorrection, a_residual, true);
  this->assign(a_residual, r);
}

/*****/
void
EBViscousTensorOp::
applyOp(LevelData<EBCellFAB>&             a_lhs,
        const LevelData<EBCellFAB>&       a_phi,
        bool                              a_homogeneous)
{

  //contains an exchange
  this->fillGrad(a_phi);

  EBLevelDataOps::setToZero(a_lhs);
  incr( a_lhs, a_phi, m_alpha);
  for (DataIterator dit = m_eblg.getDBL().dataIterator(); dit.ok(); ++dit)
    {
      for(int idir = 0; idir < SpaceDim; idir++)
        {
          a_lhs[dit()].mult((*m_acoef)[dit()], 0, idir, 1);
        }

      for(int idir = 0; idir < SpaceDim; idir++)
        {
          incrOpRegularDir(a_lhs[dit()], a_phi[dit()], a_homogeneous, idir, dit());
        }
      applyOpIrregular(a_lhs[dit()], a_phi[dit()], a_homogeneous, dit());
    }
}

/*****/
void
EBViscousTensorOp::
getFlux(FArrayBox&                    a_flux,
        const FArrayBox&              a_phi,
        const FArrayBox&              a_gradPhi,
        const Box&                    a_faceBox,
        const int&                    a_idir,
        const DataIndex&              a_datInd)
{
  ProblemDomain domain(m_eblg.getDomain());
  Real dx(m_dx);

  CH_assert(a_flux.nComp() == SpaceDim);
  CH_assert( a_phi.nComp() == SpaceDim);

  FArrayBox  faceDiv(a_faceBox, 1);
  FArrayBox  faceGrad(a_faceBox, a_gradPhi.nComp());
  ViscousTensorOp::getFaceDivAndGrad(faceDiv, faceGrad, a_phi, a_gradPhi, domain, a_faceBox, a_idir, dx);
  const FArrayBox& lamFace =(const FArrayBox&)((*m_lambda)[a_datInd][a_idir].getSingleValuedFAB());
  const FArrayBox& etaFace =(const FArrayBox&)((*m_eta   )[a_datInd][a_idir].getSingleValuedFAB());

  ViscousTensorOp::getFluxFromDivAndGrad(a_flux, faceDiv, faceGrad, etaFace, lamFace, a_faceBox, a_idir);
}
/*****/
void
EBViscousTensorOp::
incrOpRegularDir(EBCellFAB&             a_lhs,
                 const EBCellFAB&       a_phi,
                 const bool&            a_homogeneous,
                 const int&             a_dir,
                 const DataIndex&       a_datInd)
{
  const Box& grid = m_eblg.getDBL()[a_datInd];
  Box domainFaces = m_eblg.getDomain().domainBox();
  domainFaces.surroundingNodes(a_dir);
  Box interiorFaces = grid;
  interiorFaces.surroundingNodes(a_dir);
  interiorFaces.grow(a_dir, 1);
  interiorFaces &=  domainFaces;
  interiorFaces.grow( a_dir, -1);

  //do flux difference for interior points
  FArrayBox interiorFlux(interiorFaces, SpaceDim);
  const FArrayBox& grad = (FArrayBox&)(m_grad[a_datInd].getSingleValuedFAB());
  const FArrayBox& phi  = (FArrayBox&)(a_phi.getSingleValuedFAB());
  getFlux(interiorFlux, phi,  grad, interiorFaces, a_dir, a_datInd);

  Box loBox, hiBox, centerBox;
  int hasLo, hasHi;
  EBArith::loHiCenter(loBox, hasLo, hiBox, hasHi, centerBox, m_eblg.getDomain(),grid, a_dir);

  //do the low high center thing
  BaseFab<Real>& reglhs        = a_lhs.getSingleValuedFAB();
  Box dummyBox(IntVect::Zero, IntVect::Unit);
  FArrayBox domainFluxLo(dummyBox, SpaceDim);
  FArrayBox domainFluxHi(dummyBox, SpaceDim);

  RealVect origin = RealVect::Zero;
  Real time = 0.0;
  RealVect dxVect = m_dx*RealVect::Unit;
  if(hasLo==1)
    {
      Box loBoxFace = loBox;
      loBoxFace.shiftHalf(a_dir, -1);
      domainFluxLo.resize(loBoxFace, SpaceDim);
      if(!s_turnOffBCs)
        {
          m_domainBC->getFaceFlux(domainFluxLo,phi,origin,dxVect,a_dir,Side::Lo,a_datInd,time,a_homogeneous);
        }
      else
        {
          //extrapolate to domain flux if there is no bc
          for(BoxIterator boxit(loBoxFace); boxit.ok(); ++boxit)
            {
              const IntVect& iv = boxit();
              IntVect ivn= iv;
              ivn[a_dir]++;
              for(int ivar = 0; ivar < SpaceDim; ivar++)
                {
                  domainFluxLo(iv, ivar) = interiorFlux(ivn, ivar);
                }
            }
        }
    }
  if(hasHi==1)
    {
      Box hiBoxFace = hiBox;
      hiBoxFace.shiftHalf(a_dir, 1);
      domainFluxHi.resize(hiBoxFace, SpaceDim);
      if(!s_turnOffBCs)
        {
          m_domainBC->getFaceFlux(domainFluxHi,phi,origin,dxVect,a_dir,Side::Hi,a_datInd,time,a_homogeneous);
        }
      else
        {
          //extrapolate to domain flux if there is no bc
          //extrapolate to domain flux if there is no bc
          for(BoxIterator boxit(hiBoxFace); boxit.ok(); ++boxit)
            {
              const IntVect& iv = boxit();
              IntVect ivn= iv;
              ivn[a_dir]--;
              for(int ivar = 0; ivar < SpaceDim; ivar++)
                {
                  domainFluxHi(iv, ivar) = interiorFlux(ivn, ivar);
                }
            }
        }
    }
  for(int icomp =0 ; icomp < SpaceDim; icomp++)
    {
      FORT_INCRAPPLYEBVTOP(CHF_FRA1(reglhs,icomp),
                           CHF_CONST_FRA1(interiorFlux, icomp),
                           CHF_CONST_FRA1(domainFluxLo, icomp),
                           CHF_CONST_FRA1(domainFluxHi, icomp),
                           CHF_CONST_REAL(m_beta),
                           CHF_CONST_REAL(m_dx),
                           CHF_BOX(loBox),
                           CHF_BOX(hiBox),
                           CHF_BOX(centerBox),
                           CHF_CONST_INT(hasLo),
                           CHF_CONST_INT(hasHi),
                           CHF_CONST_INT(a_dir));
    }
}
/*****/
void
EBViscousTensorOp::
applyOpIrregular(EBCellFAB&             a_lhs,
                 const EBCellFAB&       a_phi,
                 const bool&            a_homogeneous,
                 const DataIndex&       a_datInd)
{
  RealVect vectDx = m_dx*RealVect::Unit;
  for(int ivar = 0; ivar < SpaceDim; ivar++)
    {
      m_stencil[ivar][a_datInd]->apply(a_lhs, a_phi, false);
    }
//  //debug put in slow apply
//  for(m_vofIter[a_datInd].reset(); m_vofIter[a_datInd].ok();  ++m_vofIter[a_datInd])
//    {
//      const VolIndex& vof = m_vofIter[a_datInd]();
//      for(int ivar = 0; ivar < SpaceDim; ivar++)
//        {
//          VoFStencil vofSten;
//          getVoFStencil(vofSten, vof, a_datInd, ivar);
//          Real lphivar = applyVoFStencil(vofSten, a_phi,  ivar);
//          a_lhs(vof, ivar) = lphivar;
//        }
//    }
//
//  //suspect faces
//  IntVect ivHi(D_DECL(6,5,0));
//  IntVect ivLo(D_DECL(6,4,0));
//
//  VolIndex vofLo(ivLo, 0);
//  VolIndex vofHi(ivHi, 0);
//  FaceIndex dFace(vofLo, vofHi, 1);
//  Real grad[SpaceDim][SpaceDim];
//  Real flux[SpaceDim];
//  for(int ivar = 0; ivar < SpaceDim; ivar++)
//    {
//      VoFStencil fluxSten;
//      getFluxStencil(fluxSten, dFace, a_datInd, ivar);
//      flux[ivar] = applyVoFStencil(fluxSten, a_phi, ivar);
//      for(int ider = 0; ider < SpaceDim; ider++)
//        {
//          VoFStencil gradSten;
//          getGradientStencil(gradSten, ivar, ider, dFace, a_datInd);
//          grad[ivar][ider] = applyVoFStencil(gradSten, a_phi, ivar);
//        }
//    }
//  //end debug

  if(!a_homogeneous)
    {
      const Real factor = 1.0;
      m_ebBC->applyEBFlux(a_lhs, a_phi, m_vofIter[a_datInd], (*m_eblg.getCFIVS()),
                          a_datInd, RealVect::Zero, vectDx, factor,
                          a_homogeneous, 0.0);
    }
  const EBISBox& ebisBox = m_eblg.getEBISL()[a_datInd];
  for(int idir = 0; idir < SpaceDim; idir++)
    {
      for(int comp = 0; comp < SpaceDim; comp++)
        {
          for(m_vofIterDomLo[idir][a_datInd].reset(); m_vofIterDomLo[idir][a_datInd].ok();  ++m_vofIterDomLo[idir][a_datInd])
            {
              Real flux;
              const VolIndex& vof = m_vofIterDomLo[idir][a_datInd]();
              m_domainBC->getFaceFlux(flux,vof,comp,a_phi,
                                      RealVect::Zero,vectDx,idir,Side::Lo, a_datInd, 0.0,
                                      a_homogeneous);

              Real area = ebisBox.sumArea(vof, idir, Side::Lo);
              a_lhs(vof,comp) -= flux*area/m_dx;
            }
          for(m_vofIterDomHi[idir][a_datInd].reset(); m_vofIterDomHi[idir][a_datInd].ok();  ++m_vofIterDomHi[idir][a_datInd])
            {
              Real flux;
              const VolIndex& vof = m_vofIterDomHi[idir][a_datInd]();
              m_domainBC->getFaceFlux(flux,vof,comp,a_phi,
                                      RealVect::Zero,vectDx,idir,Side::Hi,a_datInd,0.0,
                                      a_homogeneous);

              Real area = ebisBox.sumArea(vof, idir, Side::Hi);
              a_lhs(vof,comp) += flux*area/m_dx;
            }
        }
    }
}
/*****/
void
EBViscousTensorOp::
fillGrad(const LevelData<EBCellFAB>&       a_phi)
{
  LevelData<EBCellFAB>& phi = (LevelData<EBCellFAB>&)a_phi;
  phi.exchange(phi.interval());
  const DisjointBoxLayout& grids = a_phi.disjointBoxLayout();
  //compute gradient of phi for parts NOT in ghost cells
  for(DataIterator dit = grids.dataIterator(); dit.ok(); ++dit)
    {
      cellGrad(m_grad[dit()],
               a_phi [dit()],
               grids.get(dit()),
               dit());
    }
  m_grad.exchange();
}
/**/
void
EBViscousTensorOp::
cfinterp(const LevelData<EBCellFAB>&       a_phi,
         const LevelData<EBCellFAB>&       a_phiCoarse)
{
  if(m_hasCoar)
    {
      LevelData<EBCellFAB>& phi = (LevelData<EBCellFAB>&)a_phi;
      EBLevelDataOps::setToZero(m_grad);
      if (a_phiCoarse.isDefined())
        {
          m_interpWithCoarser->coarseFineInterp(phi, m_grad, a_phiCoarse);
        }
    }
}
/**/
void
EBViscousTensorOp::
cellGrad(EBCellFAB&             a_gradPhi,
         const  EBCellFAB&      a_phi,
         const Box&             a_grid,
         const DataIndex&       a_datInd)
{
  CH_assert(a_gradPhi.nComp() == SpaceDim*SpaceDim);
  CH_assert(a_phi.nComp() == SpaceDim);
  const EBISBox& ebisBox = m_eblg.getEBISL()[a_datInd];

  for(int derivDir = 0; derivDir < SpaceDim; derivDir++)
    {
      Box loBox, hiBox, centerBox;
      int hasLo, hasHi;
      EBArith::loHiCenter(loBox, hasLo, hiBox, hasHi, centerBox, m_eblg.getDomain(),a_grid, derivDir);
      for(int phiDir = 0; phiDir < SpaceDim; phiDir++)
        {
          int gradcomp = TensorCFInterp::gradIndex(phiDir,derivDir);
          //
          BaseFab<Real>&       regGrad = a_gradPhi.getSingleValuedFAB();
          const BaseFab<Real>& regPhi  =     a_phi.getSingleValuedFAB();
          FORT_CELLGRADEBVTOP(CHF_FRA1(regGrad, gradcomp),
                              CHF_CONST_FRA1(regPhi, phiDir),
                              CHF_CONST_REAL(m_dx),
                              CHF_BOX(loBox),
                              CHF_BOX(hiBox),
                              CHF_BOX(centerBox),
                              CHF_CONST_INT(hasLo),
                              CHF_CONST_INT(hasHi),
                              CHF_CONST_INT(derivDir));

          for(m_vofIter[a_datInd].reset(); m_vofIter[a_datInd].ok(); ++m_vofIter[a_datInd])
            {
              const VolIndex& vof =  m_vofIter[a_datInd]();
              Vector<FaceIndex> loFaces = ebisBox.getFaces(vof, derivDir, Side::Lo);
              Vector<FaceIndex> hiFaces = ebisBox.getFaces(vof, derivDir, Side::Hi);
              bool hasLoPt = ((loFaces.size() == 1)  && (!loFaces[0].isBoundary()));
              bool hasHiPt = ((hiFaces.size() == 1)  && (!hiFaces[0].isBoundary()));
              Real valLo=0, valHi=0;
              if(hasLoPt) valLo = a_phi(loFaces[0].getVoF(Side::Lo), phiDir);
              if(hasHiPt) valHi = a_phi(hiFaces[0].getVoF(Side::Hi), phiDir);

              Real valCe = a_phi(vof, phiDir);

              if(hasLoPt && hasHiPt)
                {
                  a_gradPhi(vof, gradcomp) = (valHi - valLo)/(2.*m_dx);
                }
              else if(hasHiPt)
                {
                  a_gradPhi(vof, gradcomp) = (valHi - valCe)/(m_dx);
                }
              else if(hasLoPt)
                {
                  a_gradPhi(vof, gradcomp) = (valCe - valLo)/(m_dx);
                }
              else
                {
                  a_gradPhi(vof, gradcomp) = 0.0;
                }
            }
        }
    }
}

#include "NamespaceFooter.H"

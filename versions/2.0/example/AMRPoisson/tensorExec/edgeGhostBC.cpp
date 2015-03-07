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


/*****************/
/*****************/


#include "edgeGhostBC.H"
#include "edgeGhostBCF_F.H"
#include "MayDay.H"
#include "UsingNamespace.H"

BoxEdgeBC::BoxEdgeBC() : BoxGhostBC()
{
}

BoxEdgeBC::BoxEdgeBC(int a_dir, Side::LoHiSide a_sd, const Interval& a_comps)
  : BoxGhostBC()
{
  define(a_dir, a_sd, a_comps);
}

BoxEdgeBC::~BoxEdgeBC()
{
}


BoxGhostBC*
BoxEdgeBC::new_boxghostbc() const 
{
  BoxEdgeBC* newop = new BoxEdgeBC();
  if (newop == NULL) {
    MayDay::Error("Out of memory in BoxEdgeBC::new_boxghostBC");
  } 
  // set new edge type to be the same as this one
  newop->ixType(m_ix_type);
  
  return static_cast<BoxGhostBC*> (newop);
}



void
BoxEdgeBC::define(int a_dir, Side::LoHiSide a_sd) 
{
  Interval comps(0,0);
  define(a_dir, a_sd, comps);
}

void
BoxEdgeBC::define(int a_dir, Side::LoHiSide a_sd, const Interval& a_comps)
{
  BoxGhostBC::define(a_dir, a_sd, a_comps);
}


void
BoxEdgeBC::define(int a_dir, Side::LoHiSide a_sd, const Interval& a_comps,
		  IndexType& a_ixtype)
{
  BoxGhostBC::define(a_dir, a_sd, a_comps);
  m_ix_type = a_ixtype;
}

const IndexType&
BoxEdgeBC::ixType() const
{
  return m_ix_type;
}

void
BoxEdgeBC::ixType(const IndexType& a_ixtype) 
{
  m_ix_type = a_ixtype;
}

void
BoxEdgeBC::applyInhomogeneousBCs(FArrayBox& a_state, 
				 const Box& a_domain,
				 Real a_dx) const
{
  ProblemDomain physdomain(a_domain);
  applyInhomogeneousBCs(a_state, physdomain, a_dx);
}

void
BoxEdgeBC::applyInhomogeneousBCs(FArrayBox& a_state, 
				 const ProblemDomain& a_domain,
				 Real a_dx) const
{
  // only do anything if not periodic
  if (!a_domain.isPeriodic(m_direction))
  {

    // apply edge-centered bc to edge-centered data at box boundary
    Box bx = a_state.box();
    CH_assert(bx.ixType() == m_ix_type);
    // convert box to cell-centered box if appropriate
    // (since all of our domain testing is done with CC boxes)
    if (m_ix_type.ixType(m_direction) == IndexType::NODE) {
      bx.enclosedCells(m_direction);
      // ensure that edges on the boundary get included
      bx.grow(m_direction,1);
    }
    bx.grow(-1);
    Box bc_box;
    int idir = m_direction;
    bool isbc;
    if (m_side == Side::Lo) {
      isbc = (bx.smallEnd(idir) <= a_domain.domainBox().smallEnd(idir));
      if (isbc) {
        Box testBox = a_domain.domainBox();
        testBox.surroundingNodes(m_direction);
        int ichop = testBox.smallEnd(m_direction);
        bc_box = a_state.box();
        //bc_box.chop(m_direction, ichop);
        bc_box.setBig(m_direction,ichop);
      }
    }
    else if(m_side == Side::Hi) {
      isbc = (bx.bigEnd(idir) >= a_domain.domainBox().bigEnd(idir));
      if (isbc) {
        Box testBox = a_domain.domainBox();
        testBox.surroundingNodes(m_direction);
        int ichop = testBox.bigEnd(m_direction);
        bc_box = a_state.box();
        //bc_box = chop_box.chop(m_direction,ichop);
        bc_box.setSmall(m_direction,ichop);
      }
    } else {
      cerr << "DomainGhostBC::applyghostbc: bogus side" << endl;
      abort();
    }
    
    if(isbc) {
      int ncomps = a_state.nComp();
      FArrayBox neumfac(bc_box, ncomps);
      FArrayBox dircfac(bc_box, ncomps);
      FArrayBox inhmval(bc_box, ncomps);
      fillBCValues(neumfac, dircfac, inhmval, a_dx, a_domain);
      applyBCs(bc_box, a_state, neumfac, dircfac, inhmval, a_dx);
    }
    
  }
}

void
BoxEdgeBC::applyHomogeneousBCs(FArrayBox& a_state, 
			       const Box& a_domain,
			       Real a_dx) const
{
  ProblemDomain physdomain(a_domain);
  applyHomogeneousBCs(a_state, a_domain, a_dx);
}



void
BoxEdgeBC::applyHomogeneousBCs(FArrayBox& a_state, 
			       const ProblemDomain& a_domain,
			       Real a_dx) const
{
  
  // only do anything if not periodic
  if (!a_domain.isPeriodic(m_direction))
  {
    // apply edge centered bc to edge-centered data at domain boundary
    Box bx = a_state.box();
    CH_assert(bx.ixType() == m_ix_type);
    // convert box to cell-centered box if appropriate
    // (since all of our domain testing is done with CC boxes)
    if (m_ix_type.ixType(m_direction) == IndexType::NODE) {
      bx.enclosedCells(m_direction);
    }
    // this is inappropriate for edge-centering
    //bx.grow(-1);
    Box bc_box;
    int idir = m_direction;
    bool isbc;
    if (m_side == Side::Lo) {
      isbc = (bx.smallEnd(idir) == a_domain.domainBox().smallEnd(idir));
      if (isbc) {
        // chopping doesn't work so well for edge-centered boxes --
        // do this instead
        bc_box = bdryLo(a_state.box(), m_direction);
        int stateDirHi = a_domain.domainBox().smallEnd(m_direction);
        if (stateDirHi > bc_box.bigEnd(m_direction)) {
          bc_box.setBig(stateDirHi, m_direction);
        }
      }
    }
    else if(m_side == Side::Hi) {
      isbc = (bx.bigEnd(idir) == a_domain.domainBox().bigEnd(idir));
      if (isbc) {
        // do similar operation for high end
        bc_box = bdryHi(a_state.box(),m_direction);
        int stateDirLo = a_domain.domainBox().bigEnd(m_direction)+1;
        if (stateDirLo < bc_box.smallEnd(m_direction)) {
          bc_box.setSmall(stateDirLo, m_direction);
        }
      }
    } else {
      cerr << "DomainGhostBC::applyghostbc: bogus side" << endl;
      abort();
    }
    
    if(isbc) {
      int ncomps = a_state.nComp();
      FArrayBox neumfac(bc_box, ncomps);
      FArrayBox dircfac(bc_box, ncomps);
      FArrayBox inhmval(bc_box, ncomps);
      fillBCValues(neumfac, dircfac, inhmval, a_dx, a_domain);
      inhmval.setVal(0.0);
      applyBCs(bc_box, a_state, neumfac, dircfac, inhmval, a_dx);
    }
    
  }
}


void
BoxEdgeBC::fillBCValues(FArrayBox& a_neumfac,
			FArrayBox& a_dircfac,
			FArrayBox& a_inhmval,
			Real a_dx,
			const Box& a_domain) const
{
  // shouldn't be in this function
  MayDay::Error("Should not be in BoxEdgeBC::fillBCValues!");
}



void
BoxEdgeBC::fillBCValues(FArrayBox& a_neumfac,
			FArrayBox& a_dircfac,
			FArrayBox& a_inhmval,
			Real a_dx,
			const ProblemDomain& a_domain) const
{
  // shouldn't be in this function
  MayDay::Error("Should not be in BoxEdgeBC::fillBCValues!");
}


void 
BoxEdgeBC::applyBCs(const Box& a_bcbox,
		    FArrayBox& a_state,
		    const FArrayBox& a_neumfac,
		    const FArrayBox& a_dircfac,
		    const FArrayBox& a_inhmval,
		    Real a_dx) const
{

  int iSide = sign(m_side);
  int iDir = m_direction;
  int startcomp = m_components.begin();
  int endcomp = m_components.end();

  FORT_BOXEDGEBC(CHF_FRA(a_state),
		 CHF_CONST_FRA(a_neumfac),
		 CHF_CONST_FRA(a_dircfac),
		 CHF_CONST_FRA(a_inhmval),
		 CHF_BOX(a_bcbox),
		 CHF_CONST_INT(iDir),
		 CHF_CONST_INT(iSide),
		 CHF_CONST_REAL(a_dx),
		 CHF_CONST_INT(startcomp),
		 CHF_CONST_INT(endcomp));
}


edgeDirichletBC::edgeDirichletBC() : BoxEdgeBC()
{
  m_bcVal = 0.0;
}

edgeDirichletBC::edgeDirichletBC(int a_dir, 
					 Side::LoHiSide a_sd,
					 const Interval& a_comps) 
  : BoxEdgeBC(a_dir, a_sd, a_comps)
{
  m_bcVal = 0.0;
}

edgeDirichletBC::~edgeDirichletBC() 
{
}


void 
edgeDirichletBC::setBCVal(const Real a_bcVal)
{
  m_bcVal = a_bcVal;
}

Real 
edgeDirichletBC::BCVal() const
{
  return m_bcVal;
}


void 
edgeDirichletBC::fillBCValues(FArrayBox& a_neumfac,
				  FArrayBox& a_dircfac,
				  FArrayBox& a_inhmval,
				  Real a_dx,
				  const Box& a_domain) const 
{
  ProblemDomain physdomain(a_domain);
  fillBCValues(a_neumfac, a_dircfac, a_inhmval, a_dx, physdomain);
}


void 
edgeDirichletBC::fillBCValues(FArrayBox& a_neumfac,
				  FArrayBox& a_dircfac,
				  FArrayBox& a_inhmval,
				  Real a_dx,
				  const ProblemDomain& a_domain) const 
{
  CH_assert(!a_neumfac.box().isEmpty());
  CH_assert(!a_dircfac.box().isEmpty());
  CH_assert(!a_inhmval.box().isEmpty());
  a_neumfac.setVal(0.0);
  a_dircfac.setVal(1.0);
  a_inhmval.setVal(m_bcVal);
}


BoxGhostBC*
edgeDirichletBC::new_boxghostbc() const 
{
  edgeDirichletBC* newop = new edgeDirichletBC();
  if (newop == NULL) {
    MayDay::Error("Out of memory in edgeDirichletBC::new_boxghostBC");
  } 

  // set new edge type to be the same as this one
  newop->ixType(ixType());
  newop->setBCVal(m_bcVal);

  return static_cast<BoxGhostBC*> (newop);
}


NoOpBC::NoOpBC() :BoxGhostBC()
{
}


NoOpBC::NoOpBC(int a_dir, Side::LoHiSide a_sd) :BoxGhostBC(a_dir, a_sd)
{
}

NoOpBC::NoOpBC(int a_dir, Side::LoHiSide a_sd, const Interval& a_comps)
  : BoxGhostBC(a_dir, a_sd, a_comps)
{
}

NoOpBC::~NoOpBC() 
{
}

BoxGhostBC*
NoOpBC::new_boxghostbc() const
{
  NoOpBC* newop = new NoOpBC();
  if (newop == NULL) {
    MayDay::Error("Out of memory in edgeDirichletBC::new_boxghostBC");
  } 

  return static_cast<BoxGhostBC*> (newop);
}


void 
NoOpBC::define(int a_dir, Side::LoHiSide a_sd) 
{
  BoxGhostBC::define(a_dir, a_sd);
}

void 
NoOpBC::define(int a_dir, Side::LoHiSide a_sd, const Interval& a_comps) 
{
  BoxGhostBC::define(a_dir, a_sd, a_comps);
}

void
NoOpBC::applyInhomogeneousBCs(FArrayBox& a_state,
			      const Box& a_domain,
			      Real a_dx) const
{
  // do nothing here
}

void
NoOpBC::applyHomogeneousBCs(FArrayBox& a_state, 
		      const Box& a_domain,
		      Real a_dx) const
{
  // true to its name, this class does nothing here
}


void
NoOpBC::fillBCValues(FArrayBox& a_neumfac,
	       FArrayBox& a_dircfac,
	       FArrayBox& a_inhmval,
	       Real a_dx,
	       const Box& a_domain) const
{
  // this also doesn't do anything...
}



void
NoOpBC::applyInhomogeneousBCs(FArrayBox& a_state,
			      const ProblemDomain& a_domain,
			      Real a_dx) const
{
  // do nothing here
}

void
NoOpBC::applyHomogeneousBCs(FArrayBox& a_state, 
		      const ProblemDomain& a_domain,
		      Real a_dx) const
{
  // true to its name, this class does nothing here
}


void
NoOpBC::fillBCValues(FArrayBox& a_neumfac,
	       FArrayBox& a_dircfac,
	       FArrayBox& a_inhmval,
	       Real a_dx,
	       const ProblemDomain& a_domain) const
{
  // this also doesn't do anything...
}


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

#ifndef _edgeGhostBC_H_
#define _edgeGhostBC_H_

#include "GhostBC.H"
#include "UsingNamespace.H"

/// class to do edge-centered physical boundary conditions
/** BoxEdgeBC is designed to do physical BC's on edge centered
    data.  It is derived from BoxGhostBC so that it can be used
    seamlessly with BoxGhostBCs.  */

class BoxEdgeBC: public BoxGhostBC
{
public:
  /// 
  virtual ~BoxEdgeBC();
  
  ///
  BoxEdgeBC();

  ///  sets interval to be (0,0)
  BoxEdgeBC(int a_dir, Side::LoHiSide a_sd);

  ///
  BoxEdgeBC(int a_dir, Side::LoHiSide a_sd, const Interval& a_comps);

  /// 
  BoxEdgeBC(int a_dir, Side::LoHiSide a_sd, const Interval& a_comps, 
	    IndexType& a_ixtype);

  /// returns the type of centering of this edge-centered BC
  const IndexType& ixType() const;
  
  ///
  void ixType(const IndexType& a_ixType);


protected:

  /// virtual constructor workaround
  virtual BoxGhostBC* new_boxghostbc() const;

  /// sets interval to be (0,0)
  void define(int a_dir, Side::LoHiSide a_sd);

  /// 
  void define(int a_dir, Side::LoHiSide a_sd, const Interval& a_comps);

  /// 
  void define(int a_dir, Side::LoHiSide a_sd, const Interval& a_comps,
	      IndexType& a_ixtype);

  ///
  virtual void
  applyInhomogeneousBCs(FArrayBox& a_state,
			const ProblemDomain& a_domain,
			Real a_dx) const;
  ///
  virtual void
  applyHomogeneousBCs(FArrayBox& a_state, 
		      const ProblemDomain& a_domain,
		      Real a_dx) const;

  ///
  virtual void
  fillBCValues(FArrayBox& a_neumfac,
	       FArrayBox& a_dircfac,
	       FArrayBox& a_inhmval,
	       Real a_dx,
	       const ProblemDomain& a_domain) const;


  /// deprecated interface
  virtual void
  applyInhomogeneousBCs(FArrayBox& a_state,
			const Box& a_domain,
			Real a_dx) const;
  /// deprecated interface
  virtual void
  applyHomogeneousBCs(FArrayBox& a_state, 
		      const Box& a_domain,
		      Real a_dx) const;

  /// deprecated interface
  virtual void
  fillBCValues(FArrayBox& a_neumfac,
	       FArrayBox& a_dircfac,
	       FArrayBox& a_inhmval,
	       Real a_dx,
	       const Box& a_domain) const;

  ///
  virtual void 
  applyBCs(const Box& a_bcbox, 
           FArrayBox& a_state,
           const FArrayBox& a_neumfac,
           const FArrayBox& a_dircfac,
           const FArrayBox& a_inhmval,
           Real a_dx) const;


private:
  BoxEdgeBC(const BoxEdgeBC&) {;}
  virtual void operator=(const BoxEdgeBC&) {;}
  
  IndexType m_ix_type;
  
};


/// class to do Dirichlet edge-centered BC's 
/** this class implements the edge-centered instantiation
    of the edge-centered homogeneous BC's */
class edgeDirichletBC : public BoxEdgeBC
{
public:
  friend class DomainGhostBC;

  /// 
  ~edgeDirichletBC();
  
  ///
  edgeDirichletBC();

  /// sets interval to be (0,0)
  edgeDirichletBC(int a_dir, Side::LoHiSide a_sd);

  /// 
  edgeDirichletBC(int a_dir, Side::LoHiSide a_sd, const Interval& a_comps);

  /// set edge-BC value (single-valued)
  void setBCVal(const Real a_bcVal);

  ///
  Real BCVal() const;

protected:
  virtual void 
  fillBCValues(FArrayBox& a_neumfac,
	       FArrayBox& a_dircfac,
	       FArrayBox& a_inhmval,
	       Real a_dx,
	       const ProblemDomain& a_domain) const;


  virtual void 
  fillBCValues(FArrayBox& a_neumfac,
	       FArrayBox& a_dircfac,
	       FArrayBox& a_inhmval,
	       Real a_dx,
	       const Box& a_domain) const;

  ///
  virtual  BoxGhostBC* new_boxghostbc() const;

  /// 
  Real m_bcVal;


private:
  edgeDirichletBC(const edgeDirichletBC& dirin) {;}
  void operator=(const edgeDirichletBC& dirin) {;}
};
    

/// class to do nothing on an edge
/** This class functions as a placeholder -- if you don't want
    to change the values on this edge, but it's a physical boundary,
    it's the clas you want! */
class NoOpBC : public BoxGhostBC
{
public:
  
  /// 
  NoOpBC();

  ///
  NoOpBC(int a_dir, Side::LoHiSide a_sd);

  ///
  NoOpBC(int a_dir, Side::LoHiSide a_sd, const Interval& a_comps);

  ///
  virtual ~NoOpBC();

protected:

  /// virtual constructor workaround
  virtual BoxGhostBC* new_boxghostbc() const;

  /// sets interval to be (0,0)
  void define(int a_dir, Side::LoHiSide a_sd);

  /// 
  void define(int a_dir, Side::LoHiSide a_sd, const Interval& a_comps);

  ///
  virtual void
  applyInhomogeneousBCs(FArrayBox& a_state,
			const ProblemDomain& a_domain,
			Real a_dx) const;

  /// deprecated interface using Box instead of ProblemDomain
  virtual void
  applyInhomogeneousBCs(FArrayBox& a_state,
			const Box& a_domain,
			Real a_dx) const;
  ///
  virtual void
  applyHomogeneousBCs(FArrayBox& a_state, 
		      const ProblemDomain& a_domain,
		      Real a_dx) const;

  /// deprecated interface using Box instead of ProblemDomain
  virtual void
  applyHomogeneousBCs(FArrayBox& a_state, 
		      const Box& a_domain,
		      Real a_dx) const;

  
  ///
  virtual void
  fillBCValues(FArrayBox& a_neumfac,
	       FArrayBox& a_dircfac,
	       FArrayBox& a_inhmval,
	       Real a_dx,
	       const ProblemDomain& a_domain) const;


  ///
  virtual void
  fillBCValues(FArrayBox& a_neumfac,
	       FArrayBox& a_dircfac,
	       FArrayBox& a_inhmval,
	       Real a_dx,
	       const Box& a_domain) const;

private:
  NoOpBC(const NoOpBC&) {;}
  virtual void operator =(const NoOpBC&) {;}


};
#endif




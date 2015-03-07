#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

#include "EBArith.H"
#include "InflowOutflowAdvectBC.H"
#include "DirichletPoissonEBBC.H"

#include "NamespaceHeader.H"

/****************************/
void
InflowOutflowAdvectBC::
fluxBC(EBFluxFAB&            a_primGdnv,
       const EBCellFAB&      a_primCenter,
       const EBCellFAB&      a_primExtrap,
       const Side::LoHiSide& a_side,
       const Real&           a_time,
       const EBISBox&        a_ebisBox,
       const DataIndex&      a_dit,
       const Box&            a_box,
       const Box&            a_faceBox,
       const int&            a_dir)
{
  CH_assert(m_isDefined);
  CH_assert(!m_domain.isPeriodic(a_dir));

  //gets face centered region of data
  Box FBox = a_primGdnv[a_dir].getRegion();
  Box cellBox = a_faceBox;
  CH_assert(a_primGdnv[a_dir].nComp() == 1);

  // Determine which side and thus shifting directions
  int isign = sign(a_side);
  cellBox.shiftHalf(a_dir,isign);

  //figure out if we are on a wall and if its an inflow and the inflow value
  bool isInflow   = ((a_side == Side::Lo) && (a_dir==m_flowDir));
  bool isOutflow  = ((a_side == Side::Hi) && (a_dir==m_flowDir));

  // Is there a domain boundary next to this grid
  if (!m_domain.contains(cellBox))
    {
      cellBox &= m_domain;
      // Find the strip of cells next to the domain boundary
      Box bndryBox = adjCellBox(cellBox, a_dir, a_side, 1);

      // Shift things to all line up correctly
      bndryBox.shift(a_dir,-isign);

      IntVectSet ivs(bndryBox);
      for(VoFIterator vofit(ivs, a_ebisBox.getEBGraph()); vofit.ok(); ++vofit)
        {
          const VolIndex& vof = vofit();

          Vector<FaceIndex> bndryFaces = a_ebisBox.getFaces(vof, a_dir, a_side);
          for(int iface= 0; iface < bndryFaces.size(); iface++)
            {
              //this all will work for the case of spatially varying inflow
              //for inhomogeneous other faces, this is wrong.
              const FaceIndex& face = bndryFaces[iface];
              Real velComp = 1.e99;
              if(isInflow && (m_velComp == m_flowDir))
                {
                  velComp = m_inflowVel;
                }
              else if(isInflow && (m_velComp != m_flowDir))
                {
                  velComp = 0.0;
                }
              else if(isOutflow)
                {
                  velComp = a_primExtrap(vof, 0);
                }
              else //solid wall
                {
                  if(a_dir == m_velComp)
                    {
                      velComp = 0.0;
                    }
                  else
                    {
                      velComp = a_primExtrap(vof, 0);
                    }
                }
              a_primGdnv[a_dir](face, 0) = velComp;
            }
        }
    }
}

#include "NamespaceFooter.H"

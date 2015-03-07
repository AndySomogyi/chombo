C     -*- Mode: Fortran; -*-
c-----------------------------------------------------------------
c     Ravi Samtaney
c     Copyright 2002
c     Princeton Plasma Physics Laboratory
c     All Rights Reserved
c-----------------------------------------------------------------
c     $Log: init4.F,v $
c     Revision 1.1  2005/05/12 19:22:45  ligocki
c     Moved current "{exec,src}MHD" to "{exec,src}IdealMHD"...
c
c     Revision 1.2  2004/10/06 20:52:56  samtaney
c     Intermediate checkin
c
c     Revision 1.1.1.1  2004/07/09 19:27:04  samtaney
c     Original source for Unsplit Ideal MHD based on new Chombo interface
c
c     Revision 1.1  2002/06/13 22:16:18  samtaney
c     Original unsplit 8-wave Ideal MHD source (includes B filter).
c
c-----------------------------------------------------------------
c     Subroutine to set up the initial conditions for the Orszag
c     Tang vortex. 
c---------------------------------------------------------------------------
c     
      subroutine OrszagTangInitialConditions(ux,
     &     iulo,iuhi,julo,juhi,
     &     ilo,ihi,jlo,jhi,
     &     klo,khi,
     &     adx,ady,adz)
c----------------------------------------------------------------------
c
      implicit none
c
#include "EosCommon.fh"
c
c     Incoming arguments
      integer ilo,ihi,jlo,jhi,klo,khi
      integer iulo,iuhi,julo,juhi
      real*8 adx,ady,adz
      real*8 ux(iulo:iuhi,
     &     julo:juhi,
     &     1:1,8)
c     
c     Local variables
      integer i,j,k
      real*8 xl,yl,xr,yr
      real*8 rho,press, Bi, Bj, Bk, u,v,w
      real*8 press0,rho0
      real*8 twopi
c     
      real*8 x,y
c     
      twopi=8.D0*datan(1.D0)
c     
      xl=0.D0
      xr=1.D0
      yl=0.D0
      yr=1.D0       
c     
      press0=(10.D0/3.D0)/(2.D0)
      rho0=gamma*press0
      rho=rho0
      press=press0
      bk=0.D0
      w=0.D0
c     
      do k=klo,khi,1
         do j=jlo,jhi,1
            y = yl+(j - 0.5D0)*ady;
            do i=ilo,ihi,1
               x = xl+(i - 0.5D0)*adx;
               u=-dsin(twopi*y)
               v=dsin(twopi*x)
               bi=u
               bj=dsin(2.D0*twopi*x)
               ux(i,j,k,1)=rho
               ux(i,j,k,2)=rho*u
               ux(i,j,k,3)=rho*v
               ux(i,j,k,4)=rho*w
               ux(i,j,k,5)=bi
               ux(i,j,k,6)=bj
               ux(i,j,k,7)=Bk
               ux(i,j,k,8)=press/(gamma-1.D0)+
     &              0.5*rho*(u*u+v*v+w*w)+0.5D0*(bi*bi+bj*bj+bk*bk)
            enddo
         enddo
      enddo
c     
      return
      end








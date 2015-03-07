C     -*- Mode: Fortran; -*-
c-----------------------------------------------------------------
c     Ravi Samtaney
c     Copyright 2002
c     Princeton Plasma Physics Laboratory
c     All Rights Reserved
c-----------------------------------------------------------------
c     $Log: init.F,v $
c     Revision 1.1  2005/05/12 19:33:37  ligocki
c     Copied "IdealMHD" to "ResistiveMHD" as a starting point...
c
c     Revision 1.1  2005/05/12 19:22:45  ligocki
c     Moved current "{exec,src}MHD" to "{exec,src}IdealMHD"...
c
c     Revision 1.2  2004/10/06 20:52:56  samtaney
c     Intermediate checkin
c
c     Revision 1.1.1.1  2004/07/09 19:27:04  samtaney
c     Original source for Unsplit Ideal MHD based on new Chombo interface
c
c     Revision 1.2  2002/07/09 19:00:59  samtaney
c     Intermediate commit of files towards a single 2D/3D Ideal MHD code.
c
c     Revision 1.1  2002/06/13 22:16:18  samtaney
c     Original unsplit 8-wave Ideal MHD source (includes B filter).
c
c-----------------------------------------------------------------
c     Subroutine to set up the initial conditions for the
c     "Rotor" problem. 
c---------------------------------------------------------------------------
c     
      subroutine RotorInitialConditions(ux,
     &     iulo,iuhi,julo,juhi,kulo,kuhi,
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
      integer iulo,iuhi,julo,juhi,kulo,kuhi
      real*8 adx,ady,adz
      real*8 ux(iulo:iuhi,
     &     julo:juhi,
     &     kulo:kuhi,8)
c     
c     Local variables
      integer i,j,k
      real*8 xl,yl,xr,yr,zl,zr
      real*8 rho,press, Bi, Bj, Bk, u,v,w
      real*8 rho1,rho2,rad0,rad1,rad,fac,v0
c     
      real*8 x,y,z
c     
      write(6,*) 'ROTOR INIT',ilo,ihi,jlo,jhi
      write(6,*) 'ROTOR INIT',iulo,iuhi,julo,juhi
      write(6,*) 'ROTOR INIT',adx,ady,adz,8
      write(6,*) 'ROTOR INIT',gamma,rgas,wmol
c     
c$$$      xl=-1.D0
c$$$      xr=1.D0
c$$$      yl=-1.D0
c$$$      yr=1.D0       
c$$$      zl=-1.D0
c$$$      zr=1.D0       
c
      xl=-0.5D0
      xr=0.5D0
      yl=-0.5D0
      yr=0.5D0       

      rho1=10.D0
      rho2=1.D0
      press=0.5D0
      v0=1.D0
      rad0=0.1D0
      rad1=0.115D0
      w=0.D0
      Bi=2.5D0/dsqrt(16.D0*datan(1.D0))
      bj=0.00D0
      bk=0.00D0
c     
      do k=klo,khi,1
         z = zl+(k - 0.5D0)*adz
         do j=jlo,jhi,1
            y = yl+(j - 0.5D0)*ady
            do i=ilo,ihi,1
               x = xl+(i - 0.5D0)*adx
               rad=dsqrt(x*x+y*y)
               fac=(rad1-rad)/(rad1-rad0)
               if(rad.le.rad0) then
                  rho=rho1
                  u=-v0*(y)/rad0
                  v=v0*(x)/rad0
c                  w=v0*(z)/rad0
                 else if(rad.gt.rad0.and.rad.le.rad1) then
                    rho=(rho1-rho2)*fac+rho2
                    u=-v0*fac*(y)/rad0
                    v=v0*fac*(x)/rad0
                 else 
                  rho=rho2
                  u=0.D0
                  v=0.D0
c                  w=0.D0
               endif
c               u=v0
c               v=0.D0
               ux(i,j,k,1)=rho
               ux(i,j,k,2)=rho*u
               ux(i,j,k,3)=rho*v
               ux(i,j,k,4)=rho*w
               ux(i,j,k,5)=Bi
               ux(i,j,k,6)=Bj
               ux(i,j,k,7)=Bk
               ux(i,j,k,8)=press/(gamma-1.D0)+
     &              0.5*rho*(u*u+v*v+w*w)+0.5D0*(bi*bi+bj*bj+bk*bk)
c               write(6,*) 'INIT',i,j,ux(i,j,k,1),ux(i,j,k,2),ux(i,j,k,3),
c     &              ux(i,j,k,8)
            enddo
         enddo
      enddo
      
      return
      end








C -*- Mode: Fortran; -*-
c-----------------------------------------------------------------
c       Ravi Samtaney
c       Copyright 2001
c       Princeton Plasma Physics Laboratory
c       All Rights Reserved
c-----------------------------------------------------------------
c       $Log: init2b.F,v $
c       Revision 1.1  2005/05/24 21:44:03  samtaney
c       Initial Resistive MHD implementation. Implicit handling of diffusion terms.
c
c       Revision 1.2  2004/07/09 21:05:33  samtaney
c       Fixed density and pressure.
c
c       Revision 1.1  2004/04/05 22:42:46  samtaney
c       Original source. For GEM Reconnection problem.
c
c       Revision 1.1  2002/06/13 22:16:18  samtaney
c       Original unsplit 8-wave Ideal MHD source (includes B filter).
c
c-----------------------------------------------------------------
c Subroutine to set up the initial conditions
c---------------------------------------------------------------------------
c
c----------------------------------------------------------------------
      subroutine ReconInitialConditions(ux,
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
c
      real*8 adx,ady,adz
      real*8 ux(iulo:iuhi,
     &     julo:juhi,
     &     kulo:kuhi,8)
c
	real*8 psi(ilo-1:ihi+1,
     &                        jlo-1:jhi+1,
     &                        1:1)
	integer i,j,k,l
	real*8 dtodx,dtody
        real*8 rho,press, Bi, Bj, Bk, u,v,w
        real*8 dxi,dyi
        real*8 xl,yl,xr,yr
c
        real*8 x,y,z
        real*8 pi, kx,ky,lamda
        real*8 psi0
c
        xl=-12.8D0
        xr=12.8D0
        yl=-6.4D0
        yr=6.4D0
c
	pi=4.D0*datan(1.D0)
        dxi=0.5D0/adx
        dyi=0.5D0/ady
c
        psi0=0.1D0
c        psi0=0.0D0
c
        rho=1.D0
        press=0.2D0
	u=0.D0
	v=0.D0
	w=0.D0
c
	kx=pi/(xr)
	ky=pi/(2.D0*yr)
c        
        do k=klo,khi,1
           do j=jlo,jhi,1
              y = yl+(j - 0.5D0)*ady;
              do i=ilo,ihi,1
                 x = xl+(i - 0.5D0)*adx;
c                 Bi=tanh(y/0.5D0)+psi0*ky*dcos(kx*x)*dsin(ky*y)
c                 Bj=-psi0*kx*dsin(kx*x)*dcos(ky*y)
                 Bi=tanh(y/0.5D0)-psi0*ky*dcos(kx*x)*dsin(ky*y)
                 Bj=psi0*kx*dsin(kx*x)*dcos(ky*y)
		 Bk=0.D0
                 rho=0.2D0+(1-tanh(y/0.5D0)**2)
                 press=0.5D0*rho
                 ux(i,j,k,1)=rho
                 ux(i,j,k,2)=rho*u
                 ux(i,j,k,3)=rho*v
                 ux(i,j,k,4)=rho*w
                 ux(i,j,k,5)=Bi
                 ux(i,j,k,6)=Bj
                 ux(i,j,k,7)=Bk
                 ux(i,j,k,8)=press/(gamma-1.D0)+
     &                0.5*rho*(u*u+v*v+w*w)+0.5D0*(bi*bi+bj*bj+bk*bk)
              enddo
           enddo
        enddo
c     
      return
      end








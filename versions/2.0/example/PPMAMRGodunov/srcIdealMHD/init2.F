C -*- Mode: Fortran; -*-
c-----------------------------------------------------------------
c       Ravi Samtaney
c       Copyright 2001
c       Princeton Plasma Physics Laboratory
c       All Rights Reserved
c-----------------------------------------------------------------
c       $Log: init2.F,v $
c       Revision 1.1  2005/05/12 19:22:45  ligocki
c       Moved current "{exec,src}MHD" to "{exec,src}IdealMHD"...
c
c       Revision 1.2  2004/10/06 20:52:56  samtaney
c       Intermediate checkin
c
c       Revision 1.1.1.1  2004/07/09 19:27:04  samtaney
c       Original source for Unsplit Ideal MHD based on new Chombo interface
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
c
c$$$	gamma=1.667D0
c$$$	rgas=8317.D0
c$$$	wmol=1.D0
c$$$	mu=1.D-03
c$$$	eta=1.D-03
c$$$	kappa=1.D-03
c
        xl=-1.0D0
        xr=1.D0
        yl=0.D0
c        yr=0.5D0
        yr=1.D0
c
	pi=4.D0*datan(1.D0)
        dxi=0.5D0/adx
        dyi=0.5D0/ady
c
        rho=1.D0
        press=0.2D0
	u=0.D0
	v=0.D0
	w=0.D0
c
	kx=3.D0*pi/2.D0
	ky=pi
c	kx=3.D0*pi
c	ky=4.D0*pi
	lamda=dsqrt(kx*kx+ky*ky)
c        
c        do k=klo,khi,1
c           do j=jlo-1,jhi+1,1
c              y = yl+(j - 0.5D0)*ady;
c              do i=ilo-1,ihi+1,1
c                 x = xl+(i - 0.5D0)*adx;
c                 psi(i,j,k)=cos(kx*x)*sin(ky*y)
c              enddo
c           enddo
c        enddo
c
c$$$        do k=klo,khi,1
c$$$           do j=jlo,jhi,1
c$$$              y = yl+(j - 0.5D0)*ady;
c$$$              do i=ilo,ihi,1
c$$$                 x = xl+(i - 0.5D0)*adx;
c$$$                 write(11,*) x,y,psi(i,j,k)
c$$$              enddo
c$$$              write(11,*) 
c$$$           enddo
c$$$        enddo
c
c        call flush(11)
        do k=klo,khi,1
           do j=jlo,jhi,1
              y = yl+(j - 0.5D0)*ady;
              do i=ilo,ihi,1
                 x = xl+(i - 0.5D0)*adx;
		 Bi=-ky*cos(kx*x)*cos(ky*y)
		 Bj=-kx*sin(kx*x)*sin(ky*y)
c                 Bi=-(psi(i,j+1,k)-psi(i,j-1,k))*dyi
c                 Bj=(psi(i+1,j,k)-psi(i-1,j,k))*dxi
		 Bk=-lamda*cos(kx*x)*sin(ky*y)
c		 if(cos(kx*x)*sin(ky*y).gt.0.D0) then
c                    Bi=0.1*Bi; Bj=0.1*Bj; Bk=0.1*Bk;
c                    Bi=0.D0; Bj=0.D0; Bk=0.D0
c                 endif
                 ux(i,j,k,1)=rho
                 ux(i,j,k,2)=rho*u
                 ux(i,j,k,3)=rho*v
                 ux(i,j,k,4)=rho*w
                 ux(i,j,k,5)=Bi
                 ux(i,j,k,6)=Bj
                 ux(i,j,k,7)=Bk
                 ux(i,j,k,8)=press/(gamma-1.D0)+
     &                0.5*rho*(u*u+v*v+w*w)+0.5D0*(bi*bi+bj*bj+bk*bk)
c$$$                 write(6,*) i,j,ux(i,j,k,1)
c$$$                 write(6,*) ux(i,j,k,2)
c$$$                 write(6,*) ux(i,j,k,3)
c$$$                 write(6,*) ux(i,j,k,4)
c$$$                 write(6,*) ux(i,j,k,5)
c$$$                 write(6,*) ux(i,j,k,6)
c$$$                 write(6,*) ux(i,j,k,7)
c$$$                 write(6,*) ux(i,j,k,8)
              enddo
           enddo
        enddo
c     
      return
      end








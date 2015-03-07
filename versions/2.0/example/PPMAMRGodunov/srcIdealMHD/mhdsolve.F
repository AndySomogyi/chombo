C -*- Mode: Fortran; -*-
c-----------------------------------------------------------------
c     Ravi Samtaney
c     Copyright 2002
c     Princeton Plasma Physics Laboratory
c     All Rights Reserved
c-----------------------------------------------------------------
c     $Log: mhdsolve.F,v $
c     Revision 1.2  2005/09/19 20:51:41  dmartin
c     changes to clean up output from code -- removed the code which was causing fort.21 and fort.22 files to be generated, also cleaned up output in Stone correction fortran and switched LevelSolver verbosity off if s_verbosity <= 3.
c
c     Revision 1.1  2005/05/12 19:22:45  ligocki
c     Moved current "{exec,src}MHD" to "{exec,src}IdealMHD"...
c
c     Revision 1.2  2004/10/06 20:52:57  samtaney
c     Intermediate checkin
c
c     Revision 1.1.1.1  2004/07/09 19:27:04  samtaney
c     Original source for Unsplit Ideal MHD based on new Chombo interface
c
c     Revision 1.1  2002/10/22 16:54:29  samtaney
c     Initial revision
c
c     Revision 1.2  2002/07/09 19:00:59  samtaney
c     Intermediate commit of files towards a single 2D/3D Ideal MHD code.
c
c     Revision 1.1  2002/06/13 22:16:18  samtaney
c     Original unsplit 8-wave Ideal MHD source (includes B filter).
c
c-----------------------------------------------------------------------
#define TWO_D
c#define TWO_HALF_D
c
c----------------------------------------------------------------------
       subroutine MHDSolveInit(ux,iulo,iuhi,julo,juhi,ilo,ihi,jlo,jhi,
     &                     klo,khi,
     &                     adx,ady,adz)
c----------------------------------------------------------------------
c
c     This routine is really a place holder
      implicit none
c 
#include "EosCommon.fh"
c
c
c     Incoming arguments
	integer iulo, iuhi,julo,juhi
	integer ilo,ihi,jlo,jhi,klo,khi
	real*8 ux(iulo:iuhi,
     &                        julo:juhi,
     &                        1:1,8)
        real*8 adx,ady,adz
c
c$$$	write(6,*) 'MHDSOLVE bounds=',ilo,ihi,jlo,jhi,klo,khi
c$$$	write(6,*) 'MHDSOLVE ubounds=', iulo, iuhi,julo,juhi
c
	return
	end
c     
c----------------------------------------------------------------------
       subroutine MHDSolve(ux,ilo,ihi,jlo,jhi,
     &                     klo,khi,
     &                     adx,ady,adz)
c----------------------------------------------------------------------
c     This routine is really a place holder
c
      implicit none
c
#include "EosCommon.fh"
c
c     Incoming arguments
	integer ilo,ihi,jlo,jhi,klo,khi
	real*8 ux(ilo:ihi,
     &                        jlo:jhi,
     &                        1:1,8)
        real*8 adx,ady,adz
c
c	write(6,*) 'MHDSOLVE bounds=',ilo,ihi,jlo,jhi,klo,khi
c
	return
	end
c
c----------------------------------------------------------------------
       subroutine SetMHDProperties(a_gamma)
c----------------------------------------------------------------------
c 
      implicit none
c
#include "EosCommon.fh"
c
c
c     Incoming arguments
	real*8 a_gamma
c
c       Set properties here - later worry about getting this
c       from the C++ part
	gamma=a_gamma
c	rgas=8317.D0
	rgas=1.D0
	wmol=1.D0
	mu=1.D-03
	eta=1.D-05
	kappa=1.D-02
c        etaFactor=1000.D0
        etaFactor=10000.D0
c
c        write(6,*) 'MHD PROPS', rgas,wmol,gamma,mu,eta,kappa
	return
	end
c     
c     
c-----------------------------------------------------------------------
      subroutine ConservedToPrimitive(u,iulo,iuhi,julo,juhi,kulo,kuhi,
     &     v,ivlo,ivhi,jvlo,jvhi,kvlo,kvhi,
     &     ilo,ihi,jlo,jhi,klo,khi)
c     
      implicit none
c     
#include "EosCommon.fh"
c     
c     Incoming arguments
      integer iulo, iuhi,julo,juhi, kulo,kuhi
      integer ivlo, ivhi,jvlo,jvhi,kvlo,kvhi
      integer ilo, ihi,jlo,jhi,klo,khi
c     
      real*8 u(iulo:iuhi,
     &     julo:juhi,
     &     kulo:kuhi,8)
      real*8 v(ivlo:ivhi,
     &     jvlo:jvhi,
     &     kvlo:kvhi,8)
c     
c
c     Local variables
      integer i,j,k,m
      real*8  rho_inv, vi,vj,vk,bi,bj,bk
c     
c$$$c     
c$$$      write(6,*) 'CONS2PRIM',ilo,ihi,jlo,jhi,klo,khi
c$$$      write(6,*) 'CONS2PRIM', iulo, iuhi,julo,juhi, kulo,kuhi
c$$$      write(6,*) 'CONS2PRIM', ivlo, ivhi,jvlo,jvhi,kvlo,kvhi
c
      do k=klo,khi,1
         do j=jlo,jhi,1
            do i=ilo,ihi,1
               rho_inv=1.0D0/u(i,j,k,1)
               v(i,j,k,1)=u(i,j,k,1)
               v(i,j,k,2)=u(i,j,k,2)*rho_inv
               v(i,j,k,3)=u(i,j,k,3)*rho_inv
               v(i,j,k,4)=u(i,j,k,4)*rho_inv
               v(i,j,k,5)=u(i,j,k,5)
               v(i,j,k,6)=u(i,j,k,6)
               v(i,j,k,7)=u(i,j,k,7)
               vi=v(i,j,k,2)
               vj=v(i,j,k,3)
               vk=v(i,j,k,4)
               bi=v(i,j,k,5)
               bj=v(i,j,k,6)
               bk=v(i,j,k,7)
               v(i,j,k,8)=(gamma-1.D0)*(u(i,j,k,8)-0.5D0*v(i,j,k,1)
     &              *(vi*vi+vj*vj+vk*vk)-0.5D0*(bi*bi+bj*bj+bk*bk))
               if(v(i,j,k,8).lt.0.01D0) then
C                   write(6,*) 'WARNING: Negative Pressure',i,j,k,
C      &                 v(i,j,k,8)
C                   write(6,*) u(i,j,k,1)
C                   write(6,*) u(i,j,k,2)
C                   write(6,*) u(i,j,k,3)
C                   write(6,*) u(i,j,k,4)
C                   write(6,*) u(i,j,k,5)
C                   write(6,*) u(i,j,k,6)
C                   write(6,*) u(i,j,k,7)
C                   write(6,*) u(i,j,k,8)
                  v(i,j,k,8)=0.01D0
                  u(i,j,k,8)=v(i,j,k,8)/(gamma-1.D0)+0.5D0*v(i,j,k,1)
     &              *(vi*vi+vj*vj+vk*vk)+0.5D0*(bi*bi+bj*bj+bk*bk)
               endif
            enddo
         enddo
      enddo
c
c     Corners
c$$$      if(ivlo.ne.ilo) then
c$$$      k=1
c$$$      i=ilo
c$$$      j=jlo
c$$$      do m=1, 8,1
c$$$         v(i-1,j-1,k,1)=v(i,j,k,m)
c$$$      enddo
c$$$      i=ihi
c$$$      j=jlo
c$$$      do m=1, 8,1
c$$$         v(i+1,j-1,k,1)=v(i,j,k,m)
c$$$      enddo
c$$$      i=ilo
c$$$      j=jhi
c$$$      do m=1, 8,1
c$$$         v(i-1,j+1,k,1)=v(i,j,k,m)
c$$$      enddo
c$$$      i=ihi
c$$$      j=jhi
c$$$      do m=1, 8,1
c$$$         v(i+1,j+1,k,1)=v(i,j,k,m)
c$$$      enddo
c$$$c
c$$$      endif
      return
      end
c-----------------------------------------------------------------------
      subroutine PrimitiveToConserved(u,iulo,iuhi,julo,juhi,kulo,kuhi,
     &     v,ivlo,ivhi,jvlo,jvhi,kvlo,kvhi,
     &     ilo,ihi,jlo,jhi,klo,khi)
c
      implicit none
c
#include "EosCommon.fh"
c
c     Incoming arguments
      integer iulo, iuhi,julo,juhi,kulo,kuhi
      integer ivlo, ivhi,jvlo,jvhi,kvlo,kvhi
      integer ilo, ihi,jlo,jhi,klo,khi
	real*8 u(iulo:iuhi,
     &     julo:juhi,
     &     kulo:kuhi,8)
	real*8 v(ivlo:ivhi,
     &       jvlo:jvhi,
     &       kvlo:kvhi,8)
c
c     Local variables
      integer i,j,k
      real*8  rho, vi,vj,vk,bi,bj,bk
c
c$$$	write(6,*) 'PRIM2CONS',ilo,ihi,jlo,jhi,klo,khi
c$$$	write(6,*) 'PRIM2CONS',8,gamma
c
	do k=klo,khi,1
	   do j=jlo,jhi,1
	      do i=ilo,ihi,1
               rho=v(i,j,k,1)
               u(i,j,k,1)=rho
               u(i,j,k,2)=v(i,j,k,2)*rho
               u(i,j,k,3)=v(i,j,k,3)*rho
               u(i,j,k,4)=v(i,j,k,4)*rho
               u(i,j,k,5)=v(i,j,k,5)
               u(i,j,k,6)=v(i,j,k,6)
               u(i,j,k,7)=v(i,j,k,7)
               vi=v(i,j,k,2)
               vj=v(i,j,k,3)
               vk=v(i,j,k,4)
               bi=v(i,j,k,5)
               bj=v(i,j,k,6)
               bk=v(i,j,k,7)
               u(i,j,k,8)=v(i,j,k,8)/(gamma-1.D0)+0.5D0*rho
     &              *(vi*vi+vj*vj+vk*vk)+0.5D0*(bi*bi+bj*bj+bk*bk)
            enddo
         enddo
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine MHDStabilizingSource(ux,ilo,ihi,jlo,jhi,
     &     klo,khi,
     &     adx,ady,adz,
     &     dt)
c     
      implicit none
c     
#include "EosCommon.fh"
c     
c     Incoming arguments
      integer ilo,ihi,jlo,jhi,klo,khi
      real*8 ux(ilo:ihi,
     &     jlo:jhi,
     &     1:1,8)
      real*8 adx,ady,adz
      real*8 dt
c     
c     Local variables
      real*8 divergB(ilo-1:ihi+1,jlo-1:jhi+1,1)
      real*8 vi,vj,vk,bi,bj,bk
      real*8 dxi,dyi,dzi
      integer i,j,k
c
      dxi=0.5D0/adx
      dyi=0.5D0/ady
      dzi=0.5D0/adz
c     
#ifdef TWO_D
      do k=klo,khi,1
         do j=jlo,jhi,1
            do i=ilo,ihi,1
               divergB(i,j,k)=dxi*(ux(i+1,j,k,5)-ux(i-1,j,k,5))+
     &              dyi*(ux(i,j+1,k,6)-ux(i,j-1,k,6))
            enddo
         enddo
      enddo
c
#else 
      do k=klo,khi,1
         do j=jlo,jhi,1
            do i=ilo,ihi,1
               divergB(i,j,k)=dxi*(ux(i+1,j,k,5)-ux(i-1,j,k,5))+
     &              dyi*(ux(i,j+1,k,6)-ux(i,j-1,k,6))+
     &              dzi*(ux(i,j,k+1,7)-ux(i,j,k-1,7))
            enddo
         enddo
      enddo
#endif
c     
c     
      do k=klo,khi,1
         do j=jlo,jhi,1
            do i=ilo,ihi,1
               bi=ux(i,j,k,5)
               bj=ux(i,j,k,6)
               bk=ux(i,j,k,7)
               vi=ux(i,j,k,2)/ux(i,j,k,1)
               vj=ux(i,j,k,3)/ux(i,j,k,1)
               vk=ux(i,j,k,4)/ux(i,j,k,1)
               ux(i,j,k,2)=ux(i,j,k,2)-dt*divergB(i,j,k)*bi
               ux(i,j,k,3)=ux(i,j,k,3)-dt*divergB(i,j,k)*bj
               ux(i,j,k,4)=ux(i,j,k,4)-dt*divergB(i,j,k)*bk
               ux(i,j,k,5)=ux(i,j,k,5)-dt*divergB(i,j,k)*vi
               ux(i,j,k,6)=ux(i,j,k,6)-dt*divergB(i,j,k)*vj
               ux(i,j,k,7)=ux(i,j,k,7)-dt*divergB(i,j,k)*vk
               ux(i,j,k,8)=ux(i,j,k,8)
     &              -dt*divergB(i,j,k)*(vi*bi+vj*bj+vk*bk)
            enddo
         enddo
      enddo
c     
      return
      end
c-----------------------------------------------------------------------
      subroutine PerpendicularCurrent(bpar,current,ilo,ihi,jlo,jhi,
     &     klo,khi,
     &     adx,ady,level)
c     
      implicit none
c     
#include "EosCommon.fh"
c     
c     Incoming arguments
      integer ilo,ihi,jlo,jhi,klo,khi
      real*8 bpar(ilo-1:ihi+1,
     &     jlo-1:jhi+1,
     &     1:1,2)
      real*8 current(ilo:ihi,jlo:jhi,1)
      real*8 adx, ady
      integer level
c     
c     Local variables
      real*8 dyi, dxi
      real*8 xl,yl,xr,yr,xc,yc
      integer i,j,k
c     
      xl=-1.0D0
      xr=1.D0
      yl=0.D0
      yr=1.D0
c     
c     
c     write(6,*) 'PERPCURRENT',ilo,ihi,jlo,jhi,
c     &                     klo,khi,
c     &                     adx,ady
      dxi=0.5D0/adx
      dyi=0.5D0/ady
c     
      do k=klo,khi,1
         do j=jlo+1,jhi-1,1
            do i=ilo+1,ihi-1,1
               current(i,j,k)=0.D0
            enddo
         enddo
      enddo
c     
      do k=klo,khi,1
         do j=jlo+1,jhi-1,1
            do i=ilo+1,ihi-1,1
c     
               yc = yl+(j - 0.5D0)*ady
               xc = xl+(i - 0.5D0)*adx
c     
               current(i,j,k)=dyi*(Bpar(i,j+1,k,1)-Bpar(i,j-1,k,1))-
     &              dxi*(Bpar(i+1,j,k,2)-Bpar(i-1,j,k,2))
               current(i,j,k)=current(i,j,k)/(level+1)
            enddo
         enddo
      enddo
c     
      return
      end
c
c-----------------------------------------------------------------------
c
       subroutine TotalKineticEnergy(tke,ux,ilo,ihi,jlo,jhi,
     &                     klo,khi,
     &                     adx,ady,adz,level)
c
      implicit none
c
#include "EosCommon.fh"
c
c     Incoming arguments
	integer ilo,ihi,jlo,jhi,klo,khi
	real*8 ux(ilo:ihi,
     &       jlo:jhi,
     &       1:1,8)
	real*8 tke
	real*8 adx, ady,adz
	integer level
c
c     Local variables
	integer i,j,k
c
	tke=0.D0
c
	return
	end
c
c-----------------------------------------------------------------------
c
       subroutine TotalMagneticEnergy(tme,ux,ilo,ihi,jlo,jhi,
     &                     klo,khi,
     &                     adx,ady,adz,level)

      implicit none
#include "EosCommon.fh"
c     Incoming arguments
	integer ilo,ihi,jlo,jhi,klo,khi
	real*8 ux(ilo:ihi,
     &       jlo:jhi,
     &       1:1,8)
	real*8 tme
	real*8 adx, ady,adz
	integer level
c
c     Local variables
	integer i,j,k
c
	tme=0.D0
c
c
	return
	end
c     
c-----------------------------------------------------------------------
c     
      subroutine PsiDot(dpsidt,ux,ilo,ihi,jlo,jhi,
     &     klo,khi,
     &     adx,ady,adz,level,ierr)
c     
      implicit none
c     
#include "EosCommon.fh"
c     
c     Incoming arguments
      integer ilo,ihi,jlo,jhi,klo,khi
      real*8 ux(ilo:ihi,
     &     jlo:jhi,
     &     1:1,8)
      real*8 dpsidt
      real*8 adx, ady,adz
      integer level,ierr
c     
c     Local variables
      integer i,j,k
      real*8 xl,yl,xr,yr
      real*8 xp,yp
      real*8 xlo,ylo,xhi,yhi
      real*8 etaprime,bk,u,v
      integer ip, jp

      ierr=0
      xl=-1.0D0
      xr=1.D0
      yl=0.D0
      yr=1.D0
      dpsidt=0.D0
      xp=0.D0
      yp=0.5D0
      
c     
      xlo=xl+(ilo-0.5D0)*adx
      ylo=yl+(jlo-0.5D0)*ady
      xhi=xl+(ihi-0.5D0)*adx
      yhi=yl+(jhi-0.5D0)*ady
c     write(6,*) 'PsiDot',level,xlo,ylo,xhi,yhi
c     write(6,*) 'PsiDot',level,ilo,jlo,ihi,jhi
      if((xp-xlo)*(xhi-xp).ge.0.D0 .and. 
     &     (yp-ylo)*(yhi-yp).ge.0.D0)then
         jp=(yp-yl)/ady+0.5D0
         ip=(xp-xl)/ady+0.5D0
c     write(6,*) 'Point found ',level,ip,jp,xlo,ylo,xhi,yhi
      else
         ierr=1
         return
      endif
c     
      do k=klo,khi,1
         do j=jp,jp+1,1
            do i=ip,ip+1,1
               u=ux(i,j,k,2)/ux(i,j,k,1)
               v=ux(i,j,k,3)/ux(i,j,k,1)
               bk=ux(i,j,k,7)
               if(bk.le.0.D0) then		 
                  etaPrime=eta+
     &                 etafactor*eta*(1.D0-exp(-(bk/0.075D0/7.8D0)**2))
               else
                  etaPrime=eta
               endif
               dpsidt=dpsidt
     &              -0.5D0*u*(ux(i+1,j,k,8)-ux(i-1,j,k,8))/adx
     &              -0.5D0*v*(ux(i,j+1,k,8)-ux(i,j-1,k,8))/ady
     &              +etaprime*((ux(i+1,j,k,8)+ux(i-1,j,k,8)-	 
     &              2.D0*ux(i,j,k,8))/adx/adx+
     &              (ux(i,j+1,k,8)+ux(i,j-1,k,8)-		 
     &              2.D0*ux(i,j,k,8))/ady/ady)
            enddo
         enddo
      enddo
c     
      dpsidt=0.25D0*dpsidt
c     
      return
      end
c     
c-----------------------------------------------------------------------
      subroutine FilterMagneticField(u,iulo,iuhi,julo,juhi,kulo,kuhi,
     1     f,iflo,ifhi,jflo,jfhi,kflo,kfhi,
     2     ilo,ihi,jlo,jhi,klo,khi,
     3     dt,dx,dy,dz)
c     
      implicit none
c     
#include "EosCommon.fh"
c     
c     Incoming arguments
      integer iulo, iuhi,julo,juhi,kulo,kuhi
      integer iflo, ifhi,jflo,jfhi,kflo,kfhi
      integer ilo, ihi,jlo,jhi,klo,khi
      real*8 u(iulo:iuhi,
     1     julo:juhi,
     2     kulo:kuhi,8)
      real*8 f(iflo:ifhi,
     1     jflo:jfhi,
     2     kflo:kfhi,8)
      real*8  dt,dx,dy,dz
c     
c     Local variables
      real*8 bnew(iulo:iuhi,
     1     julo:juhi,
     2     kulo:kuhi,3)
      integer i,j,k,l
c     
      real*8 lamda
      real*8 dxi,dyi,dzi, dxi2,dyi2,dzi2,dxidyi 
c     
      lamda=dx*dx/8.D0
c     
c$$$      write(6,*) 'Filter bounds',iulo,iuhi,julo,juhi
c$$$      write(6,*) 'Filter bounds',ilo,ihi,jlo,jhi
c$$$      write(6,*) 'Filter bounds',iflo,ifhi,jflo,jfhi, lamda
c$$$      write(6,*) 'Filter bounds',dx,dy,dz,dt
c     
      return
c$$$      do k=klo,khi,1
c$$$         do j=jlo,jhi,1
c$$$            do i=ilo,ihi,1
c$$$               do l=1,8,1
c$$$                  write(6,*) 'soln',i,j,u(i,j,k,l)
c$$$               enddo
c$$$            enddo
c$$$         enddo
c$$$      enddo
c$$$      return
c
      dxi=0.5D0/dx
      dyi=0.5D0/dy
      dzi=0.5D0/dz
      dxi2=1.D0/dx/dx
      dyi2=1.D0/dy/dy
      dzi2=1.D0/dz/dz
      dxidyi=0.25D0/dx/dy
c     
      do k=klo,khi,1
         do j=jlo,jhi,1
            do i=ilo,ihi,1
               bnew(i,j,k,1)=u(i,j,k,5)
               bnew(i,j,k,2)=u(i,j,k,6)
               bnew(i,j,k,3)=u(i,j,k,7)
            enddo
         enddo
      enddo
c     
      do k=klo,khi,1
         do j=jlo,jhi,1
            do i=ilo,ihi,1
               bnew(i,j,k,1)=u(i,j,k,5)+lamda*(   
     &              dxi2*(u(i+1,j,k,5)+u(i-1,j,k,5)-2.D0*u(i,j,k,5))
     &              +dxidyi*(u(i+1,j+1,k,6)+u(i-1,j-1,k,6)-
     &              u(i-1,j+1,k,6)-u(i+1,j-1,k,6)))
               bnew(i,j,k,2)=u(i,j,k,6)+lamda*(
     &              dyi2*(u(i,j+1,k,6)+u(i,j-1,k,6)-2.D0*u(i,j,k,6))
     &              +dxidyi*(u(i+1,j+1,k,5)+u(i-1,j-1,k,5)-
     &              u(i-1,j+1,k,5)-u(i+1,j-1,k,5)))
            enddo
         enddo
      enddo



      do k=klo,khi,1
         do j=jlo,jhi,1
            do i=ilo,ihi,1
               u(i,j,k,8)=u(i,j,k,8)+
     &              0.5D0*(bnew(i,j,k,1)*bnew(i,j,k,1)+bnew(i,j,k,2)*
     &              bnew(i,j,k,2)+bnew(i,j,k,3)*bnew(i,j,k,3)-
     &              (U(i,j,k,5)*U(i,j,k,5)+U(i,j,k,6)*
     &              U(i,j,k,6)+U(i,j,k,7)*U(i,j,k,7)))
               u(i,j,k,5)=bnew(i,j,k,1)
               u(i,j,k,6)=bnew(i,j,k,2)   
            enddo
         enddo
      enddo

      return
      end

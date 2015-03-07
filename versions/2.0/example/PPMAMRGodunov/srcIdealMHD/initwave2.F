C -*- Mode: Fortran; -*-
c-----------------------------------------------------------------
c       Ravi Samtaney
c       Copyright 2001
c       Princeton Plasma Physics Laboratory
c       All Rights Reserved
c-----------------------------------------------------------------
c       $Log: initwave2.F,v $
c       Revision 1.2  2005/09/19 20:51:41  dmartin
c       changes to clean up output from code -- removed the code which was causing fort.21 and fort.22 files to be generated, also cleaned up output in Stone correction fortran and switched LevelSolver verbosity off if s_verbosity <= 3.
c
c       Revision 1.1  2005/05/12 19:22:45  ligocki
c       Moved current "{exec,src}MHD" to "{exec,src}IdealMHD"...
c
c       Revision 1.3  2005/04/13 00:11:37  samtaney
c       Intermediate checkin - passes wave tests.
c
c       Revision 1.2  2004/10/06 20:52:56  samtaney
c       Intermediate checkin
c
c       Revision 1.1.1.1  2004/07/09 19:27:04  samtaney
c       Original source for Unsplit Ideal MHD based on new Chombo interface
c
c       Revision 1.2  2002/07/09 19:00:59  samtaney
c       Intermediate commit of files towards a single 2D/3D Ideal MHD code.
c
c       Revision 1.1  2002/07/09 15:00:15  samtaney
c       Original source of wave initial conditions file.
c
c-----------------------------------------------------------------
c Subroutine to set up the initial conditions
c---------------------------------------------------------------------------
c
      subroutine WaveInitialConditions(ux,
     &                     iulo,iuhi,julo,juhi,kulo,kuhi,
     &                     ilo,ihi,jlo,jhi,
     &                     klo,khi,
     &                     adx,ady,adz,alpha,
     &     pdir,kratio,waveNumber,pertAmplitude)
c----------------------------------------------------------------------
      implicit none
c
#include "EosCommon.fh"
c=======Declarations=========
	integer ilo,ihi,jlo,jhi,klo,khi
	integer iulo,iuhi,julo,juhi,kulo,kuhi
	real*8 ux(iulo:iuhi,
     &                        julo:juhi,
     &                        kulo:kuhi,8)
c     &                        klo-nghost:khi+nghost,nvar)
c
        integer  pdir,kratio,waveNumber
        real*8 pertAmplitude 
        real*8 perturbation
        real*8 adx,ady,adz
	integer  i,j,k,l,m
        real*8 dxi,dyi
        real*8 xl,yl,xr,yr,zl,zr
c
        real*8 x(3), kwave(3)
        real*8 v(8), v0(8),w0(8), w(8)
        real*8 evl(8,8), evr(8,8),lambda(8)
        real*8 theta, alpha, beta
        real*8 un,ut,bn,bt, uz,bz
        real*8 bx0,by0
        real*8 pi
c
c        write(6,*) 'WAVE INIT',ilo,ihi,jlo,jhi
c	write(6,*) 'WAVE INIT', iulo,iuhi,julo,juhi,kulo,kuhi
c        write(6,*) 'WAVE INIT',adx,ady,adz
c        write(6,*) 'WAVE INIT',gamma,rgas,wmol
c        write(6,*) 'WAVE INIT', pdir,kratio,waveNumber,
c     &       pertAmplitude
c
c     Eigenvector perturbation direction
c        pdir=1
c     Direction of propagation
c        theta=datan(1.D0)*4.D0/3.D0
c        theta=0.D0
c
        xl=-1.D0
        xr=1.D0
        yl=-1.D0
        yr=1.D0       
        zl=-1.D0
        zr=1.D0       
c
        v(1)=1.D0
        v(2)=0.D0
        v(3)=0.D0
        v(4)=0.D0
c        v(5)=0.0001D0
        v(5)=1.D0
        v(6)=0.D0
        v(7)=0.D0
        v(8)=0.1D0
c
        pi=4.D0*datan(1.D0)
        kwave(1)=pi*waveNumber
        kwave(2)=kratio*kwave(1)
        kwave(3)=0.D0
c     Wave vector is at an angle theta to x
        theta=datan(kwave(2)/kwave(1))
c
c     Kludge - wave moving in Y
c        kwave(2)=pi*waveNumber
c        kwave(1)=0.D0
c        theta=pi/2.D0
c
c     Kludge - wave moving in X
c        kwave(1)=pi*waveNumber
c        kwave(2)=0.D0
c        theta=0.D0
c
        x(3)=0.D0
c
c     B0 is aligned at an angle alpha to x
c        alpha=-pi/4.D0
c        alpha=pi/180.D0*(44.D0)
c        alpha=0.D0
        alpha=pi/180.D0*alpha
        bx0=v(5)*dcos(alpha)
        by0=v(5)*dsin(alpha)
c
        beta=alpha+theta
c
c        kwave(1)=kwave(1)*dcos(beta)
c        kwave(2)=kwave(1)*dsin(beta)
c
        v(2)=0.D0
        v(3)=0.D0
        v(4)=0.D0
        v(5)=bx0*dcos(theta)+by0*dsin(theta)
        v(6)=-bx0*dsin(theta)+by0*dcos(theta)
        v(7)=0.D0
c
c        write(6,*) 'Wave no=',kwave
c        write(6,*) 'Angles', theta*180.D0/pi, alpha*180.D0/pi,
c     &       beta*180.D0/pi
c
        call ProjectToCharacteristicSpace(v,w0,evl,evr,lambda)
c
c        write(6,*) 'Primitive Vars',v
c        write(6,*) 'Characteristic Vars',w0
c        write(6,*) 'Eigenvalues',lambda, dsqrt(gamma*v(8)/v(1))
c        do l=1,8,1
c           do k=1,8,1
c              write(6,*) 'Eigenvectors',l,k,evl(l,k),evr(l,k)
c           enddo
c        enddo
c
c     
c        write(6,*) 'ORTHO CHECK', matmul(evl,evr)
c     
c        v0=v
        do k=klo,khi,1
c           x(3) = zl+(k - 0.5D0)*adz;
           do j=jlo,jhi,1
              x(2) = yl+(j - 0.5D0)*ady;
              do i=ilo,ihi,1
                 x(1) = xl+(i - 0.5D0)*adx;
c                 perturbation=pertAmplitude*
c     &                dcos((kwave(1)*dcos(beta)*x(1)+
c     &                kwave(1)*dsin(beta)*x(2)+kwave(3)*x(3)))
                 perturbation=pertAmplitude*
     &                dcos((kwave(1)*x(1)+kwave(2)*x(2)+kwave(3)*x(3)))
c     &                dcos((kwave(1)*x(1)+kwave(2)*x(2)))
                 do l=1,8,1
                    w(l)=w0(l)
                    if(l.eq.pdir) then
                       w(l)=w0(l)+perturbation
                    endif
                 enddo
c
                 do l=1,8,1
                    v(l)=0.D0
                    do m=1,8,1
                       v(l)=v(l)+evr(l,m)*w(m)
                    enddo
                 enddo
c
c                 v(4)=0.D0
c                 v(7)=0.D0
                 un=v(2)
                 ut=v(3)
                 bn=v(5)
                 bt=v(6)
                 if(bn.eq.0.D0) then
                    write(6,*) 'ZERO',v
                 endif
                 v(2)=un*dcos(theta)-ut*dsin(theta)
                 v(3)=un*dsin(theta)+ut*dcos(theta)
                 v(5)=bn*dcos(theta)-bt*dsin(theta)
                 v(6)=bn*dsin(theta)+bt*dcos(theta)
c                 v(2)=un/dsqrt(3.D0)
c                 v(3)=un/dsqrt(3.D0)
c                 v(4)=un/dsqrt(3.D0)
c                 v(5)=bn/dsqrt(3.D0)
c                 v(6)=bn/dsqrt(3.D0)
c                 v(7)=bn/dsqrt(3.D0)
c     
                 ux(i,j,k,1)=v(1)
                 ux(i,j,k,2)=v(1)*v(2)
                 ux(i,j,k,3)=v(1)*v(3)
                 ux(i,j,k,4)=v(1)*v(4)
                 ux(i,j,k,5)=v(5)
                 ux(i,j,k,6)=v(6)
                 ux(i,j,k,7)=v(7)
                 ux(i,j,k,8)=v(8)/(gamma-1.D0)+
     &                0.5D0*v(1)*(v(2)*v(2)+v(3)*v(3)+v(4)*v(4))
     &                +0.5D0*(v(5)*v(5)+v(6)*v(6)+v(7)*v(7))
              enddo
           enddo
        enddo
                 
        i=ilo
        do k=klo,khi,1
           do j=jlo,jhi,1
              x(2) = (j - 0.5D0)*ady;
              x(1) = (i - 0.5D0)*adx;
c              write(21,*) x(2),
c     &             ux(i,j,k,1),ux(i,j,k,2),
c     &             ux(i,j,k,3),
c     &             ux(i,j,k,4),
c     &             ux(i,j,k,5),
c     &             ux(i,j,k,6),
c     &             ux(i,j,k,7),
c     &             ux(i,j,k,8)
           enddo
        enddo
c     
        j=jlo
        do k=klo,khi,1
           do i=ilo,ihi,1
              x(2) = (j - 0.5D0)*ady;
              x(1) = (i - 0.5D0)*adx;
c              write(22,*) x(1),
c     &             ux(i,j,k,1),ux(i,j,k,2),
c     &             ux(i,j,k,3),
c     &             ux(i,j,k,4),
c     &             ux(i,j,k,5),
c     &             ux(i,j,k,6),
c     &             ux(i,j,k,7),
c     &             ux(i,j,k,8)
           enddo
        enddo
        return
c     
      return
      end

c-----------------------------------------------------------------------
      subroutine ProjectToCharacteristicSpace(v,w,evl,evr,lambda)
#include "EosCommon.fh"
c     
      real*8 v(8)
      real*8 w(8)
      real*8 lambda(8)
      real*8 vi, vj, vk
      real*8 bi, bj, bk
      real*8 evl(8,8), evr(8,8), rhs(8)
      real*8 rhoInv
      real*8 Asq
      real*8 axsq
      real*8 csndsq
      real*8 cfast
      real*8 cslow
      real*8 tmp
      real*8 alphaS,alphaF
      real*8 betaY,betaZ
      real*8 betaPerpMagInv
      real*8 Sbx
      real*8 evec
      real*8 gamm
      real*8 dtodx
      real*8 pl(8), pr(8), up, um
      real*8 betaPerpMag
      integer l,k
c     
      real*8 TwoSqrtInv
      TwoSqrtInv=0.70710678D0
c     
         do k=1,8,1
            do l=1,8,1
               evr(k,l)=0.D0
               evl(k,l)=0.D0
            enddo
         enddo
         rhoInv=1.D0/v(1)
c     
         axsq=v(5)**2*rhoInv
         Asq=(v(5)**2+v(6)**2+v(7)**2)*rhoInv
         csndsq=gamma*v(8)*rhoInv
         tmp=dsqrt(dmax1(((csndsq+Asq)**2-4.D0*csndsq*axsq),0.D0))
         cfast=0.5D0*(csndsq+Asq+tmp)
         cslow=0.5D0*(csndsq+Asq-tmp)
c     
         alphaF=dsqrt(dmax1(csndsq-cslow,0.D0)/(cfast-cslow))
         alphaS=dsqrt(dmax1(cfast-csndsq,0.D0)/(cfast-cslow))
c     
         Sbx=dsign(1.D0,v(5))
c     
         betaPerpMag=dsqrt(v(6)**2+v(7)**2)
c         if(v(6).eq.0.D0.and.v(7).eq.0.D0) then
         if(betaPerpMag.eq.0.D0) then
            betaPerpMagInv=1.D0
            betaY=TwoSqrtInv
            betaZ=TwoSqrtInv
         else
c            betaPerpMagInv=1.D0/dsqrt(v(6)**2+v(7)**2)
            betaPerpMagInv=1.D0/betaPerpMag
            betaY=v(6)*betaPerpMagInv
            betaZ=v(7)*betaPerpMagInv
         endif
c     
c         if(betaPerpMagInv.gt.1.D08) then
c            betaPerpMagInv=1.D0
c            betaY=TwoSqrtInv
c            betaZ=TwoSqrtInv
c         endif
c
         cfast=dsqrt(cfast)
         cslow=dsqrt(dmax1(cslow,0.D0))
c
         axsq=dsqrt(dmax1(axsq,0.D0))
c
         lambda(1)=v(2)
         lambda(2)=v(2)
         lambda(3)=v(2)+axsq
         lambda(4)=v(2)-axsq
         lambda(5)=v(2)+cfast
         lambda(6)=v(2)-cfast
         lambda(7)=v(2)+cslow
         lambda(8)=v(2)-cslow
c
c         write(6,*) 'LAMDAs',lambda
c     Eigenvector corresponding to u
         evl(1,1)=1.D0
         evl(1,8)=-1.D0/csndsq
c     
c     Eigenvector corresponding to u
         evl(2,5)=1.D0
c     
c     Eigenvector corresponding to u+ax
         evl(3,3)=-betaZ*TwoSqrtInv
         evl(3,4)=betaY*TwoSqrtInv
         evl(3,6)=1.D0/dsqrt(v(1))*TwoSqrtInv*betaZ*SBx
         evl(3,7)=-1.D0/dsqrt(v(1))*TwoSqrtInv*betaY*SBx
c     
c     Eigenvector corresponding to u-ax
         evl(4,3)=-betaZ*TwoSqrtInv
         evl(4,4)=betaY*TwoSqrtInv
         evl(4,6)=-evl(3,6)
         evl(4,7)=-evl(3,7)
c     
c     Eigenvector corresponding to u+cfast
         evl(5,1)=0.D0
         evl(5,2)=0.5D0*alphaF*cfast/csndsq
         evl(5,3)=-0.5D0*alphaS*cslow*betaY*Sbx/csndsq
         evl(5,4)=-0.5D0*alphaS*cslow*betaZ*Sbx/csndsq
         evl(5,6)=0.5D0*alphaS/dsqrt(v(1)*csndsq)*betaY
         evl(5,7)=0.5D0*alphaS/dsqrt(v(1)*csndsq)*betaZ
         evl(5,8)=0.5D0*alphaF/v(1)/csndsq
c     
c     Eigenvector corresponding to u-cfast
         evl(6,1)=evl(5,1)
         evl(6,2)=-evl(5,2)
         evl(6,3)=-evl(5,3)
         evl(6,4)=-evl(5,4)
         evl(6,6)=evl(5,6)
         evl(6,7)=evl(5,7)
         evl(6,8)=evl(5,8)
c     
c     Eigenvector corresponding to u+cslow
         evl(7,1)=0.D0
         evl(7,2)=0.5D0*alphaS*cslow/csndsq
         evl(7,3)=0.5D0*alphaF*cfast*betaY*Sbx/csndsq
         evl(7,4)=0.5D0*alphaF*cfast*betaZ*Sbx/csndsq
         evl(7,6)=-0.5D0*alphaF/dsqrt(v(1)*csndsq)*betaY
         evl(7,7)=-0.5D0*alphaF/dsqrt(v(1)*csndsq)*betaZ
         evl(7,8)=0.5D0*alphaS/v(1)/csndsq
c     
c     
c     Eigenvector corresponding to u-cslow
         evl(8,1)=evl(7,1)
         evl(8,2)=-evl(7,2)
         evl(8,3)=-evl(7,3)
         evl(8,4)=-evl(7,4)
         evl(8,6)=evl(7,6)
         evl(8,7)=evl(7,7)
         evl(8,8)=evl(7,8)
c
c     Right eigenvectors
c     Eigenvector corresponding to u
         evr(1,1)=1.D0
c     
c     Eigenvector corresponding to u
         evr(5,2)=1.D0
c     
c     Eigenvector corresponding to u+ax
         evr(3,3)=-betaZ*TwoSqrtInv
         evr(4,3)=betaY*TwoSqrtInv
         evr(6,3)=dsqrt(v(1))*TwoSqrtInv*betaZ*SBx
         evr(7,3)=-dsqrt(v(1))*TwoSqrtInv*betaY*SBx
c     
c     Eigenvector corresponding to u-ax
         evr(3,4)=-betaZ*TwoSqrtInv
         evr(4,4)=betaY*TwoSqrtInv
         evr(6,4)=-evr(6,3)
         evr(7,4)=-evr(7,3)
c     
c     Eigenvector corresponding to u+cfast
         evr(1,5)=v(1)*alphaF
         evr(2,5)=alphaF*cfast
         evr(3,5)=-alphaS*cslow*betaY*Sbx
         evr(4,5)=-alphaS*cslow*betaZ*Sbx
         evr(6,5)=alphaS*dsqrt(v(1)*csndsq)*betaY
         evr(7,5)=alphaS*dsqrt(v(1)*csndsq)*betaZ
         evr(8,5)=v(1)*alphaF*csndsq
c     
c     Eigenvector corresponding to u-cfast
         evr(1,6)=evr(1,5)
         evr(2,6)=-evr(2,5)
         evr(3,6)=-evr(3,5)
         evr(4,6)=-evr(4,5)
         evr(6,6)=evr(6,5)
         evr(7,6)=evr(7,5)
         evr(8,6)=evr(8,5)
c     
c     Eigenvector corresponding to u+cslow
         evr(1,7)=v(1)*alphaS
         evr(2,7)=alphaS*cslow
         evr(3,7)=alphaF*cfast*betaY*Sbx
         evr(4,7)=alphaF*cfast*betaZ*Sbx
         evr(6,7)=-alphaF*dsqrt(v(1)*csndsq)*betaY
         evr(7,7)=-alphaF*dsqrt(v(1)*csndsq)*betaZ
         evr(8,7)=v(1)*alphaS*csndsq
c     
c     
c     Eigenvector corresponding to u-cslow
         evr(1,8)=evr(1,7)
         evr(2,8)=-evr(2,7)
         evr(3,8)=-evr(3,7)
         evr(4,8)=-evr(4,7)
         evr(6,8)=evr(6,7)
         evr(7,8)=evr(7,7)
         evr(8,8)=evr(8,7)
c     
         do l=1,8,1
            w(l)=0.D0
            do k=1,8,1
               w(l)=w(l)+evl(l,k)*v(k)
            enddo
         enddo
      return
      end subroutine ProjectToCharacteristicSpace
c-----------------------------------------------------------------------
      subroutine WaveSolutionError(ux,
     &                     iulo,iuhi,julo,juhi,kulo,kuhi,
     &                     ilo,ihi,jlo,jhi,
     &                     klo,khi,
     &                     adx,ady,adz,alpha,
     &     pdir,kratio,waveNumber,pertAmplitude,ttot,linferr,l2err)
c----------------------------------------------------------------------
      implicit none
c
#include "EosCommon.fh"
c=======Declarations=========
	integer ilo,ihi,jlo,jhi,klo,khi
	integer iulo,iuhi,julo,juhi,kulo,kuhi
	real*8 ux(iulo:iuhi,
     &                        julo:juhi,
     &                        kulo:kuhi,8)
	real*8 uxa(iulo:iuhi,
     &                        julo:juhi,
     &                        kulo:kuhi,8)
c     &                        klo-nghost:khi+nghost,nvar)
c
        integer  pdir,kratio,waveNumber
        real*8 pertAmplitude 
        real*8 perturbation
        real*8 adx,ady,adz
	integer  i,j,k,l,m
        real*8 dxi,dyi
        real*8 xl,yl,xr,yr,zl,zr
c
        real*8 x(3), kwave(3)
        real*8 v(8), v0(8),w0(8), w(8)
        real*8 evl(8,8), evr(8,8),lambda(8)
        real*8 theta, alpha, beta
        real*8 un,ut,bn,bt, uz,bz
        real*8 bx0,by0
        real*8 pi
        real*8 lambdax,lambday
        real*8 linferr, l2err
        real*8 ttot
        real*8 l2sum
c
        write(6,*) 'WAVE SOLUTION ERROR',ilo,ihi,jlo,jhi, ttot
c        return
c$$$        write(6,*) 'WAVE INIT',adx,ady,adz
c$$$        write(6,*) 'WAVE INIT',gamma,rgas,wmol
c$$$        write(6,*) 'WAVE INIT', pdir,kratio,waveNumber,
c$$$     &       pertAmplitude
c
c     Eigenvector perturbation direction
c        pdir=1
c     Direction of propagation
c        theta=datan(1.D0)*4.D0/3.D0
c        theta=0.D0
c
        xl=-1.D0
        xr=1.D0
        yl=-1.D0
        yr=1.D0       
        zl=-1.D0
        zr=1.D0       
c
        v(1)=1.D0
        v(2)=0.D0
        v(3)=0.D0
        v(4)=0.D0
c        v(5)=0.0001D0
        v(5)=1.0D0
        v(6)=0.D0
        v(7)=0.D0
        v(8)=0.1D0
c
        pi=4.D0*datan(1.D0)
        kwave(1)=pi*waveNumber
c        kwave(1)=4.D0*datan(1.D0)*waveNumber
c        kwave(2)=0.D0
c        kwave(2)=kwave(1)*dtan(theta)
        kwave(2)=kratio*kwave(1)
        kwave(3)=0.D0
c        kwave(2)=kwave(1)
c        kwave(3)=kwave(1)
c        theta=datan(2.D0)
c     Wave vector is at an angle theta to x
        theta=datan(kwave(2)/kwave(1))
        x(3)=0.D0
c
c     B0 is aligned at an angle alpha to x
c        alpha=-pi/4.D0
c        alpha=pi/180.D0*(44.D0)
c        alpha=0.D0
        bx0=v(5)*dcos(alpha)
        by0=v(5)*dsin(alpha)
c
        beta=alpha+theta
c
c        kwave(1)=kwave(1)*dcos(beta)
c        kwave(2)=kwave(1)*dsin(beta)
c
        v(2)=0.D0
        v(3)=0.D0
        v(4)=0.D0
        v(5)=bx0*dcos(theta)+by0*dsin(theta)
        v(6)=-bx0*dsin(theta)+by0*dcos(theta)
        v(7)=0.D0
c
c        write(6,*) 'Wave no=',kwave
c        write(6,*) 'Angles', theta*180.D0/pi, alpha*180.D0/pi,
c     &       beta*180.D0/pi
c
        call ProjectToCharacteristicSpace(v,w0,evl,evr,lambda)
c
c        write(6,*) 'Primitive Vars',v
c        write(6,*) 'Characteristic Vars',w0
c        write(6,*) 'Eigenvalues',lambda, dsqrt(gamma*v(8)/v(1))
c        do l=1,8,1
c           do k=1,8,1
c              write(6,*) 'Eigenvectors',l,k,evl(l,k),evr(l,k)
c           enddo
c        enddo
c
c     
c        write(6,*) 'ORTHO CHECK', matmul(evl,evr)
c     
c        v0=v
        l2sum=0.D0
        linferr=-1.D0
        lambdax=lambda(pdir)*cos(theta)
        lambday=lambda(pdir)*sin(theta)
        do k=klo,khi,1
c           x(3) = zl+(k - 0.5D0)*adz;
           do j=jlo,jhi,1
              x(2) = yl+(j - 0.5D0)*ady;
              do i=ilo,ihi,1
                 x(1) = xl+(i - 0.5D0)*adx;
c                 perturbation=pertAmplitude*
c     &                dcos((kwave(1)*dcos(beta)*x(1)+
c     &                kwave(1)*dsin(beta)*x(2)+kwave(3)*x(3)))
                 perturbation=pertAmplitude*
     &                dcos((kwave(1)*x(1)+kwave(2)*x(2)+kwave(3)*x(3))
     &                -kwave(1)*lambdax*ttot
     &                -kwave(2)*lambday*ttot)
c     &                -kwave(1)*lambda(pdir)*ttot
c     &                -kwave(2)*lambda(pdir)*ttot)
c     &                dcos((kwave(1)*x(1)+kwave(2)*x(2)))
                 do l=1,8,1
                    w(l)=w0(l)
                    if(l.eq.pdir) then
                       w(l)=w0(l)+perturbation
                    endif
                 enddo
c
                 do l=1,8,1
                    v(l)=0.D0
                    do m=1,8,1
                       v(l)=v(l)+evr(l,m)*w(m)
                    enddo
                 enddo
c
                 un=v(2)
                 ut=v(3)
                 bn=v(5)
                 bt=v(6)
                 if(bn.eq.0.D0) then
                    write(6,*) 'ZERO',v
                 endif
                 v(2)=un*dcos(theta)-ut*dsin(theta)
                 v(3)=un*dsin(theta)+ut*dcos(theta)
                 v(5)=bn*dcos(theta)-bt*dsin(theta)
                 v(6)=bn*dsin(theta)+bt*dcos(theta)
c     
                 uxa(i,j,k,1)=v(1)
                 uxa(i,j,k,2)=v(2)*v(1)
                 uxa(i,j,k,3)=v(3)*v(1)
                 uxa(i,j,k,4)=v(4)*v(1)
                 uxa(i,j,k,5)=v(5)
                 uxa(i,j,k,6)=v(6)
                 uxa(i,j,k,7)=v(7)
                 uxa(i,j,k,8)=v(8)/(gamma-1.D0)+
     &                0.5D0*v(1)*(v(2)*v(2)+v(3)*v(3)+v(4)*v(4))
     &                +0.5D0*(v(5)*v(5)+v(6)*v(6)+v(7)*v(7))
                 l2sum=l2sum+(ux(i,j,k,5)-v(5))**2
                 linferr=max(linferr,dabs(ux(i,j,k,5)-v(5)))
c$$$                 ux(i,j,k,1)=v(1)
c$$$                 ux(i,j,k,2)=v(1)*v(2)
c$$$                 ux(i,j,k,3)=v(1)*v(3)
c$$$                 ux(i,j,k,4)=v(1)*v(4)
c$$$                 ux(i,j,k,5)=v(5)
c$$$                 ux(i,j,k,6)=v(6)
c$$$                 ux(i,j,k,7)=v(7)
c$$$                 ux(i,j,k,8)=v(8)/(gamma-1.D0)+
c$$$     &                0.5D0*v(1)*(v(2)*v(2)+v(3)*v(3)+v(4)*v(4))
c$$$     &                +0.5D0*(v(5)*v(5)+v(6)*v(6)+v(7)*v(7))
              enddo
           enddo
        enddo
        l2err=dsqrt(l2sum/(jhi-jlo+1)/(ihi-ilo+1))
        write(6,*) 'LERROR',l2err, linferr
c                 
c$$$        do k=klo,khi,1
c$$$           do j=jlo,jhi,1
c$$$              i=j
c$$$              x(2) = (j - 0.5D0)*ady;
c$$$              x(1) = (i - 0.5D0)*adx;
c$$$              write(21,*) dsqrt(x(1)**2+x(2)**2),
c$$$     &             ux(i,j,k,1),ux(i,j,k,2),
c$$$     &             ux(i,j,k,3),
c$$$     &             ux(i,j,k,4),
c$$$     &             ux(i,j,k,5),
c$$$     &             ux(i,j,k,7),
c$$$     &             ux(i,j,k,8)
c$$$           enddo
c$$$        enddo
c     
        i=ilo
        do k=klo,khi,1
           do j=jlo,jhi,1
              x(2) = (j - 0.5D0)*ady;
              x(1) = (i - 0.5D0)*adx;
c              write(21,*) x(2),
c     &             ux(i,j,k,1),ux(i,j,k,2),
c     &             ux(i,j,k,3),
c     &             ux(i,j,k,4),
c     &             ux(i,j,k,5),
c     &             ux(i,j,k,6),
c     &             ux(i,j,k,7),
c     &             ux(i,j,k,8)
           enddo
        enddo
c     
        j=jlo
        do k=klo,khi,1
           do i=ilo,ihi,1
              x(2) = (j - 0.5D0)*ady;
              x(1) = (i - 0.5D0)*adx;
c              write(22,*) x(1),
c     &             ux(i,j,k,1),ux(i,j,k,2),
c     &             ux(i,j,k,3),
c     &             ux(i,j,k,4),
c     &             ux(i,j,k,5),
c     &             ux(i,j,k,6),
c     &             ux(i,j,k,7),
c     &             ux(i,j,k,8)
           enddo
        enddo
c     
        j=jlo
        do k=klo,khi,1
           do i=ilo,ihi,1
              x(2) = (j - 0.5D0)*ady;
              x(1) = (i - 0.5D0)*adx;
              do l=1,8,1
                 uxa(i,j,k,l)=
     &                dabs(ux(i,j,k,l)-uxa(i,j,k,l))/ttot
              enddo
c              write(23,*) x(1),
c     &             uxa(i,j,k,1),uxa(i,j,k,2),
c     &             uxa(i,j,k,3),
c     &             uxa(i,j,k,4),
c     &             uxa(i,j,k,5),
c     &             uxa(i,j,k,6),
c     &             uxa(i,j,k,7),
c     &             uxa(i,j,k,8)
           enddo
        enddo
c     
      return
      end

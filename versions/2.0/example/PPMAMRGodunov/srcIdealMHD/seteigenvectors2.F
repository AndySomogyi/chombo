C     -*- Mode: Fortran; -*-
c-----------------------------------------------------------------
c     Ravi Samtaney
c     Copyright 2002
c     Princeton Plasma Physics Laboratory
c     All Rights Reserved
c-----------------------------------------------------------------------
c     Use these eigenvectors when using the Riemann solver 
c     These are not conservative eigenvectors
c     The absolute value of the eigenvalues is not taken.
c     
      subroutine SetEigenSystem(utilde,evl,evr,alamda,ilo,ihi,direction)
c
      implicit none
c
#include "EosCommon.fh"
c
c     Incoming arguments
      integer ilo,ihi,direction
      real*8 utilde(ilo:ihi,8)
      real*8 evl(ilo:ihi,8,8)
      real*8 evr(ilo:ihi,8,8)
      real*8 alamda(ilo:ihi,8)
c
c     Local variables
      integer i,l,m
      real*8 uabs(ilo:ihi)
      real*8 rhoInv(ilo:ihi)
      real*8 Asq(ilo:ihi)
      real*8 axsq(ilo:ihi)
      real*8 csndsq(ilo:ihi)
      real*8 cfast(ilo:ihi)
      real*8 cslow(ilo:ihi)
      real*8 tmp(ilo:ihi)
      real*8 alphaS(ilo:ihi),alphaF(ilo:ihi)
      real*8 betaY(ilo:ihi),betaZ(ilo:ihi)
      real*8 betaPerpMagInv(ilo:ihi)
      real*8 Sbx(ilo:ihi)
      real*8 betaPerpMag
c     
      real*8 TwoSqrtInv,epsilon
      real*8 tolerance
c
      TwoSqrtInv=0.70710678D0
      epsilon=1.D-02
      tolerance=1.D-08
c
c     utilde(1)=rho, 
c     utilde(2)=u
c     utilde(3)=v
c     utilde(4)=w
c     utilde(5)=Bx
c     utilde(6)=By
c     utilde(7)=Bz
c     utilde(8)=p
c     
      do i=ilo,ihi,1
         cslow(i)=0.D0
         cfast(i)=0.D0
         tmp(i)=0.D0
         axsq(i)=0.D0
         Asq(i)=0.D0
         csndsq(i)=0.D0
         rhoInv(i)=0.D0
      enddo

c     
      do i=ilo,ihi,1
         uabs(i)=abs(utilde(i,2))
         rhoInv(i)=1.D0/utilde(i,1)
c     
         axsq(i)=utilde(i,5)**2*rhoInv(i)
         Asq(i)=(utilde(i,5)**2+utilde(i,6)**2+utilde(i,7)**2)*rhoInv(i)
         csndsq(i)=gamma*utilde(i,8)*rhoInv(i)
         tmp(i)=dsqrt(dmax1(((csndsq(i)+Asq(i))**2
     &        -4.D0*csndsq(i)*axsq(i)),0.D0))
         cfast(i)=0.5D0*(csndsq(i)+Asq(i)+tmp(i))
         cslow(i)=0.5D0*(csndsq(i)+Asq(i)-tmp(i))
c     
         alphaS(i)=dsqrt((dmax1(cfast(i)-csndsq(i),0.D0))
     &        /(cfast(i)-cslow(i)))
         alphaF(i)=dsqrt((dmax1(csndsq(i)-cslow(i),0.D0))
     &        /(cfast(i)-cslow(i)))
         
c     
         Sbx(i)=dsign(1.D0,utilde(i,5))
c     
      enddo
c     
      do i=ilo,ihi,1
         betaPerpMag=dsqrt(utilde(i,6)**2+utilde(i,7)**2)
c     if(utilde(i,6).eq.0.D0.and.utilde(i,7).eq.0.D0) then
c         if(betaPerpMag.eq.0.D0)then
         if(dabs(betaPerpMag).lt.tolerance) then
            betaPerpMagInv(i)=1.D0
            betaY(i)=TwoSqrtInv
            betaZ(i)=TwoSqrtInv
         else
c     betaPerpMagInv(i)=1.D0/dsqrt(utilde(i,6)**2+utilde(i,7)**2)
            betaPerpMagInv(i)=1.D0/betaPerpMag
            betaY(i)=utilde(i,6)*betaPerpMagInv(i)
            betaZ(i)=utilde(i,7)*betaPerpMagInv(i)
         endif
      enddo
c     
      do i=ilo,ihi,1
         do l=1,8,1
            do m=1,8,1
               evl(i,l,m)=0.D0
               evr(i,l,m)=0.D0
            enddo
         enddo
c     
         cfast(i)=dsqrt(cfast(i))
         cslow(i)=dsqrt(dmax1(cslow(i),0.D0))
         axsq(i)=dsqrt(axsq(i))
c     
         alamda(i,1)=utilde(i,2)
         alamda(i,2)=utilde(i,2)
         alamda(i,3)=(utilde(i,2)+axsq(i))
         alamda(i,4)=(utilde(i,2)-axsq(i))
         alamda(i,5)=(utilde(i,2)+cfast(i))
         alamda(i,6)=(utilde(i,2)-cfast(i))
         alamda(i,7)=(utilde(i,2)+cslow(i))
         alamda(i,8)=(utilde(i,2)-cslow(i))
c     
c     Eigenvector corresponding to u
         evr(i,1,1)=1.D0
c     
c     Eigenvector corresponding to u
         evr(i,5,2)=1.D0
c     
c     Eigenvector corresponding to u+ax
         evr(i,3,3)=-betaZ(i)*TwoSqrtInv
         evr(i,4,3)=betaY(i)*TwoSqrtInv
         evr(i,6,3)=dsqrt(utilde(i,1))*TwoSqrtInv*betaZ(i)
     &        *SBx(i)
         evr(i,7,3)=-dsqrt(utilde(i,1))*TwoSqrtInv*betaY(i)
     &        *SBx(i)
c     
c     Eigenvector corresponding to u-ax
         evr(i,3,4)=-betaZ(i)*TwoSqrtInv
         evr(i,4,4)=betaY(i)*TwoSqrtInv
         evr(i,6,4)=-evr(i,6,3)
         evr(i,7,4)=-evr(i,7,3)
c     
c     Eigenvector corresponding to u+cfast
         evr(i,1,5)=utilde(i,1)*alphaF(i)
         evr(i,2,5)=alphaF(i)*cfast(i)
         evr(i,3,5)=-alphaS(i)*cslow(i)*betaY(i)*Sbx(i)
         evr(i,4,5)=-alphaS(i)*cslow(i)*betaZ(i)*Sbx(i)
         evr(i,6,5)=alphaS(i)*dsqrt(utilde(i,1)*csndsq(i))
     &        *betaY(i)
         evr(i,7,5)=alphaS(i)*dsqrt(utilde(i,1)*csndsq(i))
     &        *betaZ(i)
         evr(i,8,5)=utilde(i,1)*alphaF(i)*csndsq(i)
c     
c     Eigenvector corresponding to u-cfast
         evr(i,1,6)=evr(i,1,5)
         evr(i,2,6)=-evr(i,2,5)
         evr(i,3,6)=-evr(i,3,5)
         evr(i,4,6)=-evr(i,4,5)
         evr(i,6,6)=evr(i,6,5)
         evr(i,7,6)=evr(i,7,5)
         evr(i,8,6)=evr(i,8,5)
c     
c     Eigenvector corresponding to u+cslow
         evr(i,1,7)=utilde(i,1)*alphaS(i)
         evr(i,2,7)=alphaS(i)*cslow(i)
         evr(i,3,7)=alphaF(i)*cfast(i)*betaY(i)*Sbx(i)
         evr(i,4,7)=alphaF(i)*cfast(i)*betaZ(i)*Sbx(i)
         evr(i,6,7)=-alphaF(i)*dsqrt(utilde(i,1)*csndsq(i))
     &        *betaY(i)
         evr(i,7,7)=-alphaF(i)*dsqrt(utilde(i,1)*csndsq(i))
     &        *betaZ(i)
         evr(i,8,7)=utilde(i,1)*alphaS(i)*csndsq(i)
c     
c     
c     Eigenvector corresponding to u-cslow
         evr(i,1,8)=evr(i,1,7)
         evr(i,2,8)=-evr(i,2,7)
         evr(i,3,8)=-evr(i,3,7)
         evr(i,4,8)=-evr(i,4,7)
         evr(i,6,8)=evr(i,6,7)
         evr(i,7,8)=evr(i,7,7)
         evr(i,8,8)=evr(i,8,7)

c     
c     Eigenvector corresponding to u
         evl(i,1,1)=1.D0
         evl(i,1,8)=-1.D0/csndsq(i)
c     
c     Eigenvector corresponding to u
         evl(i,2,5)=1.D0
c     
c     Eigenvector corresponding to u+ax
         evl(i,3,3)=-betaZ(i)*TwoSqrtInv
         evl(i,3,4)=betaY(i)*TwoSqrtInv
         evl(i,3,6)=1.D0/dsqrt(utilde(i,1))*TwoSqrtInv*betaZ(i)
     &        *SBx(i)
         evl(i,3,7)=-1.D0/dsqrt(utilde(i,1))*TwoSqrtInv*betaY(i)
     &        *SBx(i)
c     
c     Eigenvector corresponding to u-ax
         evl(i,4,3)=-betaZ(i)*TwoSqrtInv
         evl(i,4,4)=betaY(i)*TwoSqrtInv
         evl(i,4,6)=-evl(i,3,6)
         evl(i,4,7)=-evl(i,3,7)
c     
c     Eigenvector corresponding to u+cfast
         evl(i,5,1)=0.D0
         evl(i,5,2)=0.5D0*alphaF(i)*cfast(i)/csndsq(i)
         evl(i,5,3)=-0.5D0*alphaS(i)*cslow(i)*betaY(i)*Sbx(i)
     &        /csndsq(i)
         evl(i,5,4)=-0.5D0*alphaS(i)*cslow(i)*betaZ(i)*Sbx(i)
     &        /csndsq(i)
         evl(i,5,6)=0.5D0*alphaS(i)/dsqrt(utilde(i,1)*csndsq(i))
     &        *betaY(i)
         evl(i,5,7)=0.5D0*alphaS(i)/dsqrt(utilde(i,1)*csndsq(i))
     &        *betaZ(i)
         evl(i,5,8)=0.5D0*alphaF(i)/utilde(i,1)/csndsq(i)
c     
c     Eigenvector corresponding to u-cfast
         evl(i,6,1)=evl(i,5,1)
         evl(i,6,2)=-evl(i,5,2)
         evl(i,6,3)=-evl(i,5,3)
         evl(i,6,4)=-evl(i,5,4)
         evl(i,6,6)=evl(i,5,6)
         evl(i,6,7)=evl(i,5,7)
         evl(i,6,8)=evl(i,5,8)
c     
c     Eigenvector corresponding to u+cslow
         evl(i,7,1)=0.D0
         evl(i,7,2)=0.5D0*alphaS(i)*cslow(i)/csndsq(i)
         evl(i,7,3)=0.5D0*alphaF(i)*cfast(i)*betaY(i)*Sbx(i)
     &        /csndsq(i)
         evl(i,7,4)=0.5D0*alphaF(i)*cfast(i)*betaZ(i)*Sbx(i)
     &        /csndsq(i)
         evl(i,7,6)=-0.5D0*alphaF(i)/dsqrt(utilde(i,1)
     &        *csndsq(i))*betaY(i)
         evl(i,7,7)=-0.5D0*alphaF(i)/dsqrt(utilde(i,1)*csndsq(i))
     &        *betaZ(i)
         evl(i,7,8)=0.5D0*alphaS(i)/utilde(i,1)/csndsq(i)
c     
c     
c     Eigenvector corresponding to u-cslow
         evl(i,8,1)=evl(i,7,1)
         evl(i,8,2)=-evl(i,7,2)
         evl(i,8,3)=-evl(i,7,3)
         evl(i,8,4)=-evl(i,7,4)
         evl(i,8,6)=evl(i,7,6)
         evl(i,8,7)=evl(i,7,7)
         evl(i,8,8)=evl(i,7,8)
c     
      enddo
c     
      return
c     Harten's entropy fix
c     HARTEN'S correction for transonic case.
      do i=ilo,ihi,1
         do l=1,8,1
            if(dabs(alamda(i,l)).lt.epsilon) then
               alamda(i,l)=dsign(1.D0,alamda(i,l))*
     &              0.5D0*(alamda(i,l)**2
     &              +epsilon**2)/epsilon
            endif
         enddo
      enddo
c     
c     
      return
      end 








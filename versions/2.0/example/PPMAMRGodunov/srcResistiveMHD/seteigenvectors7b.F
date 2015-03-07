C     -*- Mode: Fortran; -*-
c-----------------------------------------------------------------
c     Ravi Samtaney
c     Copyright 2005
c     Princeton Plasma Physics Laboratory
c     All Rights Reserved
c-----------------------------------------------------------------------
c     Use these eigenvectors when using the Riemann solver 
c     These are not conservative eigenvectors
c     The absolute value of the eigenvalues is not taken.
c     These are for a 7 wave formulation.
c     Eigenvectors from Powell et al eliminating the one corresponding
c     to second u.
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
      real*8 evl(ilo:ihi,7,7)
      real*8 evr(ilo:ihi,7,7)
      real*8 alamda(ilo:ihi,7)
c
c     Local variables
      integer i,k,l,m
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
      real*8 evlp(7,7)
      real*8 evrp(7,7)
      real*8 dwdwt(7,7)
      real*8 dwtdw(7,7)
      real*8 press
c     
c     
      real*8 TwoSqrtInv,epsilon
      real*8 tolerance
c
      TwoSqrtInv=1.D0/dsqrt(2.D0)
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
         press=utilde(i,8)-0.5D0*(
     &        utilde(i,5)**2+utilde(i,6)**2+utilde(i,7)**2)
c         csndsq(i)=gamma*(utilde(i,8)*rhoInv(i)-0.5D0*Asq(i))
         csndsq(i)=gamma*press*rhoInv(i)
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
      do l=1,7,1
         do m=1,7,1
            evlp(l,m)=0.D0
            evrp(l,m)=0.D0
         enddo
      enddo
c
      do i=ilo,ihi,1
         do l=1,7,1
            do m=1,7,1
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
         alamda(i,2)=(utilde(i,2)+axsq(i))
         alamda(i,3)=(utilde(i,2)-axsq(i))
         alamda(i,4)=(utilde(i,2)+cfast(i))
         alamda(i,5)=(utilde(i,2)-cfast(i))
         alamda(i,6)=(utilde(i,2)+cslow(i))
         alamda(i,7)=(utilde(i,2)-cslow(i))
c     
c     Eigenvector corresponding to u
         evrp(1,1)=1.D0
c     
c     Eigenvector corresponding to u+ax
         evrp(3,2)=-betaZ(i)*TwoSqrtInv
         evrp(4,2)=betaY(i)*TwoSqrtInv
         evrp(5,2)=dsqrt(utilde(i,1))*TwoSqrtInv*betaZ(i)
     &        *SBx(i)
         evrp(6,2)=-dsqrt(utilde(i,1))*TwoSqrtInv*betaY(i)
     &        *SBx(i)
c     
c     Eigenvector corresponding to u-ax
         evrp(3,3)=-betaZ(i)*TwoSqrtInv
         evrp(4,3)=betaY(i)*TwoSqrtInv
         evrp(5,3)=-evrp(5,2)
         evrp(6,3)=-evrp(6,2)
c     
c     Eigenvector corresponding to u+cfast
         evrp(1,4)=utilde(i,1)*alphaF(i)
         evrp(2,4)=alphaF(i)*cfast(i)
         evrp(3,4)=-alphaS(i)*cslow(i)*betaY(i)*Sbx(i)
         evrp(4,4)=-alphaS(i)*cslow(i)*betaZ(i)*Sbx(i)
         evrp(5,4)=alphaS(i)*dsqrt(utilde(i,1)*csndsq(i))
     &        *betaY(i)
         evrp(6,4)=alphaS(i)*dsqrt(utilde(i,1)*csndsq(i))
     &        *betaZ(i)
         evrp(7,4)=utilde(i,1)*alphaF(i)*csndsq(i)
c     
c     Eigenvector corresponding to u-cfast
         evrp(1,5)=evrp(1,4)
         evrp(2,5)=-evrp(2,4)
         evrp(3,5)=-evrp(3,4)
         evrp(4,5)=-evrp(4,4)
         evrp(5,5)=evrp(5,4)
         evrp(6,5)=evrp(6,4)
         evrp(7,5)=evrp(7,4)
c     
c     Eigenvector corresponding to u+cslow
         evrp(1,6)=utilde(i,1)*alphaS(i)
         evrp(2,6)=alphaS(i)*cslow(i)
         evrp(3,6)=alphaF(i)*cfast(i)*betaY(i)*Sbx(i)
         evrp(4,6)=alphaF(i)*cfast(i)*betaZ(i)*Sbx(i)
         evrp(5,6)=-alphaF(i)*dsqrt(utilde(i,1)*csndsq(i))
     &        *betaY(i)
         evrp(6,6)=-alphaF(i)*dsqrt(utilde(i,1)*csndsq(i))
     &        *betaZ(i)
         evrp(7,6)=utilde(i,1)*alphaS(i)*csndsq(i)
c     
c     
c     Eigenvector corresponding to u-cslow
         evrp(1,7)=evrp(1,6)
         evrp(2,7)=-evrp(2,6)
         evrp(3,7)=-evrp(3,6)
         evrp(4,7)=-evrp(4,6)
         evrp(5,7)=evrp(5,6)
         evrp(6,7)=evrp(6,6)
         evrp(7,7)=evrp(7,6)

c     
c     Eigenvector corresponding to u
         evlp(1,1)=1.D0
         evlp(1,7)=-1.D0/csndsq(i)
c     
c     
c     Eigenvector corresponding to u+ax
         evlp(2,3)=-betaZ(i)*TwoSqrtInv
         evlp(2,4)=betaY(i)*TwoSqrtInv
         evlp(2,5)=1.D0/dsqrt(utilde(i,1))*TwoSqrtInv*betaZ(i)
     &        *SBx(i)
         evlp(2,6)=-1.D0/dsqrt(utilde(i,1))*TwoSqrtInv*betaY(i)
     &        *SBx(i)
c     
c     Eigenvector corresponding to u-ax
         evlp(3,3)=-betaZ(i)*TwoSqrtInv
         evlp(3,4)=betaY(i)*TwoSqrtInv
         evlp(3,5)=-evlp(2,5)
         evlp(3,6)=-evlp(2,6)
c     
c     Eigenvector corresponding to u+cfast
         evlp(4,1)=0.D0
         evlp(4,2)=0.5D0*alphaF(i)*cfast(i)/csndsq(i)
         evlp(4,3)=-0.5D0*alphaS(i)*cslow(i)*betaY(i)*Sbx(i)
     &        /csndsq(i)
         evlp(4,4)=-0.5D0*alphaS(i)*cslow(i)*betaZ(i)*Sbx(i)
     &        /csndsq(i)
         evlp(4,5)=0.5D0*alphaS(i)/dsqrt(utilde(i,1)*csndsq(i))
     &        *betaY(i)
         evlp(4,6)=0.5D0*alphaS(i)/dsqrt(utilde(i,1)*csndsq(i))
     &        *betaZ(i)
         evlp(4,7)=0.5D0*alphaF(i)/utilde(i,1)/csndsq(i)
c     
c     Eigenvector corresponding to u-cfast
         evlp(5,1)=evlp(4,1)
         evlp(5,2)=-evlp(4,2)
         evlp(5,3)=-evlp(4,3)
         evlp(5,4)=-evlp(4,4)
         evlp(5,5)=evlp(4,5)
         evlp(5,6)=evlp(4,6)
         evlp(5,7)=evlp(4,7)
c     
c     Eigenvector corresponding to u+cslow
         evlp(6,1)=0.D0
         evlp(6,2)=0.5D0*alphaS(i)*cslow(i)/csndsq(i)
         evlp(6,3)=0.5D0*alphaF(i)*cfast(i)*betaY(i)*Sbx(i)
     &        /csndsq(i)
         evlp(6,4)=0.5D0*alphaF(i)*cfast(i)*betaZ(i)*Sbx(i)
     &        /csndsq(i)
         evlp(6,5)=-0.5D0*alphaF(i)/dsqrt(utilde(i,1)
     &        *csndsq(i))*betaY(i)
         evlp(6,6)=-0.5D0*alphaF(i)/dsqrt(utilde(i,1)*csndsq(i))
     &        *betaZ(i)
         evlp(6,7)=0.5D0*alphaS(i)/utilde(i,1)/csndsq(i)
c     
c     
c     Eigenvector corresponding to u-cslow
         evlp(7,1)=evlp(6,1)
         evlp(7,2)=-evlp(6,2)
         evlp(7,3)=-evlp(6,3)
         evlp(7,4)=-evlp(6,4)
         evlp(7,5)=evlp(6,5)
         evlp(7,6)=evlp(6,6)
         evlp(7,7)=evlp(6,7)
c     
         do k=1,7,1
            do l=1,7,1
               dwtdw(k,l)=0.D0
               dwdwt(k,l)=0.D0
            enddo
         enddo
c         dwtdw=0.D0
c         dwdwt=0.D0
         do l=1,7,1
            dwtdw(l,l)=1.D0
            dwdwt(l,l)=1.D0
         enddo
         dwtdw(7,5)=utilde(i,6)
         dwtdw(7,6)=utilde(i,7)
         dwdwt(7,5)=-utilde(i,6)
         dwdwt(7,6)=-utilde(i,7)
         do k=1,7,1
            do l=1,7,1
               do m=1,7,1
                  evr(i,k,l)=evr(i,k,l)+dwtdw(k,m)*evrp(m,l)
                  evl(i,k,l)=evl(i,k,l)+evlp(k,m)*dwdwt(m,l)
               enddo
            enddo
         enddo   
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








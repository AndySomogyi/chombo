C     -*- Mode: Fortran; -*-
c-----------------------------------------------------------------
c     Ravi Samtaney
c     Copyright 2002
c     Princeton Plasma Physics Laboratory
c     All Rights Reserved
c-----------------------------------------------------------------
c     $Log: limitslopech.F,v $
c     Revision 1.1  2005/05/12 19:22:45  ligocki
c     Moved current "{exec,src}MHD" to "{exec,src}IdealMHD"...
c
c     Revision 1.2  2004/10/06 20:52:57  samtaney
c     Intermediate checkin
c
c     Revision 1.1.1.1  2004/07/09 19:27:04  samtaney
c     Original source for Unsplit Ideal MHD based on new Chombo interface
c
c     Revision 1.1  2002/06/13 22:16:18  samtaney
c     Original unsplit 8-wave Ideal MHD source (includes B filter).
c
c-----------------------------------------------------------------
      subroutine LimitSlopeCh(dvdx,dvl,dvr,v,ilo,ihi,direction)
c     
      implicit none
c     
#include "EosCommon.fh"
c     
c     Incoming arguments
      integer ilo,ihi,direction
      real*8 v(ilo:ihi,8)
      real*8 dvl(ilo:ihi,8)
      real*8 dvr(ilo:ihi,8)
      real*8 dvdx(ilo:ihi,8)
      real*8 dwl(8)
      real*8 dwr(8)
      real*8 dwlim(8)
      real*8 dwdx(8)
      real*8 dtodx
c     
c     Local variables
      integer i,k,l,m
      real*8 lambda(8)
      real*8 evl(8,8), evr(8,8)
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
      real*8 sum1
      real*8 pl(8), pr(8), up, um
c     
      real*8 TwoSqrtInv
      real*8 betaPerpMag
      real*8 tolerance
      real*8 epsilon
c     
c      TwoSqrtInv=0.70710678D0
      TwoSqrtInv=1.D0/dsqrt(2.D0)
      tolerance=1.D-08
      epsilon=1.D-02
c     
      do i=ilo,ihi,1
         do l=1,8,1
            dwl(l)=dvl(i,l)
            dwr(l)=dvr(i,l)
            dwdx(l)=dvdx(i,l)
            dwlim(l) = min(2.D0*dabs(dWl(l)),2.D0*dabs(dWr(l)))
            dWlim(l) = min(dWlim(l),dabs(dwdx(l)))

            if (dWl(l)*dWr(l) .lt. 0.D0) then
              dWlim(l) = 0.D0
            else
              dWlim(l) = dWlim(l)*sign(1.D0,dWl(l))
            endif
            
         enddo
         do l=1,8,1
            dWdx(l) = dWlim(l)
            dvdx(i,l)=dwdx(l)
         enddo
      enddo
      return
c     
      do i=ilo,ihi,1
c
         do l=1,8,1
            do m=1,8,1
               evl(l,m)=0.D0
               evr(l,m)=0.D0
            enddo
         enddo
c     
         rhoInv=1.D0/v(i,1)
         axsq=v(i,5)**2*rhoInv
         Asq=(v(i,5)**2+v(i,6)**2+v(i,7)**2)*rhoInv
         csndsq=gamma*v(i,8)*rhoInv
         tmp=dsqrt(dmax1(((csndsq+Asq)**2-4.D0*csndsq*axsq),0.D0))
         cfast=0.5D0*(csndsq+Asq+tmp)
         cslow=0.5D0*(csndsq+Asq-tmp)
c     
         alphaF=dsqrt(dmax1(csndsq-cslow,0.D0)/(cfast-cslow))
         alphaS=dsqrt(dmax1(cfast-csndsq,0.D0)/(cfast-cslow))
c     
         Sbx=dsign(1.D0,v(i,5))
c     
         betaPerpMag=dsqrt(v(i,6)**2+v(i,7)**2)
c     if(v(i,6).eq.0.D0.and.v(i,7).eq.0.D0) then
         if(dabs(betaPerpMag).lt.tolerance) then
c            write(6,*) 'BETAPERPMAG=',betaPerpMag
            betaPerpMagInv=1.D0
            betaY=TwoSqrtInv
            betaZ=TwoSqrtInv
         else
c     betaPerpMagInv=1.D0/dsqrt(v(i,6)**2+v(i,7)**2)
            betaPerpMagInv=1.D0/betaPerpMag
            betaY=v(i,6)*betaPerpMagInv
            betaZ=v(i,7)*betaPerpMagInv
         endif
c     
         cfast=dsqrt(cfast)
         cslow=dsqrt(dmax1(cslow,0.D0))
c     
         axsq=dsqrt(dmax1(axsq,0.D0))
c     
         lambda(1)=v(i,2)
         lambda(2)=v(i,2)
         lambda(3)=v(i,2)+axsq
         lambda(4)=v(i,2)-axsq
         lambda(5)=v(i,2)+cfast
         lambda(6)=v(i,2)-cfast
         lambda(7)=v(i,2)+cslow
         lambda(8)=v(i,2)-cslow
c         do l=1,8,1
c            if(dabs(lambda(l)).lt.0.01D0) then
c               lambda(l)=sign(1.D0,lambda(l))*
c     &              0.5D0*(lambda(l)**2
c     &              +epsilon**2)/epsilon
c            endif
c         enddo
c         write(6,*) 'LAMDAS= ',i,lambda
c     
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
         evl(3,6)=1.D0/dsqrt(v(i,1))*TwoSqrtInv*betaZ*SBx
         evl(3,7)=-1.D0/dsqrt(v(i,1))*TwoSqrtInv*betaY*SBx
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
         evl(5,6)=0.5D0*alphaS/dsqrt(v(i,1)*csndsq)*betaY
         evl(5,7)=0.5D0*alphaS/dsqrt(v(i,1)*csndsq)*betaZ
         evl(5,8)=0.5D0*alphaF/v(i,1)/csndsq
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
         evl(7,6)=-0.5D0*alphaF/dsqrt(v(i,1)*csndsq)*betaY
         evl(7,7)=-0.5D0*alphaF/dsqrt(v(i,1)*csndsq)*betaZ
         evl(7,8)=0.5D0*alphaS/v(i,1)/csndsq
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
         evr(6,3)=dsqrt(v(i,1))*TwoSqrtInv*betaZ*SBx
         evr(7,3)=-dsqrt(v(i,1))*TwoSqrtInv*betaY*SBx
c     
c     Eigenvector corresponding to u-ax
         evr(3,4)=-betaZ*TwoSqrtInv
         evr(4,4)=betaY*TwoSqrtInv
         evr(6,4)=-evr(6,3)
         evr(7,4)=-evr(7,3)
c     
c     Eigenvector corresponding to u+cfast
         evr(1,5)=v(i,1)*alphaF
         evr(2,5)=alphaF*cfast
         evr(3,5)=-alphaS*cslow*betaY*Sbx
         evr(4,5)=-alphaS*cslow*betaZ*Sbx
         evr(6,5)=alphaS*dsqrt(v(i,1)*csndsq)*betaY
         evr(7,5)=alphaS*dsqrt(v(i,1)*csndsq)*betaZ
         evr(8,5)=v(i,1)*alphaF*csndsq
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
         evr(1,7)=v(i,1)*alphaS
         evr(2,7)=alphaS*cslow
         evr(3,7)=alphaF*cfast*betaY*Sbx
         evr(4,7)=alphaF*cfast*betaZ*Sbx
         evr(6,7)=-alphaF*dsqrt(v(i,1)*csndsq)*betaY
         evr(7,7)=-alphaF*dsqrt(v(i,1)*csndsq)*betaZ
         evr(8,7)=v(i,1)*alphaS*csndsq
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
c     Project the slopes to the eigenspace
         do l=1,8,1
            dwl(l)=0.D0
            dwr(l)=0.D0
            dwdx(l)=0.D0
            do m=1,8,1
c     uncomment line below to scale slope by 1/2
c                  sum1=sum1+evl(k,m)*dvdx(i,m)*0.5D0
               dwl(l)=dwl(l)+evl(l,m)*dvl(i,m)
               dwr(l)=dwr(l)+evl(l,m)*dvr(i,m)
               dwdx(l)=dwdx(l)+evl(l,m)*dvdx(i,m)
            enddo
            dwlim(l) = min(2.D0*dabs(dWl(l)),2.D0*dabs(dWr(l)))
            dWlim(l) = min(dWlim(l),dabs(dwdx(l)))

            if (dWl(l)*dWr(l) .lt. 0.D0) then
              dWlim(l) = 0.D0
            else
              dWlim(l) = dWlim(l)*sign(1.D0,dWl(l))
            endif
            
         enddo
         do l=1,8,1
            dWdx(l) = 0.D0
            do m=1,8,1
               dWdx(l)=dWdx(l)+dWlim(m)*evr(l,m)
            enddo
            dvdx(i,l)=dwdx(l)
         enddo
c     
      enddo
c     
      return
      end 
c-----------------------------------------------------------------------


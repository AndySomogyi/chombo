C -*- Mode: Fortran; -*-
c-----------------------------------------------------------------------
c     Ravi Samtaney
c     Copyright 2002
c     Princeton Plasma Physics Laboratory
c     All Rights Reserved
c-----------------------------------------------------------------------
c     $Log: setalphas.F,v $
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
c-----------------------------------------------------------------------
      subroutine SetAlphas(utilde,alpha,durl,evl,evr,ilo,ihi,direction)
c
      implicit none
c
c     Incoming arguments
      integer ilo,ihi,direction
      real*8 utilde(ilo:ihi,8)
      real*8 alpha(ilo:ihi,8)
      real*8 durl(ilo:ihi,8)
      real*8 evl(ilo:ihi,8,8)
      real*8 evr(ilo:ihi,8,8)
c
c     Local variables
      integer i,k,l

c     
      do i=ilo,ihi,1
         do l=1,8,1
            alpha(i,l)=0.D0
            do k=1,8,1
               alpha(i,l)=alpha(i,l)+evl(i,l,k)*durl(i,k)
            enddo
         enddo
      enddo
c
c      do i=ilo,ihi,1
c            do k=1,8,1
c               write(6,*) 'CHECK ',sum(alpha(i,:)*evr(i,k,:)),durl(i,k)
c         enddo
c      enddo
      return
      end 


c-----------------------------------------------------------------------
      subroutine SetAlphas2(dw,idwlo,idwhi,jdwlo,jdwhi,
     &     kdwlo,kdwhi,
     &     w,iWlo,iWhi,jWlo,jWhi,
     &     kwlo,kwhi,
     &     ilo,ihi,
     &     jlo,jhi,klo,khi,idir)
c     
      implicit none
c     
#include "EosCommon.fh"
c     
c     Incoming arguments
      integer idwlo, idwhi,jdwlo,jdwhi, kdwlo,kdwhi
      integer iwlo, iwhi,jwlo,jwhi,kwlo,kwhi
      integer ilo, ihi,jlo,jhi,klo,khi
      integer idir, inorm,itan, itanb
c     
      real*8 dw(idwlo:idwhi,
     &     jdwlo:jdwhi,
     &     kdwlo:kdwhi,8)
c     
      real*8 w(iwlo:iwhi,
     &     jwlo:jwhi,
     &     kwlo:kwhi,8)
c     
c     Local variables
      integer i,j,k,m
      real*8  rho, rho_inv,vi,vj,vk,bi,bj,bk, p
      real*8  drho, dvi,dvj,dvk,dbi,dbj,dbk, dp
c     
      real*8 evl(8,8)
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
      real*8 pllambda(8), prlambda(8)
      real*8 dwlocal(8)
c     
      real*8 TwoSqrtInv
      real*8 betaPerpMag
      real*8 tolerance
      real*8 epsilon
c     
c     Sanity check
      if(idwlo.lt.ilo.or. idwlo.gt.ihi) then
         write(6,*) 'I -BOUNDS ERROR on DW',idwlo,ilo,idwhi,ihi
      endif
      if(jdwlo.lt.jlo.or. jdwlo.gt.jhi) then
         write(6,*) 'J -BOUNDS ERROR on DW',jdwlo,jlo,jdwhi,jhi
      endif
      if(kdwlo.lt.klo.or. kdwlo.gt.khi) then
         write(6,*) 'K -BOUNDS ERROR on DW',kdwlo,klo,kdwhi,khi
      endif
c     
      TwoSqrtInv=1.D0/dsqrt(2.D0)
      tolerance=1.D-08
      epsilon=1.D-02
c     
      inorm = 2 +     idir              
      itan = 2 + mod(idir + 1,3)
      itanb = 2 + mod(idir + 2,3)
      do k=klo,khi,1
         do j=jlo,jhi,1
            do i=ilo,ihi,1
               rho=w(i,j,k,1)
               rhoInv=1.D0/rho
               vi=w(i,j,k,inorm)
               vj=w(i,j,k,itan)
               vk=w(i,j,k,itanb)
               bi=w(i,j,k,inorm+3)
               bj=w(i,j,k,itan+3)
               bk=w(i,j,k,itanb+3)
               p=w(i,j,k,8)
c     
               drho=dw(i,j,k,1)
               dvi=dw(i,j,k,inorm)
               dvj=dw(i,j,k,itan)
               dvk=dw(i,j,k,itanb)
               dbi=dw(i,j,k,inorm+3)
               dbj=dw(i,j,k,itan+3)
               dbk=dw(i,j,k,itanb+3)
               dp=dw(i,j,k,8)
               dwlocal(1)=dw(i,j,k,1)
               dwlocal(2)=dw(i,j,k,inorm)
               dwlocal(3)=dw(i,j,k,itan)
               dwlocal(4)=dw(i,j,k,itanb)
               dwlocal(5)=dw(i,j,k,inorm+3)
               dwlocal(6)=dw(i,j,k,itan+3)
               dwlocal(7)=dw(i,j,k,itanb+3)
               dwlocal(8)=dw(i,j,k,8)
               axsq=bi**2*rhoInv
               Asq=(bi**2+bj**2+bk**2)*rhoInv
               csndsq=gamma*p*rhoInv
               tmp=dsqrt(dmax1(((csndsq+Asq)**2-4.D0*csndsq*axsq),0.D0))
               cfast=0.5D0*(csndsq+Asq+tmp)
               cslow=0.5D0*(csndsq+Asq-tmp)
c     
               alphaF=dsqrt(dmax1(csndsq-cslow,0.D0)/(cfast-cslow))
               alphaS=dsqrt(dmax1(cfast-csndsq,0.D0)/(cfast-cslow))
c     
               Sbx=dsign(1.D0,bi)
c     
               betaPerpMag=dsqrt(bj**2+bk**2)
               if(dabs(betaPerpMag).lt.tolerance) then
                  betaPerpMagInv=1.D0
                  betaY=TwoSqrtInv
                  betaZ=TwoSqrtInv
               else
                  betaPerpMagInv=1.D0/betaPerpMag
                  betaY=bj*betaPerpMagInv
                  betaZ=bk*betaPerpMagInv
               endif
c     
               cfast=dsqrt(cfast)
               cslow=dsqrt(dmax1(cslow,0.D0))
c     
               axsq=dsqrt(dmax1(axsq,0.D0))
c     
c$$$  lambda(1)=w(i,j,k,2)
c$$$  lambda(2)=w(i,j,k,2)
c$$$  lambda(3)=w(i,j,k,2)+axsq
c$$$  lambda(4)=w(i,j,k,2)-axsq
c$$$  lambda(5)=w(i,j,k,2)+cfast
c$$$  lambda(6)=w(i,j,k,2)-cfast
c$$$  lambda(7)=w(i,j,k,2)+cslow
c$$$  lambda(8)=w(i,j,k,2)-cslow

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
               evl(3,6)=1.D0/dsqrt(w(i,j,k,1))*TwoSqrtInv*betaZ*SBx
               evl(3,7)=-1.D0/dsqrt(w(i,j,k,1))*TwoSqrtInv*betaY*SBx
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
               evl(5,6)=0.5D0*alphaS/dsqrt(w(i,j,k,1)*csndsq)*betaY
               evl(5,7)=0.5D0*alphaS/dsqrt(w(i,j,k,1)*csndsq)*betaZ
               evl(5,8)=0.5D0*alphaF/w(i,j,k,1)/csndsq
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
               evl(7,6)=-0.5D0*alphaF/dsqrt(w(i,j,k,1)*csndsq)*betaY
               evl(7,7)=-0.5D0*alphaF/dsqrt(w(i,j,k,1)*csndsq)*betaZ
               evl(7,8)=0.5D0*alphaS/w(i,j,k,1)/csndsq
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
               dw(i,j,k,1)=evl(1,k)*dwlocal(k)
            enddo
         enddo
      enddo
      return
      end

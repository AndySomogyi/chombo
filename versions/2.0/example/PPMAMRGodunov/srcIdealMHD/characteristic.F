C -*- Mode: Fortran; -*-
c-----------------------------------------------------------------------
c     Ravi Samtaney
c     Copyright 2004
c     Princeton Plasma Physics Laboratory
c     All Rights Reserved
c-----------------------------------------------------------------------
c     $Log: characteristic.F,v $
c     Revision 1.1  2005/05/12 19:22:45  ligocki
c     Moved current "{exec,src}MHD" to "{exec,src}IdealMHD"...
c
c     Revision 1.2  2004/10/06 20:52:56  samtaney
c     Intermediate checkin
c
c     Revision 1.1.1.1  2004/07/09 19:27:04  samtaney
c     Original source for Unsplit Ideal MHD based on new Chombo interface
c
c-----------------------------------------------------------------------
      subroutine MHDAlphas(dw,idwlo,idwhi,jdwlo,jdwhi,
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
      integer i,j,k,l,m
      real*8  rho, rho_inv,vi,vj,vk,bi,bj,bk, p
      real*8  drho, dvi,dvj,dvk,dbi,dbj,dbk, dp
      real*8  dwlocal(8)
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
      write(6,*) 'SETALPHAS',idir,inorm,itan,itanb
      do l=1,8,1
         do m=1,8,1
            evl(l,m)=0.D0
         enddo
      enddo
c      evl=0.D0
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
               dwlocal(1)=dw(i,j,k,1)
               dwlocal(2)=dw(i,j,k,inorm)
               dwlocal(3)=dw(i,j,k,itan)
               dwlocal(4)=dw(i,j,k,itanb)
               dwlocal(5)=dw(i,j,k,inorm+3)
               dwlocal(6)=dw(i,j,k,itan+3)
               dwlocal(7)=dw(i,j,k,itanb+3)
               dwlocal(8)=dw(i,j,k,8)
c
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
               do l=1,8,1
                  dw(i,j,k,l)=0.D0
               enddo
               do m=1,8,1
                  dw(i,j,k,1)=dw(i,j,k,1)+evl(1,m)*dwlocal(m)
                  dw(i,j,k,inorm)=dw(i,j,k,inorm)+evl(2,m)*dwlocal(m)
                  dw(i,j,k,itan)=dw(i,j,k,itan)+evl(3,m)*dwlocal(m)
                  dw(i,j,k,itanb)=dw(i,j,k,itanb)+evl(4,m)*dwlocal(m)
                  dw(i,j,k,inorm+3)=dw(i,j,k,inorm+3)
     &                 +evl(5,m)*dwlocal(m)
                  dw(i,j,k,itan+3)=dw(i,j,k,itan+3)
     &                 +evl(6,m)*dwlocal(m)
                  dw(i,j,k,itanb+3)=dw(i,j,k,itanb+3)
     &                 +evl(7,m)*dwlocal(m)
                  dw(i,j,k,8)=dw(i,j,k,8)
     &                 +evl(8,m)*dwlocal(m)
c     dw(i,j,k,l)=dw(i,j,k,l)+evl(l,m)*dwlocal(m)
               enddo
c     enddo
            enddo
         enddo
      enddo
      return
      end
c     
c-----------------------------------------------------------------------
      subroutine MHDCharacteristicSynthesis(dw,idwlo,idwhi,jdwlo,jdwhi,
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
c     
c     Local variables
      integer i,j,k,l,m
      real*8  rho, rho_inv,vi,vj,vk,bi,bj,bk, p
      real*8  drho, dvi,dvj,dvk,dbi,dbj,dbk, dp
c     
      real*8 dwlocal(8)
      real*8 evr(8,8)
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
      do l=1,8,1
         do m=1,8,1
            evr(l,m)=0.D0
         enddo
      enddo
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
               dwlocal(1)=dw(i,j,k,1)
               dwlocal(2)=dw(i,j,k,inorm)
               dwlocal(3)=dw(i,j,k,itan)
               dwlocal(4)=dw(i,j,k,itanb)
               dwlocal(5)=dw(i,j,k,inorm+3)
               dwlocal(6)=dw(i,j,k,itan+3)
               dwlocal(7)=dw(i,j,k,itanb+3)
               dwlocal(8)=dw(i,j,k,8)
c
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
               evr(6,3)=dsqrt(rho)*TwoSqrtInv*betaZ*SBx
               evr(7,3)=-dsqrt(rho)*TwoSqrtInv*betaY*SBx
c     
c     Eigenvector corresponding to u-ax
               evr(3,4)=-betaZ*TwoSqrtInv
               evr(4,4)=betaY*TwoSqrtInv
               evr(6,4)=-evr(6,3)
               evr(7,4)=-evr(7,3)
c     
c     Eigenvector corresponding to u+cfast
               evr(1,5)=rho*alphaF
               evr(2,5)=alphaF*cfast
               evr(3,5)=-alphaS*cslow*betaY*Sbx
               evr(4,5)=-alphaS*cslow*betaZ*Sbx
               evr(6,5)=alphaS*dsqrt(rho*csndsq)*betaY
               evr(7,5)=alphaS*dsqrt(rho*csndsq)*betaZ
               evr(8,5)=rho*alphaF*csndsq
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
               evr(1,7)=rho*alphaS
               evr(2,7)=alphaS*cslow
               evr(3,7)=alphaF*cfast*betaY*Sbx
               evr(4,7)=alphaF*cfast*betaZ*Sbx
               evr(6,7)=-alphaF*dsqrt(rho*csndsq)*betaY
               evr(7,7)=-alphaF*dsqrt(rho*csndsq)*betaZ
               evr(8,7)=rho*alphaS*csndsq
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
                  dw(i,j,k,l)=0.D0
               enddo
               do m=1,8,1
c                  dw(i,j,k,l)=dw(i,j,k,l)+evr(l,m)*dwlocal(m)
                  dw(i,j,k,1)=dw(i,j,k,1)+evr(1,m)*dwlocal(m)
                  dw(i,j,k,inorm)=dw(i,j,k,inorm)+evr(2,m)*dwlocal(m)
                  dw(i,j,k,itan)=dw(i,j,k,itan)+evr(3,m)*dwlocal(m)
                  dw(i,j,k,itanb)=dw(i,j,k,itanb)+evr(4,m)*dwlocal(m)
                  dw(i,j,k,inorm+3)=dw(i,j,k,inorm+3)
     &                 +evr(5,m)*dwlocal(m)
                  dw(i,j,k,itan+3)=dw(i,j,k,itan+3)
     &                 +evr(6,m)*dwlocal(m)
                  dw(i,j,k,itanb+3)=dw(i,j,k,itanb+3)
     &                 +evr(7,m)*dwlocal(m)
                  dw(i,j,k,8)=dw(i,j,k,8)+evr(8,m)*dwlocal(m)
               enddo
c            enddo
            enddo
         enddo
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine MHDEigenvalues(lambda,
     &     ilambdalo,ilambdahi,jlambdalo,jlambdahi,
     &     klambdalo,klambdahi,
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
      integer ilambdalo, ilambdahi,jlambdalo,jlambdahi, 
     &klambdalo,klambdahi
      integer iwlo, iwhi,jwlo,jwhi,kwlo,kwhi
      integer ilo, ihi,jlo,jhi,klo,khi
      integer idir
c     
      real*8 lambda(ilambdalo:ilambdahi,
     &     jlambdalo:jlambdahi,
     &     klambdalo:klambdahi,8)
c     
      real*8 w(iwlo:iwhi,
     &     jwlo:jwhi,
     &     kwlo:kwhi,8)
c     
c     
c     Local variables
      integer i,j,k,l,m
      integer inorm,itan, itanb
      real*8  rho, rho_inv,vi,vj,vk,bi,bj,bk, p
      real*8  drho, dvi,dvj,dvk,dbi,dbj,dbk, dp
c     
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
c     
      real*8 TwoSqrtInv
      real*8 betaPerpMag
      real*8 tolerance
      real*8 epsilon
c     
c     Sanity check
      if(ilambdalo.lt.ilo.or. ilambdalo.gt.ihi) then
         write(6,*) 'I -BOUNDS ERROR on LAMBDA',ilambdalo,ilo,
     &        ilambdahi,ihi
      endif
      if(jlambdalo.lt.jlo.or. jlambdalo.gt.jhi) then
         write(6,*) 'J -BOUNDS ERROR on LAMBDA',jlambdalo,jlo,
     &        jlambdahi,jhi
      endif
      if(klambdalo.lt.klo.or. klambdalo.gt.khi) then
         write(6,*) 'K -BOUNDS ERROR on LAMBDA',klambdalo,klo,
     &        klambdahi,khi
      endif
c     
      TwoSqrtInv=1.D0/dsqrt(2.D0)
      tolerance=1.D-08
      epsilon=1.D-02
c     
      inorm = 2 +     idir              
      itan = 2 + mod(idir + 1,3)
      itanb = 2 + mod(idir + 2,3)
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
c$$$               lambda(i,j,k,1)=w(i,j,k,2)
c$$$               lambda(i,j,k,2)=w(i,j,k,2)
c$$$               lambda(i,j,k,3)=w(i,j,k,2)+axsq
c$$$               lambda(i,j,k,4)=w(i,j,k,2)-axsq
c$$$               lambda(i,j,k,5)=w(i,j,k,2)+cfast
c$$$               lambda(i,j,k,6)=w(i,j,k,2)-cfast
c$$$               lambda(i,j,k,7)=w(i,j,k,2)+cslow
c$$$               lambda(i,j,k,8)=w(i,j,k,2)-cslow
               lambda(i,j,k,1)=vi
               lambda(i,j,k,2)=vi
               lambda(i,j,k,3)=vi+axsq
               lambda(i,j,k,4)=vi-axsq
               lambda(i,j,k,5)=vi+cfast
               lambda(i,j,k,6)=vi-cfast
               lambda(i,j,k,7)=vi+cslow
               lambda(i,j,k,8)=vi-cslow
            enddo
         enddo
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine MHDQuasiLinearUpdate(adwdx,
     &     iadwdxlo,iadwdxhi,jadwdxlo,jadwdxhi,
     &     kadwdxlo,kadwdxhi,
     &     whalf,iWhalflo,iWhalfhi,jWhalflo,jWhalfhi,
     &     kwhalflo,kwhalfhi,
     &     w,iWlo,iWhi,jWlo,jWhi,
     &     kwlo,kwhi,
     &     ilo,ihi,
     &     jlo,jhi,klo,khi,idir,dtodx)
c
      integer iadwdxlo, iadwdxhi,jadwdxlo,jadwdxhi, 
     &kadwdxlo,kadwdxhi
      integer iwhalflo, iwhalfhi,jwhalflo,jwhalfhi, 
     &kwhalflo,kwhalfhi
      
      integer iwlo, iwhi,jwlo,jwhi,kwlo,kwhi
      integer ilo, ihi,jlo,jhi,klo,khi
      integer idir

      real*8 dtodx
      real*8 adwdx(iadwdxlo:iadwdxhi,
     &     jadwdxlo:jadwdxhi,
     &     kadwdxlo:kadwdxhi,8)
      real*8 whalf(iwhalflo:iwhalfhi,
     &     jwhalflo:jwhalfhi,
     &     kwhalflo:kwhalfhi,8)
      real*8 w(iwlo:iwhi,
     &     jwlo:jwhi,
     &     kwlo:kwhi,8)
c
#include "EosCommon.fh"
c     Local variables
      integer i,j,k,l,m
      real*8  rho, rho_inv,vi,vj,vk,bi,bj,bk, p
      real*8  drho, dvi,dvj,dvk,dbi,dbj,dbk, dp, divb
      integer inorm,itan,itanb
      integer ioff, joff, koff
      real*8 uc(8), fl(8), fr(8), sdivb(8)
c

c
      inorm = 2 +     idir              
      itan = 2 + mod(idir + 1,3)
      itanb = 2 + mod(idir + 2,3)
c
      ioff = 0
      joff = 0
      koff = 0
      if(idir.eq.0) ioff=1
      if(idir.eq.1) joff=1
      if(idir.eq.2) koff=1
c
      do k=klo,khi,1
         do j=jlo,jhi,1
            do i=ilo,ihi,1
               rho=whalf(i,j,k,1)
               vi=whalf(i,j,k,inorm)
               vj=whalf(i,j,k,itan)
               vk=whalf(i,j,k,itanb)
               bi=whalf(i,j,k,inorm+3)
               bj=whalf(i,j,k,itan+3)
               bk=whalf(i,j,k,itanb+3)
               p=whalf(i,j,k,8)
               divb=whalf(i,j,k,inorm+3)
               sdivb(2)=bi
               sdivb(3)=bj
               sdivb(4)=bk
               sdivb(5)=vi
               sdivb(6)=vj
               sdivb(7)=vk
               sdivb(8)=(vi*bi+vj*bj+vk*bk)
               fl(1)=rho*vi
               fl(2)=rho*vi*vi +p+0.5D0*(bj*bj+bk*bk-bi*bi)
               fl(3)=rho*vi*vj-bi*bj
               fl(4)=rho*vi*vk-bi*bk
               fl(5)=0.D0
               fl(6)=vi*bj-vj*bi
               fl(7)=vi*bk-vk*bi
               fl(8)=(0.5D0*rho*(vi*vi+vj*vj+vk*vk)+p/(gamma-1.D0)+
     &           p+(bi*bi+bj*bj+bk*bk))*vi-
     &           bi*(vi*bi+vj*bj+vk*bk)
               rho=whalf(i+ioff,j+joff,k+koff,1)
               vi=whalf(i+ioff,j+joff,k+koff,inorm)
               vj=whalf(i+ioff,j+joff,k+koff,itan)
               vk=whalf(i+ioff,j+joff,k+koff,itanb)
               bi=whalf(i+ioff,j+joff,k+koff,inorm+3)
               bj=whalf(i+ioff,j+joff,k+koff,itan+3)
               bk=whalf(i+ioff,j+joff,k+koff,itanb+3)
               p=whalf(i+ioff,j+joff,k+koff,8)
               fr(1)=rho*vi
               fr(2)=rho*vi*vi +p+0.5D0*(bj*bj+bk*bk-bi*bi)
               fr(3)=rho*vi*vj-bi*bj
               fr(4)=rho*vi*vk-bi*bk
               fr(5)=0.D0
               fr(6)=vi*bj-vj*bi
               fr(7)=vi*bk-vk*bi
               fr(8)=(0.5D0*rho*(vi*vi+vj*vj+vk*vk)+p/(gamma-1.D0)+
     &           p+(bi*bi+bj*bj+bk*bk))*vi-
     &           bi*(vi*bi+vj*bj+vk*bk)
               divb=bi-divb
               sdivb(2)=0.5D0*(sdivb(2)+bi)
               sdivb(3)=0.5D0*(sdivb(3)+bj)
               sdivb(4)=0.5D0*(sdivb(4)+bk)
               sdivb(5)=0.5D0*(sdivb(5)+vi)
               sdivb(6)=0.5D0*(sdivb(6)+vj)
               sdivb(7)=0.5D0*(sdivb(7)+vk)
               sdivb(8)=0.5D0*(sdivb(8)+(vi*bi+vj*bj+vk*bk))
               rho=w(i,j,k,1)
               vi=w(i,j,k,inorm)
               vj=w(i,j,k,itan)
               vk=w(i,j,k,itanb)
               bi=w(i,j,k,inorm+3)
               bj=w(i,j,k,itan+3)
               bk=w(i,j,k,itanb+3)
               uc(1)=rho
               uc(2)=rho*vi
               uc(3)=rho*vj
               uc(4)=rho*vk
               uc(5)=bi
               uc(6)=bj
               uc(7)=bk
               uc(8)=0.5D0*rho*(vi*vi+vj*vj+vk*vk)+p/(gamma-1.D0)+
     &           0.5D0*(bi*bi+bj*bj+bk*bk)
c
               do l=1,8,1
                  uc(l)=uc(l)-dtodx*(fr(l)-fl(l))
               enddo
c     8-wave non-cons source
               do l=2,8,1
                  uc(l)=uc(l)-dtodx*divb*sdivb(l)
               enddo
               rho=uc(1)
               vi=uc(2)/uc(1)
               vj=uc(3)/uc(1)
               vk=uc(4)/uc(1)
               bi=uc(5)
               bj=uc(6)
               bk=uc(7)
               p=(gamma-1.D0)*(uc(8)-0.5D0*rho*(vi*vi+vj*vj+vk*vk)
     &              -0.5D0*(bi*bi+bj*bj+bk*bk))
               adwdx(i,j,k,1)=rho-w(i,j,k,1)
               adwdx(i,j,k,inorm)=vi-w(i,j,k,inorm)
               adwdx(i,j,k,itan)=vj-w(i,j,k,itan)
               adwdx(i,j,k,itanb)=vk-w(i,j,k,itanb)
               adwdx(i,j,k,inorm+3)=bi-w(i,j,k,inorm+3)
               adwdx(i,j,k,itan+3)=bj-w(i,j,k,itan+3)
               adwdx(i,j,k,itanb+3)=bk-w(i,j,k,itanb+3)
               adwdx(i,j,k,8)=p-w(i,j,k,8)
            enddo
         enddo
      enddo
c
      return
      end
c-----------------------------------------------------------------------

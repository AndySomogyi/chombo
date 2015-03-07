C -*- Mode: Fortran; -*-
c-----------------------------------------------------------------
c       Ravi Samtaney
c       Copyright 2002
c       Princeton Plasma Physics Laboratory
c       All Rights Reserved
c-----------------------------------------------------------------
c       $Log: newdt.F,v $
c       Revision 1.1  2005/05/12 19:22:45  ligocki
c       Moved current "{exec,src}MHD" to "{exec,src}IdealMHD"...
c
c       Revision 1.3  2005/04/13 00:11:37  samtaney
c       Intermediate checkin - passes wave tests.
c
c       Revision 1.2  2004/10/06 20:52:57  samtaney
c       Intermediate checkin
c
c       Revision 1.1.1.1  2004/07/09 19:27:04  samtaney
c       Original source for Unsplit Ideal MHD based on new Chombo interface
c
c       Revision 1.2  2002/07/09 19:00:59  samtaney
c       Intermediate commit of files towards a single 2D/3D Ideal MHD code.
c
c       Revision 1.1  2002/06/13 22:16:18  samtaney
c       Original unsplit 8-wave Ideal MHD source (includes B filter).
c
c-----------------------------------------------------------------
       subroutine maxWaveSpeedMHD(cmax,ux,iulo,iuhi,julo,juhi,
     &                    kulo,kuhi,
     &                    ilo,ihi,jlo,jhi,
     &                     klo,khi)
c----------------------------------------------------------------------
c
      implicit none
c
#include  "EosCommon.fh"
c
c     Incoming arguments
        integer iulo,iuhi,julo,juhi,kulo,kuhi
	integer ilo,ihi,jlo,jhi,klo,khi
        real*8 cmax,dx
	real*8 ux(iulo:iuhi,
     &            julo:juhi,
     &            kulo:kuhi,8)
c     &            1:1,8)
c       
c     Local variables
	real*8 rho
	real*8 u,v,w,press
	real*8 bi,bj,bk
	real*8 rhoInv
	real*8 Asq
	real*8 axsq,aysq,azsq
	real*8 csndsq
	real*8 cfastx,cfasty,cfastz
c	real*8 cslowx,cslowy,cslowz
	real*8 csnd,wavespeedx,wavespeedy,wavespeedz
	real*8 tmp,tmpx,tmpy,tmpz
        real*8 diffmax,etaPrime
	integer i,j,k

c#if CH_SPACEDIM==2
c        write(6,*) 'MAXWAVE SPEED in 2D'
c#elif CH_SPACEDIM==3
c        write(6,*) 'MAXWAVE SPEED in 3D'
c#endif
c       compute maximum velocity and adjust time step to satisfy
c       cfl criterion
c        write(6,*) 'TIME STEP',gamma,rgas,eta,kappa,mu
        cmax=0.D0
        diffmax=0.D0
	do k = klo,khi
           do j = jlo,jhi
              do i = ilo,ihi
                 rho = ux(i,j,k,1)
                 u = ux(i,j,k,2)/rho
                 v = ux(i,j,k,3)/rho
                 w = ux(i,j,k,4)/rho
		 bi=ux(i,j,k,5)
		 bj=ux(i,j,k,6)
		 bk=ux(i,j,k,7)
		 rhoInv=1.D0/ux(i,j,k,1)
		 press=(gamma-1.D0)*(ux(i,j,k,8)-0.5D0*rho
     &		      *(u*u+v*v+w*w)-0.5D0*(bi*bi+bj*bj+bk*bk))
c       
		 axsq=bi**2*rhoInv
		 aysq=bj**2*rhoInv
		 azsq=bk**2*rhoInv
		 Asq=(bi**2+bj**2+bk**2)*rhoInv
		 csndsq=gamma*press*rhoInv
		 tmpx=dsqrt(dmax1(((csndsq+Asq)**2-4.0*csndsq*axsq),
     &		      0.D0))
		 tmpy=dsqrt(dmax1(((csndsq+Asq)**2-4.0*csndsq*aysq),
     &		      0.D0))
		 tmpz=dsqrt(dmax1(((csndsq+Asq)**2-4.0*csndsq*azsq),
     &		      0.D0))
		 tmp=csndsq+Asq
		 cfastx=dsqrt(0.5D0*(tmp+tmpx))
c		 cslowx=dsqrt(0.5D0*(tmp-tmpx))
		 cfasty=dsqrt(0.5D0*(tmp+tmpy))
c		 cslowy=dsqrt(0.5D0*(tmp-tmpy))
		 cfastz=dsqrt(0.5D0*(tmp+tmpz))
c		 cslowz=dsqrt(0.5D0*(tmp-tmpz))
		 csnd=dsqrt(csndsq)
		 waveSpeedx=dabs(u)+dmax1(cfastx,dsqrt(axsq))
		 waveSpeedy=dabs(v)+dmax1(cfasty,dsqrt(aysq))
		 waveSpeedz=dabs(w)+dmax1(cfastz,dsqrt(azsq))
c		 waveSpeedx=dabs(u)+dsqrt(csndsq+asq)
c		 waveSpeedy=dabs(v)+dsqrt(csndsq+asq)
c		 waveSpeedz=dabs(w)+dsqrt(csndsq+asq)
c$$$                 xspmax = dmin1(xspmax,dx/(dabs(u)
c$$$     &                +waveSpeedx))
c$$$                 yspmax = dmin1(yspmax,dy/(dabs(v)
c$$$     &                +waveSpeedy))
c                 zspmax = dmin1(zspmax,dz/(dabs(w)
c     &                +waveSpeed))
c                 dt_local = dmin1(xspmax,yspmax,zspmax)*cfl
                 cmax = dmax1(cmax,waveSpeedx,waveSpeedy,waveSpeedz)
c                 cmax = 
c$$$                 if(bk.le.0.D0) then
c$$$                    etaPrime=eta+
c$$$     &                etafactor*eta*(1.D0-dexp(-(bk/0.075D0/7.8D0)**2))
c$$$c                    write(6,*) 'etaprime=', i,j,etaprime
c$$$                 else
c$$$                   etaPrime=eta
c$$$                 endif
c$$$                  diffmax = dmax1(2.D0*etaPrime/dx,diffmax)
c$$$                  diffmax = dmax1(2.D0*mu/rho/dx,diffmax)
c$$$                  diffmax = dmax1(2.D0*kappa/dx,diffmax)
c$$$                  if(cmax.gt.diffmax) then
c$$$                     write(6,*) 'Convection dominated',cmax,diffmax
c$$$                  else
c$$$                     write(6,*) 'Diffusion dominated',cmax,diffmax
c$$$                  endif
c                  cmax=dmax1(cmax,diffmax)
c
              enddo
           enddo
        enddo
c$$$        if(cmax.gt.diffmax) then
c$$$           write(6,*) 'Convection dominated',cmax,diffmax
c$$$        else
c$$$           write(6,*) 'Diffusion dominated',cmax,diffmax
c$$$        endif
c     For some reason code is more robust when cfl no. is small
c     when we are diffusion dominated. Hence the factor 5 below - it
c     is a rotten kludge which needs to be addressed later.
c        cmax=dmax1(cmax,5.D0*diffmax)
c        write(6,*) 'MAX Wave speed',cmax
        return
	end




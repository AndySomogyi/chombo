C     -*- Mode: Fortran; -*-
c-----------------------------------------------------------------------
c     Ravi Samtaney
c     Copyright 2002
c     Princeton Plasma Physics Laboratory
c     All Rights Reserved
c-----------------------------------------------------------------------
c     $Log: inviscidfluxrp.F,v $
c     Revision 1.1  2005/05/12 19:22:45  ligocki
c     Moved current "{exec,src}MHD" to "{exec,src}IdealMHD"...
c
c     Revision 1.2  2004/10/06 20:52:56  samtaney
c     Intermediate checkin
c
c     Revision 1.1.1.1  2004/07/09 19:27:04  samtaney
c     Original source for Unsplit Ideal MHD based on new Chombo interface
c
c     Revision 1.2  2002/07/09 19:00:59  samtaney
c     Intermediate commit of files towards a single 2D/3D Ideal MHD code.
c
c     Revision 1.1  2002/06/13 22:16:18  samtaney
c     Original unsplit 8-wave Ideal MHD source (includes B filter).
c
c-----------------------------------------------------------------------
      subroutine RiemannSolverMHD(um,ul,ur,
     &     iumlo,iumhi,jumlo,jumhi,kumlo,kumhi,
     &     iullo,iulhi,jullo,julhi,kullo,kulhi,
     &     iurlo,iurhi,jurlo,jurhi,kurlo,kurhi,
     &     ilo,ihi,jlo,jhi,klo,khi,
     &     direction,llo,lhi)
c
      implicit none
c
#include "EosCommon.fh"
c     
c     incoming arguments
      integer iumlo,iumhi,jumlo,jumhi,kumlo,kumhi,
     &     iullo,iulhi,jullo,julhi,kullo,kulhi,
     &     iurlo,iurhi,jurlo,jurhi,kurlo,kurhi,
     &     ilo,ihi,jlo,jhi,klo,khi,
     &     direction,llo,lhi
c
      real*8 ul(iullo:iulhi,jullo:julhi,kullo:kulhi,8)
      real*8 ur(iurlo:iurhi,jurlo:jurhi,kurlo:kurhi,8)
      real*8 um(iumlo:iumhi,jumlo:jumhi,kumlo:kumhi,8)
c
c     Local variables
      integer i,j,k,l,m
c     
      real*8 vm(llo:lhi,8)
      real*8 vl(llo:lhi,8)
      real*8 vr(llo:lhi,8)
c
c     This subroutine determines the solution to the Riemann problem
c     at the cell interfaces. The solution is based upon linearization
c     and no guarantee can be made that it is accurate.
c     v: Vector of variables - rho, u,v,w,Bx,By,Bz,p
c     rho=density, u=x-velocity, v=y-velocity
c     w=z-velocity
c     Bx,By,Bz=x,y,z components of the magnetic field
c     p=thermal pressure
c     
c     vl(i) = Left state at the i-th interface
c     vr(i) = Right state at the i-th interface
c     
c     
c$$$      write(6,*) 'RiemannSolverMHD', iumlo,iumhi,jumlo,jumhi,
c$$$     &     iullo,iulhi,jullo,julhi,
c$$$     &     iurlo,iurhi,jurlo,jurhi,
c$$$     &     ilo,ihi,jlo,jhi,klo,khi,
c$$$     &     direction,llo,lhi
      if(direction.eq.0) then
         do k=klo,khi,1
            do j=jlo,jhi,1
               do i=ilo,ihi,1
                  do m=1,8,1
                     vl(i,m)=ul(i,j,k,m)
                     vr(i,m)=ur(i,j,k,m)
                  enddo
               enddo
               call RiemannProblem(vm,vl,vr,ilo,ihi,direction)
c               if(j.eq.jlo) then
c               do i=ilo,ihi-1,1
c                  write(6,*) 'RPXU ',i,vm(i,8)+0.5D0*(vm(i,6)**2+
c     &                 vm(i,7)**2-vm(i,5)**2),
c     &             vm(i+1,8)+0.5D0*(vm(i+1,6)**2+
c     &                 vm(i+1,7)**2-vm(i+1,5)**2),
c                  write(6,*) 'RPXU ',i,vm(i+1,8)-vm(i,8),
c     &                 vm(i+1,5)**2-vm(i,5)**2,
c     &             vm(i+1,6)**2-vm(i,6)**2,
c     &                 vm(i+1,7)**2-vm(i,7)**2
c                  write(6,*) 'RPXV ',i,vm(i,5)*vm(i,6),
c     &                 vm(i+1,5)*vm(i+1,6)
c                  write(6,*) 'RPBX ',i,vm(i,5),
c     &                 vl(i,5),vr(i,5),vm(i+1,5)
c                  write(6,*) 'RPUX ',i,vm(i,2),
c     &                 vl(i,2),vr(i,2),vm(i+1,2)
c               enddo
c               endif
               do i=ilo,ihi,1
                  do m=1,8,1
                     um(i,j,k,m)=vm(i,m)
                  enddo
               enddo
c$$$               if(j.eq.16) then
c$$$                  do i=ilo,ihi
c$$$                     do l=1,8,1
c$$$                     write(6,*) 'RP1:',i,l,um(i,j,k,l),
c$$$     &                       ul(i,j,k,l), ur(i,j,k,l)
c$$$                     enddo
c$$$                  enddo
c$$$               endif
            enddo
         enddo
      elseif (direction.eq.1) then
         do k=klo,khi,1
            do i=ilo,ihi,1
               do j=jlo,jhi,1
                  do m=1,8,1
                     vl(j,m)=ul(i,j,k,m)
                     vr(j,m)=ur(i,j,k,m)
                  enddo
               enddo
               call RemapToLocalCoords(vl,jlo,jhi,2)
               call RemapToLocalCoords(vr,jlo,jhi,2)
               call RiemannProblem(vm,vl,vr,jlo,jhi,direction)
c$$$               if(i.eq.ilo) then
c$$$               do j=jlo,jhi-1,1
c$$$               write(6,*) 'RPY ',j,vm(j,8)+0.5D0*(vm(j,6)**2+
c$$$     &                 vm(j,7)**2-vm(j,5)**2),
c$$$     &             vm(j+1,8)+0.5D0*(vm(j+1,6)**2+
c$$$     &                 vm(j+1,7)**2-vm(j+1,5)**2)
c$$$               enddo
c$$$               endif
               call RemapToGlobalCoords(vm,jlo,jhi,2)
               do j=jlo,jhi,1
                  do m=1,8,1
                     um(i,j,k,m)=vm(j,m)
                  enddo
               enddo
c$$$               if(i.eq.200) then
c$$$                  do j=jlo,jhi
c$$$                     do l=1,8,1
c$$$                     write(6,*) 'RP2:',j,l,um(i,j,k,l),
c$$$     &                       ul(i,j,k,l), ur(i,j,k,l)
c$$$                     enddo
c$$$                  enddo
c$$$               endif
            enddo
         enddo
      elseif (direction.eq.2) then
         do j=jlo,jhi,1
            do i=ilo,ihi,1
               do k=klo,khi,1
                  do m=1,8,1
                     vl(k,m)=ul(i,j,k,m)
                     vr(k,m)=ur(i,j,k,m)
                  enddo
               enddo
               call RemapToLocalCoords(vl,klo,khi,3)
               call RemapToLocalCoords(vr,klo,khi,3)
               call RiemannProblem(vm,vl,vr,klo,khi,direction)
               call RemapToGlobalCoords(vm,klo,khi,3)
               do k=klo,khi,1
                  do m=1,8,1
                     um(i,j,k,m)=vm(k,m)
                  enddo
               enddo
            enddo
         enddo
      else
         write(6,*) 'FATAL ERROR: direcion !=0,1,2'
c     
      endif
c     
      return
      end
c     
c-----------------------------------------------------------------------
      subroutine RiemannProblem(vm,vl,vr,ilo,ihi,direction)
c     
      implicit none
c     
#include "EosCommon.fh"
c     
c     Incoming arguments
      integer ilo, ihi, direction
      real*8 vl(ilo:ihi,8)
      real*8 vr(ilo:ihi,8)
      real*8 vm(ilo:ihi,8)
c     
c     Local variables
      real*8 alamda(ILO:IHI,8)
      real*8 evl(ILO:IHI,8,8)
      real*8 evr(ILO:IHI,8,8)
      real*8 durl(ILO:IHI,8)
      real*8 alpha(ILO:IHI,8)
      real*8 utilde(ilo:ihi,8)
      integer i,m
c     
c     v: Vector of variables - rho, u,v,w,Bx,By,Bz,p
c     rho=density, u=x-velocity, v=y-velocity
c     w=z-velocity
c     Bx,By,Bz= Magnetic Field
c     p=pressure
c     
c     utilde = vector of Roe averaged variables.
c     utilde(1,2,3,4,5,6,7,8) = {rhotilde, utilde, vtilde, wtilde, 
c     Bxtilde,Bytilde,Bztilde,ptilde}
c
c     Set the Roe averaged variables depending  on vl,vr
      call SetRoeVariables(utilde,vl,vr,ilo,ihi,direction)
c     
      do i=ilo,ihi,1
         do m=1,8,1
            alamda(i,m)=0.D0
            durl(i,m)=vr(i,m)-vl(i,m)
         enddo
      enddo
c     
      call SetEigenSystem(utilde,evl,evr,alamda,ilo,ihi,direction)
      call SetAlphas(utilde,alpha,durl,evl,evr,ilo,ihi,direction)
      call RPSolve(evl,evr,alamda,alpha,vl,vr,vm,ilo,ihi,direction)
c     
c$$$      if(direction.eq.1) then
c$$$      do i=ilo,ihi,1
c$$$         do m=1,8,1
c$$$            write(6,*) 'ALPHA=',m,i,alpha(i,m)
c$$$            write(6,*) 'LAMDA=',m,i,alamda(i,m)
c$$$            write(6,*) 'EVR5=',i,m,evr(i,5,m)
c$$$         enddo
c$$$      enddo
c$$$      endif
c$$$      if(ilo.eq.ihi) then
c$$$         write(6,*) 'RPSOLVE Utilde: ',utilde
c$$$         write(6,*) 'RPSOLVE durl: ',durl
c$$$         write(6,*) 'RPSOLVE vl: ',vl
c$$$         write(6,*) 'RPSOLVE vr: ',vr
c$$$         write(6,*) 'RPSOLVE vm: ',vm
c$$$         write(6,*) 'RPSOLVE alamda: ',alamda
c$$$         write(6,*) 'RPSOLVE alpha: ',alpha
c$$$         do i=1,8,1
c$$$            do m=1,8,1
c$$$               write(6,*) 'RPSOLVE evr: ',i,m,evr(1,m,i),
c$$$            enddo
c$$$         enddo
c$$$      endif
      return
      end

c-----------------------------------------------------------------------
c     Routine to map velocities to local coordinates.
      subroutine RemapToLocalCoords(vx,ilo,ihi,direction)
c     
      implicit none
c     
#include "EosCommon.fh"
c     
c     
c     Incoming arguments
      integer ilo, ihi
      integer direction
      real*8 vx(ilo:ihi,8)
c     
c     Local variables
      real*8 tmp(ilo:ihi,3)
      integer i
c     
      if(direction.eq.2) then
         do i=ilo,ihi,1
            tmp(i,1)=vx(i,3)
            tmp(i,2)=vx(i,4)
            tmp(i,3)=vx(i,2)
            vx(i,2)=tmp(i,1)
            vx(i,3)=tmp(i,2)
            vx(i,4)=tmp(i,3)
            tmp(i,1)=vx(i,6)
            tmp(i,2)=vx(i,7)
            tmp(i,3)=vx(i,5)
            vx(i,5)=tmp(i,1)
            vx(i,6)=tmp(i,2)
            vx(i,7)=tmp(i,3)
         enddo
      elseif(direction.eq.3) then
         do i=ilo,ihi,1
            tmp(i,1)=vx(i,4)
            tmp(i,2)=vx(i,2)
            tmp(i,3)=vx(i,3)
            vx(i,2)=tmp(i,1)
            vx(i,3)=tmp(i,2)
            vx(i,4)=tmp(i,3)
            tmp(i,1)=vx(i,7)
            tmp(i,2)=vx(i,5)
            tmp(i,3)=vx(i,6)
            vx(i,5)=tmp(i,1)
            vx(i,6)=tmp(i,2)
            vx(i,7)=tmp(i,3)
         enddo
      endif
c     
      return
      end
c     
c-----------------------------------------------------------------------
      subroutine RemapToGlobalCoords(vx,ilo,ihi,direction)
      implicit none
c
#include "EosCommon.fh"
c     
c     Routine to map velocities to global coordinates.
c     Incoming arguments
      integer ilo, ihi
      integer direction
      real*8 vx(ilo:ihi,8)
c     
c     Local variables
      real*8 tmp(ilo:ihi,3)
      integer i
c     
      if(direction.eq.2) then
         do i=ilo,ihi,1
            tmp(i,2)=vx(i,2)
            tmp(i,3)=vx(i,3)
            tmp(i,1)=vx(i,4)
            vx(i,2)=tmp(i,1)
            vx(i,3)=tmp(i,2)
            vx(i,4)=tmp(i,3)
            tmp(i,2)=vx(i,5)
            tmp(i,3)=vx(i,6)
            tmp(i,1)=vx(i,7)
            vx(i,5)=tmp(i,1)
            vx(i,6)=tmp(i,2)
            vx(i,7)=tmp(i,3)
         enddo
      elseif(direction.eq.3) then
         do i=ilo,ihi,1
            tmp(i,3)=vx(i,2)
            tmp(i,1)=vx(i,3)
            tmp(i,2)=vx(i,4)
            vx(i,2)=tmp(i,1)
            vx(i,3)=tmp(i,2)
            vx(i,4)=tmp(i,3)
            tmp(i,3)=vx(i,5)
            tmp(i,1)=vx(i,6)
            tmp(i,2)=vx(i,7)
            vx(i,5)=tmp(i,1)
            vx(i,6)=tmp(i,2)
            vx(i,7)=tmp(i,3)
         enddo
      endif
c     
      return
      end
c     
c-----------------------------------------------------------------------
c     The actual Riemann solver
      subroutine RPSolve(evl,evr,alamda,alpha,vl,vr,vm,ilo,ihi,idir)
c
      implicit none
c
#include "EosCommon.fh"
c
c     Incoming arguments
      integer ilo,ihi
      integer idir
      real*8 vl(ilo:ihi,8)
      real*8 vr(ilo:ihi,8)
      real*8 vm(ilo:ihi,8)
      real*8 alamda(ILO:IHI,8)
      real*8 evl(ILO:IHI,8,8)
      real*8 evr(ILO:IHI,8,8)
      real*8 alpha(ILO:IHI,8)
c
c     Local variables
      integer i,k,l,m
      real*8 vml(ilo:ihi,8)
      real*8 vmr(ilo:ihi,8)
c     
c     
      do i=ilo,ihi,1
         do m=1,8,1
            vmr(i,m)=vr(i,m)
            vml(i,m)=vl(i,m)
         enddo
      enddo
      do i=ilo,ihi,1
c     do l=1,8,1
c     Original RP
c$$$         do k=1,8,1
c$$$            if(alamda(i,k).ge.0.D0) then
c$$$               do l=1,8,1
c$$$                  vmr(i,l)=vmr(i,l)-alpha(i,k)*evr(i,l,k)
c$$$               enddo
c$$$            else 
c$$$               do l=1,8,1
c$$$                  vml(i,l)=vml(i,l)+alpha(i,k)*evr(i,l,k)
c$$$               enddo
c$$$            endif
c$$$         enddo
c     Modified RP solve on 09/07/02
         do k=1,8,1
            if(alamda(i,k).gt.0.D0) then
               do l=1,8,1
                  vmr(i,l)=vmr(i,l)-alpha(i,k)*evr(i,l,k)
               enddo
            else if(alamda(i,k).lt.0.D0) then
               do l=1,8,1
                  vml(i,l)=vml(i,l)+alpha(i,k)*evr(i,l,k)
               enddo
c             else 
c                do l=1,8,1
c                   vml(i,l)=vml(i,l)+alpha(i,k)*evr(i,l,k)
c                   vmr(i,l)=vmr(i,l)-alpha(i,k)*evr(i,l,k)
c                enddo
            endif
         enddo
c     enddo
      enddo
c     
c     vm=0.5D0*(vml+vmr)
      do i=ilo,ihi,1
         do m=1,8,1
c     Original
c            vm(i,m)=vml(i,m)
c     Modified on 09/07/02
            vm(i,m)=0.5D0*(vmr(i,m)+vml(i,m))
         enddo
         vm(i,5)=0.5D0*(vl(i,5)+vr(i,5))
      enddo
      return
      end

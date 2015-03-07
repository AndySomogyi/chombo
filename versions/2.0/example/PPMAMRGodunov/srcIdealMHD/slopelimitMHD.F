C     -*- Mode: Fortran; -*-
c-----------------------------------------------------------------------
c     Ravi Samtaney
c     Copyright 2002
c     Princeton Plasma Physics Laboratory
c     All Rights Reserved
c-----------------------------------------------------------------------
c     $Log: slopelimitMHD.F,v $
c     Revision 1.1  2005/05/12 19:22:45  ligocki
c     Moved current "{exec,src}MHD" to "{exec,src}IdealMHD"...
c
c     Revision 1.2  2004/10/06 20:52:57  samtaney
c     Intermediate checkin
c
c     Revision 1.1.1.1  2004/07/09 19:27:04  samtaney
c     Original source for Unsplit Ideal MHD based on new Chombo interface
c
c     Revision 1.4  2002/07/19 13:41:23  samtaney
c     Fixed order of arguments.
c
c     Revision 1.3  2002/07/10 13:18:14  samtaney
c     Fixed typos - jlo,jhi were used instead of klo,khi in a call to constructlrstates.
c
c     Revision 1.2  2002/07/09 19:00:59  samtaney
c     Intermediate commit of files towards a single 2D/3D Ideal MHD code.
c
c     Revision 1.1  2002/06/13 22:16:18  samtaney
c     Original unsplit 8-wave Ideal MHD source (includes B filter).
c
c-----------------------------------------------------------------------
c      subroutine NormalPredMHD(ul,ur,ux,dudx,
c     &     iullo,iulhi,jullo,julhi,kullo,kulhi,
c     &     iurlo,iurhi,jurlo,jurhi,kurlo,kurhi,
      subroutine slopeLimitMHD(dul,dur,dudx,ux,
     &     iullo,iulhi,jullo,julhi,kullo,kulhi,
     &     iurlo,iurhi,jurlo,jurhi,kurlo,kurhi,
     &     idudxlo,idudxhi,jdudxlo,jdudxhi,kdudxlo,kdudxhi,
     &     iuxlo,iuxhi,juxlo,juxhi,kuxlo,kuxhi,
     &     ilo,ihi,jlo,jhi,klo,khi,
     &     direction,llo,lhi)
c
      implicit none
c
#include "EosCommon.fh"
c     
c     Incoming arguments
      integer iullo,iulhi,jullo,julhi,kullo,kulhi,
     &     iurlo,iurhi,jurlo,jurhi,kurlo,kurhi,
     &     iuxlo,iuxhi,juxlo,juxhi,kuxlo,kuxhi,
     &     idudxlo,idudxhi,jdudxlo,jdudxhi,kdudxlo,kdudxhi,
     &     ilo,ihi,jlo,jhi,klo,khi
      integer llo,lhi,direction
c     
      real*8 dul(iullo:iulhi,jullo:julhi,kullo:kulhi,8)
      real*8 dur(iurlo:iurhi,jurlo:jurhi,kurlo:kurhi,8)
      real*8 ux(iuxlo:iuxhi,juxlo:juxhi,kuxlo:kuxhi,8)
      real*8 dudx(idudxlo:idudxhi,jdudxlo:jdudxhi,
     &     kdudxlo:kdudxhi,8)
c     
c     Local variables
      integer i,j,k,l,m
      real*8 vx(llo:lhi,8)
      real*8 dvdx(llo:lhi,8)
      real*8 dvl(llo:lhi,8)
      real*8 dvr(llo:lhi,8)
c     
c$$$  write(6,*) 'Normalpred', iullo,iulhi,jullo,julhi,
c$$$  &                         iurlo,iurhi,jurlo,jurhi,
c$$$  &                         iuxlo,iuxhi,juxlo,juxhi,
c$$$  &                         idudxlo,idudxhi,jdudxlo,jdudxhi,
c$$$  &                         ilo,ihi,jlo,jhi,klo,khi
c     
c      return
      if(direction.eq.0) then
         do k=klo,khi,1
            do j=jlo,jhi,1
               do i=ilo,ihi,1
                  do m=1,8,1
                     vx(i,m)=ux(i,j,k,m)
                     dvL(i,m)=duL(i,j,k,m)
                     dvR(i,m)=duR(i,j,k,m)
                  enddo
               enddo
c               call LimitSlopeCh(dvdx,dvL,dvR,vx,ilo,ihi,direction)
               do i=ilo,ihi,1
                  do m=1,8,1
                     dudx(i,j,k,m)=dvdx(i,m)
                  enddo
              enddo
            enddo
         enddo
      elseif (direction.eq.1) then
         do k=klo,khi,1
            do i=ilo,ihi,1
               do j=jlo,jhi,1
                  do m=1,8,1
                     vx(j,m)=ux(i,j,k,m)
                     dvL(i,m)=duL(i,j,k,m)
                     dvR(i,m)=duR(i,j,k,m)
                  enddo
               enddo
               call RemapToLocalCoords(vx,jlo,jhi,2)
               call RemapToLocalCoords(dvl,jlo,jhi,2)
               call RemapToLocalCoords(dvr,jlo,jhi,2)
c               call LimitSlopeCh(dvdx,dvL,dvR,vx,jlo,jhi,direction)
               call RemapToGlobalCoords(dvdx,jlo,jhi,2)
               do j=jlo,jhi,1
                  do m=1,8,1
                     dudx(i,j,k,m)=dvdx(j,m)
                  enddo
               enddo
            enddo
         enddo
      elseif (direction.eq.2) then
         do j=jlo,jhi,1
            do i=ilo,ihi,1
               do k=klo,khi,1
                  do m=1,8,1
                     vx(k,m)=ux(i,j,k,m)
                     dvl(k,m)=dul(i,j,k,m)
                     dvr(k,m)=dur(i,j,k,m)
                  enddo
               enddo
               call RemapToLocalCoords(vx,klo,khi,3)
               call RemapToLocalCoords(dvl,klo,khi,3)
               call RemapToLocalCoords(dvr,klo,khi,3)
               call LimitSlopeCh(dvdx,dvL,dvR,vx,klo,khi,direction)
               call RemapToGlobalCoords(dvdx,klo,khi,3)
               do k=klo,khi,1
                  do m=1,8,1
                     dudx(i,j,k,m)=dvdx(k,m)
                  enddo
               enddo
            enddo
         enddo
c     
      else
         write(6,*) 'FATAL ERROR: direcion !=0,1,2'
c     
      endif
c     
      return
      end
c
c-----------------------------------------------------------------------

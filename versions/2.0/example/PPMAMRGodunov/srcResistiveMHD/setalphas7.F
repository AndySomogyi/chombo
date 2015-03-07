C -*- Mode: Fortran; -*-
c-----------------------------------------------------------------------
c     Ravi Samtaney
c     Copyright 2005
c     Princeton Plasma Physics Laboratory
c     All Rights Reserved
c-----------------------------------------------------------------------
c     $Log: setalphas7.F,v $
c     Revision 1.1  2005/05/26 20:28:07  samtaney
c     These are used in the 7 wave formulation; this actually uses 7 eigenvectors
c     while previously we were using 8 and ignoring one. Essentially the same but
c     leads to less computations.
c
c-----------------------------------------------------------------------
      subroutine SetAlphas(utilde,alpha,durl,evl,evr,ilo,ihi,direction)
c
      implicit none
c
c     Incoming arguments
      integer ilo,ihi,direction
      real*8 utilde(ilo:ihi,8)
      real*8 alpha(ilo:ihi,7)
      real*8 durl(ilo:ihi,7)
      real*8 evl(ilo:ihi,7,7)
      real*8 evr(ilo:ihi,7,7)
c
c     Local variables
      integer i,k,l

c     
      do i=ilo,ihi,1
         do l=1,7,1
            alpha(i,l)=0.D0
            do k=1,7,1
               alpha(i,l)=alpha(i,l)+evl(i,l,k)*durl(i,k)
            enddo
         enddo
      enddo
c
      return
      end 





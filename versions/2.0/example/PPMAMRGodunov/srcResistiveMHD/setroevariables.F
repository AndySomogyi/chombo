C -*- Mode: Fortran; -*-
c-----------------------------------------------------------------------
c     Ravi Samtaney
c     Copyright 2002
c     Princeton Plasma Physics Laboratory
c     All Rights Reserved
c-----------------------------------------------------------------------
c     $Log: setroevariables.F,v $
c     Revision 1.1  2005/05/12 19:33:37  ligocki
c     Copied "IdealMHD" to "ResistiveMHD" as a starting point...
c
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
        subroutine SetRoeVariables(utilde,vl,vr,ilo,ihi,direction)
c
        implicit none
c  
c     Incoming arguments
        integer ilo,ihi
        integer direction
        real*8 vl(ilo:ihi,8)
        real*8 vr(ilo:ihi,8)
        real*8 utilde(ilo:ihi,8)
c
c     Local variables
        integer i,l
        real*8 rhol, rhor

c     Simple arithmetic average
        do i=ilo,ihi,1
           do l=1,8,1
c              rhol=sqrt(vl(i,1))
c              rhor=sqrt(vr(i,1))
               utilde(i,l)=0.5D0*(vl(i,l)+ vr(i,l))
c              utilde(i,l)=(vl(i,l)*rhol+ vr(i,l)*rhor)/(rhol+rhor)
c     utilde(i,l)=vl(i,l)
           enddo
        enddo	
	return
	end

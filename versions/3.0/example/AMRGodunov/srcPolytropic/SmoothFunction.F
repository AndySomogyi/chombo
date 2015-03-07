#include "REAL.H"
#include "CONSTANTS.H"
      real_t function smoothened(x)
      implicit none
      real_t x

      if (abs(x) .ge. half) then
         smoothened = zero
      else
         smoothened = (sin(Pi * (x + half)))**6
      endif

      return
      end

#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/

// Constants used between C++ and Fortran to index the conserved and
// primitive variables and to know the number of each.  Also, the number
// of slopes and fluxes are defined here.
#endif

#ifndef _LGINTEGRATOR_H_
#define _LGINTEGRATOR_H_

#if (CH_SPACEDIM == 2)

#define WRHO     0
#define WVELX    1
#define WVELY    2
#define WVELZ    3
#define WBX      4
#define WBY      5
#define WBZ      6
#define WPRES    7
#define WTEMP    7
#define WC       8

#define WNUM     8
#define WSLOPE   8

#define URHO     0
#define UMOMX    1
#define UMOMY    2
#define UMOMZ    3
#define UBX      4
#define UBY      5
#define UBZ      6
#define UENG     7
#define DIVB     8
#define UFACE    9
#define VFACE   10
#define WFACE   11
#define BXFACE  12
#define BYFACE  13
#define BZFACE  14

#define UNUM     8
#define UFLUX   15

#endif

#if (CH_SPACEDIM == 3)

#define WRHO     0
#define WVELX    1
#define WVELY    2
#define WVELZ    3
#define WBX      4
#define WBY      5
#define WBZ      6
#define WPRES    7
#define WTEMP    7
#define WC       8

#define WNUM     8
#define WSLOPE   8

#define URHO     0
#define UMOMX    1
#define UMOMY    2
#define UMOMZ    3
#define UBX      4
#define UBY      5
#define UBZ      6
#define UENG     7
#define DIVB     8
#define UFACE    9
#define VFACE   10
#define WFACE   11
#define BXFACE  12
#define BYFACE  13
#define BZFACE  14

#define UNUM     8
#define UFLUX    8

#endif

#endif

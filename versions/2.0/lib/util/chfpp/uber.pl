#!/usr/local/gnu/bin/perl -w

####  _______              __
#### / ___/ /  ___  __ _  / /  ___
####/ /__/ _ \/ _ \/  ' \/ _ \/ _ \
####\___/_//_/\___/_/_/_/_.__/\___/
####
####
#### This software is copyright (C) by the Lawrence Berkeley
#### National Laboratory.  Permission is granted to reproduce
#### this software for non-commercial purposes provided that
#### this notice is left intact.
####
#### It is acknowledged that the U.S. Government has rights to
#### this software under Contract DE-AC03-765F00098 between
#### the U.S.  Department of Energy and the University of
#### California.
####
#### This software is provided as a professional and academic
#### contribution for joint exchange. Thus it is experimental,
#### is provided ``as is'', with no warranties of any kind
#### whatsoever, no support, no promise of updates, or printed
#### documentation. By using this software, you acknowledge
#### that the Lawrence Berkeley National Laboratory and
#### Regents of the University of California shall have no
#### liability with respect to the infringement of other
#### copyrights by any part of this software.
####
#### perl uber.pl  -f input  -p fortranout -c cout -D $(DIM) -d

#########################################
#### This is an augmented version of ChomboFortran.
#### It also removes the thrice-accursed D_TERM
####  DTGraves July 2000
#########################################
use strict;

###parse the command line
###and generate appropriate file names
use ChF;
    &ChF::parseInputs();

####    #install atexit-style handler so that when we exit or die,
####    #we automatically delete this temporary file

    END{unlink($ChF::T1File)
            or die "cannot unlink $ChF::T1File \n";}
    END{unlink($ChF::T2File)
            or die "cannot unlink $ChF::T2File \n";}

    ###need two temp files because final output can be /dev/null
###Strip comments because they are evil
    require "stripcompp.pm";
    StripComProc->import();
    unless(my $flag = &StripComProc::StripComments($ChF::FINFile, $ChF::T2File,
                                                   $ChF::SpaceDim, $ChF::debug))
    {
        die "problem $flag in commentstripper\n";
    }

###parse all the chombo chacha and put fortran
###into T1File and C header into COUT
###this also substitutes CHF_DDECL for CHF_IX
    use File::Basename ;
# extract leading directory and file basename.  Throw away the trailing suffix.
    ($ChF::baseName,$ChF::basePath) = fileparse( $ChF::FINFile ,'\..*' ) ;

    require "subrout.pm";
    SubroutProc->import();
    unless(my $flag =
           &SubroutProc::procSubrout($ChF::T2File, $ChF::COUTFile, $ChF::T1File,
                                     $ChF::SpaceDim, $ChF::debug, $ChF::baseName))
    {
        die "problem $flag in chombofortran\n";
    }

###process multido statements
    require "multidopp.pm";
    MultiDoProc->import();
    unless(my $flag = &MultiDoProc::procMultiDoMacros($ChF::T1File, $ChF::T2File,
                                                      $ChF::SpaceDim, $ChF::debug))
    {
        die "problem $flag in multido processor\n";
    }

###process LBOUND macros.
    require "lboundpp.pm";
    LBoundProc->import();
    unless(my $flag = &LBoundProc::procLBoundMacros($ChF::T2File, $ChF::T1File,
                                                    $ChF::SpaceDim, $ChF::debug))
    {
        die "problem $flag in lbound processor\n";
    }

###process UBOUND macros.
    require "uboundpp.pm";
    UboundProc->import();
    unless(my $flag = &UboundProc::procUboundMacros($ChF::T1File, $ChF::T2File,
                                                    $ChF::SpaceDim, $ChF::debug))
    {
        die "problem $flag in ubound processor\n";
    }

###process NCOMP macros.
    require "ncomppp.pm";
    NcompProc->import();
    unless(my $flag = &NcompProc::procNcompMacros($ChF::T2File, $ChF::T1File,
                                                  $ChF::SpaceDim, $ChF::debug))
    {
        die "problem $flag in ncomp processor\n";
    }

###process CHF_IX macros.
    require "chfixpp.pm";
    ChfixProc->import();
    unless(my $flag = &ChfixProc::procChfixMacros($ChF::T1File, $ChF::T2File,
                                                  $ChF::SpaceDim, $ChF::debug))
    {
        die "problem $flag in chfix processor\n";
    }


###process DTERM macros.
    require "dtermpp.pm";
    DTERMProc->import();
    unless(my $flag = &DTermProc::procDTermMacros($ChF::T2File, $ChF::T1File,
                                                  $ChF::SpaceDim, $ChF::debug))
    {
        die "problem $flag in dterm processor\n";
    }


###process DDECL macros.
    require "ddeclpp.pm";
    DDeclProc->import();
    unless(my $flag = &DDeclProc::procDDeclMacros($ChF::T1File, $ChF::T2File,
                                                  $ChF::SpaceDim, $ChF::debug))
    {
        die "problem $flag in ddecl processor\n";
    }


###process DSELECT macros.
    require "dselectpp.pm";
    DSelectProc->import();
    unless(my $flag = &DSelectProc::procDSelectMacros($ChF::T2File, $ChF::FOUTFile,
                                                      $ChF::SpaceDim, $ChF::debug))
    {
        die "problem $flag in dselect processor\n";
    }



###blow town
    exit 0 ;

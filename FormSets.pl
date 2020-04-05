#!/usr/bin/perl
use strict;
{
    my $target = "TB254";
    chdir "data/$target/complexes/";
    my (@entry) = `ls *Complex.pdb`;
    chomp @entry;
    
    print "PDBid C_RESIDUE INERTIA_3 SMALLEST_SIZE SURFACE_HULL VOLUME_HULL hydrophobic_kyte hydrophobicity_pocket p_Ccoo p_N_atom p_Ooh p_aliphatic_residue p_aromatic_residue p_negative_residue\n";
    foreach(@entry){
        
        my $pdb = substr $_, 0,-12;
        
        my $filename = $pdb."_prox4_pock.pdb";
   
        if(-e $filename && -f _){
            my $f = "Features_".$target."_".$pdb.".rpt";
                if ( -e $f ) {
                    
                    my $awkProp = `awk 'NR==5{print \$0}' $f`;
                    chomp $awkProp;
                    print "$awkProp\n";
                    #last;
                }
        }
    }
}

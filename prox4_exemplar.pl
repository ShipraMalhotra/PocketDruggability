#!/usr/bin/perl
use strict;
{
    chdir "data/TB254/complexes/"; # Target Protein Chain: B and Target residue number: 254
    my (@entry) = `ls *Complex.pdb`;
    chomp @entry;
   
    foreach(@entry){
        get_lig_pock($_);
    }
}

sub get_lig_pock{
    
    my $het = "TMP";
    my $chain  = "A";
    my ($f) = @_;
    my @title = split(/\_/,$f);

    my $pdb = substr $f, 0,4;
    
    my $out = $pdb."_".$title[5]."_".$title[9]."_prox4_pock.pdb";
    my $ligout = $pdb."_".$title[5]."_".$title[9]."_ligand.pdb";

    my @filepocket_lines;
    my @finalpocket;

    my @fdata = `awk '\$1 ~ "HETATM" {print \$0}' $f`;
    chomp @fdata;
    my $oldresnu;
    for my $ligline (@fdata){
        if ($ligline =~ /^HETATM/ && substr($ligline,17,3) eq "$het" && substr ($ligline,21,1) eq $chain){
            chomp $ligline;
            $oldresnu = substr($ligline,22,4); #for multiple copies, consider the first one
            $oldresnu =~ s/^\s*//;
            last;
        }
    }
    open LIG, ">", $ligout;
    for my $ligline (@fdata){
        if ($ligline =~ /^HETATM/ && substr($ligline,17,3) eq "$het" && substr ($ligline,21,1) eq $chain){
            chomp $ligline;
            my $newresnu = substr($ligline,22,4);
            $newresnu =~ s/^\s*//;
            if ($newresnu == $oldresnu){
            
                print LIG $ligline, "\n";
            }
        }
    }
    close LIG;

    my @ligdata = `awk '\$1 ~ "HETATM" {print \$0}' $ligout`;
    chomp @ligdata;

    if (!-z $ligout){

        my @data = `awk '\$1 == "ATOM" {print \$0}' $f`;
        for my $line (@data){
            chomp $line;
            if (substr ($line, 77,1) ne "H" ){
                for my $fline (@ligdata){
                    if ($fline =~ /^HETATM/ && substr ($fline, 77,1) ne "H" && substr($fline,17,3) eq "$het" && substr ($fline,21,1) eq $chain){
                        chomp $fline;
                    
                        my $dist = distanceAB($line, $fline);
                        
                        if ($dist < 4.0000001){
                            push (@filepocket_lines, $line);
                            last;
                        }
                    }
                }
            }
        }

        my @finalpocket = do { my %seen; grep { !$seen{$_}++ } @filepocket_lines };
        if (scalar @finalpocket > 5){
            open OUT, ">", $out;
            foreach (@finalpocket){
                chomp $_;
                print OUT "$_\n";
            }
            print OUT "TER\nEND";
            close OUT;
        }
    
    }
}
sub distanceAB
{
    my $distance = 0;
    
    #PDB file format
    my $Ax = substr($_[0], 30, 8) + 0;
    my $Ay = substr($_[0], 38, 8) + 0;
    my $Az = substr($_[0], 46, 8) + 0;

   #PDB file format
   my $Bx = substr($_[1], 30, 8) + 0;
   my $By = substr($_[1], 38, 8) + 0;
   my $Bz = substr($_[1], 46, 8) + 0;

    $distance = sqrt(($Ax - $Bx)**2 + ($Ay - $By)**2 + ($Az - $Bz)**2);
    return sprintf("%4.2f", $distance);
}

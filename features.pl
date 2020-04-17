#!/usr/bin/perl
use strict;

my ($PDBname, $target) = @ARGV;

my $file = $PDBname."_prox4_pock.pdb";
my $pdbfile = "$PDBname.pdb";

my $pockname = substr ($file, 0,-4);

my $ext = "_";

open RPT, ">", "Features_".$target."_".$PDBname.".rpt";
print RPT "$PDBname $pockname\n";

my $seq_ref = get_sql($file);
my @seq = @$seq_ref;

##Step 1: Get hydrophobic_kyte value
my @hydro;
print RPT "Seq: ";
for my $seq (@seq){
    print RPT $seq;
    my ($hydropathy) = get_hydropathy_values($seq);
    push (@hydro, $hydropathy);
}
print RPT "\n";
my ($final_value_r, $finalvalue) = hydrophobicity_kyte(\@hydro);
print RPT "\n";

#Step 2: Get values from RADI software
my $RADI_path = "~/bin.RADI.4.0.1.macosx/";
my $in = "Input_$pockname.txt";
my $out = "Output_$pockname.txt";

open I, ">", $in;
print I "PDB\n$file\nC\n";
close I;

`$RADI_path < $in >$out`;
#Get values for six properties from RADI's output

my ($surf_hull,$vol_hull) = volume_and_surface_hull($pockname);
my ($small_size) = smallest_size($pockname);
my ($inertia_3) = inertia($pockname);


#Step 3: Residues Properties 
my $num_res = @seq; #C_RESIDUE
my ($p_aro_res, $p_ali_res, $p_negative_res, $p_positive_res) = p_res(\@seq, $num_res);

#Step 4: Milletti Paper properties and other atom properties
my ($p_Ccoo, $p_N_atom, $p_Ooh) = milletti($file, $num_res);

#Step 5: Naccess
my $asa_h = naccess("~/PocketDruggability/data/$target/apo/$pdbfile", $file);

my %Properties = (
    "VOLUME_HULL" => $vol_hull,
    "INERTIA_3" => $inertia_3,
    "SMALLEST_SIZE" => $small_size,
    "hydrophobic_kyte" => $finalvalue,
    "C_RESIDUE" => $num_res,
    "p_Ooh" => $p_Ooh,
    "p_N_atom" => $p_N_atom,
    "SURFACE_HULL" => $surf_hull,
    "p_aliphatic_residue" => $p_ali_res,
    "p_negative_residue" => $p_negative_res,
    "p_Ccoo" => $p_Ccoo,
    "hydrophobicity_pocket" => $asa_h,
    "p_aromatic_residue" => $p_aro_res,
);

print RPT "PDBid ";
print RPT "$_ " for (sort keys %Properties);
print RPT "\n";
print RPT "$PDBname ";
printf RPT ("%.3f ","$Properties{$_}") for (sort keys %Properties);
print RPT"\n";
close RPT;

`rm Input* Output* *asa* *log*`;

sub naccess{
    my ($f,$pock) = @_;
    my $cmd = "~/Naccess/naccess -a $f";
    print `pwd`;
    print $cmd, "\n";
    system($cmd);
    
    my @pock = `awk \'\$1==\"ATOM\" {print substr(\$0, 1, 55)}\' $pock`;
    chomp @pock;
    #foreach (@pock){print $_,"\n\n";}
    
    my $name = substr($pock, 0, -15);
    my $asafile = "$name.asa";
    
    my $asa;
    my $hydrophobic_asa = 0;
    open PDB, "<", $asafile;
    while (my $line = <PDB>){
        
        if ($line =~ /^ATOM/ && (substr($line, 13, 1) eq "C" || substr($line, 13, 1) eq "S")){
            chomp $line;
            my $var = substr($line, 0, 55);
            if ( grep { $_ eq $var} @pock ){
                
                substr($line, 54, 8) =~ /(\S+)/;
                $asa = $1;
                
                $hydrophobic_asa = $hydrophobic_asa + $asa;
            }
        }
    }
    return ($hydrophobic_asa);
}
sub p_res{
    my @sql = @{$_[0]};
    my $num_res = $_[1];
    my $aromatic_res = 0;
    my $aliphatic = 0;
    my $negative = 0;
    
    foreach (@sql){
        if ($_ eq "F" || $_ eq "Y" || $_ eq "H" || $_ eq "W"){
            $aromatic_res++;
        }
        if ($_ eq "I" || $_ eq "L" || $_ eq "V"){
            $aliphatic++;
        }
        if ($_ eq "D" || $_ eq "E" ){
            $negative++;
        }
    }
    my $p_aromatic_res = ($aromatic_res/$num_res);
    my $p_aliphatic_res = ($aliphatic/$num_res);
    my $p_negative = ($negative/$num_res);

    return ($p_aromatic_res, $p_aliphatic_res, $p_negative);
}
sub milletti{
    my ($pdbfile, $nres) = @_;
    my $num_atoms = `grep ^ATOM $pdbfile |wc -l`;#C_ATOM
    chomp $num_atoms;
    
    my $atomname;
    
    my $Ccoo = 0;
    my $Ooh = 0;
    my $N = 0;

    open PDB, "<", $pdbfile;
    while (my $line = <PDB>){
        if ($line =~ /^ATOM/){
            chomp $line;
            
            substr($line, 12, 4) =~ /(\S+)/;
            $atomname = $1;
            
            ## Milletti paper
            if ($atomname eq "N" || (substr($line, 17, 3) eq "GLN" && $atomname eq "NE2") || (substr($line, 17, 3) eq "ASN" && $atomname eq "ND1")){
                $N++;
            }
            if ((substr($line, 17, 3) eq "ASP" && $atomname eq "CB" ) || (substr($line, 17, 3) eq "GLU" && $atomname eq "CG")){
                $Ccoo++;
            }
            if ((substr($line, 17, 3) eq "SER" && $atomname eq "CB" ) || (substr($line, 17, 3) eq "THR" && $atomname eq "CB")){
                $Ooh++;
            }
    	}
    }
    close PDB;
    
    my $p_Ccoo = ($Ccoo/$num_atoms);
    my $p_N_atom = ($N/$num_atoms);
    my $p_Ooh = ($Ooh/$num_atoms);
    
    return ($p_Ccoo, $p_N_atom, $p_Ooh);
}
sub volume_and_surface_hull{
    my ($pock) = @_;
    my $output = "Output_$pock.txt";
    my $output_lines = `grep "VOLUME" $output |grep -v "VOLUME =        0.000000"`;
    my @output_segments = split (/\:|,|\n/,$output_lines);
    my $surface_hull_line = @output_segments[1];
    my $volume_hull_line = @output_segments[$#output_segments];
    my ($rh, $buf, $surface_hull) = split (/[=\s+]+/,$surface_hull_line);
    my ($rc, $bf, $volume_hull) = split (/[=\s+]+/,$volume_hull_line);
    return($surface_hull,$volume_hull);
}
sub smallest_size{
    my ($pock) = @_;
    my $output = "Output_$pock.txt";
    my $output_lines = `grep "SMALLEST" $output |grep -v "SMALLEST SIZE :        0.000000 ;"`;
    my @output_segments = split (/\;|\n/,$output_lines);
    my $smallest_size_line = @output_segments[0];
    my ($rh, $buf,$sz, $smallest_size) = split (/[:\s+]+/,$smallest_size_line);
    return($smallest_size);
}
sub inertia{
   my ($pock) = @_;
   my $output = "Output_$pock.txt";
   my $output_lines = `grep "EIGENVALUES OF THE INERTIA    MATRIX:" $output |grep -v "0.000000        0.000000        0.000000"`;
   my @output_segments = split (/\:/, $output_lines);
   my ($space,$inertia_1, $inertia_2, $inertia_3) = split (/\s+/,$output_segments[1]);
   return($inertia_3);
}
sub get_sql{
my ($pdb) = @_;
my %aa=qw(
	ALA 	A
	CYS 	C
	ASP 	D
	GLU 	E 
	PHE 	F
	GLY	G
	HIS	H
	ILE	I
	LYS	K
	LEU	L
	MET	M
	ASN	N
	PRO	P
	GLN	Q
	ARG	R
	SER	S
	THR	T
	VAL	V
	TRP	W
	TYR	Y);
open IN, "<", $pdb;
my $resnu;
my $resname;
my %seen;
my @sql;
while (my $line = <IN>){
    if ($line =~ /^ATOM/){
	substr($line, 22, 4) =~ /(\S+)/;
        $resnu = $1;
        $resname = substr($line, 17, 3);	
        push (@sql, $aa{$resname}) if ! $seen{$resnu}++;
    }
}
close IN;
return (\@sql);
}
sub hydrophobicity_kyte{
    my @value_array = @{$_[0]};
    my $size = @value_array;
    #print $size, "\n";

    my $sum = 0;
    for ( @value_array ) {
    	$sum += $_;
    }
    #print $sum, "\n";
    my $avg = $sum/$size;
    #print $avg, "\n";
    my $rounded = sprintf "%.2f", $avg;
    #return ($rounded);
    return ($rounded,$avg);
}

sub get_hydropathy_values{
my ($input_aa)= @_;

my %hydropathy = (
    "I" => 4.5,
    "V" => 4.2,
    "L" => 3.8,
    "F" => 2.8,
    "C" => 2.5,
    "M" => 1.9,
    "A" => 1.8,
    "G" => -0.4,
    "T" => -0.7,
    "W" => -0.9,
    "S" => -0.8,
    "Y" => -1.3,
    "P" => -1.6,
    "H" => -3.2,
    "E" => -3.5,
    "Q" => -3.5,
    "D" => -3.5,
    "N" => -3.5,
    "K" => -3.9,
    "R" => -4.5,
);

my $value = $hydropathy{$input_aa};
return ($value);
}

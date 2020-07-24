#!/usr/bin/perl

# dihenmr.pl       June 29, 2017

# Roland Dunbrack
# Fox Chase Cancer Center
# 333 Cottman Avenue
# Philadelphia PA 19111

# Roland.Dunbrack@fccc.edu

# usage 
# dihenmr.pl filename [chain]
#
# dihenmr.pl 2hla.pdb
#     will print out all dihedrals from file pdb2hla.ent
# dihenmr.pl 2hla.pdb A
#     will print out dihedrals only from chain A

# res resno chain   phi      psi     omega    chi1    chi2     chi3     chi4     model
# GLY      1 A    999.00   155.90   999.00   999.00   999.00   999.00   999.00     1
# SER      2 A   -153.26  -172.14  -178.86    43.76   999.00   999.00   999.00     1
# HIS      3 A   -119.53   158.24  -177.09   -66.08   -84.84  -176.47   999.00     1
# SER      4 A   -131.70   146.26  -173.82   -55.49   999.00   999.00   999.00     1
# MET      5 A   -120.93   133.32   167.75  -164.81   176.03    65.75   999.00     1
# ARG      6 A   -146.79   139.82  -172.87  -126.80    93.81    62.32    99.93     1
# TYR      7 A   -121.45   142.38   177.50   -85.23    95.87  -177.02    -0.19     1
# PHE      8 A   -125.64   129.46  -176.42   -69.90    87.68  -178.18    -0.70     1
# TYR      9 A   -135.75   138.87  -179.16   -65.28   167.47   178.19    -0.44     1
# THR     10 A   -129.48   139.96  -179.56   -52.99   999.00   999.00   999.00     1

# Dihedrals that are missing because of missing or undefined 
# atoms are listed as 999.00.

use POSIX;
our %field=();

# convert one-letter code to three-letter code
our %one23= (
	 A => "ALA",	 C => "CYS",	 D => "ASP",	 E => "GLU",
	 F => "PHE",	 G => "GLY",	 H => "HIS",	 I => "ILE",
	 K => "LYS",	 L => "LEU",	 M => "MET",	 N => "ASN",
	 P => "PRO",	 Q => "GLN",	 R => "ARG",	 S => "SER",
	 T => "THR",	 V => "VAL",	 W => "TRP",	 Y => "TYR", X => "UNK"
    );

# convert three-letter code to one letter code, including modified amino acids
%three21 = (
	    ALA => "A", CYS => "C", ASP => "D", GLU => "E", PHE => "F", GLY => "G",
	    HIS => "H", ILE => "I", LYS => "K", LEU => "L", MET => "M", ASN => "N",
	    PRO => "P", GLN => "Q", ARG => "R", SER => "S", THR => "T", VAL => "V", 
	    TRP => "W", TYR => "Y", ASX => "N", GLX => "Q", UNK => "X", INI => "K",
	    AAR => "R", ACE => "X", ACY => "G", AEI => "T", AGM => "R", ASQ => "D",
	    AYA => "A", BHD => "D", CAS => "C", CAY => "C", CEA => "C", CGU => "E",
	    CME => "C", CMT => "C", CSB => "C", CSD => "C", CSE => "C", CSO => "C",
	    CSP => "C", CSS => "C", CSW => "C", CSX => "C", CXM => "M", CYG => "C",
	    CYM => "C", DOH => "D", EHP => "F", FME => "M", FTR => "W", GL3 => "G",
	    H2P => "H", HIC => "H", HIP => "H", HTR => "W", HYP => "P", KCX => "K",
	    LLP => "K", LLY => "K", LYZ => "K", M3L => "K", MEN => "N", MGN => "Q",
	    MHO => "M", MHS => "H", MIS => "S", MLY => "K", MLZ => "K", MSE => "M",
	    NEP => "H", NPH => "C", OCS => "C", OCY => "C", OMT => "M", OPR => "R",
	    PAQ => "Y", PCA => "Q", PHD => "D", PRS => "P", PTH => "Y", PYX => "C",
	    SEP => "S", SMC => "C", SME => "M", SNC => "C", SNN => "D", SVA => "S",
	    TPO => "T", TPQ => "Y", TRF => "W", TRN => "W", TRO => "W", TYI => "Y",
	    TYN => "Y", TYQ => "Y", TYS => "Y", TYY => "Y", YOF => "Y", FOR => "X"
	   );
$three21{"5HP"}="Q";


#constants
$pi=acos(-1);
$r2d=180.0/$pi;
$d2r=$pi/180.0;

# process arguments
$file=$ARGV[0];
if ($file =~ /\.gz$/) {open FILE,"gunzip -c $file |";}
else {open FILE,"<$file";}

$file =~ s/\.gz//;
$file =~ s/\.pdb//;
$file =~ s/\.cif//;

if ($#ARGV==1) {$chain=$ARGV[1];}
else {$chain=".";}
$nres=-1;


@file=split(/\//,$file);

$filename=$file[$#file];

# read coordinates from file using &pdbread from coor.pm
# put hash from pdbread into @pdbres


#loop_
#_atom_site.group_PDB 
#_atom_site.id 
#_atom_site.type_symbol 
#_atom_site.label_atom_id 
#_atom_site.label_alt_id 
#_atom_site.label_comp_id 
#_atom_site.label_asym_id 
#_atom_site.label_entity_id 
#_atom_site.label_seq_id 
#_atom_site.pdbx_PDB_ins_code 
#_atom_site.Cartn_x 
#_atom_site.Cartn_y 
#_atom_site.Cartn_z 
#_atom_site.occupancy 
#_atom_site.B_iso_or_equiv 
#_atom_site.pdbx_formal_charge 
#_atom_site.auth_seq_id 
#_atom_site.auth_comp_id 
#_atom_site.auth_asym_id 
#_atom_site.auth_atom_id 
#_atom_site.pdbx_PDB_model_num 
#ATOM   1    N N   . GLU A 1 2   ? 40.331 112.546 89.262  1.00 44.35 ? 2    GLU A N   1 
#ATOM   2    C CA  . GLU A 1 2   ? 40.678 111.124 89.548  1.00 45.01 ? 2    GLU A CA  1 
#ATOM   3    C C   . GLU A 1 2   ? 39.412 110.329 89.861  1.00 45.27 ? 2    GLU A C   1 
#ATOM   4    O O   . GLU A 1 2   ? 38.805 110.488 90.922  1.00 45.54 ? 2    GLU A O   1 
#ATOM   5    C CB  . GLU A 1 2   ? 41.656 111.042 90.724  1.00 44.52 ? 2    GLU A CB  1 
#ATOM   6    N N   . MET A 1 3   ? 39.028 109.473 88.919  1.00 44.97 ? 3    MET A N   1 
#ATOM   7    C CA  . MET A 1 3   ? 37.831 108.642 89.028  1.00 44.48 ? 3    MET A CA  1 

$modelnum=1;
$nfield=0;
$natom=0;
while(<FILE>) {
    if (/^MODEL/) {
	($model, $modelnum)=split;
    }
    if (/^ATOM/ || /^HETATM/) {
	%pdbcoor=();
	chomp;
	$line=$_;
	%pdbcoor=pdbread($line);
	if (!%pdbcoor) {next;}
	$pdbcoor{modelnum}=$modelnum;
	$asymid=$pdbcoor{pdb_chainid};
	$ncoor=$pdbcoor{pdb_ncoor};
	$altloc=$pdbcoor{altloc};
	if ($altloc eq ".") {$altloc = "A";}
	$atomtypeno=$pdbcoor{atomtypeno};
	$natom++;
	$pdbres{$modelnum}{$asymid}{$ncoor}{$altloc}{$atomtypeno}=$natom;
	$natoms{$modelnum}{$asymid}{$ncoor}{$altloc}++;
	$occ=$pdbcoor{occ};
	if (!defined($occ{$modelnum}{$asymid}{$ncoor}{$altloc})) {
	    $occ{$modelnum}{$asymid}{$ncoor}{$altloc}=$occ;
	}
	else {
	    if ($occ<$occ{$modelnum}{$asymid}{$ncoor}{$altloc}) { #save the lowest value of occ
		$occ{$modelnum}{$asymid}{$ncoor}{$altloc}=$occ;
	    }
	}
	%{$atoms[$natom]}=%pdbcoor;
#	print "$modelnum $asymid natom $natom ncoor $ncoor altloc $altloc $pdbcoor{pdbres}\n";
    }
}
close FILE;
    
#    $pdbcoor{atomnum}  = $array[$rfield->{'id'}]; # atom serial number
#    $pdbcoor{atomtype} = $array[$rfield->{'label_atom_id'}]; # atom type
#    $pdbcoor{element}  = $array[$rfield->{'type_symbol'}]; # element N,C,O,S,H
#    $pdbcoor{altloc}   = $array[$rfield->{'label_alt_id'}]; # alternate atom location
#    $pdbcoor{pdbres}   = $array[$rfield->{'label_comp_id'}]; # residue
#    $pdbcoor{auth_chainid} =  $array[$rfield->{'auth_asym_id'}]; # chain ID
#    $pdbcoor{pdb_chainid} =  $array[$rfield->{'label_asym_id'}]; # chain ID
#    $pdbcoor{auth_ncoor}    = $array[$rfield->{'label_seq_id'}]; # ATOM coord number 
#    $pdbcoor{auth_ncoor}    = $array[$rfield->{'auth_seq_id'}]; # ATOM coord number
#    $pdbcoor{pdb_ncoor}    = $array[$rfield->{'label_seq_id'}]; # ATOM coord number
#    $pdbcoor{icode}    = $array[$rfield->{'pdbx_PDB_ins_code'}]; # insertion code
    
#    $pdbcoor{coor}[0]  = $array[$rfield->{'Cartn_x'}]; # x coord
#    $pdbcoor{coor}[1]  = $array[$rfield->{'Cartn_y'}]; # y coord
#    $pdbcoor{coor}[2]  = $array[$rfield->{'Cartn_z'}]; # z coord
#    $pdbcoor{occ}     =  $array[$rfield->{'occupancy'}]; # occupancy
#    $pdbcoor{bfac}    =  $array[$rfield->{'B_iso_or_equiv'}]; # b-factor

#    $pdbcoor{auth_ncoor} .= $pdbcoor{icode};   # include insertion code
#    $pdbcoor{pdb_ncoor} .= $pdbcoor{icode};   # include insertion code

#    $pdbcoor{pdbres} = $one23{$three21{$pdbcoor{pdbres}}};
#    $pdbcoor{modelnum} = $array[$rfield->{'pdbx_PDB_model_num'}]; #model number


foreach $modelnum (sort {$a<=>$b} keys %pdbres) {
    foreach $asymid (sort keys %{$pdbres{$modelnum}}) {
	foreach $ncoor (sort {$a<=>$b} keys %{$pdbres{$modelnum}{$asymid}}) {
	    foreach $altloc (sort keys %{$pdbres{$modelnum}{$asymid}{$ncoor}}) {
#		print "loop $modelnum $asymid $ncoor $altloc $natoms{$modelnum}{$asymid}{$ncoor}{$altloc}\n";
		if ($natoms{$modelnum}{$asymid}{$ncoor}{$altloc}<1) {next;}

		undef(%A0); undef(%C0); undef(%N1); undef(%A1); undef(%C1); undef(%O1); undef(%N2); undef(%A2);
		undef(%B1); undef(%G1); undef(%D1); undef(%E1); undef(%Z1); undef(%H1);
		
#		$natom= $pdbres{$modelnum}{$asymid}{$ncoor}{$altloc}{0};  print "natom 0 $natom\n"; if ($natom) {%N1=%{$atoms[$natom]};}
#		$natom= $pdbres{$modelnum}{$asymid}{$ncoor}{$altloc}{1};  print "natom 1 $natom\n"; if ($natom) {%A1=%{$atoms[$natom]};}
#		$natom= $pdbres{$modelnum}{$asymid}{$ncoor}{$altloc}{2};  print "natom 2 $natom\n"; if ($natom) {%C1=%{$atoms[$natom]};}
#		$natom= $pdbres{$modelnum}{$asymid}{$ncoor}{$altloc}{3};  print "natom 3 $natom\n"; if ($natom) {%O1=%{$atoms[$natom]};}

		$natom= $pdbres{$modelnum}{$asymid}{$ncoor}{$altloc}{0};  if ($natom) {%N1=%{$atoms[$natom]};}
		$natom= $pdbres{$modelnum}{$asymid}{$ncoor}{$altloc}{1};  if ($natom) {%A1=%{$atoms[$natom]};}
		$natom= $pdbres{$modelnum}{$asymid}{$ncoor}{$altloc}{2};  if ($natom) {%C1=%{$atoms[$natom]};}
		$natom= $pdbres{$modelnum}{$asymid}{$ncoor}{$altloc}{3};  if ($natom) {%O1=%{$atoms[$natom]};}
		if (defined($pdbres{$modelnum}{$asymid}{$ncoor-1})) {
#		    $natom= $pdbres{$modelnum}{$asymid}{$ncoor-1}{$altloc}{1}; print "natom -1 $natom\n"; if ($natom) {%A0=%{$atoms[$natom]};}
#		    $natom= $pdbres{$modelnum}{$asymid}{$ncoor-1}{$altloc}{2}; print "natom -2 $natom\n"; if ($natom) {%C0=%{$atoms[$natom]};} 
		    $natom= $pdbres{$modelnum}{$asymid}{$ncoor-1}{$altloc}{1}; if ($natom) {%A0=%{$atoms[$natom]};}
		    $natom= $pdbres{$modelnum}{$asymid}{$ncoor-1}{$altloc}{2}; if ($natom) {%C0=%{$atoms[$natom]};} 
		}
		if (defined($pdbres{$modelnum}{$asymid}{$ncoor+1})) {
#		    $natom= $pdbres{$modelnum}{$asymid}{$ncoor+1}{$altloc}{0}; print "natom +0 $natom\n"; if ($natom) {%N2=%{$atoms[$natom]};} 
#		    $natom= $pdbres{$modelnum}{$asymid}{$ncoor+1}{$altloc}{1}; print "natom +1 $natom\n"; if ($natom) {%A2=%{$atoms[$natom]};}
		    $natom= $pdbres{$modelnum}{$asymid}{$ncoor+1}{$altloc}{0}; if ($natom) {%N2=%{$atoms[$natom]};} 
		    $natom= $pdbres{$modelnum}{$asymid}{$ncoor+1}{$altloc}{1}; if ($natom) {%A2=%{$atoms[$natom]};}
		}
		$natom= $pdbres{$modelnum}{$asymid}{$ncoor}{$altloc}{4}; if ($natom) {%B1=%{$atoms[$natom]};} 
		$natom= $pdbres{$modelnum}{$asymid}{$ncoor}{$altloc}{5}; if ($natom) {%G1=%{$atoms[$natom]};} 
		$natom= $pdbres{$modelnum}{$asymid}{$ncoor}{$altloc}{6}; if ($natom) {%D1=%{$atoms[$natom]};} 
		$natom= $pdbres{$modelnum}{$asymid}{$ncoor}{$altloc}{7}; if ($natom) {%E1=%{$atoms[$natom]};} 
		$natom= $pdbres{$modelnum}{$asymid}{$ncoor}{$altloc}{8}; if ($natom) {%Z1=%{$atoms[$natom]};} 
		$natom= $pdbres{$modelnum}{$asymid}{$ncoor}{$altloc}{9}; if ($natom) {%H1=%{$atoms[$natom]};} 
		
		if ($altloc ne "A") {
		    if (!%N1) {$natom= $pdbres{$modelnum}{$asymid}{$ncoor}{"A"}{0}; if ($natom) {%N1=%{$atoms[$natom]};}}
		    if (!%A1) {$natom= $pdbres{$modelnum}{$asymid}{$ncoor}{"A"}{1}; if ($natom) {%A1=%{$atoms[$natom]};}}
		    if (!%C1) {$natom= $pdbres{$modelnum}{$asymid}{$ncoor}{"A"}{2}; if ($natom) {%C1=%{$atoms[$natom]};}}
		    if (!%O1) {$natom= $pdbres{$modelnum}{$asymid}{$ncoor}{"A"}{3}; if ($natom) {%O1=%{$atoms[$natom]};}}
		    if (!%A0) {$natom= $pdbres{$modelnum}{$asymid}{$ncoor-1}{"A"}{1}; if ($natom) {%A0=%{$atoms[$natom]};}}
		    if (!%C0) {$natom= $pdbres{$modelnum}{$asymid}{$ncoor-1}{"A"}{2}; if ($natom) {%C0=%{$atoms[$natom]};} }
		    if (!%N2) {$natom= $pdbres{$modelnum}{$asymid}{$ncoor+1}{"A"}{0}; if ($natom) {%N2=%{$atoms[$natom]};} }
		    if (!%A2) {$natom= $pdbres{$modelnum}{$asymid}{$ncoor+1}{"A"}{1}; if ($natom) {%A2=%{$atoms[$natom]};}}
		    if (!%B1) {$natom= $pdbres{$modelnum}{$asymid}{$ncoor}{"A"}{4}; if ($natom) {%B1=%{$atoms[$natom]};} }
		    if (!%G1) {$natom= $pdbres{$modelnum}{$asymid}{$ncoor}{"A"}{5}; if ($natom) {%G1=%{$atoms[$natom]};} }
		    if (!%D1) {$natom= $pdbres{$modelnum}{$asymid}{$ncoor}{"A"}{6}; if ($natom) {%D1=%{$atoms[$natom]};} }
		    if (!%E1) {$natom= $pdbres{$modelnum}{$asymid}{$ncoor}{"A"}{7}; if ($natom) {%E1=%{$atoms[$natom]};} }
		    if (!%Z1) {$natom= $pdbres{$modelnum}{$asymid}{$ncoor}{"A"}{8}; if ($natom) {%Z1=%{$atoms[$natom]};} }
		    if (!%H1) {$natom= $pdbres{$modelnum}{$asymid}{$ncoor}{"A"}{9}; if ($natom) {%H1=%{$atoms[$natom]};} }
		}

		$phi=999.0;
		$psi=999.0;
		$omega=999.0;
		$chi1=999.0;
		$chi2=999.0;
		$chi3=999.0;
		$chi4=999.0;
		$chi5=999.0;
	    
# OMEGA		
		if (%A0 && %C0 && %N1 && %A1) {
		    $omega=dihedral( @{$A0{coor}},@{$C0{coor}}, @{$N1{coor}}, @{$A1{coor}} );
		}
		
# PHI
		if (%C0 && %N1 && %A1 && %C1) {
		    $phi=dihedral( @{$C0{coor}},@{$N1{coor}}, @{$A1{coor}}, @{$C1{coor}} );
		}
		
# PSI
		if (%N1 && %A1 && %C1 && %N2) {
		    $psi=dihedral( @{$N1{coor}},@{$A1{coor}}, @{$C1{coor}}, @{$N2{coor}} );
		}
		
# CHI1
		if (%N1 && %A1 && %B1 && %G1) {
		    $chi1=dihedral( @{$N1{coor}},@{$A1{coor}}, @{$B1{coor}}, @{$G1{coor}} );
		}
		
# CHI2
		if (%A1 && %B1 && %G1 && %D1) {
		    $chi2=dihedral( @{$A1{coor}},@{$B1{coor}}, @{$G1{coor}}, @{$D1{coor}} );
		}
		
# CHI3
		if (%B1 && %G1 && %D1 && %E1) {
		    $chi3=dihedral( @{$B1{coor}},@{$G1{coor}}, @{$D1{coor}}, @{$E1{coor}} );
		}
		
# CHI4:
		if (%G1 && %D1 && %E1 && %Z1) {
		    $chi4=dihedral( @{$G1{coor}},@{$D1{coor}}, @{$E1{coor}}, @{$Z1{coor}} );
		}
		
# CHI5:
		if (%D1 && %E1 && %Z1 && %H1) {
		    $chi5=dihedral( @{$D1{coor}},@{$E1{coor}}, @{$Z1{coor}}, @{$H1{coor}} );
		}
		

		$altloc2=$altloc;
		if ($natoms{$modelnum}{$asymid}{$ncoor}{"B"}<1) {$altloc2=".";}
		printf("%3s %s %5.2f %6s %4s %9.2f%9.2f%9.2f%9.2f%9.2f%9.2f%9.2f%9.2f %3d  %6s %4s  %s\n",
		       $A1{pdbres},$altloc2, $occ{$modelnum}{$asymid}{$ncoor}{$altloc},
		       $A1{auth_ncoor},
		       $A1{auth_chainid},
		       $phi,$psi,$omega,$chi1,$chi2,$chi3,$chi4, $chi5, $modelnum, 
		       $A1{pdb_ncoor},
		       $A1{pdb_chainid},$filename
		    );
	    }
	}
    }   
}

# constants

# convert one-letter code to three-letter code
%one23= (
	 A => "ALA",	 C => "CYS",	 D => "ASP",	 E => "GLU",
	 F => "PHE",	 G => "GLY",	 H => "HIS",	 I => "ILE",
	 K => "LYS",	 L => "LEU",	 M => "MET",	 N => "ASN",
	 P => "PRO",	 Q => "GLN",	 R => "ARG",	 S => "SER",
	 T => "THR",	 V => "VAL",	 W => "TRP",	 Y => "TYR", X => "UNK");

# convert three-letter code to one letter code, including modified amino acids
%three21 = (
	    ALA => "A", CYS => "C", ASP => "D", GLU => "E", PHE => "F", GLY => "G",
	    HIS => "H", ILE => "I", LYS => "K", LEU => "L", MET => "M", ASN => "N",
	    PRO => "P", GLN => "Q", ARG => "R", SER => "S", THR => "T", VAL => "V", 
	    TRP => "W", TYR => "Y", ASX => "N", GLX => "Q", UNK => "X", INI => "K",
	    AAR => "R", ACE => "X", ACY => "G", AEI => "T", AGM => "R", ASQ => "D",
	    AYA => "A", BHD => "D", CAS => "C", CAY => "C", CEA => "C", CGU => "E",
	    CME => "C", CMT => "C", CSB => "C", CSD => "C", CSE => "C", CSO => "C",
	    CSP => "C", CSS => "C", CSW => "C", CSX => "C", CXM => "M", CYG => "C",
	    CYM => "C", DOH => "D", EHP => "F", FME => "M", FTR => "W", GL3 => "G",
	    H2P => "H", HIC => "H", HIP => "H", HTR => "W", HYP => "P", KCX => "K",
	    LLP => "K", LLY => "K", LYZ => "K", M3L => "K", MEN => "N", MGN => "Q",
	    MHO => "M", MHS => "H", MIS => "S", MLY => "K", MLZ => "K", MSE => "M",
	    NEP => "H", NPH => "C", OCS => "C", OCY => "C", OMT => "M", OPR => "R",
	    PAQ => "Y", PCA => "Q", PHD => "D", PRS => "P", PTH => "Y", PYX => "C",
	    SEP => "S", SMC => "C", SME => "M", SNC => "C", SNN => "D", SVA => "S",
	    TPO => "T", TPQ => "Y", TRF => "W", TRN => "W", TRO => "W", TYI => "Y",
	    TYN => "Y", TYQ => "Y", TYS => "Y", TYY => "Y", YOF => "Y", FOR => "X"
	   );
$three21{"5HP"}="Q";

%codon2aa = (
    ATG => "M",                                                             # Met
    GCT => "A", GCC => "A", GCA => "A", GCG => "A",                         # Ala
    TGT => "C", TGC => "C",                                                 # Cys
    GAT => "D", GAC => "D",                                                 # Asp
    GAA => "E", GAG => "E",                                                 # Gln
    TTT => "F", TTC => "F",                                                 # Phe
    GGT => "G", GGC => "G", GGA => "G", GGG => "G",                         # Gly
    CAT => "H", CAC => "H",                                                 # His
    ATT => "I", ATC => "I", ATA => "I",                                     # Ile
    AAA => "K", AAG => "K",                                                 # Lys
    TTA => "L", TTG => "L", CTT => "L", CTC => "L", CTA => "L", CTG => "L", # Leu
    AAT => "N", AAC => "N",                                                 # Asn
    CCT => "P", CCC => "P", CCA => "P", CCG => "P",                         # Pro
    CAA => "Q", CAG => "Q",                                                 # Gln
    TCT => "S", TCC => "S", TCA => "S", TCG => "S", AGT => "S", AGC => "S", # Ser
    ACT => "T", ACC => "T", ACA => "T", ACG => "T",                         # Thr
    CGT => "R", CGC => "R", CGA => "R", CGG => "R", AGA => "R", AGG => "R", # Arg
    GTT => "V", GTC => "V", GTA => "V", GTG => "V",                         # Val
    TGG => "W",                                                             # Trp
    TAT => "Y", TAC => "Y",                                                 # Tyr
    TGA => "*", TAA => "*", TAG => "*"                                      # Ter/stop
    );

sub getyymmdd {
    my $month=(localtime)[4]+1;
    my $year=(localtime)[5];
    my $day=(localtime)[3];

    if ($day<10) {$day="0$day";}
    else {$day="$day";}
    if ($month<10) {$month="0$month";}
    else {$month="$month";}
    if (length($year)>2) {$year=substr($year,length($year)-2,2);}
    my $date=$year . $month . $day;
    return $date;
}

sub getmmddyy {
    my @months=(January, February, March, April, May, June, July, August, September, October, November, December);
    my $month=$months[(localtime)[4]];
    my $year=(localtime)[5];
    my $day=(localtime)[3];
    $year=$year+1900;
    my $date="$month $day\, $year";
    return $date;
}

# returns rounded integer; e.g., rint(1.5) is 2; rint(1.4) is 1
sub rint {
  $x=$_[0];
  my $ix=int($x);
  if ($x-$ix>=0.5) {$ix+=1;}
  if ($x-$ix<=-0.5) {$ix-=1;}
  return $ix;
}

%month2num = (
	      Jan => "01", Feb => "02", Mar => "03", Apr => "04", May => "05", Jun => "06",
	      Jul => "07", Aug => "08", Sep => "09", Oct => "10", Nov => "11", Dec => "12");


# atom types for each amino acid by one-letter code
%iupac=();
%{$iupac{A}}=(N=>1,CA=>1,C=>1,O=>1,CB=>1);
%{$iupac{C}}=(N=>1,CA=>1,C=>1,O=>1,CB=>1,SG =>1);
%{$iupac{D}}=(N=>1,CA=>1,C=>1,O=>1,CB=>1,CG =>1,OD1=>1,OD2=>1);

%{$iupac{E}}=(N=>1,CA=>1,C=>1,O=>1,CB=>1,CG =>1,CD =>1,OE1=>1,OE2=>1);
%{$iupac{Q}}=(N=>1,CA=>1,C=>1,O=>1,CB=>1,CG =>1,CD =>1,OE1=>1,NE2=>1);

%{$iupac{F}}=(N=>1,CA=>1,C=>1,O=>1,CB=>1,CG =>1,CD1=>1,CE1=>1,CD2=>1,CE2=>1,CZ=>1);
%{$iupac{G}}=(N=>1,CA=>1,C=>1,O=>1);
%{$iupac{H}}=(N=>1,CA=>1,C=>1,O=>1,CB=>1,CG =>1,ND1=>1,CE1=>1,CD1=>1,NE2=>1);
%{$iupac{I}}=(N=>1,CA=>1,C=>1,O=>1,CB=>1,CG1=>1,CD1=>1,CG2=>1);
%{$iupac{K}}=(N=>1,CA=>1,C=>1,O=>1,CB=>1,CG =>1,CD =>1,CE =>1,NZ =>1);
%{$iupac{L}}=(N=>1,CA=>1,C=>1,O=>1,CB=>1,CG =>1,CD1=>1,CD2=>1);
%{$iupac{M}}=(N=>1,CA=>1,C=>1,O=>1,CB=>1,CG =>1,SD =>1,CE =>1);
%{$iupac{N}}=(N=>1,CA=>1,C=>1,O=>1,CB=>1,CG =>1,OD1=>1,ND2=>1);
%{$iupac{P}}=(N=>1,CA=>1,C=>1,O=>1,CB=>1,CG =>1,CD =>1);
%{$iupac{R}}=(N=>1,CA=>1,C=>1,O=>1,CB=>1,CG =>1,CD =>1,NE =>1,CZ =>1,NH1=>1,NH2=>1);
%{$iupac{S}}=(N=>1,CA=>1,C=>1,O=>1,CB=>1,OG =>1);
%{$iupac{T}}=(N=>1,CA=>1,C=>1,O=>1,CB=>1,OG1=>1,CG2=>1);
%{$iupac{V}}=(N=>1,CA=>1,C=>1,O=>1,CB=>1,CG1=>1,CG2=>1);
%{$iupac{W}}=(N=>1,CA=>1,C=>1,O=>1,CB=>1,CG =>1,CD1=>1,CD2=>1,NE1=>1,CE2=>1,CE3=>1,CZ2=>1,CZ3=>1,CH2=>1);
%{$iupac{Y}}=(N=>1,CA=>1,C=>1,O=>1,CB=>1,CG =>1,CD1=>1,CE1=>1,CD2=>1,CE2=>1,CZ =>1,OH =>1);





=pod
    $pdbcoor{atomnum}  = $array[$rfield->{'id'}]; # atom serial number
    $pdbcoor{atomtype} = $array[$rfield->{'label_atom_id'}]; # atom type
    $pdbcoor{element}  = $array[$rfield->{'type_symbol'}]; # element N,C,O,S,H
    $pdbcoor{altloc}   = $array[$rfield->{'label_alt_id'}]; # alternate atom location
    $pdbcoor{pdbres}   = $array[$rfield->{'label_comp_id'}]; # residue
    if (!defined($three21{$pdbcoor{pdbres}})) {%pdbcoor=(); return %pdbcoor;}
    $pdbcoor{auth_chainid} =  $array[$rfield->{'auth_asym_id'}]; # chain ID
    $pdbcoor{pdb_chainid} =  $array[$rfield->{'label_asym_id'}]; # chain ID
    $pdbcoor{auth_ncoor}    = $array[$rfield->{'label_seq_id'}]; # ATOM coord number 
    if (defined($rfield->{'auth_seq_id'})) {
	$pdbcoor{auth_ncoor}    = $array[$rfield->{'auth_seq_id'}]; # ATOM coord number
    }
    $pdbcoor{pdb_ncoor}    = $array[$rfield->{'label_seq_id'}]; # ATOM coord number
    $pdbcoor{icode}    = $array[$rfield->{'pdbx_PDB_ins_code'}]; # insertion code
    $pdbcoor{coor}[0]  = $array[$rfield->{'Cartn_x'}]; # x coord
    $pdbcoor{coor}[1]  = $array[$rfield->{'Cartn_y'}]; # y coord
    $pdbcoor{coor}[2]  = $array[$rfield->{'Cartn_z'}]; # z coord
    $pdbcoor{occ}     =  $array[$rfield->{'occupancy'}]; # occupancy
    $pdbcoor{bfac}    =  $array[$rfield->{'B_iso_or_equiv'}]; # b-factor
    if ($pdbcoor{icode} ne "?") {
	$pdbcoor{auth_ncoor} .= $pdbcoor{icode};   # include insertion code in author numbering
    }
    # Translate any modified residues
    $pdbcoor{pdbres} = $one23{$three21{$pdbcoor{pdbres}}};
    $pdbcoor{modelnum} = $array[$rfield->{'pdbx_PDB_model_num'}]; #model number    
=cut

# subroutine pdbread - read lines from PDB format files and return hash
# with data from line
sub pdbread {
  my %pdbcoor=();
  my $line=$_[0];
  chomp $line;
  $pdbcoor{atomnum} =  unpack("x6 a5",  $line); # atom serial number
  $pdbcoor{atomtype} = unpack("x12 a4", $line); # atom type
  $pdbcoor{element} = unpack("x13 a1",$line); # element N,C,O,S,H
  $pdbcoor{altloc} =   unpack("x16 a1", $line); # alternate atom location
  if ($pdbcoor{altloc} eq " ") {$pdbcoor{altloc}=".";}
  $pdbcoor{pdbres} =   unpack("x17 a3", $line); # residue
  if (!defined($three21{$pdbcoor{pdbres}})) {%pdbcoor=(); return %pdbcoor;}
  $pdbcoor{auth_chainid} =  unpack("x21 a1", $line); # chain ID
  $pdbcoor{pdb_chainid} =  unpack("x21 a1", $line); # chain ID
  $pdbcoor{auth_ncoor} =    unpack("x22 a4", $line); # ATOM coord number
  $pdbcoor{pdb_ncoor} =    unpack("x22 a4", $line); # ATOM coord number
  $pdbcoor{icode} =    unpack("x26 a1", $line); # insertion code
#  if ($pdbcoor{icode} eq " ") {$pdbcoor{icode} = ".";}
  $pdbcoor{coor}[0] =  unpack("x30 a8", $line); # x coord
  $pdbcoor{coor}[1] =  unpack("x38 a8", $line); # y coord
  $pdbcoor{coor}[2] =  unpack("x46 a8", $line); # z coord
  $pdbcoor{occ}=1;
  $pdbcoor{bfac}=0.0;
  if (length($line)>55) {
    $pdbcoor{occ}     =  unpack("x54 a6", $line); # occupancy
  }
  if (length($line)>61) {
    $pdbcoor{bfac}    =  unpack("x60 a6", $line); # b-factor
  }
  $pdbcoor{bfac} =~ s/ //g;
  $pdbcoor{occ} =~ s/ //g;
  $pdbcoor{auth_ncoor} .= $pdbcoor{icode};   # include insertion code
  $pdbcoor{auth_ncoor} =~ s/ //g;   # Strip spaces
  $pdbcoor{pdb_ncoor} .= $pdbcoor{icode};   # include insertion code
  $pdbcoor{pdb_ncoor} =~ s/ //g;   # Strip spaces
  $pdbcoor{coor}[0] =~ s/ //g;
  $pdbcoor{coor}[1] =~ s/ //g;
  $pdbcoor{coor}[2] =~ s/ //g;
  
  # Translate any modified residues
  $pdbcoor{pdbres} = $one23{$three21{$pdbcoor{pdbres}}};
  $pdbcoor{atom}=unpack("x13 a1",$line); 
  $pdbcoor{atomtype}=~s/\s+//g;

  # assign atomtypeno according to following scheme:
  # N            0
  # CA           1
  # C            2
  # O            3
  # CB           4
  # CG/CG1/OG1   5
  # CD/CD1/OD1   6
  # CE/CE1/OE1   7
  # CZ/CZ1/OZ1   8
  # CH/CH1/OH1   9
  # CG2/OG2/NG2 10
  # CD2/OD2/ND2 11
  # CE2/OE2/NE2 12
  # CZ2/OZ2/NZ2 13
  # CH2/OH2/NH2 14
  # CE3/OE3/NE3 15
  # CZ3/OZ3/NZ3 16
  # OXT         17
  # H*          20

  $pdbcoor{atomtypeno}=100;
  my $atomtype=$pdbcoor{atomtype};
  if ($pdbcoor{element} eq "H") {$pdbcoor{atomtypeno}=20; return %pdbcoor;}
  if ($atomtype eq "N"   ) {$pdbcoor{atomtypeno}=0; return %pdbcoor;}
  if ($atomtype eq "CA"  ) {$pdbcoor{atomtypeno}=1; return %pdbcoor;}
  if ($atomtype eq "C"   ) {$pdbcoor{atomtypeno}=2; return %pdbcoor;}
  if ($atomtype eq "O"   ) {$pdbcoor{atomtypeno}=3; return %pdbcoor;}
  if ($atomtype eq "OXT" ) {$pdbcoor{atomtypeno}=17; return %pdbcoor;}
  if ($atomtype eq "CB"  ) {$pdbcoor{atomtypeno}=4; return %pdbcoor;}
  if ($atomtype eq "SE" && $pdbcoor{pdbres} eq "MSE" ) {
    $pdbcoor{atomtypeno}=6; return %pdbcoor;}

  if ($atomtype =~ /[COSN]G1|[COSN]G$/) {
    $pdbcoor{atomtypeno}=5; return %pdbcoor;}
  if ($atomtype =~ /[COSN]D1|[COSN]D$/) {
    $pdbcoor{atomtypeno}=6; return %pdbcoor;}
  if ($atomtype =~ /[COSN]E1|[COSN]E$/) {
    $pdbcoor{atomtypeno}=7; return %pdbcoor;}
  if ($atomtype =~ /[COSN]Z1|[COSN]Z$/) {
    $pdbcoor{atomtypeno}=8; return %pdbcoor;}
  if ($atomtype =~ /[COSN]H1|[COSN]H$/) {
    $pdbcoor{atomtypeno}=9; return %pdbcoor;}
  if ($atomtype =~ /[COSN]G2/) {
    $pdbcoor{atomtypeno}=10; return %pdbcoor;}
  if ($atomtype =~ /[COSN]D2/) {
    $pdbcoor{atomtypeno}=11; return %pdbcoor;}
  if ($atomtype =~ /[COSN]E2/) {
    $pdbcoor{atomtypeno}=12; return %pdbcoor;}
  if ($atomtype =~ /[COSN]Z2/) {
    $pdbcoor{atomtypeno}=13; return %pdbcoor;}
  if ($atomtype =~ /[COSN]H2/) {
    $pdbcoor{atomtypeno}=14; return %pdbcoor;}
  if ($atomtype =~ /[COSN]E3/) {
    $pdbcoor{atomtypeno}=15; return %pdbcoor;}
  if ($atomtype =~ /[COSN]Z3/) {
    $pdbcoor{atomtypeno}=16; return %pdbcoor;}
  return %pdbcoor;
}




# subroutine cifread - read lines from PDB format files and return hash
# with data from line

sub cifread {

    my $rfield= \%field;
    my %pdbcoor=();
    my $line=$_[0];
    chomp $line;
    my @array=split(/\s+/,$line);
    
    my $nfield=0;

#    print $rfield->{'label_seq_id'}, " ", $array[$rfield->{'label_seq_id'}],"\n";

    $pdbcoor{atomnum}  = $array[$rfield->{'id'}]; # atom serial number
    $pdbcoor{atomtype} = $array[$rfield->{'label_atom_id'}]; # atom type
    $pdbcoor{element}  = $array[$rfield->{'type_symbol'}]; # element N,C,O,S,H
    $pdbcoor{altloc}   = $array[$rfield->{'label_alt_id'}]; # alternate atom location
    $pdbcoor{pdbres}   = $array[$rfield->{'label_comp_id'}]; # residue
    if (!defined($three21{$pdbcoor{pdbres}})) {%pdbcoor=(); return %pdbcoor;}
    $pdbcoor{auth_chainid} =  $array[$rfield->{'auth_asym_id'}]; # chain ID
    $pdbcoor{pdb_chainid} =  $array[$rfield->{'label_asym_id'}]; # chain ID
    $pdbcoor{auth_ncoor}    = $array[$rfield->{'label_seq_id'}]; # ATOM coord number 
    if (defined($rfield->{'auth_seq_id'})) {
	$pdbcoor{auth_ncoor}    = $array[$rfield->{'auth_seq_id'}]; # ATOM coord number
    }
    $pdbcoor{pdb_ncoor}    = $array[$rfield->{'label_seq_id'}]; # ATOM coord number
    $pdbcoor{icode}    = $array[$rfield->{'pdbx_PDB_ins_code'}]; # insertion code
    $pdbcoor{coor}[0]  = $array[$rfield->{'Cartn_x'}]; # x coord
    $pdbcoor{coor}[1]  = $array[$rfield->{'Cartn_y'}]; # y coord
    $pdbcoor{coor}[2]  = $array[$rfield->{'Cartn_z'}]; # z coord
    $pdbcoor{occ}     =  $array[$rfield->{'occupancy'}]; # occupancy
    $pdbcoor{bfac}    =  $array[$rfield->{'B_iso_or_equiv'}]; # b-factor

    if ($pdbcoor{icode} ne "?") {
	$pdbcoor{auth_ncoor} .= $pdbcoor{icode};   # include insertion code in author numbering
    }
    
    # Translate any modified residues
    $pdbcoor{pdbres} = $one23{$three21{$pdbcoor{pdbres}}};
    $pdbcoor{modelnum} = $array[$rfield->{'pdbx_PDB_model_num'}]; #model number
    
    # assign atomtypeno according to following scheme:
    # N            0
    # CA           1
    # C            2
    # O            3
    # CB           4
    # CG/CG1/OG1   5
    # CD/CD1/OD1   6
    # CE/CE1/OE1   7
    # CZ/CZ1/OZ1   8
    # CH/CH1/OH1   9
    # CG2/OG2/NG2 10
    # CD2/OD2/ND2 11
    # CE2/OE2/NE2 12
    # CZ2/OZ2/NZ2 13
    # CH2/OH2/NH2 14
    # CE3/OE3/NE3 15
    # CZ3/OZ3/NZ3 16
    # OXT         17
    # H*          20
    
    $pdbcoor{atomtypeno}=100;
    my $atomtype=$pdbcoor{atomtype};
    if ($pdbcoor{element} eq "H") {$pdbcoor{atomtypeno}=20; return %pdbcoor;}
    if ($atomtype eq "N"   ) {$pdbcoor{atomtypeno}=0; return %pdbcoor;}
    if ($atomtype eq "CA"  ) {$pdbcoor{atomtypeno}=1; return %pdbcoor;}
    if ($atomtype eq "C"   ) {$pdbcoor{atomtypeno}=2; return %pdbcoor;}
    if ($atomtype eq "O"   ) {$pdbcoor{atomtypeno}=3; return %pdbcoor;}
    if ($atomtype eq "OXT" ) {$pdbcoor{atomtypeno}=17; return %pdbcoor;}
    if ($atomtype eq "CB"  ) {$pdbcoor{atomtypeno}=4; return %pdbcoor;}
    if ($atomtype eq "SE" && $pdbcoor{pdbres} eq "MSE" ) {
	$pdbcoor{atomtypeno}=6; return %pdbcoor;}
    
    if ($atomtype =~ /[COSN]G1|[COSN]G$/) {
	$pdbcoor{atomtypeno}=5; return %pdbcoor;}
    if ($atomtype =~ /[COSN]D1|[COSN]D$/) {
	$pdbcoor{atomtypeno}=6; return %pdbcoor;}
    if ($atomtype =~ /[COSN]E1|[COSN]E$/) {
	$pdbcoor{atomtypeno}=7; return %pdbcoor;}
    if ($atomtype =~ /[COSN]Z1|[COSN]Z$/) {
	$pdbcoor{atomtypeno}=8; return %pdbcoor;}
    if ($atomtype =~ /[COSN]H1|[COSN]H$/) {
	$pdbcoor{atomtypeno}=9; return %pdbcoor;}
    if ($atomtype =~ /[COSN]G2/) {
	$pdbcoor{atomtypeno}=10; return %pdbcoor;}
    if ($atomtype =~ /[COSN]D2/) {
	$pdbcoor{atomtypeno}=11; return %pdbcoor;}
    if ($atomtype =~ /[COSN]E2/) {
	$pdbcoor{atomtypeno}=12; return %pdbcoor;}
    if ($atomtype =~ /[COSN]Z2/) {
	$pdbcoor{atomtypeno}=13; return %pdbcoor;}
    if ($atomtype =~ /[COSN]H2/) {
	$pdbcoor{atomtypeno}=14; return %pdbcoor;}
    if ($atomtype =~ /[COSN]E3/) {
	$pdbcoor{atomtypeno}=15; return %pdbcoor;}
    if ($atomtype =~ /[COSN]Z3/) {
	$pdbcoor{atomtypeno}=16; return %pdbcoor;}
    return %pdbcoor;
}

sub vecprint {
  return sprintf("%8.3f %8.3f %8.3f\n",$_[0],$_[1],$_[2]);
}

sub transpose {
  my $a11=$_[0];
  my $a12=$_[1];
  my $a13=$_[2];
  my $a21=$_[3];
  my $a22=$_[4];
  my $a23=$_[5];
  my $a31=$_[6];
  my $a32=$_[7];
  my $a33=$_[8];
  
  return ($a11,$a21,$a31,$a12,$a22,$a32,$a13,$a23,$a33);
}

sub cvec {
  my $a1=$_[1];
  my $a2=$_[2];
  my $a3=$_[3];
  my $c=$_[0];
  return ($c*$a1,$c*$a2,$c*$a3);
}

sub vecadd {
  my $a1=$_[0];
  my $a2=$_[1];
  my $a3=$_[2];
  my $b1=$_[3];
  my $b2=$_[4];
  my $b3=$_[5];
  return ($a1+$b1,$a2+$b2,$a3+$b3);
}

sub vecsubtract {
  my $a1=$_[0];
  my $a2=$_[1];
  my $a3=$_[2];
  my $b1=$_[3];
  my $b2=$_[4];
  my $b3=$_[5];
  return ($a1-$b1,$a2-$b2,$a3-$b3);
}

sub matrixmult {
  my $a11=$_[0];
  my $a12=$_[1];
  my $a13=$_[2];
  my $a21=$_[3];
  my $a22=$_[4];
  my $a23=$_[5];
  my $a31=$_[6];
  my $a32=$_[7];
  my $a33=$_[8];
  my $b11=$_[9];
  my $b12=$_[10];
  my $b13=$_[11];
  my $b21=$_[12];
  my $b22=$_[13];
  my $b23=$_[14];
  my $b31=$_[15];
  my $b32=$_[16];
  my $b33=$_[17];
  my $c11=$a11*$b11 + $a12*$b21 + $a13*$b31;
  my $c12=$a11*$b12 + $a12*$b22 + $a13*$b32;
  my $c13=$a11*$b13 + $a12*$b23 + $a13*$b33;
  my $c21=$a21*$b11 + $a22*$b21 + $a23*$b31;
  my $c22=$a21*$b12 + $a22*$b22 + $a23*$b32;
  my $c23=$a21*$b13 + $a22*$b23 + $a23*$b33;
  my $c31=$a31*$b11 + $a32*$b21 + $a33*$b31;
  my $c32=$a31*$b12 + $a32*$b22 + $a33*$b32;
  my $c33=$a31*$b13 + $a32*$b23 + $a33*$b33;
  return ($c11,$c12,$c13,$c21,$c22,$c23,$c31,$c32,$c33);
}

sub matrixvecmult {
  my $a11=$_[0];
  my $a12=$_[1];
  my $a13=$_[2];
  my $a21=$_[3];
  my $a22=$_[4];
  my $a23=$_[5];
  my $a31=$_[6];
  my $a32=$_[7];
  my $a33=$_[8];
  my $b11=$_[9];
  my $b21=$_[10];
  my $b31=$_[11];
  my $c11=$a11*$b11 + $a12*$b21 + $a13*$b31;
  my $c21=$a21*$b11 + $a22*$b21 + $a23*$b31;
  my $c31=$a31*$b11 + $a32*$b21 + $a33*$b31;
  return ($c11,$c21,$c31);
}

sub dot {
  my($x1,$y1,$z1,$x2,$y2,$z2)=@_;
  return ($x1*$x2+$y1*$y2+$z1*$z2);
}

sub dist {
  my($x1,$y1,$z1,$x2,$y2,$z2)=@_;
  my($dist);
  $dist=($x1-$x2)**2 + ($y1-$y2)**2 + ($z1-$z2)**2 ;
  $dist=sqrt($dist);
  return $dist;
}

sub dist2 {
  my($x1,$y1,$z1,$x2,$y2,$z2)=@_;
  my($dist2);
  $dist2=($x1-$x2)**2 + ($y1-$y2)**2 + ($z1-$z2)**2 ;
  return $dist2;
}

sub angl {
  my @r1=($_[0],$_[1], $_[2]);
  my @r2=($_[3],$_[4], $_[5]);
  my @r3=($_[6],$_[7], $_[8]);
  my ($angl);
  
  my @v1=pts2normvec(@r1,@r2);
  my @v2=pts2normvec(@r3,@r2);

  my $dot=dot(@v1,@v2);
  $angl=acos($dot);
  return $r2d*$angl;
}

#pts2normvec -- calculates unit vector between two points
sub pts2normvec {
  my($x1,$y1,$z1,$x2,$y2,$z2)=@_;
  my($dist);
  $dist=sqrt(($x1-$x2)**2 + ($y1-$y2)**2 + ($z1-$z2)**2);
  
  if ($dist>0) {
    return (($x2-$x1)/$dist,($y2-$y1)/$dist,($z2-$z1)/$dist);
  }
  else {return (0,0,0);}
}

# cross  -- return unit normal vector to 2 vectors

sub cross {
  my($x1,$y1,$z1,$x2,$y2,$z2)=@_;
  my $dist ;
  my @vec;

  @vec=($y1*$z2 - $z1*$y2,
	$z1*$x2 - $x1*$z2,
	$x1*$y2 - $y1*$x2);
  $dist=sqrt($vec[0]*$vec[0] + $vec[1]*$vec[1] + $vec[2]*$vec[2]);

  if ($dist>1.0e-8) {
    return ($vec[0]/$dist,$vec[1]/$dist,$vec[2]/$dist);
  }
  else {return (0,0,0);}
}

sub euler {
  my $a11=$_[0];
  my $a12=$_[1];
  my $a13=$_[2];
  my $a21=$_[3];
  my $a22=$_[4];
  my $a23=$_[5];
  my $a31=$_[6];
  my $a32=$_[7];
  my $a33=$_[8];

  my $phi=   0;
  my $theta= 0;
  my $psi=   0;

  my $pi=acos(-1.0);

  my $cosphi;
  my $cospsi;
  my $sinphi;
  my $sinpsi;

# Calculate the rotation matrix

  $theta=acos($a33);

  if ($theta!=0) {
    $sinpsi= $a13/sin($theta);
    $cospsi= $a23/sin($theta);
    $sinphi= $a31/sin($theta);
    $cosphi=-$a32/sin($theta);
    if ($cosphi> 1.0) {$cosphi= 1.0;}
    if ($cosphi<-1.0) {$cosphi=-1.0;}
    if ($cospsi> 1.0) {$cospsi= 1.0;}
    if ($cospsi<-1.0) {$cospsi=-1.0;}
    if ($sinphi> 1.0) {$sinphi= 1.0;}
    if ($sinphi<-1.0) {$sinphi=-1.0;}
    if ($sinpsi> 1.0) {$sinpsi= 1.0;}
    if ($sinpsi<-1.0) {$sinpsi=-1.0;}

    if ($sinphi>=0) {$phi=acos($cosphi);}
    else {$phi=2*$pi-acos($cosphi);}

    if ($sinpsi>=0) {$psi=acos($cospsi);}
    else {$psi=2*$pi-acos($cospsi);}
    
  }
  
  else  {
    $sinphi=$a12;
    $cosphi=$a11;
    if ($cosphi> 1.0) {$cosphi= 1.0;}
    if ($cosphi<-1.0) {$cosphi=-1.0;}
    if ($sinphi> 1.0) {$sinphi= 1.0;}
    if ($sinphi<-1.0) {$sinphi=-1.0;}
    if ($sinphi>=0) {$phi=acos($cosphi);}
    else {$phi=2*$pi-acos($cosphi);}
    $psi=0;
  }
  
  return ($phi,$theta,$psi);
}

sub dihedral {
  my @r1=($_[0],$_[1], $_[2]);
  my @r2=($_[3],$_[4], $_[5]);
  my @r3=($_[6],$_[7], $_[8]);
  my @r4=($_[9],$_[10],$_[11]);
  
  my @v1=pts2normvec(@r1,@r2);
  my @v2=pts2normvec(@r3,@r2);
  my @n1=cross(@v1,@v2);
  my @v3=pts2normvec(@r2,@r3);
  my @v4=pts2normvec(@r4,@r3);
  my @n2=cross(@v3,@v4);
  my @n3=cross(@n1,@n2);

  my $alpha=0;
  if (dot(@n3,@v3)>0) {$alpha=acos(dot(@n1,@n2)); }
  elsif (dot(@n3,@v3)<0) {$alpha=-acos(dot(@n1,@n2));}
  else {$alpha=$pi;}
  return $r2d*$alpha;
}

sub build {
  my @r1=($_[0],$_[1], $_[2]);
  my @r2=($_[3],$_[4], $_[5]);
  my @r3=($_[6],$_[7], $_[8]);
  my $b=$_[9];
  my $a=$d2r*$_[10];
  my $d=$d2r*$_[11];
  
  my @j=pts2normvec(@r2,@r3);
  my @k=cross(pts2normvec(@r2,@r1),@j);
  my @i=cross(@j,@k);
  my $cj=-$b*cos($a);
  my $ci= $b*sin($a)*cos($d);
  my $ck=-$b*sin($a)*sin($d);
  my @r4=(0,0,0);
  @r4=vecadd(@r3,cvec($cj,@j));
  @r4=vecadd(@r4,cvec($ci,@i));
  @r4=vecadd(@r4,cvec($ck,@k));
  return ($r4[0],$r4[1],$r4[2]);
}

sub cart2sph {
  my $x=$_[0];
  my $y=$_[1];
  my $z=$_[2];
  my $r=sqrt($x*$x + $y*$y + $z*$z);
  my $phi=0.0;
  if ($x>0) {$phi=$r2d*atan($y/$x);}
  elsif ($x==0) {$phi=$r2d*$pi/2.0;}
  else {$phi=$r2d*(atan($y/$x)+$pi);}
  my $theta=0.0;
  if ($r>0) {$theta=$r2d*acos($z/$r);}
  if ($phi<0) {$phi+=360.0;}
  return ($r,$phi,$theta);
}

sub transrot {
  my @c1=( $_[0], $_[1], $_[2]);
  my @n1=( $_[3], $_[4], $_[5]);
  my @a1=( $_[6], $_[7], $_[8]);
  my @c2=( $_[9], $_[10],$_[11]);
  my @n2=( $_[12],$_[13],$_[14]);
  my @a2=( $_[15],$_[16],$_[17]);

  my @i1= pts2normvec(@n1,@c1);
  my @i2= pts2normvec(@n2,@c2);
  my @na1=pts2normvec(@n1,@a1);
  my @na2=pts2normvec(@n2,@a2);
  my @k1= cross(@i1,@na1);
  my @k2= cross(@i2,@na2);
  my @j1= cross(@k1,@i1);
  my @j2= cross(@k2,@i2);

  my @M01=(@i1,@j1,@k1);
  my @M02=(@i2,@j2,@k2);
  my @M20=transpose(@M02);

  my @Mst=matrixmult(@M01,@M20);
  
  my @rst=matrixvecmult(@M01,vecsubtract(@n2,@n1));
  my @est=euler(@Mst);
  
  @est=cvec($r2d,@est);
  @rst=cart2sph(@rst);
  return (@rst,@est);

}


#!/usr/bin/perl -w

### The master perl script of the IS-score program.
### For documentation, type "./isscore.pl -m" in the command line.
###
### Mu Gao
### E-mail: <mu.gao@gatech.edu>

use strict;
use File::Basename;
use File::Path;
use Getopt::Long;
use Pod::Usage;
use Time::Local;

use FindBin;
use lib "$FindBin::Bin";
use Constants;


############ Global Parameters #####################################
#### Do not modify these parameters unless you know what you are doing
my $CONT_CUTOFF  = 4.5;    #### 4.5 Angstrom distance cutoff for contact calculations
my $MIN_NUM_RES  = 25;     #### minimum number of protein residues, ignore short peptide in a PDB file
my $MIN_INT_RES  = 15;     #### minimum number of protein residues of a interface

my $VERSION = '1.1b5';  ### version number of this package


###########  Constants for Sequence Alignment  #####################
my @NOS = (1,2,3);
my ($UP,$LEFTUP,$LEFT) = @NOS;
my ($d,$e) = (-11,-1);  #### gap opening and extension penalty
my $matrix = \@Constants::blosum62;
####################################################################


###################### Input & options #############################
my $pdb_file1;
my $pdb_file2;
my $pchains1 = '';
my $pchains2 = '';

my $work_path = '.';
my $pdb_path;

my $out_file;
my $set_file1;

my $help;
my $manual;
my $version;

my $aln_flag = 2;
my $tra_flag;
my $lib_flag;
my $sel_flag;
my $seq_aln_once_flag;
my $ca_flag;
my $ow_flag;

my $measure  = 'is';
my $vmd;



### get the options
  GetOptions(
	     #### input files
	     'l|liblst1|l1=s'=> \$set_file1,    ### a list of entries
	     'p1|pdb1=s'     => \$pdb_file1,    ### pdb file1
	     'p2|pdb2=s'     => \$pdb_file2,    ### pdb file2
	     'c1|chains1=s'  => \$pchains1,     ### list of protein chains of pdb1
	     'c2|chains2=s'  => \$pchains2,     ### list of protein chains of pdb2
	     #### paths
	     'p|pdbpath=s'   => \$pdb_path,     ### path to original pdb files
	     'w|workpath=s'  => \$work_path,    ### path to input files of interfaces
	     #### output file
             'o|output=s'    => \$out_file,     ### save results to an output file
	     #### options
             'a|aln=i'       => \$aln_flag,     ### choose a format of alignment printout
             'ow|=s'         => \$ow_flag,      ### choose overwrite options
             't|trans'       => \$tra_flag,     ### the first structure aligned to the second
             'e|measure=s'   => \$measure,      ### similarity measure: TM-score or IS-score
             'vmd=s'         => \$vmd,          ### generate tcl script for VMD
	     #### advanced options
             'ca'            => \$ca_flag,      ### Calpha models
	     'dc=f'          => \$CONT_CUTOFF,  ### distance cutoff for extracting interface
	     'minp=i'        => \$MIN_NUM_RES,  ### minimum number of protein residues
	     'mini=i'        => \$MIN_INT_RES,  ### minimum number of interface residues
             'seqonce'       => \$seq_aln_once_flag,     ### only perform sequence alignment once
	     #### help
             'h|help'        => \$help,         ### print the help message
             'm|manual'      => \$manual,       ### print the manual page
             'v|version'     => \$version       ### print the version number;
            ) or pod2usage(0);


  pod2usage(-verbose => 0) if $help;
  pod2usage(-verbose => 2) if $manual;

  if ($version) {
    print "IS-score version $VERSION\n" ;
    exit;
  }



  my $num_arg = $#ARGV + 1;
  if( $num_arg == 4 ) {
    $pdb_file1 = $ARGV[0];
    $pchains1  = $ARGV[1];
    $pdb_file2 = $ARGV[2];
    $pchains2  = $ARGV[3];
  }
  elsif( $num_arg == 2 and not defined $set_file1 ) {
    $pdb_file1 = $ARGV[0];
    $pdb_file2 = $ARGV[1];
  }
  elsif( $num_arg == 2 and defined $set_file1 ) {
    $pdb_file2 = $ARGV[0];   ### template is the first, target/query is the second
    $pchains2  = $ARGV[1];
  }
  elsif( $num_arg == 1 and defined $set_file1 ) {
    $pdb_file2 = $ARGV[0];
  }
  elsif( $num_arg > 0 ) {
    pod2usage(-verbose =>0, -msg=>'Error: invalid argument.');
  }


  unless( (defined $pdb_file1 and defined $pdb_file2) or
	  (defined $set_file1 and defined $pdb_file2) ) {
    pod2usage(-verbose =>0, -msg=>'Error: invalid argument.');
  }


  if( defined $vmd ) {
    $aln_flag = 2;
    if( $vmd ne 'cartoon' and $vmd ne 'ribbon' ) {
      $vmd = 'cartoon';
    }
  }

  if( defined $ow_flag ) {
    if( $ow_flag eq 'y' or $ow_flag eq '1' ) {
      $ow_flag = 'a';   ### always overwrite interface files if they already exists
    }
    elsif( $ow_flag eq 'n' or $ow_flag eq '0' ) {
      $ow_flag = 'n';
    }
    else {
      print "Warning: $ow_flag is invalid, setting the overwrite option as always overwrite\n";
      $ow_flag = 'a';
    }
  }


#### creating two sets of structures for IS-align runs  ####
my @pdb_set1 = ();
my @pdb_set2 = ();

if( defined $pdb_file1 ) { $pdb_set1[0] = {'pdb_file'=>$pdb_file1,'pchains'=>$pchains1}; }
if( defined $pdb_file2 ) { $pdb_set2[0] = {'pdb_file'=>$pdb_file2,'pchains'=>$pchains2}; }
if( defined $set_file1 ) { readSetFile( $set_file1, \@pdb_set1 ); }



### create the work_path if it does not exist
unless( -d $work_path ) {
  mkpath( $work_path ) or die "Error: could not create the directory $work_path\n";
}

### get the path of executables
my ($bin_base, $bin_path, $bin_ext) = fileparse( $0, qr/\.pl/ );

### options passing to IS-score
my $switch = '';
if( $measure eq 'tm' )  { $switch .= ' -t '; }
elsif( $measure ne 'is' ) {
  die "Error: the measure option can only be \"tm\" or \"is\".\n";
}


if( $aln_flag == 0 or $aln_flag == 1 or $aln_flag == 2 ) {
  $switch .= " -v $aln_flag ";
}
else { die "Error: invalid option for alignment printout format\n"; }





##-------------------------- end of input & options   -------------------------##


#### determine where to ouput results ####
my $out;
if( defined $out_file ) {
  open $out, ">$out_file" or die "Error: could not open $out_file\n";
}
else {
  open $out, ">&STDOUT";
}
print_banner( $out );


print $out "IS-score starts at ", timeStamp(), "\n";

######################################################################
######### Stage 1: processing raw PDB files, sanity check   ##########
######################################################################
print $out "\n+++++ Stage 1: Processing PDB files ...\n";
print "+++++ Stage 1: Processing PDB files ...\n" if(defined $out_file);

my %parsed_model = ();   ### save information of processed PDB files
for (@pdb_set1) {
  my $name = $_->{pdb_file} . '_' . $_->{pchains};
  $parsed_model{$name} = $_;
}

my %parsed_target = ();
for (@pdb_set2) {
  my $name = $_->{pdb_file} . '_' . $_->{pchains};
  $parsed_target{$name} = $_;
}

my $target_name;
for my $name (sort keys %parsed_target) {
  my ($parsed_pdb_file, $pchains) = parsePDBFile( $parsed_target{$name}, $pdb_path, $work_path, $bin_path, $out );
  $target_name = $name;
  $parsed_target{$name}->{'parsed_file'} = $parsed_pdb_file;
  $parsed_target{$name}->{'chains'}      = $pchains;
}



######################################################################
######### Stage 2:  extracting protein-protein interfaces  ###########
######################################################################
print $out "\n+++++ Stage 2: Extracting protein-protein interfaces ...\n";
print $out "Mininum number of interface residue required: $MIN_INT_RES AAs\n";
print "+++++ Stage 2: Extracting protein-protein interfaces ...\n" if (defined $out_file);


my $target_pdb_file = $parsed_target{$target_name}->{parsed_file};
my $target_pchains  = $parsed_target{$target_name}->{chains};
$parsed_target{$target_name}->{ppi} = extractProtProtInt( $bin_path, $target_pdb_file, $target_pchains, $out );
print "Target: $target_pdb_file\n";

my $target_ppi      = $parsed_target{$target_name}->{ppi};
my $target_con_file = $target_ppi->[0]->{confile};

if( not defined $target_con_file ) {
  die "No valid ppi found for the target, stopping IS-score.\n";
}

my ($tar_int_lst, $tar_int_res) = readContFile( $target_con_file );
my ($tar_struct, $tar_seq_of, $tar_ind2resid_of) = readPDBStructALL( $target_pdb_file );

my $tar_resid_to_mod_resid;
my $counter = 0;
for my $mod_name (sort keys %parsed_model) {
  my ($parsed_pdb_file, $pchains) = parsePDBFile( $parsed_model{$mod_name}, $pdb_path, $work_path, $bin_path, $out );

  $parsed_model{$mod_name}->{'parsed_file'} = $parsed_pdb_file;
  $parsed_model{$mod_name}->{'chains'}      = $pchains;

  print "Model: $parsed_pdb_file\n";
  my ($mod_struct, $mod_seq_of, $mod_ind2resid_of) = readPDBStructALL( $parsed_pdb_file );

  if( not defined $seq_aln_once_flag or $counter == 0 ) {
    $tar_resid_to_mod_resid = alignModTar( $mod_seq_of, $pchains, $mod_ind2resid_of,
					    $tar_seq_of, $target_pchains, $tar_ind2resid_of );
  }

  $parsed_model{$mod_name}->{'tar_resid_to_mod_resid'} = $tar_resid_to_mod_resid;

  $parsed_model{$mod_name}->{ppi} = writeModFile( $tar_struct, $tar_int_lst, $tar_int_res, $target_pchains,
						  $mod_struct, $bin_path, $parsed_pdb_file, $pchains, $out, $tar_resid_to_mod_resid );

  $counter++;
}






#####################################################################
############ Stage 3:  performing interface alignment  ##############
#####################################################################
print $out "\n+++++ Stage 3: Calculating IS-score ...\n";
print "+++++ Stage 3: Calculating IS-score ...\n" if(defined $out_file);

my $set1 = \@pdb_set1;
my $set2 = \@pdb_set2;


my $num1 = scalar @$set1;
my $num2 = scalar @$set2;

$counter  = 1;
for(my $m=0; $m<$num1; $m++) {
  my $name1    = $$set1[$m]->{pdb_file} . '_' . $$set1[$m]->{pchains};
  my $ppi1     = $parsed_model{$name1}->{ppi};
  my $num_ppi1 = scalar @$ppi1;
  my $tar_resid_to_mod_resid = $parsed_model{$name1}->{tar_resid_to_mod_resid};

  for(my $i=0; $i<$num_ppi1; $i++) {

    for(my $k=0; $k<$num2; $k++) {
      my $name2     = $$set2[$k]->{pdb_file} . '_' . $$set2[$k]->{pchains};
      my $ppi2      = $parsed_target{$name2}->{ppi};
      my $num_ppi2  = scalar @$ppi2;

      for(my $j=0; $j<$num_ppi2; $j++) {

	my ($result, $trans_vec, $rot_mat, $parsed_out) = callISscore( $$ppi1[$i], $$ppi2[$j], $switch, $out, $tra_flag );
	if( defined $out_file ) { print ">>>$counter $result\n"; }
	if( defined $vmd ) {
	  my $pair        = $$ppi1[$i]->{name} . '_' . $$ppi2[$j]->{name};
	  my $script_file = "$pair.vmd";
	  my $pdbname1    = $$set1[$m]->{pdb_file};
	  my $pdbname2    = $$set2[$k]->{pdb_file};

	  mk_vmdscript( $script_file, $pdbname1, $pdbname2, $trans_vec, $rot_mat, $parsed_out, $vmd, $tar_resid_to_mod_resid );
	}

	$counter++;
      } ### j interface 2
    } ### k protein 2

  } ### i interface 1
} ### m protein 1



print $out "IS-score ends at ", timeStamp(), ".\n\n";

close $out;

####------------------------ End of main routine --------------------------------####





#####################################################################################
############################         Subroutine        ##############################
#####################################################################################




#### perform interface alignment with IS-align
sub callISscore {
  my ($ppi1, $ppi2, $switch, $out, $tra_flag) = @_;

  my $name1     = $ppi1->{name};
  my $int_file1 = $ppi1->{intfile};
  my $con_file1 = $ppi1->{confile};
  my $mat_file  = $ppi1->{matfile};

  my $name2     = $ppi2->{name};
  my $int_file2 = $ppi2->{intfile};
  my $con_file2 = $ppi2->{confile};

  my $pair = "$name1 vs $name2";
  print $out ">>>$name1 vs $name2\n";

  ### call IS-score
  my $alnout = `$bin_path/IS-score $switch -m $mat_file $int_file1 $int_file2 $con_file1 $con_file2`;

  ### parsing output of IS-score
  my ($trans_vec, $rot_mat, $score_ln, $parsed_out) = parseISscoreOut( $alnout );

  ### outputing results
  print $out "\n$parsed_out\n\n";

  ### transform original pdbfile1 to give the best inteface alignment to pdbfile2
  if( defined $tra_flag ) {
    my $org_pdbfile = $ppi1->{pdbfile};

    my $len1 = length( $name1 );
    my $org_name = substr( $name1, 0, $len1-2 );
    my $tra_pdbfile = "$name1\_to\_$name2.pdb";

    if( transform_pdb( $org_pdbfile, $tra_pdbfile, $rot_mat, $trans_vec ) ) {
      print $out "Transformed $org_name was saved to $tra_pdbfile\n";
    }
  }

  my $result = "$pair: $score_ln";
  return ($result, $trans_vec, $rot_mat, $parsed_out);
}



#### read a list of pdb files ####
sub readSetFile {
  my ($set_file, $pdb_set) = @_;

  open SET, "<$set_file" or die "Error: could not open $set_file\n";
  while(<SET>){
    next if /^#/;
    chomp;
    my @fields   = split(' ',$_);
    my $pdb_file = $fields[0];
    my $pchains  = '';
    if( defined $fields[1] and not ($fields[1] =~ /\W/) ) { $pchains = $fields[1]; }

    push( @$pdb_set, { 'pdb_file'=>$pdb_file, 'pchains'=>$pchains } );
  }
  close SET;
}



######### parse a raw PDB file  ######
sub parsePDBFile {
  my ($parsed_rec, $pdb_path, $work_path, $bin_path, $outfile) = @_;

  my $pdb_file   = $parsed_rec->{pdb_file};
  my $pchain_lst = $parsed_rec->{pchains};

  my ($base, $path, $ext) = fileparse( $pdb_file, qr/\.\w*$/ );
  my $parsed_pdb_file = "$work_path/$base.parsed";

  my @pchains = ();
  unless ( -s $parsed_pdb_file ) {
    if( not (-e $pdb_file) and defined $pdb_path ) {
      $pdb_file = "$pdb_path/$pdb_file";
    }

    unless( -s $pdb_file ) {
      if( -s "$pdb_file.pdb" )    { $pdb_file = "$pdb_file.pdb"; }
      elsif( -s "$pdb_file.ent" ) { $pdb_file = "$pdb_file.ent"; }
      else {
	die "Error: could not find the pdb file $pdb_file\n";
      }
    }

    my $out = `perl $bin_path/process_pdb.pl -pdb $pdb_file -o $parsed_pdb_file -nowarn`;
  }

  @pchains = split(//,$pchain_lst);
  if( scalar @pchains == 0 ) { @pchains = getProtChains( $parsed_pdb_file ); }

  print $outfile "$pdb_file:  found ", scalar @pchains, " protein chains ", join(' ',@pchains), "\n";

  return( $parsed_pdb_file, \@pchains );
}






########### get protein chains from  Calpha atoms ##############
sub getProtChains {
  my $pdbfile = shift @_;

  open PDB, "<$pdbfile" or die "Error: could open file $pdbfile\n";
  my %pchain = ();
  while (my $line = <PDB>) {
    next if ( $line =~ /^REMARK/ );
    last if ( $line =~ /^ENDMDL/ );  #### only consider the first model

    #### PDB format: http://www.wwpdb.org/documentation/format23/sect9.html ####
    if ( $line =~ /^(ATOM  |HETATM)/ ) {
      (my $atomname = substr($line, 12, 4))=~ s/ //g;
      (my $altloc = substr($line, 16, 1));
      (my $resname = substr($line,17,3)) =~ s/ //g;
      (my $chain = substr($line, 21, 1));
      (my $resid = substr($line, 22, 4)) =~ s/ //g;
      (my $icode = substr($line, 26, 1)) =~ s/ //g;

      if( $atomname eq 'CA' && convAA( $resname ) ne 'X' ) {

	if ( $chain eq '' ) { $chain = '_'; }

	unless( exists $pchain{$chain} ) {
	  $pchain{$chain} = 0;
	}

	$pchain{$chain}++;
      }
    }

  }

  close PDB;

  my @protein_chains = ();
  foreach my $chain (sort keys %pchain) {
    if( $pchain{$chain} >= $MIN_NUM_RES ) {
      #push( @protein_chains, { 'id'=>$chain, 'len'=>$pchain{$chain} } );
      push( @protein_chains, $chain );
    }
  }

  return @protein_chains;
}


### read protein Calpha atoms or DNA P atoms from a pdb file
sub readPDBStructALL {
  my $pdbfile = shift @_;

  my %struct = (); my %seq_of = (); my %ind2resid_of = ();
  my $res_index = 1; my $atm_index = 1;
  open PDB, "<$pdbfile" or die "Eorror: could not open $pdbfile\n";
  while (my $line = <PDB>) {
    next if ( $line =~ /^REMARK/ );
    last if ( $line =~ /^ENDMDL/ );  #### only consider the first model

    #### PDB format: http://www.wwpdb.org/documentation/format23/sect9.html ####
    if ( $line =~ /^(ATOM  |HETATM)/ ) {
      (my $atomnum = substr($line, 6, 5)) =~ s/ //g;  #remove spaces
      (my $atomname = substr($line, 12, 4)) =~ s/ //g;
      (my $altloc = substr($line, 16, 1));
      (my $resname = substr($line,17,3)) =~ s/ //g;
      (my $chain = substr($line, 21, 1));
      (my $resid = substr($line, 22, 4)) =~ s/ //g;
      (my $icode = substr($line, 26, 1)) =~ s/ //g;
      (my $x = substr($line,30,8)) =~ s/ //g;
      (my $y = substr($line,38,8)) =~ s/ //g;
      (my $z = substr($line,46,8)) =~ s/ //g;
#      (my $occupancy = substr($line, 54, 6)) =~ s/ //g;
#      (my $bfactor = substr($line, 60, 6)) =~ s/ //g;
#      (my $element = substr($line, 76, 2)) =~ s/ //g;

      ### save atom information to the corresponding residue records ###
      my $residue = "$chain.$resid$icode";

      unless (exists $struct{$residue}) {
	$struct{$residue} = { 'resind'=>$res_index, 'atom'=>{} };
	$res_index++;
      }

      unless (exists $seq_of{$chain}) {
	$seq_of{$chain} = '';
	$ind2resid_of{$chain} = {};
      }

      if( exists $struct{$residue}->{atom}->{$atomname} ) {
	print "Warning: residue $residue atom $atomname appears more than once!\n";
	next;
      }

      if( $atomname eq 'CA' ) {
	my $resnm = convAA($resname);
	if( $resnm ne 'X' ) {
	  my $index          = length( $seq_of{$chain} );
	  $seq_of{$chain}   .= $resnm;
	  $ind2resid_of{$chain}->{$index} = $residue;
	}
      }

      $struct{$residue}->{atom}->{$atomname} = { 'coor'=>[$x, $y, $z], 'atmind'=>$atm_index, 'line'=>$line };
      $atm_index++;

    }
  }
  close PDB;

  return (\%struct, \%seq_of, \%ind2resid_of);
}


######## convert residues from three-letter to one-letter and vice visa ############
sub convAA {
  my $res = shift @_;

  my %ts=(
     'GLY'=>'G',     'ALA'=>'A',     'VAL'=>'V',     'LEU'=>'L',
     'ILE'=>'I',     'SER'=>'S',     'THR'=>'T',     'CYS'=>'C',
     'MET'=>'M',     'PRO'=>'P',     'ASP'=>'D',     'ASN'=>'N',
     'GLU'=>'E',     'GLN'=>'Q',     'LYS'=>'K',     'ARG'=>'R',
     'HIS'=>'H',     'PHE'=>'F',     'TYR'=>'Y',     'TRP'=>'W',
     'HSD'=>'H',     'HSP'=>'H',
  );

  if (exists $ts{$res}) {
    return $ts{$res};
  } else {
    return "X";
  }
}


### read interface contact file ####
sub readContFile{
  my $file = shift @_;

  my %contlist = ();
  my %contres  = ();
  open INP, "<$file" or die "Error: could not open $file\n";
  my $flag = 0;
  while(<INP>) {
    next if /^#/;

    if( /^Total residue-residue contacts found .* (\d+)/ ) {
      my $num = $1;
      if( $num != scalar keys %contlist ) {
	die "Error: expect $num contact in $file, found ", scalar keys %contlist, "\n";
      }
      last;
    }

    if( $flag ) {
      my @fields = split(' ',$_);
      next unless (scalar @fields == 9);
      my $chain1 = $fields[0];
      my $resid1 = $fields[1];
      my $resnm1 = $fields[2];
      my $chain2 = $fields[4];
      my $resid2 = $fields[5];
      my $resnm2 = $fields[6];

      my $res1 = "$chain1.$resid1";
      my $res2 = "$chain2.$resid2";
      my $contact = "$res1-$res2";

      unless(exists $contlist{$contact}) {
	$contlist{$contact} = 1;
      }

      unless(exists $contres{$res1}) {
	$contres{$res1} = 1;
      }

      unless(exists $contres{$res2}) {
	$contres{$res2} = 1;
      }
    }

    $flag = 1 if( /^R_ch/ );

  }
  close INP;


  return (\%contlist, \%contres);
}

######## extract coordinates and contacts of native interface residues from models  #######
sub writeModFile {
  my ($nat_struct, $nat_int_lst, $nat_int_res, $nat_pchains,
      $mod_struct, $bin_path, $mod_pdb_file, $mod_pchains, $outfile, $nat_resid_to_mod_resid ) = @_;

  my ($base, $path, $ext) = fileparse( $mod_pdb_file, qr/\.\w*$/ );

  my $mod_re_ch = $$mod_pchains[0];
  my $mod_li_ch = $$mod_pchains[1];

  my $nat_re_ch = $$nat_pchains[0];
  my $nat_li_ch = $$nat_pchains[1];



  ###### find locations of the native interfacial residues in the model structure  ########
  my $struct = "";
  my $matchlst = "";
  my $prev_chain = ""; my $curr_chain = "";
  my $counter=0;
  for my $nat_resid (sort {$nat_struct->{$a}->{resind} <=> $nat_struct->{$b}->{resind}} keys %$nat_struct) {
    my $curr_chain = substr($nat_resid,0,1);
    if( $prev_chain ne "" and $curr_chain ne $prev_chain ) {
      $struct .= "TER\n";
    }
    $prev_chain = $curr_chain;

    my $mod_resid = $nat_resid_to_mod_resid->{$nat_resid};
    next unless( defined $mod_resid );  ### only consider sequence aligned regions

    if(exists $nat_int_res->{$nat_resid}) {
      unless(exists $mod_struct->{$mod_resid}) {
	print "Warning: could not locate the model residue corresponding to the target residue $nat_resid\n";
	next;
      }

      my $nat_line = $nat_struct->{$nat_resid}->{atom}->{CA}->{line};
      my $nat_ch   = substr($nat_line, 21, 1);
      my $nat_resid= substr($nat_line, 22, 5);
      my $nat_resnm= substr($nat_line, 17, 3);

      my $res = $mod_struct->{$mod_resid}->{atom};

      for my $atomname (sort {$res->{$a}->{atmind} <=> $res->{$b}->{atmind}} keys %$res) {
	unless(exists $mod_struct->{$mod_resid}->{atom}->{$atomname}) {
	   print "Warning: could not locate atom corresponding to $nat_resid $atomname of the target (2nd) structure\n";
	   next;
        }
	my $line = $mod_struct->{$mod_resid}->{atom}->{$atomname}->{line};
	#$struct .= $line;

	if($atomname eq 'CA') {
	  $counter++;
	  my $mod_ch    = substr($line,21,1);
	  my $mod_resid = substr($line,22,5);
	  my $mod_resnm = substr($line,17,3);
	  $matchlst .= sprintf "%5d %1s %5s %3s %1s %5s %3s\n",$counter,$mod_ch,$mod_resid,$mod_resnm,
	    $nat_ch,$nat_resid,$nat_resnm;
	}

      }
    }

  }
  ###############################################################################################


  ########## generate input files for IS-score ##########


  my $pair         = $mod_re_ch . $mod_li_ch;
  my $mod_name     = $base  . $pair;
  my $mod_int_file = $path  . "$mod_name\_int.pdb";
  my $mod_con_file = $path  . "$mod_name\_con.lst";
  my $mod_mat_file = $path  . "$mod_name\_mat.lst";

  my $exists_flag = 0;
  if( -s $mod_mat_file ) {
    if( (defined $ow_flag and $ow_flag eq 'y') or not defined $ow_flag ) {
      print "$mod_mat_file detected, overwrite the file? [(a)ll/(y)es/(n)o]: ";
      my $choice = <STDIN>;
      chomp $choice;
      if   ( $choice eq 'all' or $choice eq 'a' ) { $ow_flag = 'a'; }
      elsif( $choice eq 'yes' or $choice eq 'y' ) { $ow_flag = 'y'; }
      elsif( $choice eq 'no'  or $choice eq 'n' ) { $ow_flag = 'n'; }
      else { $ow_flag = 'n'; print "Invalid choice, setting the option to 'n'\n"; }
    }
    $exists_flag = 1;
  }

  my $num_int_res = 0;

  if( not $exists_flag or not defined $ow_flag or $ow_flag ne 'n' ) {
    open my $mat_fh, ">$mod_mat_file" or die "Error: could not write $mod_mat_file\n";
    print $mat_fh "### Matched model residues to target interfacial residues according to sequence alignment\n";
    print $mat_fh "Index Ch1 Res1 AA1 Ch2 Res2 AA2\n";
    print $mat_fh $matchlst;
    close $mat_fh;

    my $output;
    if( not defined $ca_flag ) {
      $output = `$bin_path/extint -s $mod_pdb_file -r $mod_re_ch -l $mod_li_ch -c $mod_con_file -i $mod_int_file -res $mod_mat_file -wa -d $CONT_CUTOFF`;
    }
    else {
      $output = `$bin_path/extint -s $mod_pdb_file -r $mod_re_ch -l $mod_li_ch -c $mod_con_file -i $mod_int_file -res $mod_mat_file -wa -scom`;
    }

    if( $output =~ / (\d+) interface residues extracted/ ) {
      $num_int_res = $1;
    }
    else { die "Error: failed in extracting the interfaces\n"; }
  }

  else {
    my $output = `grep "^Total interacting residues" $mod_con_file`;
    if( $output =~ /Total interacting residues.* (\d+)/ ) {
      $num_int_res = $1;
    }
    else { die "Error: file $mod_con_file seems corrupted\n"; }
  }




  printf $outfile "$mod_pdb_file: found 1 PPI(s) $pair $num_int_res\n";

  my @ppi = ( { 'name'=>$mod_name, 'size'=>$num_int_res, 'pair'=>$pair, 'matfile'=>$mod_mat_file,
	     'intfile'=>$mod_int_file, 'confile'=>$mod_con_file, 'pdbfile'=>$mod_pdb_file } );

  return \@ppi;
}


################

sub extractProtProtInt {
  my ($bin_path, $parsed_pdb_file, $pchains, $outfile) = @_;

  my ($base, $path, $ext) = fileparse( $parsed_pdb_file, qr/\.\w*$/ );

  my @ppi = ();
  my $num_pchains = scalar @$pchains;

  if($num_pchains > 2) {
    print "Warning: more than two protein chains detected, only considers first two chains $$pchains[0] $$pchains[1]\n";
  }

  for(my $i=0; $i<1; $i++) {
    my $pchain_a = $$pchains[$i];
    for(my $j=$i+1; $j<2; $j++) {
      my $pchain_b = $$pchains[$j];
      my $pair     = $pchain_a . $pchain_b;
      my $name     = $base . $pair;
      my $int_file = "$path$name\_int.pdb";
      my $con_file = "$path$name\_con.lst";

      my $num_int_res = 0;

      my $exists_flag = 0;
      if( -s $int_file and -s $con_file ) {
	if( (defined $ow_flag and $ow_flag eq 'y') or not defined $ow_flag ) {
	  print "$int_file detected, overwrite the file? [(a)ll/(y)es/(n)o]: ";
	  my $choice = <STDIN>;
	  chomp $choice;
	  if   ( $choice eq 'all' or $choice eq 'a' ) { $ow_flag = 'a'; }
	  elsif( $choice eq 'yes' or $choice eq 'y' ) { $ow_flag = 'y'; }
	  elsif( $choice eq 'no'  or $choice eq 'n' ) { $ow_flag = 'n'; }
	  else { $ow_flag = 'n'; print "Invalid choice, setting the option to 'n'\n"; }
	}
	$exists_flag = 1;
      }

      if( $exists_flag and defined $ow_flag and $ow_flag eq 'n' ) {
	my $output = `grep "^Total interacting residues" $con_file`;
	if( $output =~ /Total interacting residues.* (\d+)/ ) {
	  $num_int_res = $1;
	}
	else { die "Error: file $con_file seems corrupted\n"; }
      }
      else {
	my $output;
	if( not defined $ca_flag ) {
	  $output = `$bin_path/extint -s $parsed_pdb_file -r $pchain_a -l $pchain_b -d $CONT_CUTOFF -i $int_file -c $con_file -e`;
	}
	else {
	  $output = `$bin_path/extint -s $parsed_pdb_file -r $pchain_a -l $pchain_b -i $int_file -c $con_file -e -scom`;
	}



	if( $output =~ / (\d+) interface residues extracted/ ) {
	  $num_int_res = $1;
	}
	else { die "Error: extraction of interfaces failed\n"; }
      }

      if( $num_int_res >= $MIN_INT_RES ) {
        push( @ppi, { 'name'=>$name, 'size'=>$num_int_res, 'pair'=>$pair,
		      'intfile'=>$int_file, 'confile'=>$con_file, 'pdbfile'=>$parsed_pdb_file } );
      }
    }
  }

  printf $outfile "$parsed_pdb_file: found %d valid PPI(s)", scalar @ppi;
  for(my $i=0; $i<scalar @ppi; $i++) {
    print $outfile " $ppi[$i]->{pair} $ppi[$i]->{size}";
  }
  print $outfile "\n";

  return \@ppi;
}




##### get translation vector and rotation matrix from IS-align output
sub parseISscoreOut {
  my $alnout = shift @_;

  my @output = split(/\n/, $alnout);

  ### remove useless information for end user
  my $num = scalar @output;
  for(my $i=0; $i < $num; $i++) {
    last if ( $output[0] =~ /^Structure 1/ );
    shift @output;
  }

  ### find measurement of the alignment
  my ($alnlen, $rmsd, $col, $score, $pvalue);
  my $score_ln;
  my @trans_vec = ();
  my @rot_mat = ();
  for (my $i=0; $i < scalar @output; $i++)  {
    my $line = $output[$i];
    #if ( $line =~ /^Structure1\:.*\s(\d+)\s+AAs/ ) { $templsize = $1; }

    if ( $line =~ /^\w{2}-score\s+=\s*(\d+\.\d+)\s+P-value =\s*(\S+)/ ) {
      $score = $1;
      $pvalue= $2;
      $score_ln = $line;
    }
    elsif ( $line =~ /^\w{2}-score\s+=\s*(\d+\.\d+)/ ) {
      $score = $1;
      $score_ln = $line;
	#print "$mod_ind $score $col\n";
    }
    elsif ( $line =~ /^Number of matched residues =\s*(\d+)/ ) {
      $alnlen = $1;
    }
    elsif ( $line =~ /^Number of matched contacts =\s*(\d+)/ ) {
      $col = $2;
    }
    elsif ( $line =~ /^RMSD\s+=\s*(\d+\.\d+)/ ) {
      $rmsd = $1;
    }


    elsif( $line =~ /rotation matrix/ ) {
      for(my $j=0; $j<3; $j++) {
        my @fields = split(' ',$output[$i+$j+2]);
        push( @trans_vec, $fields[1] );
        push( @rot_mat, [ ($fields[2], $fields[3], $fields[4]) ] );
      }
      last;
    }
  }

  my $out = join("\n", @output);
  return( \@trans_vec, \@rot_mat, $score_ln, $out );
}








########### transform a pdb structure given a transformation ##############
sub transform_pdb {
  my ($inp_file, $out_file, $rot_mat, $trans_vec) = @_;


  open PDB, "<$inp_file" or die "Error: could not open $inp_file\n";
  open OUT, ">$out_file" or die "Error: could not write $out_file\n";
  print OUT "REMARK  Generated by IS-score.\n";
  while (my $line = <PDB>) {
    next if ( $line =~ /^REMARK/ );
    last if ( $line =~ /^ENDMDL/ );  #### only consider the first model

    #### PDB format: http://www.wwpdb.org/documentation/format23/sect9.html ####
    if ( $line =~ /^(ATOM  |HETATM)/ ) {
      (my $x = substr($line,30,8)) =~ s/ //g;
      (my $y = substr($line,38,8)) =~ s/ //g;
      (my $z = substr($line,46,8)) =~ s/ //g;

      unless( defined $x and defined $y and defined $z ) {
        die "Error: cuold not find coordinates\n";
      }

      my $tx = $$trans_vec[0] + $x*$$rot_mat[0][0] + $y*$$rot_mat[0][1] + $z*$$rot_mat[0][2];
      my $ty = $$trans_vec[1] + $x*$$rot_mat[1][0] + $y*$$rot_mat[1][1] + $z*$$rot_mat[1][2];
      my $tz = $$trans_vec[2] + $x*$$rot_mat[2][0] + $y*$$rot_mat[2][1] + $z*$$rot_mat[2][2];

      substr($line,30,8) = sprintf("%8.3f", $tx);
      substr($line,38,8) = sprintf("%8.3f", $ty);
      substr($line,46,8) = sprintf("%8.3f", $tz);

      print OUT $line;
    }
    elsif( $line =~ /^TER/ ) {
      print OUT $line;
    }
  }
  close PDB;
  close OUT;

  return 1;
}


sub print_banner {
  my $outfile = shift @_;
  print $outfile "IS-score version $VERSION\n";
}


##### Time stamp
sub timeStamp {
  my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time);
  my $curr_time = sprintf "%4d-%02d-%02d %02d:%02d:%02d",
    $year+1900,$mon+1,$mday,$hour,$min,$sec;

  return $curr_time;
}


##### generate a VMD script for visualization #####
sub mk_vmdscript {
  my ($script_file, $pdb_file1, $pdb_file2, $trans_vec, $rot_mat, $result, $rep, $tar_resid_to_mod_resid ) = @_;

  my @lines = split(/\n/, $result);

  ######## decide the location of pdb files  #######
  if( not (-e $pdb_file1) and defined $pdb_path ) {
    $pdb_file1 = "$pdb_path/$pdb_file1";
  }

  unless( -s $pdb_file1 ) {
    if( -s "$pdb_file1.pdb" )    { $pdb_file1 = "$pdb_file1.pdb"; }
    elsif( -s "$pdb_file1.ent" ) { $pdb_file1 = "$pdb_file1.ent"; }
    else {
      print "Warning: could not find the pdb file $pdb_file1\n";
    }
  }


  if( not (-e $pdb_file2) and defined $pdb_path ) {
    $pdb_file2 = "$pdb_path/$pdb_file2";
  }

  unless( -s $pdb_file2 ) {
    if( -s "$pdb_file2.pdb" )    { $pdb_file2 = "$pdb_file2.pdb"; }
    elsif( -s "$pdb_file2.ent" ) { $pdb_file2 = "$pdb_file2.ent"; }
    else {
      print "Warning: could not find the pdb file $pdb_file1\n";
    }
  }
  ###--------------------------------------------###



  my @int1 = ( {'chain'=>'-'} );
  my @int2 = ( {'chain'=>'-'} );
  my @mapped_int1 = ( {'chain'=>'-'} );
  my @mapped_int2 = ( {'chain'=>'-'} );
  my $flag = 0;

  ####### find aligned interfacial residues
  foreach (@lines) {
    if( not /^\s*\d+ / ) { $flag = 0; }
    if( /^\s*Index Ch1 Resid1/ ) { $flag = 1; next; }

    next unless ( $flag );

    (my $ch1  = substr( $_, 9,  1)) =~ s/ //g;
    (my $ch2  = substr( $_, 27, 1)) =~ s/ //g;

    (my $res1 = substr( $_, 11, 5)) =~ s/ //g;
    (my $res2 = substr( $_, 29, 5)) =~ s/ //g;

    ##### get a list of all native interfacial residues in the target, and their
    ##### mapped (through sequence alignment) residues in the model
    if( $ch2 ne '' and $res2 ne '' ) {
      if ( $ch2 ne $mapped_int2[-1]->{chain} ) {
	push( @mapped_int2, { 'chain'=>$ch2, 'lst'=>[], } );
      }
      push( @{$mapped_int2[-1]->{lst}}, $res2 );

      my $tar_resid = "$ch2.$res2";
      my $mod_resid = $tar_resid_to_mod_resid->{$tar_resid};
      if( defined $mod_resid ) {
	my @subfields  = split( /\./, $mod_resid );
	my $mapped_ch1 = $subfields[0];
	my $mapped_res1= $subfields[1];
	if ( $mapped_ch1 ne $mapped_int1[-1]->{chain} ) {
	  push( @mapped_int1, { 'chain'=>$mapped_ch1, 'lst'=>[], } );
	}
	push( @{$mapped_int1[-1]->{lst}}, $mapped_res1 );
      }
    }


    next unless( $ch1 ne '' and $ch2 ne '' and $res1 ne '' and $res2 ne '' );


    if ( $ch1 ne $int1[-1]->{chain} ) {
      push( @int1, { 'chain'=>$ch1, 'lst'=>[], } );
    }
    if ( $ch2 ne $int2[-1]->{chain} ) {
      push( @int2, { 'chain'=>$ch2, 'lst'=>[], } );
    }

    push( @{$int1[-1]->{lst}}, $res1 );
    push( @{$int2[-1]->{lst}}, $res2 );
  }

  my $ch1a = $mapped_int1[1]->{chain};
  my $ch1b = $mapped_int1[2]->{chain};
  my $ch2a = $mapped_int2[1]->{chain};
  my $ch2b = $mapped_int2[2]->{chain};

  my $aln1a = join(' ',@{$int1[1]->{lst}});
  my $aln1b = join(' ',@{$int1[2]->{lst}});
  my $aln2a = join(' ',@{$int2[1]->{lst}});
  my $aln2b = join(' ',@{$int2[2]->{lst}});

  my $map1a = join(' ',@{$mapped_int1[1]->{lst}});
  my $map1b = join(' ',@{$mapped_int1[2]->{lst}});
  my $map2a = join(' ',@{$mapped_int2[1]->{lst}});
  my $map2b = join(' ',@{$mapped_int2[2]->{lst}});

  my $vmd_rep;
  if( $rep eq 'cartoon' ) {
    $vmd_rep = 'NewCartoon 0.3 6.0 4.5 0';
  }
  else {
    $vmd_rep = 'NewRibbons 0.3 6.0 3.0 0';
  }

  ###############################################################
  open TCL, ">$script_file" or die "Error: could not write $script_file\n";
  print TCL "#!/usr/local/bin/vmd\n";
  print TCL "### VMD script for visualizing the alignment generated by IS-score\n\n";

  ##### script for loading the first molecule
  print TCL "mol new $pdb_file1 type pdb\n";
  print TCL "mol delrep 0 top\n";

  print TCL "mol representation $vmd_rep\n";
  print TCL "mol color ColorID 10\n";
  print TCL "mol selection {chain $ch1a and same residue as protein within 4.5 of chain $ch1b}\n";
  print TCL "mol material BrushedMetal\n";
  print TCL "mol addrep top\n";
  print TCL "mol showrep top 0 off\n";

  print TCL "mol representation $vmd_rep\n";
  print TCL "mol color ColorID 10\n";
  print TCL "mol selection {chain $ch1a and not same residue as protein within 4.5 of chain $ch1b}\n";
  print TCL "mol material BrushedMetal\n";
  print TCL "mol addrep top\n";
  print TCL "mol showrep top 1 off\n";

  print TCL "mol representation $vmd_rep\n";
  print TCL "mol color ColorID 3\n";
  print TCL "mol selection {chain $ch1b and same residue as protein within 4.5 of chain $ch1a}\n";
  print TCL "mol material BrushedMetal\n";
  print TCL "mol addrep top\n";
  print TCL "mol showrep top 2 off\n";

  print TCL "mol representation $vmd_rep\n";
  print TCL "mol color ColorID 3\n";
  print TCL "mol selection {chain $ch1b and not same residue as protein within 4.5 of chain $ch1a}\n";
  print TCL "mol material BrushedMetal\n";
  print TCL "mol addrep top\n";
  print TCL "mol showrep top 3 off\n";

  print TCL "mol representation $vmd_rep\n";
  print TCL "mol color ColorID 10\n";
  print TCL "mol selection {chain $ch1a and resid $map1a}\n";
  print TCL "mol material BrushedMetal\n";
  print TCL "mol addrep top\n";

  print TCL "mol representation $vmd_rep\n";
  print TCL "mol color ColorID 10\n";
  print TCL "mol selection {chain $ch1a and not resid $map1a}\n";
  print TCL "mol material Transparent\n";
  print TCL "mol addrep top\n";

  print TCL "mol representation $vmd_rep\n";
  print TCL "mol color ColorID 3\n";
  print TCL "mol selection {chain $ch1b and resid $map1b}\n";
  print TCL "mol material BrushedMetal\n";
  print TCL "mol addrep top\n";

  print TCL "mol representation $vmd_rep\n";
  print TCL "mol color ColorID 3\n";
  print TCL "mol selection {chain $ch1b and not resid $map1b}\n";
  print TCL "mol material Transparent\n";
  print TCL "mol addrep top\n";


  #### show aligned Ca in VDW rep
  print TCL "mol representation VDW 0.6 8.0\n";
  print TCL "mol color ColorID 10\n";
  print TCL "mol selection {chain $ch1a and resid $aln1a and name CA}\n";
  print TCL "mol material BrushedMetal\n";
  print TCL "mol addrep top\n";

  print TCL "mol representation VDW 0.6 8.0\n";
  print TCL "mol color ColorID 3\n";
  print TCL "mol selection {chain $ch1b and resid $aln1b and name CA}\n";
  print TCL "mol material BrushedMetal\n";
  print TCL "mol addrep top\n";


  ##### script for transformation
  print TCL "\nset transMat {}\n";
  print TCL "lappend transMat {$$rot_mat[0][0] $$rot_mat[0][1] $$rot_mat[0][2] $$trans_vec[0]}\n";
  print TCL "lappend transMat {$$rot_mat[1][0] $$rot_mat[1][1] $$rot_mat[1][2] $$trans_vec[1]}\n";
  print TCL "lappend transMat {$$rot_mat[2][0] $$rot_mat[2][1] $$rot_mat[2][2] $$trans_vec[2]}\n";
  print TCL "lappend transMat {0.0 0.0 0.0 1.0}\n";
  print TCL "set myMol [atomselect top all]\n";
  print TCL "\$myMol move \$transMat\n\n";


  ##### script for loading the second molecule
  print TCL "mol new $pdb_file2 type pdb\n";
  print TCL "mol delrep 0 top\n";

  print TCL "mol representation $vmd_rep\n";
  print TCL "mol color ColorID 0\n";
  print TCL "mol selection {chain $ch2a and same residue as protein within 4.5 of chain $ch2b}\n";
  print TCL "mol material BrushedMetal\n";
  print TCL "mol addrep top\n";
  print TCL "mol showrep top 0 off\n";

  print TCL "mol representation $vmd_rep\n";
  print TCL "mol color ColorID 0\n";
  print TCL "mol selection {chain $ch2a and not same residue as protein within 4.5 of chain $ch2b}\n";
  print TCL "mol material BrushedMetal\n";
  print TCL "mol addrep top\n";
  print TCL "mol showrep top 1 off\n";

  print TCL "mol representation $vmd_rep\n";
  print TCL "mol color ColorID 1\n";
  print TCL "mol selection {chain $ch2b and same residue as protein within 4.5 of chain $ch2a}\n";
  print TCL "mol material BrushedMetal\n";
  print TCL "mol addrep top\n";
  print TCL "mol showrep top 2 off\n";

  print TCL "mol representation $vmd_rep\n";
  print TCL "mol color ColorID 1\n";
  print TCL "mol selection {chain $ch2b and not same residue as protein within 4.5 of chain $ch2a}\n";
  print TCL "mol material BrushedMetal\n";
  print TCL "mol addrep top\n";
  print TCL "mol showrep top 3 off\n";


  print TCL "mol representation $vmd_rep\n";
  print TCL "mol color ColorID 0\n";
  print TCL "mol selection {chain $ch2a and resid $map2a}\n";
  print TCL "mol material BrushedMetal\n";
  print TCL "mol addrep top\n";

  print TCL "mol representation $vmd_rep\n";
  print TCL "mol color ColorID 0\n";
  print TCL "mol selection {chain $ch2a and not resid $map2a}\n";
  print TCL "mol material Transparent\n";
  print TCL "mol addrep top\n";

  print TCL "mol representation $vmd_rep\n";
  print TCL "mol color ColorID 1\n";
  print TCL "mol selection {chain $ch2b and resid $map2b}\n";
  print TCL "mol material BrushedMetal\n";
  print TCL "mol addrep top\n";

  print TCL "mol representation $vmd_rep\n";
  print TCL "mol color ColorID 1\n";
  print TCL "mol selection {chain $ch2b and not resid $map2b}\n";
  print TCL "mol material Transparent\n";
  print TCL "mol addrep top\n";


  #### show aligned Ca in VDW rep
  print TCL "mol representation VDW 0.6 8.0\n";
  print TCL "mol color ColorID 0\n";
  print TCL "mol selection {chain $ch2a and resid $aln2a and name CA}\n";
  print TCL "mol material BrushedMetal\n";
  print TCL "mol addrep top\n";

  print TCL "mol representation VDW 0.6 8.0\n";
  print TCL "mol color ColorID 1\n";
  print TCL "mol selection {chain $ch2b and resid $aln2b and name CA}\n";
  print TCL "mol material BrushedMetal\n";
  print TCL "mol addrep top\n";

  close TCL;
}
################################################################



######## To obtain the sequence alignment between model and target (two pairs of chains) #########
sub alignModTar{
  my ($mod_seq_of, $mod_pchains, $mod_ind2resid_of, $tar_seq_of, $tar_pchains, $tar_ind2resid_of ) = @_;

  my $mod_seq; my $tar_seq;
  my $mod_ch;  my $tar_ch;
  my $mod_ind2resid;
  my $tar_ind2resid;

  ####### first pair of chains
  $mod_ch  = $$mod_pchains[0];
  $tar_ch  = $$tar_pchains[0];
  $mod_seq = $$mod_seq_of{$mod_ch};
  $tar_seq = $$tar_seq_of{$tar_ch};

  $mod_ind2resid = $$mod_ind2resid_of{$mod_ch};
  $tar_ind2resid = $$tar_ind2resid_of{$tar_ch};

  print "Model $mod_ch vs Target $tar_ch sequence alignment\n";
  my %map1 = getTwoSeqResidMap( $mod_seq, $mod_ind2resid, $tar_seq, $tar_ind2resid );

  ###### second pair of chains
  $mod_ch  = $$mod_pchains[1];
  $tar_ch  = $$tar_pchains[1];
  $mod_seq = $$mod_seq_of{$mod_ch};
  $tar_seq = $$tar_seq_of{$tar_ch};

  $mod_ind2resid = $$mod_ind2resid_of{$mod_ch};
  $tar_ind2resid = $$tar_ind2resid_of{$tar_ch};

  print "Model $mod_ch vs Target $tar_ch sequence alignment\n";
  my %map2 = getTwoSeqResidMap( $mod_seq, $mod_ind2resid, $tar_seq, $tar_ind2resid );

  my %map = ( %map1, %map2 );
  return \%map;
}


##### global sequence alignment  #####
sub getTwoSeqResidMap{
  my ($mod_seq, $mod_ind2resid, $tar_seq, $tar_ind2resid ) = @_;

  my %map = ();

  my ( $aln_mod_seq, $aln_tar_seq, $num_id, $alnlen ) = global($matrix, $mod_seq, $tar_seq);
  print "$aln_mod_seq\n$aln_tar_seq\n";
  my $len = length($aln_mod_seq);
  my $mod_ind = -1;
  my $tar_ind = -1;
  for(my $i=0;$i<$len;$i++) {
    my $mod_aa = substr($aln_mod_seq,$i,1);
    my $tar_aa = substr($aln_tar_seq,$i,1);

    $mod_ind++ if($mod_aa ne '-');
    $tar_ind++ if($tar_aa ne '-');

    next unless( $mod_aa ne '-' and $tar_aa ne '-' );  ### only consider aligned regions

    my $mod_resid = $$mod_ind2resid{$mod_ind};
    my $tar_resid = $$tar_ind2resid{$tar_ind};
    $map{$tar_resid} = $mod_resid;
  }

  #for my $resid (keys %map) { print "$resid $map{$resid}\n"; }

  return %map;
}






################### global sequence alignment ##################
sub global {

 my ($matrix,$x,$y)=@_;
 my ($m, $n)= (length($x),length($y));

 my @M;
  for (my $i=0; $i<= $m; $i++) {
    for (my $j=0; $j <= $n; $j++) {
    $M[$i][$j]= 0;
   }
  }

 my @B;
   for (my $j=0; $j<=$n; $j++) {
     for (my $i=0; $i<= $m; $i++) {
    $B[$i][$j]= 0;
   }
  }

 ## Initialization of M and B matrix
   $M[0][0]=0;
   for (my $j=1; $j<= $n; $j++) {
    $M[0][$j]= $d+$e * ($j-1);
    $B[0][$j]=$LEFT;
   }

   for (my $i=1; $i<=$m; $i++) {
    $M[$i][0]= $d+$e * ($i-1);
    $B[$i][0]=$UP;
   }

  for (my $i=1; $i<= $m; $i++){
    for (my $j=1; $j <= $n; $j++) {

    my $h= $M[$i-1][$j];
       if ($B[$i-1][$j] != $UP) {
         $h=$h+$d; } else {
         $h=$h+$e;
         }
     my $v=$M[$i][$j-1];
       if ($B[$i][$j-1] != $LEFT) {
         $v=$v+$d; } else {
         $v=$v+$e;
         }

    my $s = &score($matrix, substr ($x,$i-1,1), substr ($y,$j-1,1) );
    my $d=$M[$i-1][$j-1] + $s;

    if ($d> $h && $d > $v) {
        $M[$i][$j]=$d;
        $B[$i][$j]=$LEFTUP;
       } elsif ($h >= $d && $h >= $v) {
        $M[$i][$j]=$h; $B[$i][$j]=$UP;
       } else {
         $M[$i][$j]=$v; $B[$i][$j] =$LEFT;
       }
  }
 }

  return (&traceback($x,$y,\@B,$m,$n) );
}

sub traceback {
  my ($x, $y, $B, $i, $j) = @_;
  my ($xAlign, $yAlign) = ("", "");
  my ($m, $n)=(0,0);

    while ($$B[$i][$j] != 0 ) {
    if ($$B[$i][$j] == $UP) {
      $xAlign .= substr($x, $i-1, 1);
      $yAlign .= "-";
      $i--;
    } elsif ($$B[$i][$j] == $LEFT) {
      $yAlign .= substr($y, $j-1, 1);
      $xAlign .= "-";
      $j--;
    } elsif ($$B[$i][$j] == $LEFTUP) {
      $xAlign .= substr($x, $i-1, 1);
      $yAlign .= substr($y, $j-1, 1);
      $n++;  
     if ( substr ($x,$i-1,1) eq substr($y, $j-1, 1) )
      {
        $m++;
      }
     $i--; $j--;
   }
 }

  $xAlign = reverse $xAlign;
  $yAlign = reverse $yAlign;
  return ($xAlign,$yAlign,$m,$n);
}



















 __END__


###################### PODIATRISTS ##############################

=head1 NAME

IS-score - Metrics for assessing protein-protein complex models.

=head1 DESCRIPTION


IS-score is a tool for benchmarking the quality of structural models built for protein-protein complexes.
One can perform pairwise comparison, or compare a set of models against their native target structure. 
Two scoring functions, the Interface Similarity score (IS-score) and the interfacial Template Modeling 
score (iTM-score), are implemented to evaluate the interface similarity between a docking model and its
target. Note that the calculation of the scores requires a sequence correspondence between the model
and the target structure. Usually the two structures have an identical sequence, though gaps are allowed
in the model structure. If you want to compare structures with low/no sequence similarity, consider to
use iAlign.

Details of the IS-score program can be found in the following references:

Mu Gao and Jeffrey Skolnick, Metrics for the automated evaluation of protein-protein docking
models. Proteins, in press, 2011.

Mu Gao and Jeffrey Skolnick. iAlign: a method for the structural comparison of protein-protein
interfaces. Bioinformatics, 26:2259, 2010.

=head1 SYNOPSIS

=over 4

=item
isscore.pl [options] F<model_pdb_file> model_chain_list F<target_pdb_file> target_chain_list

=item
isscore.pl [options] F<model_pdb_file> F<target_pdb_file>

=item
isscore.pl [options] -l F<model_list_file> F<target_pdb_file> target_chain_list

=back

Use the -m option to see the full documentation with examples.

Basic Options:

=over 4

=item
-a  <0/1/2>  0 - no alignment printout, 1 - concise, 2 - detailed

=item
-c1 <string> chain list for a model PDB file

=item
-c2 <string> chain list for a target PDB file


=item
-e  <tm/is>  score selection: tm - TM-score, is - IS-score (default)


=item
-h           print this help message

=item
-l  <file>   PDB list of models


=item
-m           print the manual with examples


=item
-o  <file>   output file for saving the result

=item
-ow <y/n>    always/never overwrite intermediate interface files

=item
-p1 <file>   Model PDB file

=item
-p2 <file>   Target PDB file

=item
-p  <path>   path to PDB files

=item
-t           transform PDB1 onto PDB2 according to the optimal alignment

=item
-w  <path>   workpath or path to parsed PDB files

=item
-vmd <style> generate VMD script in ribbon or cartoon style

=item
-v           print the version number

=back



Advanced Options:

=over 4

=item
-ca          calculate score using only Calpha atoms (for low resolution models)

=item
-dc   <value> distance cutoff for an interfacial contact, default 4.5 A

=item
-minp <value> minimum number of residues for a protein chain

=item
-mini <value> minimum number of residues for an interface

=item
-seqonce    only perform sequence alignment once

=back



=head1 OPTIONS

=over 4

=item B<-a, -aln>

Choose a printout format for alignment. 0 - no printout, 1 - concise (sequential mode only), 2 -
detailed (default).

=item B<-c1, -chains1>

Specify two protein chains in a model structure, e.g., "AB", without any space between chain IDs.
If more than two chains are specified, the program will consider the first two chains.
If the list of chains is not provided, the first two chains found in the PDB record will be considered .

=item B<-c2, -chains2>

Specify a list of protein chains in a target structure. Usage is the same as described in c1.

=item B<-ca>

This option is useful for low resolution models with only Calpha atoms, which are used to estimate
the coordinates of sidechain center of mass. The resulting Ca-SCOM models are used to calculate
interfacial contacts according to optimized residue-residue-specific contact cutoffs. Note that
the option dc is ineffective when the option ca is enabled.


=item B<-dc>

Specify the distance cutoff (in Angstrom) for interfacial contact. If a pair of heavy atoms of two
residues from two separate proteins has a distance less than the cutoff, the two residues form
an interfacial contact. The default is 4.5 Angstroms. Note that the statistical significance is
estimated based on the 4.5 A cutoff value. This option is ineffective if the option -ca is enabled.


=item B<-e, -measure>

Select a scoring function for measuring interface similarity: "is" IS-score (default) and "tm" TM-score.


=item B<-h, -help>

show the brief help information.

=item B<-l, -liblst1>

Specify a file that contains a list of models. The format of the list is:

=over 12

=item pdbfile1 <chainIDs>

=item pdbfile2 <chainIDs>

=item ...

=back

ChainIDs are optional. If not specified, the first two chains found will be examined. One can give full path
in "pdbfile". If not, ialign will search specified pdb files under the directory specified by the
option -pdbpath, or the current directory otherwise.


=item B<-m, -manual>

Show the manual with examples.

=item B<-minp>

Specify the minimum number of amino acids for a peptide chain to be considered as a valid protein.
Peptide chains with AAs lower than this number will be ignored. The default is 25 AAs.

=item B<-mini>

Specify the minimum number of amino acids for a protein-protein interface to be considered as a valid.
Interfaces with AAs lower than this number will be ignored. The default is 20 AAs.

=item B<-nr <cutoff>>

Remove redundant interfaces within a PDB record. The criterion for redundancy is defined by the
cutoff value for the TM-score.


=item B<-o, -output>

Specify the output file for saving alignment results.

=item B<-ow>

Always/never overwrite intermediate interface files, the default is to ask if the file exists.

=item B<-p1, -pdb1>

Specify the PDB file of a model.

=item B<-p2, -pdb2>

Specify the PDB file of a target.

=item B<-p, -pdbpath>

Specify the Path to unprocessed PDB files.

=item B<-seqonce>

Only perform sequence alignment once. When comparing many models with identical sequence, use
the sequence alignment between the first model and the target for all models.

=item B<-t, -trans>

Superpose pdb1 onto pdb2 using the transformation matrix that generates the optimal interface
alignment. The output file is named as "pdb1chain1achain1b_to_pdb2chain2achain2b.pdb".

=item B<-w, -workpath>

Specify the path to processed PDB files and input interface files.

=item B<-vmd <style>>

Generate TCL script for visualization with VMD. Two styles are available: ribbon and cartoon.
Aligned residues are shown in VDW representations.

=item B<-v, -version>

Show the version number of this release.




=back




=head1 EXAMPLES

=over 4

=item Example 1

Compare a model (chains A and B of "model.pdb") against its target (chain C and D of "target.pdb").
The second command also transforms the model according to the optimal alignment.

isscore.pl model.pdb AB target.pdb CD

isscore.pl -t -p1 foo1.pdb -c1 AB -p2 foo2.pdb -c2 CD



=item Example 2

Scan the target interface "target.pdb" against all model interfaces listed in goo.lst.
Parsed pdb and extracted interface files are saved to the working directory "/scratch/isscore".

isscore.pl -w /scratch/isscore  -l goo.lst target.pdb


=back


=head1 AUTHOR

Mu Gao  <mu.gao@gatech.edu>.

=head1 DATE

15-Feb-2011

=cut

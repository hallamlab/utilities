#!/usr/bin/perl -wT
use strict;

# Written by Charles Howes 2008-10-15

use DB_File;
use Getopt::Long;

# Where am I?
my $HOME=$0;
$HOME=~s/\/[^\/]*$//;

my @LISTINPUT;
my @FASTAINPUT;
my @BLASTINPUT;
my $OUTPUTDIR=".";
my $SEQNAME;
my $REGEXP;
my $DUPLICATES;

# Process command line arguments
my $result=GetOptions(
  "outputdir:s"=>\$OUTPUTDIR,
  "expression:s"=>\$SEQNAME,
  "regexp"=>\$REGEXP,
  "duplicates"=>\$DUPLICATES,
  );

# Usage, or bad arguments:
if (!defined $ARGV[0] or !$result) {
  die("Usage: $0 [-d] [-o outputdir] [-r] [-e seq] *.list *.fastainput *.blastinput
  This program expands one or more sequences (.list files or
  command-line sequences) into fasta files or
  blast files containing only those sequences.

  The input filenames are not important, as the contents of the files are
  examined to determine their type.

  A list file is a text file containing lines of the format:
    >sequence1
    >sequence2
    >sequence42_life_the_universe_and_everything

  A fasta file is a text file containing lines of the format:
    >sequence1
    ASDFASDFASDFASDFASDFELVISLIVESASDFASDFASDF
    ASDFASDFASDFASDFASDFELVISLIVESASDFASDFASDF
    >sequence2
    WERTWERTWERTWERTWERTWERTWERTWERTWERTWERTW

  A blast file is a text file starting with 'BLASTP', 'BLASTN', 'Query= ' or any string
  that matches the regular expression '^T?BLAST[PXN] '.

  The output files will be named the same as the input list files, with
  '.fas' or '.blast' on the end.  The blast header from the first input blast file
  will be used as the header for all output blast files.

  -e seq: find and print this sequence to stdout.
  -r: The sequence from -e is a regular expression (in quotes).
  -d: print duplicate occurrances of a sequence, not just the first
  -o outputdir: write the files to outputdir (default: $OUTPUTDIR)
"
  );
}

if ($OUTPUTDIR!~m/^([a-z0-9.\/_+-]+)$/i) {die("The outputdir option had bad characters in it.");}
$OUTPUTDIR=$1;
if (!-e $OUTPUTDIR) {mkdir($OUTPUTDIR);}
if (!-e $OUTPUTDIR) {die("$OUTPUTDIR doesn't exist and can't be created\n");}
if (!-r $OUTPUTDIR) {die("$OUTPUTDIR isn't readable\n");}
if (!-d $OUTPUTDIR) {die("$OUTPUTDIR isn't a directory\n");}
if (!-x $OUTPUTDIR) {die("$OUTPUTDIR isn't cd-able\n");}

# Classification of input files:
foreach my $IN (@ARGV) {
  if ($IN!~m/^([\/a-z0-9._+ -]+)$/i) {die("$IN has bad characters\n");}
  $IN=$1;
  if (!-e $IN) {die("$IN doesn't exist\n");}
  if (!-f $IN) {die("$IN isn't a file\n");}
  if (!-r $IN) {die("$IN isn't readable\n");}

  # Read 2 lines from $IN:
  my @c;
  open(IN,"<$IN") or die("$IN: $!");
  while (<IN>) { push @c,$_; if (@c>1) {last;} }
  close(IN);

  # Check for short files:
  if (!defined $c[0]) {next;} #Skip empty files
  if (!defined $c[1]) {$c[1]=$c[0];}  # Turn a list file with 1 line into 2 lines

  # Recognize the files:
  if ($c[0]=~m/^(T?BLAST[PXN] |Query= )/) { push @BLASTINPUT,$IN; next; }
  if ($c[0]=~m/^>/ and $c[1]!~m/^>/) { push @FASTAINPUT,$IN; next; }
  if ($c[0]=~m/^>/ and $c[1]=~m/^>/) { push @LISTINPUT,$IN; next; }

  die("$IN was not recognized as a list file, a fasta file, or a blast output file!\n");
}

my %BLASTINPUT=map {$_=>1} @BLASTINPUT;
my %FASTAINPUT=map {$_=>1} @FASTAINPUT;
#print "Lists: @LISTINPUT\n";
#print "BLAST files: @BLASTINPUT\n";
#print "FASTA files: @FASTAINPUT\n";
#exit 0;

if (defined $SEQNAME) { push @LISTINPUT,""; }

# Error checking:
if (@LISTINPUT==0) {die("No list input files were found!");}
if (@BLASTINPUT==0 and @FASTAINPUT==0) {die("No blast or fasta input files were found!");}

# Indexing fasta files:
my %fasta;
#tie %fasta, 'DB_File', "fasta-index.db";
foreach my $ff (@FASTAINPUT) { 
  if (defined $fasta{"--$ff--loaded"}) {next;}
  findfasta($ff);
  $fasta{"--$ff--loaded"}=1;
}

# Indexing blast files:
my %blast;
#tie %blast, 'DB_File', "blast-index.db";
foreach my $bf (@BLASTINPUT) {
  if (defined $blast{"--$bf--loaded"}) {next;}
  findblast($bf);
  $blast{"--$bf--loaded"}=1; }

# Trim the blast header entry: if (defined $blast{"header"}) { $blast{"header"}=~s/.*[|]/|/; } # Cut after first

# Expanding lists
foreach (@LISTINPUT) { dolist($_); }

exit 0;

# Find individual fasta lines
sub findfasta {
  my $FILE=$_[0];
  open(IN,"<$FILE") or die("$FILE: $!");
  my $pos=0;
  my @seq;
  my $t=0;
  while (<IN>) {
    if (m/^>(\S+)/) {

      # Record the previous sequence under every name given:
      if (@seq>0) { foreach my $s (@seq) {$fasta{$s}.="|$FILE $pos $t";} };

      # Start the next sequence:
      @seq=();
      my @f=split(/[>\001]/,$_);
      foreach my $s (@f) { if ($s=~m/^(\S+)/) {push @seq,$1;} }
      $pos=$t;
    }
    $t=tell(IN); #The start of the next record
  }

  # Take care of the last sequence:
  if (@seq>0) { foreach my $s (@seq) {$fasta{$s}.="|$FILE $pos $t";} };
  close(IN);
}

# Find individual blast lines
sub findblast {
  my $FILE=$_[0];
  open(IN,"<$FILE") or die("$FILE: $!");
  my $pos=0;
  my $seq;
  my $t=0;
  while (<IN>) {
    # Sequences to keep:
    if (m/^Query= (\S+)/) { if (defined $seq) { $blast{$seq}.="|$FILE $pos $t"; }; $seq=$1; $pos=$t;}

    # Record one header (thus =, not .=):
    if (m/^T?BLAST[PXN] /) { if (defined $seq) { $blast{$seq}.="|$FILE $pos $t"; }; $seq="header"; $pos=$t;}

    # Stuff to delete:
    if (m/^  Database:/) { if (defined $seq) { $blast{$seq}.="|$FILE $pos $t"; }; $seq=undef; $pos=$t;}

    $t=tell(IN); #The start of the next record
  }
  if (defined $seq) { $blast{$seq}.="|$FILE $pos $t"; };
  close(IN);
}

# Process a list file
sub dolist {
  my $FILE=$_[0];
  my $FASTAFILE;
  my $BLASTFILE;

  # Take care of -e and -r from command line, and writing to stdout:
  if ($FILE eq "") {
    my $good=0;

    if (defined $REGEXP) {
      foreach my $k (sort keys %fasta) {
        if ($k=~m/$SEQNAME/) { $good++; processfasta($k,$fasta{$k},*STDOUT); }
      }
      my $head=0;
      foreach my $k (sort keys %blast) {
        if ($k=~m/$SEQNAME/) {
	  $good++;
	  if (!$head) {$head=1;processblast($blast{"header"},*STDOUT);}
	  processblast($blast{$k},*STDOUT);
	}
      }
      if (!$good) {die("Sequence $SEQNAME was not found in any file.\n");}
      return;
    }

    if (defined $fasta{$SEQNAME}) {
      processfasta($SEQNAME,$fasta{$SEQNAME},*STDOUT);
      $good++;
    }
    if (defined $blast{$SEQNAME}) {
      processblast($blast{"header"},*STDOUT);
      processblast($blast{$SEQNAME},*STDOUT);
      $good++;
    }
    if (!$good) {die("Sequence $SEQNAME was not found in any file.\n");}
    return;
  }

  # Normal processing of input list files:
  print "Processing $FILE\n";
  open(IN,"<$FILE") or die("$FILE: $!");
  while(<IN>) {
    if (!m/^>(\S+)/) {die("$FILE contains a bad line: $_");}
    my $seq=$1;

    my $good=0;

    if (defined $fasta{$seq}) {
      $good=1;
      if (!defined $FASTAFILE) {
        $FASTAFILE="$FILE.fas";
	$FASTAFILE=~s/.*\//$OUTPUTDIR\//;
	open(FASTA,">$FASTAFILE") or die("$FASTAFILE: $!");
      }
      processfasta($seq,$fasta{$seq},*FASTA);
    }

    if (defined $blast{$seq}) {
      $good=1;
      if (!defined $BLASTFILE) {
        $BLASTFILE="$FILE.blast";
	$BLASTFILE=~s/.*\//$OUTPUTDIR\//;
	open(BLAST,">$BLASTFILE") or die("$BLASTFILE: $!");

	# Copy a blast header from the first source file:
	processblast($blast{"header"},*BLAST);
      }
      processblast($blast{$seq},*BLAST);
    }
    if (!$good) {warn("Sequence $seq in file $FILE was not found.\n");}
  }
  close(IN);
  close(BLAST);
  close(FASTA);
}

# Process a result string.
# A result string is a series of |-separated locations
# each location consists of a filename, a start and an end position
# Read the location's data and write it to the file handle.
sub processblast {
  my $in=$_[0]; # Input filename, start and end
  my $fh=$_[1]; # Output filehandle
  $in=~s/^[|]//;
  if (!defined $DUPLICATES) {$in=~s/\|.*//;}
  my @f=split(/\|/,$in);
  foreach my $m (@f) {
    if ($m!~m/^(.*) (\d+) (\d+)$/) {die("$m doesn't match pattern.\n");}
    my ($file,$start,$end)=($1,$2,$3);
    if (!defined $BLASTINPUT{$file}) {next;} # This file was not requested
    open(IN2,"<$file") or die("$file: $!");
    my $r;
    seek(IN2,$start,0);
    my $result=read(IN2,$r,($end-$start));
    close(IN2);
    print $fh $r;
  }
}

# Process a fasta result string.
# A sequence is the first word from the fasta header.
# A result string is a series of |-separated locations
# each location consists of a filename, a start and an end position
# Read the location's data and write it to the file handle.
sub processfasta {
  my $sq=$_[0]; # Sequence we're looking for
  my $in=$_[1]; # locations
  my $fh=$_[2]; # Output filehandle
  $in=~s/^[|]//;
  if (!defined $DUPLICATES) {$in=~s/\|.*//;}
  my @f=split(/\|/,$in); # Split the locations

  foreach my $m (@f) { # For each location:

    if ($m!~m/^(.*) (\d+) (\d+)$/) {die("$m doesn't match pattern.\n");}
    my ($file,$start,$end)=($1,$2,$3);
    if (!defined $FASTAINPUT{$file}) {next;} # This file was not requested
    open(IN2,"<$file") or die("sequence $sq, $in: $!");
    my $r;
    seek(IN2,$start,0);
    my $result=read(IN2,$r,($end-$start));
    if (!defined $result or $result == 0) {die("$file: $!");}
    close(IN2);

    # print only the >fasta header for the sequence we asked for
    my @lines=split(/\n/,$r); # Split into lines
    my @headers=split(/[>\001]/,$lines[0]); # Split the first line
    my $count=0;
    my $myseq=quotemeta($sq);
    foreach my $h (@headers) {
      if ($h=~m/^$myseq([^0-9a-z]|$)/i) { # Found a match
        $count++;
        print $fh ">$h\n".join("",@lines[1..$#lines])."\n";
      }
    }
    if ($count==0) {die("Sequence $sq: $m did not have a matching header?");}

  }
}

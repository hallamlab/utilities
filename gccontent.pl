#!/usr/bin/perl 
#===============================================================================
#
#         FILE:  gccontent.pl
#
#        USAGE:  ./gccontent.pl inputfile.fas > outputfile.csv
#
#  DESCRIPTION:  Calculate the gc content of a fasta file
#
#      OPTIONS:  ---
# REQUIREMENTS:  ---
#         BUGS:  ---
#        NOTES:  ---
#       AUTHOR:  Charles Howes (), chowes+vimperlsupport@interchange.ubc.ca
#      COMPANY:  Hallam Lab, University of British Columbia
#      VERSION:  1.0
#      CREATED:  11/16/2009 09:55:05 AM
#     REVISION:  ---
#===============================================================================

use strict;
use warnings;

if (!@ARGV) {
  die("Usage: $0 input.fas [inputs.fas...] > output.csv

  This program will read each input file, and for each sequence will
  calculate the GC ratio.
");
}
foreach my $file (@ARGV) { doit($file); }

sub doit {
  my $INFILE=$_[0];
  open(IN,"<$INFILE") or die("$INFILE: $!");
  my ($name,$seq)=("","");
  while(<IN>) {
    chomp;
    if (m/^>(\S+)/) {process($name,$seq,$INFILE);$name=$1;$seq="";next;}
    $seq.=$_;
  }
  process($name,$seq);
}

sub process {
  my ($n,$s,$f)=@_;
  if ($s ne "" and $n eq "") {die("Got sequence but no name; is $f a fasta file?");}
  if ($n eq "") {return;}  # Normal start
  if ($s eq "") {die("Got a name but no sequence; is $f a fasta file?")}
  
  my %count;
  $s=lc($s);
  map {$count{substr($s,$_,1)}++} 0..(length($s)-1);
  my $seen=join("",sort keys %count);
  if ($seen!~m/^\s*a?c?g?t?$/) {
    foreach (grep {!/[actg]/} sort keys %count) {
      warn("File $f, sequence ${n}: found $count{$_} '$_' characters\n");
    }
  }
  my $gc=0;
  my $total=0;
  map {$gc+=$count{$_}} qw(c g);
  map {$total+=$count{$_}} qw(a c g t);
  if ($total eq 0) {warn("File $f, sequence ${n}: no valid nucleotides found?");}
  my $gcratio=sprintf("%.3f",$gc/$total);
  print "$n\t$gcratio\n";
}

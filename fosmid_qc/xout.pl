#!/usr/bin/perl 
#===============================================================================
#
#         FILE:  xout.pl
#
#        USAGE:  ./xout.pl [-i identity] -o output.fas -d deleteme.fas input.fas
#
#  DESCRIPTION:  Put x's over top of matches in the given fasta file
#
#      OPTIONS:  ---
# REQUIREMENTS:  ---
#         BUGS:  ---
#        NOTES:  ---
#       AUTHOR:  Charles Howes (), chowes+vimperlsupport@interchange.ubc.ca
#      COMPANY:  Hallam Lab, University of British Columbia
#      VERSION:  1.0
#      CREATED:  2010-11-05 14:04:01
#     REVISION:  ---
#===============================================================================

use strict;
use warnings;
use Getopt::Long;

my $CSVFILE;
my $OUTFILE;
my $CONTIGFILE;
my $EXCLUDE;
my $PICTURE=1;
my $DELFILE;
my $IDENTITY=100;

my $result=GetOptions(
  "deleteme=s"=>\$DELFILE,  # File to find matches
  "identity:i"=>\$IDENTITY,  # Percent identity to match at, default 100%
  "output=s"=>\$OUTFILE,  # Output fasta file
);
my $INFILE=$ARGV[0];

if (!defined $OUTFILE or !defined $DELFILE or @ARGV!=1) {
  die("Usage: $0 [-i identity] -o output.fas -d deleteme.fasta input.fasta
  This program will blast the two fasta files together and replace the
  hit nucleotides with 'X's, much like crossmatch.
  
  Options: 
    -o output.fas: the output file
    -d deleteme.fas: the sequences to X-out
    -i identity: What percent identity to match (default 100)

");
}


if(!-f $INFILE || !-r _ ) { die("$INFILE isn't a readable file!");}

my $db=formatdb($DELFILE);
my $res=blast($INFILE,$db);
createfasta($OUTFILE,$res);

sub formatdb { # formatdb(input file)=output file
  my $f=$_[0];
  my $s=`md5sum $f`;
  if ($s!~m/^([0-9a-f]{32})\s.*/) {die("md5sum reported:\n$s\n ");}
  $s="cache/blastdb.$1";
  if (! -e "$s.nhr") {
    mkdir("cache");
    my $r=`makeblastdb -in $f -out $s -dbtype nucl -title $f 2>&1`;
    print "formatdb said:\n$r\n";
  }
  return $s;
}


sub blast { # blast(Input file, database file prefix)=output file
  my ($f,$d)=@_;
  my $s=`cat $f $d.* | md5sum`;
  if ($s!~m/^([0-9a-f]{32})\s.*/) {die("md5sum reported:\n$s\n ");}
  $s="cache/blastout.$1";
  if (! -e $s) {
    mkdir("cache");
    my $r=`blastn -evalue 1e-6 -perc_identity $IDENTITY -outfmt 7 -db $d -query $f -out $s 2>&1`;
    print "blastn said:\n$r\n";
  }
  return $s;
}

sub createfasta {
  my %tocut;
  my ($file,$result)=@_;
  open(IN2,"<$result") or die("$result: $!\n");
  while(<IN2>) {
    if (m/^#/){next;}
    my @f=split(/\t/,$_); 

    # 0=query id, 1=subject id, 2=% identity, 3=alignment length,
    # 4=mismatches, 5=gap opens, 6=q. start, 7=q. end, 8=s. start, 9=s. end,
    # 10=evalue, 11=bit score

    # F07-k70-contigs.fa.cap-Contig1    pCC1Fos 100.00  2035    0       0       8204    10238   361     2395    0.0     3759
    $tocut{$f[0]}.=" $f[1]:$f[6]-$f[7]";
  }
  close(IN2);

  open(OUT,">$file") or die("$file: $!\n");
  open(IN1,"<:raw:eol(LF)",$INFILE) or die("$INFILE: $!\n");
  do {
    local $/="\n>";
    while(<IN1>) {
      s/^>?/>/s;  s/>$//;  s/\n$//s;
      if (m/^(>(\S+)[^\n]*\n)(.*)$/is) {
        my ($fullname,$name,$seq)=($1,$2,$3);
        $seq=cutup($seq,$tocut{$name});
        if ($seq ne "") { print OUT $fullname,$seq,"\n"; } # If it wasn't all vector
        next;
      }
      die("Word!?");
    }
  };
  close(IN1);
  close(OUT);
}

sub cutup {
  my ($seq,$cuts)=@_;
  if (!defined $cuts) {return $seq;}
  $seq=~s/\s+//sg; # No whitespace
  $cuts=~s/^\s*//; # Delete leading space so split works right
  foreach my $x (split(/ /,$cuts)) {
    my ($vname,$start,$end)=$x=~m/(.+):(\d+)-(\d+)$/si;
    substr($seq,$start-1,$end-$start+1,"X"x($end-$start+1));
  };
  $seq=~s/^[NX]*//s; $seq=~s/[NX]*$//s;
  return $seq;
}

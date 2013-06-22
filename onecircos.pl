#!/usr/bin/perl 
#===============================================================================
#
#         FILE:  tocircos.pl
#
#        USAGE:  ./onecircos.pl  -o output.svg fasta1
#
#  DESCRIPTION:  blasts fasta1 to fasta2 and draws links with circos
#
#      OPTIONS:  ---
# REQUIREMENTS:  ---
#         BUGS:  ---
#        NOTES:  ---
#       AUTHOR:  Charles Howes (), chowes+vimperlsupport@interchange.ubc.ca
#      COMPANY:  Hallam Lab, University of British Columbia
#      VERSION:  1.0
#      CREATED:  2010-11-02 08:44:48
#     REVISION:  ---
#===============================================================================

use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;

my $CSVFILE;
my $OUTFILE;
my $CONTIGFILE;
my $EXCLUDE;
my $PICTURE=1;
my $BSR;
my $SELF;
my $TASK="megablast";
my $PERC_IDENTITY=100;
my $EVALUE=1e-3;
my @TO_DELETE;

my $result=GetOptions(
  "bsr:f"=>\$BSR,
  "evalue:f"=>\$EVALUE,
  "self"=>\$SELF,
  "task:s"=>\$TASK,
  "perc_identity:f"=>\$PERC_IDENTITY,
  "output=s"=>\$OUTFILE,  # Output file for circos input
);
my $INFILE1=$ARGV[0];

if (!defined $OUTFILE or @ARGV!=1) {
  die("Usage: $0 [-s] [-b bsr] -o circos.svg input.fasta
  This program will blast the input fasta file against itself and draw
  a circos image showing matches.
    -s: show self-blast hits
    -p 100: set the percent_identity filter for blastn
    -e 1e-3: set the evalue filter for blastn
    -t megablast: set the task (blastn, blastn-short, megablast, dc-megablast, vecscreen)
    -b bsr: show only links with bsr values greater than this percentage (default: show all)
    -o output.svg: write the svg output to this file
");
}

my $MAKEBLASTDB=@{[map({"$_/makeblastdb"} grep({-e "$_/makeblastdb"} (split(/:/,$ENV{'PATH'}),".")))]}[0];
if (!defined $MAKEBLASTDB) {die("The program makeblastdb was not found in your path\n");}
if(!-f $MAKEBLASTDB || !-x _) { die("$MAKEBLASTDB isn't an executable file!");}

my $BLASTN=@{[map({"$_/blastn"} grep({-e "$_/blastn"} (split(/:/,$ENV{'PATH'}),".")))]}[0];
if (!defined $BLASTN) {die("The program blastn was not found in your path\n");}
if(!-f $BLASTN || !-x _) { die("$BLASTN isn't an executable file!");}

my $MD5SUM=@{[map({"$_/md5sum"} grep({-e "$_/md5sum"} (split(/:/,$ENV{'PATH'}),".")))]}[0];
if (!defined $MD5SUM) {die("The program md5sum was not found in your path\n");}
if(!-f $MD5SUM || !-x _) { die("$MD5SUM isn't an executable file!");}

my $CIRCOS=@{[map({"$_/circos"} grep({-e "$_/circos" && -f _ && -r _ && -x _} (split(/:/,$ENV{'PATH'}),".")))]}[0];
if (!defined $CIRCOS) {die("The program circos was not found in your path\n");}

if(!-f $INFILE1 || !-r _ ) { die("$INFILE1 isn't a readable file!");}

if (defined $BSR) {
  if ($BSR<0 or $BSR>100) {die("-b bsr: the value should be in the range 0-100");}
  if ($BSR>1) {$BSR/=100;} # If not in the range 0-1, make it so.
}

my $db=formatdb($INFILE1);
my $res=blast($INFILE1,$db);
my %bsr;
if (defined $BSR) {%bsr=map {bl2seqfiles($_)} @{[$INFILE1]};}
createcircos($OUTFILE,$res,\%bsr);
unlink(@TO_DELETE);

sub formatdb { # formatdb(input file)=output file
  my $f=$_[0];
  my $s=`$MD5SUM $f`;
  if ($s!~m/^([0-9a-f]{32})\s.*/) {die("md5sum reported:\n$s\n ");}
  $s="cache/blastdb.$1";
  if (! -e "$s.nhr") {
    if (! -e "cache") {mkdir("cache");}
    my $cmd="$MAKEBLASTDB -in $f -out $s -dbtype nucl -title $f";
    my $r=`$cmd 2>&1`;
	if (!defined $r) {die("The program $MAKEBLASTDB was not found in
	your path\n");}
	if (! -e "$s.nin") { die("$MAKEBLASTDB failed to create formatted db $s.n*:\ncommand: $cmd\nresult: $r\n"); }
    #print "formatdb said:\n$r\n";
  }
  return $s;
}

sub blast { # blast(Input file, database file prefix)=output file
  my ($f,$d)=@_;
  my $cmd="$BLASTN -task $TASK -evalue $EVALUE -perc_identity $PERC_IDENTITY -outfmt 7 -db $d -query $f -out ";
  my $s=`( echo \Q$cmd\E; cat $f $d.* ) | $MD5SUM`;
  if ($s!~m/^([0-9a-f]{32})\s.*/) {die("md5sum reported:\n$s\n ");}
  $s="cache/blastout.$1";
  if (! -e $s) {
    if (! -e "cache") {mkdir("cache");}
    $cmd="$cmd $s";
    print "$cmd\n";
    my $r=`$cmd 2>&1`;
	if (! defined $r) {$r="The program $BLASTN was not found in your path\n";}
	if (! -e $s) { die("$BLASTN failed to create output file $s:\ncommand: $cmd\nresult: $r\n"); }
  }
  return $s;
}

sub createcircos {
  my ($outfile,$result,$bsrref)=@_;
  open(OUT,">temp.circos.$$") or die("temp.circos.$$: $!\n");
  push @TO_DELETE,"temp.circos.$$";
  print OUT "<colors>\n";
  print OUT "<<include etc/colors.conf>>\n";
  print OUT "blackweak = 0,0,0,100\n";
  print OUT "</colors>\n";
  print OUT "\n";
  print OUT "<fonts>\n";
  print OUT "<<include etc/fonts.conf>>\n";
  print OUT "</fonts>\n";
  print OUT "\n";
  print OUT "<<include ideogram.conf>>\n";
  print OUT "<<include ticks.conf>>\n";
  print OUT "\n";
  print OUT "karyotype = temp.karyotype.$$\n";
  print OUT "\n";
  print OUT "<image>\n";
  print OUT "dir = .\n";
  print OUT "file = $outfile\n";
  print OUT "24bit = yes\n";
  print OUT "png  = no\n";
  print OUT "svg  = yes\n";
  print OUT "# radius of inscribed circle in image\n";
  print OUT "radius         = 500p\n";
  print OUT "background     = white\n";
  print OUT "# by default angle=0 is at 3 o'clock position\n";
  print OUT "angle_offset   = -90\n";
  print OUT "\n";
  print OUT "auto_alpha_colors = yes\n";
  print OUT "auto_alpha_steps  = 5\n";
  print OUT "\n";
  print OUT "</image>\n";
  print OUT "\n";
  print OUT "chromosomes_units  = 1000000\n";
  print OUT "\n";
  print OUT "# show all chromosomes\n";
  print OUT "chromosomes_display_default = yes\n";
  print OUT "\n";
  print OUT "<links>\n";
  print OUT "\n";
  print OUT "z      = 0\n";
  print OUT "radius = 0.98r\n";
  print OUT "crest  = 1\n";
  print OUT "ribbon = yes\n";
  print OUT "color  = black\n";
  print OUT "bezier_radius        = 0.2r\n";
  print OUT "bezier_radius_purity = 0.5\n";
  print OUT "\n";
  print OUT "<link segdup>\n";
  print OUT "thickness        = 5\n";
  print OUT "stroke_color     = vdgrey\n";
  print OUT "stroke_thickness = 0\n";
  print OUT "file             = temp.links.$$\n";
  print OUT "\n";
  print OUT "<rules>\n";
  print OUT "\n";
  print OUT "# set z-depth based on size\n";
  print OUT "#<rule>\n";
  print OUT "#importance = 100\n";
  print OUT "#condition  = 1\n";
  print OUT "#z = eval( scalar min(_SIZE1_,_SIZE2_) )\n";
  print OUT "#</rule>\n";
  print OUT "\n";
  print OUT "# add transparency to color by suffixing\n";
  print OUT "# color value with _a4\n";
  print OUT "<rule>\n";
  print OUT "importance = 100\n";
  print OUT "condition  = 1\n";
  print OUT "color = eval( _color_ .\"_a4\")\n";
  print OUT "</rule>\n";
  print OUT "\n";
  print OUT "</rules>\n";
  print OUT "\n";
  print OUT "</link>\n";
  print OUT "\n";
  print OUT "</links>\n";
  print OUT "\n";
  print OUT "anglestep       = 0.1\n";
  print OUT "minslicestep    = 2\n";
  print OUT "beziersamples   = 40\n";
  print OUT "debug           = no\n";
  print OUT "warnings        = no\n";
  print OUT "imagemap        = no\n";
  print OUT "\n";
  print OUT "# don't touch!\n";
  print OUT "units_ok        = bupr\n";
  print OUT "units_nounit    = n\n";
  close(OUT);

  if (!-e "ideogram.conf") {
	open(OUT,">ideogram.conf") or die("ideogram.conf: $!");
    print OUT "<ideogram>\n";
    print OUT "\n";
    print OUT "<spacing>\n";
    print OUT "\n";
    print OUT "#default = 10u\n";
    print OUT "default = 0.01u\n";
    print OUT "break   = 0.01u\n";
    print OUT "\n";
    print OUT "#axis_break_at_edge = yes\n";
    print OUT "#axis_break         = yes\n";
    print OUT "#axis_break_style   = 2\n";
    print OUT "\n";
    print OUT "<break_style 1>\n";
    print OUT "stroke_color = black\n";
    print OUT "fill_color   = blue\n";
    print OUT "thickness    = 0.25r\n";
    print OUT "stroke_thickness = 2\n";
    print OUT "</break>\n";
    print OUT "\n";
    print OUT "<break_style 2>\n";
    print OUT "stroke_color     = black\n";
    print OUT "stroke_thickness = 3\n";
    print OUT "thickness        = 1.5r\n";
    print OUT "</break>\n";
    print OUT "\n";
    print OUT "</spacing>\n";
    print OUT "\n";
    print OUT "# thickness (px) of chromosome ideogram\n";
    print OUT "#thickness        = 100p\n";
    print OUT "thickness        = 10p\n";
    print OUT "stroke_thickness = 0\n";
    print OUT "# ideogram border color\n";
    print OUT "stroke_color     = black\n";
    print OUT "fill             = yes\n";
    print OUT "# the default chromosome color is set here and any value\n";
    print OUT "# defined in the karyotype file overrides it\n";
    print OUT "fill_color       = black\n";
    print OUT "\n";
    print OUT "# fractional radius position of chromosome ideogram within image\n";
    print OUT "radius         = 0.85r\n";
    print OUT "show_label     = yes\n";
    print OUT "label_font     = default\n";
    print OUT "label_radius   = dims(ideogram,radius) + 0.01r\n";
    print OUT "#label_size     = 72\n";
    print OUT "label_size     = 12\n";
    print OUT "\n";
    print OUT "# cytogenetic bands\n";
    print OUT "band_stroke_thickness = 1\n";
    print OUT "\n";
    print OUT "# show_bands determines whether the outline of cytogenetic bands\n";
    print OUT "# will be seen\n";
    print OUT "show_bands            = yes\n";
    print OUT "# in order to fill the bands with the color defined in the karyotype\n";
    print OUT "# file you must set fill_bands\n";
    print OUT "fill_bands            = yes\n";
    print OUT "\n";
    print OUT "</ideogram>\n";
    print OUT "\n";
	close(OUT);
    print STDERR "ideogram.conf created\n";
  }

  if (!-e "ticks.conf") {
	open(OUT,">ticks.conf") or die("ticks.conf: $!");
    print OUT "show_ticks          = yes\n";
    print OUT "show_tick_labels    = yes\n";
    print OUT "\n";
    print OUT "grid_start         = dims(ideogram,radius_inner)-0.5r\n";
    print OUT "grid_end           = dims(ideogram,radius_outer)+100\n";
    print OUT "\n";
    print OUT "<ticks>\n";
    print OUT "radius               = dims(ideogram,radius_outer)\n";
    print OUT "tick_separation      = 1p\n";
    print OUT "min_label_distance_to_edge = 0p\n";
    print OUT "label_separation = 1p\n";
    print OUT "label_offset     = 2p\n";
    print OUT "label_size = 48p\n";
    print OUT "multiplier = 1e-6\n";
    print OUT "color = black\n";
    print OUT "\n";
    print OUT "<tick>\n";
    print OUT "spacing        = 5u\n";
    print OUT "size           = 5p\n";
    print OUT "thickness      = 2p\n";
    print OUT "color          = black\n";
    print OUT "show_label     = no\n";
    print OUT "label_size     = 10p\n";
    print OUT "label_offset   = 0p\n";
    print OUT "format         = %d\n";
    print OUT "grid           = yes\n";
    print OUT "grid_color     = grey\n";
    print OUT "grid_thickness = 1p\n";
    print OUT "</tick>\n";
    print OUT "\n";
    print OUT "<tick>\n";
    print OUT "spacing        = 20u\n";
    print OUT "size           = 8p\n";
    print OUT "thickness      = 2p\n";
    print OUT "color          = black\n";
    print OUT "#show_label     = yes\n";
    print OUT "show_label     = no\n";
    print OUT "label_offset   = 0p\n";
    print OUT "format         = %d\n";
    print OUT "grid           = yes\n";
    print OUT "grid_color     = dgrey\n";
    print OUT "grid_thickness = 1p\n";
    print OUT "</tick>\n";
    print OUT "</ticks>\n";
    print OUT "\n";
	close(OUT);
    print STDERR "ticks.conf created\n";
  }

  my %len;
  my %start;
  my %end;
  my %chrom;
  my $first=0;
  my %data;
  open(IN,"<$INFILE1") or die("$INFILE1: $!\n");
  do {
    local $/="\n>";
    while(<IN>) {
      if (!$first++) {if (!m/^>\S+[^\n]*\n[a-z]*\n/is) {die("$INFILE1 is not a fasta file!\n");}}
      s/^>?/>/s; s/>$//;
      if (m/^>(\S+)[^\n]*\n(.*)$/is) { # If well information is not present.
        my ($name,$seq)=($1,$2);
        $seq=~s/\s+//gs;
        if (length($seq)<1) {next;} # Empty sequence?
        $data{$name}=length($seq);
        next;
      }
      die("Fasta file $INFILE1 not readable:\n$_\n"); # What?
    }
  };
  close(IN);

  # Sort the data by well and size, or just name if no well
  my %kdata;
  foreach my $n (keys %data) {
    my $k=$n;
    if ($n=~m/^([A-H]\d{2})-/) {$k=sprintf("%3s-%08d",$1,$data{$n});}
    $kdata{$n}=$k;
  }
  my @order=sort {$kdata{$a} cmp $kdata{$b}} keys %data;

  # Assign start and end positions in well
  foreach my $l (@order) {
    my $name=$l;
    (my $well=$l)=~s/^([A-H]\d\d)-.*/$1/is;

    $start{$name}=($len{$well}||0)+1; # Start at 1 past previous end
    $end{$name}=$start{$name}+$data{$l}; # End at length(seq)+1
    $len{$well}=$end{$name};

    $chrom{$name}=$well;
    next;
  }

  # Write out the ring:
  open(OUT,">temp.karyotype.$$") or die("temp.karyotype.$$: $!\n");
  push @TO_DELETE,"temp.karyotype.$$";
  print OUT map {"chr - $_ $_ 1 $len{$_} black\n"} sort keys %len;
  print OUT map {"band $chrom{$_} $_ $_ $start{$_} $end{$_} red\n"}
    sort {sprintf("%s.%08d",$chrom{$a},$end{$a}-$start{$a}) cmp sprintf("%s.%08d",$chrom{$b},$end{$b}-$start{$b})} keys %chrom;
  close(OUT);

  
  # Write links.
  open(IN,"<$result") or die("$result: $!\n");
  open(OUT,">temp.links.$$") or die("temp.links.$$: $!\n");
  push @TO_DELETE,"temp.links.$$";
  my $count=0;
  while(<IN>) {
    if (m/^#/) {next;}
    chomp;
    my @f=split(/\t/,$_);
    
    if (!defined $SELF and $f[0] eq $f[1] and $f[6] eq $f[8] and $f[7] eq $f[9]) {next;} # Self-hit from a self-blast

    my $well1=$f[0];
    my $name1=$f[0];
    my $well2=$f[1];
    my $name2=$f[1];
    $well1=~s/^([A-H]\d\d)-.*/$1/is;
    $well2=~s/^([A-H]\d\d)-.*/$1/is;

    if (defined $BSR) {
      my $bkey="$f[0]|$f[6]|$f[7]";
      my $bsr=${$bsrref}{$bkey};
      if (!defined $bsr) {next;} # Not a full-length match, no way to compensate
      if ($bsr/$f[11]<$BSR) {next;} # Failed the bsr cutoff
    }

    # Write the new link:
    $count++;

    print OUT sprintf("link%06d ",$count);
    print OUT "$well1 ",$start{$name1}+$f[6]," ",$start{$name1}+$f[7],"\n";
    print OUT sprintf("link%06d ",$count);
    print OUT "$well2 ",$start{$name2}+$f[8]," ",$start{$name2}+$f[9],"\n";
    # 0=query id, 1=subject id, 2=% identity, 3=alignment length,
    # 4=mismatches, 5=gap opens, 6=q. start, 7=q. end, 8=s. start, 9=s. end,
    # 10=evalue, 11=bit score
  }
  close(IN);
  close(OUT);

  # Run the circos command:
  my $cmd="$CIRCOS -conf temp.circos.$$";
  my $r=`$cmd 2>&1`;
  if (!defined $r){$r="$CIRCOS not found\n";}
  if ($r !~m/^created image at/m) {print "$CIRCOS output:\n$r";}
  if (! -e $OUTFILE) {die("$CIRCOS did not create $OUTFILE\n");}
  exit 3;
}

sub bl2seqfiles {
  my ($IN,$DB)=@_;
  print STDERR "Generating bl2seq data for $IN\n";
  my ($BLAST)=blast($IN,$DB);
    # Fields: 
    #  0: query id,
    #  1: subject id,
    #  2: % identity,
    #  3: alignment length,
    #  4: mismatches,
    #  5: gap opens,
    #  6: q. start,
    #  7: q. end,
    #  8: s. start,
    #  9: s. end,
    #  10: evalue,
    #  11: bit score

  open(IN,"<$BLAST") or die("$BLAST: $!");
  my %data;
  while(<IN>) {
    next if !m/^(\S+)\t(\1)\t/s; # Skip non-self hits (and comments)
    my @f=split(/\t/,$_);
    $f[11]=~s/\s//g;
    $data{"$f[0]|$f[6]|$f[7]"}=$f[11];
  }
  close(IN);
  return %data;
}

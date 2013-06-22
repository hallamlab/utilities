#!/usr/bin/perl -w
use strict;

# Bubble.pl by Charles Howes, 2009-04-14

# If using this script, please cite:
#-------------------------------------------------------------------------
# Zaikova, E., D. A. Walsh, C. P. Stilwell, W.W. Mohn, P. D. Tortell,
# and S. J. Hallam, 2009. Microbial community dynamics in a seasonally
# anoxic fjord: Sannich Inlet British Columbia, Environmental Microbiology,
# 12(1):172-91
#-------------------------------------------------------------------------
# Thank you; it will help us fund more script development.

eval "use SVG";
if ($@) { die("Please install the 'SVG' perl module before running this again.\n");}

# These come with perl:
use Getopt::Long;
use POSIX qw(ceil);

# Rough estimate of units per normal-width letter with this font, unfortunately machine-dependent:
my $FONTWIDTH=7;
my $FONTHEIGHT=9;

my $OUTPUT;
my $INPUT;
my $XSIZE=40;
my $YSIZE=20;
my $RADIUS;
my $LOGARITHM;
my $SCALE=0;
my $result=0;
my $ITALIC1;
my $ITALIC2;

# If there are command line arguments:
if (defined $ARGV[0]) {

  # Process command line arguments
  $result=GetOptions(
    "xsize:i"=>\$XSIZE,     # XSIZE per data column
    "ysize:i"=>\$YSIZE,     # YSIZE per data row
    "radius"=>\$RADIUS,     # Area versus radius = data point
    "i1"=>\$ITALIC1,        # Make header italic
    "i2"=>\$ITALIC2,        # Make row labels italic
    "log"=>\$LOGARITHM,     # Log scale
    "scale"=>\$SCALE,       # Add a scale to the bottom
    "output=s"=>\$OUTPUT);  # One output file

  # The input files without flags, so you can use wildcards.
  $INPUT=$ARGV[0];
}

# Usage, or bad arguments:
if (!$result or !defined $OUTPUT or !defined $INPUT) {
  die("Usage: $0 [-r] [-s] [-l] [-x xsize] [-y ysize] -o output.svg input.csv
  
    This program will read the input csv file and produce an SVG image
    showing the data as a bubble chart.  The first value in the first
    row becomes the graph title, and the remainder of the first row
    becomes column labels.  The first value in each of the remaining rows
    becomes the row label, with the remaining values being interpreted
    as data values.

    -x xsize: The width of the columns ($XSIZE)
	-y ysize: The height of the rows ($YSIZE)
	-s : Add a scale to the bottom
	-r : Make the data value proportional to the radius, not the area. 
	-l : Scale the data by logarithm
	-i1 : Make column labels italic 
	-i2 : Make row labels italic

    Sample input file:

  name,libA,libB,libC
  x,3,23,0
  y,0,2,5
  z,11,12,0
  xx,3,3,27
"
  );
}

my @data;
my @header;

# Read in the data:
if(!-f $INPUT || !-r _) { die("$INPUT isn't a readable file!"); }
open(IN,"<$INPUT") or die("$INPUT: $!");
my $eolchecked=0;
while (<IN>) {
  if (!$eolchecked) {  # Check for unusual end-of-line characters:
	  if (m/\r\n./) {$/="\r\n";seek(IN,0,0);$eolchecked=1;next;} # Windows
	  if (m/\r[^\n]/) {$/="\r";seek(IN,0,0);$eolchecked=1;next;} # Mac
	  if (m/\n[^\n\r]/) {$/="\r";seek(IN,0,0);$eolchecked=1;next;} # Unix (shouldn't happen)
  }
  s/[\r\n]+$//;
  my @f=map {s/"//g;$_} split(/[\t,]/);
  if (@header==0) {@header=@f;next;}
  $f[0]=~s/_/ /g;
  push @data,[@f];
}
close(IN);

my $maxheaderlen=0; # Maximum header length

# Find maximum header length:
for (my $x=1;$x<@header;$x++) {
  if (length($header[$x])>$maxheaderlen) {$maxheaderlen=length($header[$x]);}
}

my $MAXX=0;      # Maximum number of columns (including names)
my $MAXY=@data;  # Maximum rows, easily calculated
my $MAXVAL=0;    # Maximum data value, for scaling purposes
my $MAXWIDTH=0;  # Maximum name width, for first column

# Find maximum values:
foreach my $m (@data) {
  if ($#{$m}>$MAXX) {$MAXX=$#{$m};}  # number of columns (including name column)
  if (length($m->[0])>$MAXWIDTH) {$MAXWIDTH=length($m->[0]);} # Max row label width
  for (my $n=1;$n<@$m;$n++) {
    if ($m->[$n] !~ m/^[-\d.e]+$/i) {$m->[$n]=0;} # Non-numeric data = 0
    if ($m->[$n]>$MAXVAL) {$MAXVAL=$m->[$n];}  # Highest data value
  }
}

# This error checking also avoids taking log(0) next:
if ($MAXVAL==0) {die("The maximum value found in this file was zero.  Is there an end-of-line character mismatch?\n");}

# Scale:
if ($SCALE) {
  $MAXVAL=2**ceil(log($MAXVAL)/log(2)); # Round up the maximum value to a power of 2.
  push @data,["Scale:",sprintf("%.2g",$MAXVAL/8),sprintf("%.2g",$MAXVAL/4),sprintf("%.2g",$MAXVAL/2),sprintf("%.2g",$MAXVAL)];
  if ($MAXWIDTH<6) {$MAXWIDTH=6;}
}

# X and Y offsets from top left of image:
my $XOFFSET=$FONTWIDTH*$MAXWIDTH;
my $YOFFSET=$YSIZE*2+$FONTWIDTH*$maxheaderlen*sin(30/180*3.1415926535);

# What is the maximum value after sqrt and/or log transformation:
my $mxv=$MAXVAL;
if (!defined $RADIUS) {$mxv=sqrt($mxv);}
if (defined $LOGARITHM) {$mxv=log($mxv);}
if ($mxv==0) {$mxv=1;}  # log(1)=0, which would lead to division by zero later

# Total image size:
my $WIDTH=tx($MAXX+3);
my $HEIGHT=ty($MAXY+1+$SCALE);

# Start creating an SVG file:
my $svg=SVG->new(width=>$WIDTH,height=>$HEIGHT);

# This group is for dashed lines:
my $lines=$svg->group(id => 'lines',
	style => { stroke=>'black',
		'stroke-dasharray'=>"2 2",
		'stroke-width'=>0.1,
		fill=>'none'
    }
);

# This group is for data points:
my $positive=$svg->group(id => 'positive_data' , style => { stroke=>'black' });

# This group is for left-justified text (column labels):
my %cstyle=(
  'font-size'=>10,
  'font-weight'=>'lighter',
  'font-family'=>'sans',
  'text-anchor'=>'start'
  );
if (defined $ITALIC1) { $cstyle{'font-style'}="italic"; }
my $ctext=$svg->group(id=>'left_text', style=>\%cstyle);

# This group is for right-justified text (row labels):
my %style=(
  'font-size'=>10,
  'font-weight'=>'lighter',
  'font-family'=>'sans',
  'text-anchor'=>'end'
);
if (defined $ITALIC2) { $style{'font-style'}="italic"; }
my $text=$svg->group(id=>'text', style=>\%style);

# This group does the centered title and scale labels:
my $centertext=$svg->group(id=>'center_text',
  style=>{ 'font-size'=>10,
  'font-weight'=>'lighter',
  'font-family'=>'sans',
  'text-anchor'=>'middle'
  });

# ID numbers for SVG objects, to make them unique:
my $id=0;

# Title:
$centertext->text(id=>"centerid_".($id++),x=>tx($MAXX/2),y=>$FONTHEIGHT,-cdata=>$header[0]);

# Vertical lines and headers:
for (my $x=1;$x<=$MAXX;$x++) {
  $lines->line(x1=>tx($x), x2=>tx($x), y1=>ty(-0.5), y2=>ty($MAXY-1), id=>"id_".($id++));
  my $m="translate(".tx($x)." ".ty(-0.8).") rotate(-30)";
  $ctext->text(id=>"vertid_".($id++),transform=>$m,-cdata=>$header[$x]);
}

# Data rows and horizontal lines:
for (my $y=0;$y<=$MAXY;$y++) {
  my $DS=0;
  if ($y==$MAXY and $SCALE) {$DS=0.5;}  # Draw the scale, we're at last line
  if ($y==$MAXY and !$SCALE) {last;}   # We're at the the last line

  my $thisx=$MAXX;  #  How many columns to draw?
  if ($DS) { # if draw scale is set:
    $thisx=4;  # The scale has only 4 points
    $ctext->text(id=>"nameid_".($id++),x=>tx(-0.5),y=>ty($y+$DS),-cdata=>$data[$y][0]);
    $lines->line(x1=>tx(0), x2=>tx($thisx), y1=>ty($y+$DS), y2=>ty($y+$DS), id=>"id_".($id++));
  } else {
    $text->text(id=>"nameid_".($id++),x=>tx(0.1),y=>ty($y),-cdata=>$data[$y][0]);
    $lines->line(x1=>tx(0), x2=>tx($MAXX), y1=>ty($y), y2=>ty($y), id=>"id_".($id++));
  }

  # Draw the data points:

  for (my $x=1;$x<=$thisx;$x++) {
    my $v=$data[$y][$x];
    if (!defined $v or $v==0) {next;}
    if (!defined $RADIUS) {$v=sqrt($v);}
    if (defined $LOGARITHM) {$v=log($v);}

	$v=$v/$mxv;  # $v = value between 0 and 1, even if $mxv<0

    $v=$v*min($XSIZE,$YSIZE)/2; # $v = value between 0 and halfway to the next data point

    $positive->circle(cx=>tx($x-$DS), cy=>ty($y+$DS), r=>$v, id=>"circle_${x},${y}_($v)");

    if ($DS) { # Write labels for the scale:
      $centertext->text(x=>tx($x-$DS),y=>ty($y+1.5), id=>"scale_${x},${y}", -cdata=>$data[$y][$x]);
    }
  }
}

# Write the file:
open(OUT,">$OUTPUT") or die("$OUTPUT: $!");
print OUT $svg->xmlify;
close(OUT);


# Helper functions for translating x,y values to svg coordinates:
sub tx { return $XOFFSET+$_[0]*$XSIZE; }
sub ty { return $YOFFSET+$_[0]*$YSIZE; }

sub min { return ($_[0]<$_[1])?$_[0]:$_[1]; }
sub max { return ($_[0]>$_[1])?$_[0]:$_[1]; }

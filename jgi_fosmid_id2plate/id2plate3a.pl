#!/usr/bin/perl 
#
#Program name: id2plate3.pl
#
# 010226 S. Pitluck -- Created program cosmid_origin.pl
# 010320 S. Pitluck -- Modified for id2plate.pl
# 011126 S. Pitluck -- Modified for reading a list containing ABC1234
# 051020 S. Pitluck -- Modified to produce 1 line output for each input line
#						Assuming type 1 and PGF plates
# Given a JGI read nr, calculate plate and well location
#
$numargs = @ARGV;
if ($numargs < 1) {
 print "Usage: id2plate3.pl  file_name [type]\n";
 print "       default type is 1 other possible value is currently 2\n";
 print " type = 1 for production A C\n";
 print "                         B D\n";
 print "\n";
 print " type = 2 for rnd        A B\n";
 print "                         C D\n";
 exit;
}
$input_file  = shift;
open (IN,"$input_file") or die "Cannot open input file: $input_file\n";
@list = <IN>;
close (IN);
$type = 1;
if ($numargs == 2) { $type = shift; }
if ($type < 1 || $type > 2) {
	print "Resetting type to 1\n";
	$type = 1;
}
foreach $read (@list) {
 chomp $read;
 $read =~ /(^\D+)(\d+)/;
 $lib = $1;
 $id = $2;
$plate       = int(($id - 1)/96) + 1;	
##print "\n";
$quad = ( ($plate-1) % 4) + 1;	# quadrant of 384 well plate
# 1 3	
# 2 4
# Note that this is the order for production
# sequencing at the PGF
#
##print "For PGF id: $lib$id\n";
#print "Quad: $quad of PGF 96-well plate: $tplate\n";
##print "Quad: $quad\n";
# plate_offset is location within the 96 well plate
$plate_offset = $id - ($plate-1)*96;
# column is column number within  the lanl plate
$column = int(($plate_offset-1)/8) + 1;
$row    = $plate_offset % 8;
if ($row == 0) { $row = 8; }
#
if ($quad == 1 && ($type == 1 || $type == 2)) {
	$psf_col = ($column-1)*2 + 1;
	$psf_row = ($row-1)*2 + 1;
} elsif ( ($quad == 3 && $type ==1) || ($quad == 2 && $type ==2)) {
	$psf_col = ($column-1)*2 + 2;
	$psf_row = ($row-1)*2 + 1;
} elsif ( ($quad == 2 && $type ==1) || ($quad == 3 && $type ==2)) {
	$psf_col = ($column-1)*2 + 1;
	$psf_row = ($row-1)*2 + 2;
} elsif ($quad == 4 && ($type ==1 || $type == 2)) {
	$psf_col = ($column-1)*2 + 2;
	$psf_row = ($row-1)*2 + 2;
}
#print "This is row $psf_row and col $psf_col\n";
$psf_row = chr( $psf_row + 64);
$Plate =  int(($plate-1)/4) + 1;
$row96 = chr($row + 64);
$Plate4 = ($Plate-1)*4 + 1;
$plate_end = $Plate4+3;
$plate384 = int(($id-1)/384) + 1;
$jgiplate = ($plate384-1)*4+1;
$plate384 = $plate384 -1;
print "$read Well $psf_row$psf_col for sequential 384-well Plate $plate384 (PGF Plate $jgiplate)\n";
##print "This would be 96-well plate $plate: $row96$column\n";
}

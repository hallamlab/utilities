#!/usr/bin/perl


use strict;
use warnings;
use Error;
use Getopt::Long;
sub trim;
sub standard_dev;
sub average;
sub sum;
sub max;
sub read_fasta_sequences;
sub get_file_data ;
sub open_file ;
sub binary_string ;
sub array_unique;
sub min;
sub Similarity;

use Data::Dumper;
use CGI ':standard';
use LWP 5.64; # Loads all important LWP classes, and makes
                #  sure your version is reasonably recent.
use List::Util qw(shuffle);
no strict; # TODO fix this
use lib "/Users/kishori/PROJECTS/LIBRARIES/PERL/";
use Statistics::Distributions;


use constant VERSION => 1.0;
use constant FROM => "FROM";
use constant TO => "TO";
use constant FROMID => "FROMID";
use constant TOID => "TOID";
use constant CORREL => "CORREL";

use constant ID => "ID";
use constant NAME => "NAME";
use constant LABEL => "LABEL";
use constant TYPE => "TYPE";
use constant SIZE => "SIZE";
use constant TAX => "TAX";
use constant MAX_DOUBLE => 1000000000;
use constant MIN_DOUBLE => -1000000000;
use constant DEBUG=> 0;


sub get_ec_numbers;
sub binary_string;
sub plot_graph;
sub newhash_from_array;
sub addhash_from_array;
sub fillhash_from_array;
sub tax_level_id;
sub tax_level;
sub has_sub_strings;
sub create_normalizer_vector;
sub create_metadata_vector;
sub create_label_limit_hash;
sub create_label_maps;
sub p_value_t_statistics;
sub t_statistic;
sub trim_otu_name;

sub scaling_product;
sub dot_product;
sub create_label;
sub get_label_from_limits;

#use constant PROCESSED_LIST => "processed_list.txt";

{

my $TAX_LEVEL;
my $INPUT_FILE;
my $OUTPUT_FILE;
my $PARAM_FILE;
my $FORCE_SUBSTRINGS=0;
my $SAMPLE_PARAMETERS_N_LABEL_FILE="";
my $MIN_CORREL=-0.4;
my $MAX_CORREL=0.4;
my $MAX_PVALUE=0.005;
my @METADATA_LABELS;
my $METADATA_LABEL_LIMIT_FILE=undef;
my $IDENTIFY_OTU_BY_SAMPLE_ATTRIBUTE = undef;
my $MIN_PRESENT=25;
my $MIN_THRESHOLD=0.0;
my $LABEL_MAP_FILE =undef;
my $STATISTIC="Pearson";
my $SHUFFLE;



my $result;
# If there are command line arguments:
if (defined $ARGV[0]) {
  # Process command line arguments
   $result=GetOptions(
    "i=s"=>\$INPUT_FILE,
    "o=s"=>\$OUTPUT_FILE,
    "label_file=s"=>\$LABEL_FILE,
    "mlabel=s"=>\@METADATA_LABELS,
    "max_negative_correl=f"=>\$MIN_CORREL,
    "min_positive_correl=f"=>\$MAX_CORREL,
    "max_p_value=f"=>\$MAX_PVALUE,
    "min_count=f"=>\$MIN_PRESENT,
    "min_threshold=f"=>\$MIN_THRESHOLD,
    "statistic=s"=>\$STATISTIC,
    "shuffle=i"=>\$SHUFFLE
  );
}

# Usage, or bad arguments:
if (!$result or  !defined($INPUT_FILE) or !defined($OUTPUT_FILE)  ) {
  die("Usage: $0 -i inputfile   -o outputfile   
                  [-min_count c (default 25%) ]  [-mlabel label ]  [ -max_negative_correl 0.x (default -0.4) ] 
                  [ -min_positive_correl 0.y (default 0.4) ] [-max_p_value max_pvalue (default 0.05] 
                  [ -min_threshold x] \n
                  This program will prepare the properties of the node 
                  --inputfile the otu table file
                  --outputfile is the name of the output file which produces two files: an edge file and a node file 
                  --labelfile  a label file that has the name of the sub samples, one in each row
                  --min_count  min number (in %)  of non-zero columns for an otu to be used, default is 25% of the total columns
                  --min_threshold  is the minimum normalized value of the read counts to count as non-zero

                 This script will create the edge and node files for the network
                  --statistic  [Pearson | Kendall | Spearman ] (Default Pearson)
                  --shuffle [ shuffles the columns of values row-by-row (default 10)  ] 
                  \n"
     );  
}


my @nodelines = get_file_data($INPUT_FILE);

open my $OUTFILE_NODE, ">", $OUTPUT_FILE."_nodes.txt" or die{"Failed to open output file"};
open my $OUTFILE_EDGE, ">", $OUTPUT_FILE."_edges.txt" or die{"Failed to open output file"};


## Creating the normalized data vector
my %metadata_vector;
my %data_table;
%data_table = ();
my $otu_count=0;

my @normalizer_vector = create_normalizer_vector(\@nodelines);
my @sample_names = get_sample_names($INPUT_FILE);
#print Dumper(\@sample_names);
#exit;

my %sample_parameters_n_labels;
if( defined($LABEL_FILE) ) {
    %sample_parameters_n_labels = read_column_file_data($LABEL_FILE, "\t");
}
#print Dumper(\@sample_names);
#print Dumper(\%sample_parameters_n_labels);
#print Dumper(\@normalizer_vector);

print "Sample Names \t Normalizing Constant\n";
for(my $i = 0; $i < scalar(@sample_names); $i++) {
    print $sample_names[$i]."\t".$normalizer_vector[$i]."\n";
}


#print Dumper(\@METADATA_LABELS);
my $min_num_samples = undef;

foreach my $line (@nodelines) {
   $line = trim($line);
   if( $line =~ /^#/m ) { 
      my @fields = split("\t",$line);
      my $i=0;
      map{ $fields[$i] = trim($_); $i++; } @fields;
      if( scalar(@fields) < 3 ){ next; }
      create_metadata_vector($line, \@sample_names, \%sample_parameters_n_labels, \@METADATA_LABELS, \%metadata_vector);
      next; 
   }
 
   my @fields = split("\t",$line);
   my  $i=0;
   map{ $fields[$i] = trim($_); $i++; } @fields;
   my $otuname =$fields[ scalar(@fields)-1];
   $otuname = $fields[0]."__".trim_otu_name($otuname); 
           #$otuname = $fields[0]$otu_count."__".trim_otu_name($otuname); 
           #print $otuname."  ".$fields[0]."\n";

    shift @fields;
    splice @fields, scalar(@fields)-1, 1;
    my $num_non_zero=0;
    map{ if( $_ > $MIN_THRESHOLD) { $num_non_zero++} } @fields;
   
    if( !defined($min_num_samples ) ) {
      $min_num_samples = int($MIN_PRESENT*scalar(@fields)/100);  
      print "Number of required samples to be non-zero = ".$min_num_samples."\n";
    }
    if( $num_non_zero < $min_num_samples ) { next;}
    $otu_count++; 
    #if( DEBUG==0) { print $num_non_zero."\n"; }
        
    if( sum(\@fields) <=0 ) { next; }
    @{$data_table{$otuname}} = @fields;
    eval {
        scaling_product($data_table{$otuname}, \@normalizer_vector );
    };
    if($@) {
        print $@."Error in dot product\n";
        exit();
    }
}

print "Number of candidate otus INITIALLY filtering ".scalar(@nodelines)."\n";
print "Number of candidate otus AFTER filtering ".scalar(keys %data_table)."\n";
print "Using statistic : ".$STATISTIC."\n";
if( defined($SHUFFLE) ) {
   print "Number of random shuffles : ".$SHUFFLE."\n";
}
# print Dumper(\%metadata_vector); exit;
# print Dumper(\@normalizer_vector);

#print Dumper(\%data_table); exit;


my @otus = keys %data_table;

#print Dumper(\@otus);exit;
my $count=0;
my %otus_id;
%otus_id = map { $_=>$count++; } @otus;
$count=0;
my %ids_otus = map { $count++=>$_; } @otus;



# re-enumerate the OTUs
my %vertices=();


#print Dumper(\@sample_names);


for(my $i=0; $i < scalar( keys %otus_id); $i++)  {
   %{$vertices{$ids_otus{$i}}} = (ID=>$i, TAX=>$ids_otus{$i});
   $vertices{$ids_otus{$i}}{'READ_COUNT'}=sum($data_table{$ids_otus{$i}});
}


#now compute and add the metadata fields

my @headers = qw( ID LABEL TAX READ_COUNT );
foreach my $metadata ( @METADATA_LABELS ) {
   push(@headers, $metadata);
   for(my $i=0; $i < scalar( keys %otus_id); $i++)  {
      my $metadata_value = dot_product($data_table{$ids_otus{$i}}, $metadata_vector{$metadata})/
                           $vertices{$ids_otus{$i}}{'READ_COUNT'};
      $vertices{$ids_otus{$i}}{$metadata}=$metadata_value;
   }
}


# now compute all possible edges
if( defined($SHUFFLE) ) {
  print $OUTFILE_EDGE  "FROM\tTO\tCORREL\tRANDOM_CORREL\tCORREL_TYPE\tPVALUE\n";
}
else{
  print $OUTFILE_EDGE  "FROM\tTO\tCORREL\tCORREL_TYPE\tPVALUE\n";
}
my $tot_edges = 0;
my $otu_num= scalar(keys %otus_id);

my $stat_length = 80;
for(my $i=0; $i< $stat_length; $i++) {
   print "-";
}

print "\n";
my $full_length = $otu_num*($otu_num-1)/2;
my $current_length=0;
my $actual_length=0;
my $num_positive_correl=0;
my $num_negative_correl=0; 

my %RV = ();
my $seconds=time();
my $N = scalar(@sample_names);

my $statiticMatch = 0;
my @statisticNames = qw();

my %edges;
my $numStatistics = 0;
if ($STATISTIC =~ /Pearson/ ) {
    $edges{Pearson} = [] ;
    $numStatistics++;
    push(@statisticNames, "Pearson");
}
if ($STATISTIC =~ /Kendall/ ) {
    $edges{Kendall} = [];
    $numStatistics++;
    push(@statisticNames, "Kendall");
}
if ($STATISTIC =~ /Spearman/ ) {
    $edges{Spearman} = [];
    $numStatistics++;
}



print $otu_num . " otu num \n";
# compute the pair-wise correlation measures
for(my $i=0; $i < $otu_num; $i++) {
   for(my $j=0; $j < $i; $j++) {
       my $pvalue;

       if ($STATISTIC =~ /Pearson/ ) {
           my %edge=();
           $edge{FROM} = $ids_otus{$i};
           $edge{TO} = $ids_otus{$j};
           $edge{FROMID} = $i;
           $edge{TOID} = $j;

           $edge{CORREL} = correl($data_table{$ids_otus{$i}}, $data_table{$ids_otus{$j}});

           #now compute the random correlations
           if( defined($SHUFFLE) ) {
               $edge{RANDOM_CORREL} = 0;
               for( my $j =0; $j < $SHUFFLE; $j++) {
                  my @XX = shuffle(@{$data_table{$ids_otus{$i}}});
                  my @YY = shuffle(@{$data_table{$ids_otus{$j}}});
                  $edge{RANDOM_CORREL} += correl($data_table{$ids_otus{$i}}, $data_table{$ids_otus{$j}});
               }
               $edge{RANDOM_CORREL} =$edge{RANDOM_CORREL}/$SHUFFLE; 
           }

           $pvalue= p_value_t_statistics(\%edge, \%data_table, \%ids_otus); 
           #print $edge{CORREL}."  ".$pvalue."\n"; 
           $statiticMatch = 1;
           ($edge{CORREL} > 0) ? ($edge{CORREL_TYPE}="Positive") : ($edge{CORREL_TYPE} = "Negative"); 
           $edge{PVALUE} = $pvalue;
           push(@{$edges{Pearson}}, \%edge);

#           if( ( $pvalue < $MAX_PVALUE) and ( $edge{CORREL} < $MIN_CORREL || $edge{CORREL} > $MAX_CORREL) ) {
#               $edge{ID} = $count++;
#               ($edge{CORREL} > 0) ? ($num_positive_correl++): ($num_negative_correl++); 
#               push(@{$edges{Pearson}}, \%edge);
#           }

       }

       if ($STATISTIC =~ /Kendall/ ) {
           my %edge=();
           $edge{FROM} = $ids_otus{$i};
           $edge{TO} = $ids_otus{$j};
           $edge{FROMID} = $i;
           $edge{TOID} = $j;
           $edge{CORREL} = kendall_correl($data_table{$ids_otus{$i}}, $data_table{$ids_otus{$j}});

           #now compute the random correlations
           if( defined($SHUFFLE) ) {
               $edge{RANDOM_CORREL} = 0;
               for( my $j =0; $j < $SHUFFLE; $j++) {
                  my @XX = shuffle(@{$data_table{$ids_otus{$i}}});
                  my @YY = shuffle(@{$data_table{$ids_otus{$j}}});
                  $edge{RANDOM_CORREL} += kendall_correl($data_table{$ids_otus{$i}}, $data_table{$ids_otus{$j}});
               }
               $edge{RANDOM_CORREL} =$edge{RANDOM_CORREL}/$SHUFFLE; 
           }
           $pvalue= p_value_kendall_statistic($edge{CORREL}, $N); 
           #print $edge{CORREL}."  ".$pvalue."\n"; 
           $statiticMatch = 1;
           ($edge{CORREL} > 0) ? ($edge{CORREL_TYPE}="Positive") : ($edge{CORREL_TYPE} = "Negative"); 
           $edge{PVALUE} = $pvalue;
           push(@{$edges{Kendall}}, \%edge);

#           if( ( $pvalue < $MAX_PVALUE) and ( $edge{CORREL} < $MIN_CORREL || $edge{CORREL} > $MAX_CORREL) ) {
#               $edge{ID} = $count++;
#               ($edge{CORREL} > 0) ? ($num_positive_correl++): ($num_negative_correl++); 
#               push(@{$edges{Kendall}}, \%edge);
#           }

       }

       if ($STATISTIC =~ /Spearman/ ) {
           $edge{CORREL} = spearman_correl(\%data_table, $ids_otus{$i}, $ids_otus{$j}, \%RV);
           $pvalue= p_value_spearman_statistic($edge{CORREL}, $N); 
           #now compute the random correlations
           $statiticMatch = 1;
       }
 
       if( $statiticMatch == 0) {
           print "ERROR: ".$STATISTIC." is not a supported statistic!\n\n";
           exit(0);
       }

       $tot_edges++;
       $actual_length++;
     }
     if( $seconds < time() ) {
        printf STDERR "%4.2f%%\r", 100*$actual_length/$full_length;
        $seconds = time();
     }
  }


 my @cutoffLevels = qw( 0.2 0.4 0.6 0.8);
 if( $numStatistics > 0 ) {
     foreach my $cutoff1 (@cutoffLevels) {
       foreach my $cutoff2 (@cutoffLevels) {
         my $similarity = Similarity(\%edges, $statisticNames[0], $statisticNames[1], $cutoff1, $cutoff2);
          print $statisticNames[0] . "   ". $cutoff1 . "  ". $statisticNames[1]. "   ". $cutoff2 . " = ". $similarity."\n";
       }
     }
     foreach my $statisticName (@statisticNames) {
        print "Number of edges in ". $statisticName . " : ". scalar(@{$edges{$statisticName}})."\n";
     }

 #    print Dumper(\%edges);
 }


  if( $seconds < time() ) {
     printf STDERR "%4.2f%%\r", 100*$actual_length/$full_length;
      $seconds = time();
  }




foreach my $edge  (@edges) {
  if( defined($SHUFFLE) ) {
    print $OUTFILE_EDGE "$edge->{FROMID}\t$edge->{TOID}\t$edge->{CORREL}\t$edge->{RANDOM_CORREL}\t$edge->{CORREL_TYPE}\t$edge->{PVALUE}\n";
  }
  else {
    print $OUTFILE_EDGE "$edge->{FROMID}\t$edge->{TOID}\t$edge->{CORREL}\t$edge->{CORREL_TYPE}\t$edge->{PVALUE}\n";
  }
}

#my @headers=   keys %{$vertices{$ids_otus{0}}} ; 


#write the nodes out
my $flag=0;
foreach my $header ( @headers ) {
    if( $flag==0) {
       print $OUTFILE_NODE  $header;
       $flag=1;
    }
    else {
       print $OUTFILE_NODE "\t".$header;
    }
}
print $OUTFILE_NODE  "\n";


my $barekey;
my $newkey;
while( my ($key, $value) = each %vertices ) {
    $flag=0;
    $barekey = $key;
    $barekey =~s/[0-9]*_//ig;

   # print $barekey."\n";
    $newkey = $barekey;
    if(defined($label_maps->{$barekey}))  {
       $newkey = $label_maps->{$barekey};
    }
    $vertices{$key}{LABEL} = $newkey;

    foreach my $header ( @headers ) {
        #print  $vertices{$key}{$header}."\t";
        if( $flag==0) {
           print $OUTFILE_NODE  $vertices{$key}{$header};
           $flag=1;
        }
        else {
           print $OUTFILE_NODE  "\t".$vertices{$key}{$header};
        }
    }
    print $OUTFILE_NODE "\n";
}


print "Number of nodes =".scalar( keys %vertices)."\n";
print "Number of edges = ".scalar(@edges)."/".$tot_edges."\n";
print "Number of positive edges = ".$num_positive_correl."\n";
print "Number of negative edges = ".$num_negative_correl."\n";
     

close($OUTFILE_NODE);
close($OUTFILE_EDGE);
exit;

}  #end of main


#####################################################################################
##############  SUBROUTINES #################################################
sub Similarity {
    my ($edges, $statistic1, $statistic2, $cutoff1, $cutoff2) = @_ ;

    my $key; 
    my %sets =();
    foreach my $edge (@{$edges->{$statistic1}}) {
       $key = $edge->{TOID}.$edge->{FROMID};
       if( $edge->{CORREL} < -1*$cutoff1 || $edge->{CORREL} > $cutoff1)  {
         if( !defined($sets{$key}) ) {
            $sets{$key} = 0
          }
          $sets{$key}++;
       }
    }

    foreach my $edge (@{$edges->{$statistic2}}) {
       $key = $edge->{TOID}.$edge->{FROMID};
       if( $edge->{CORREL} < -1*$cutoff2 || $edge->{CORREL} > $cutoff2)  {
          if( !defined($sets{$key}) ) {
            $sets{$key} = 0
          }
          $sets{$key}++;
       }
    }

    my $intersectionCount = 0;
    my $unionCount = 0;
    foreach my $key (keys %sets) {
       if( $sets{$key} > 1 ) {
         $intersectionCount++;
       }
       $unionCount++;
    }
    
    return $intersectionCount/$unionCount;

}

sub p_value_t_statistics{
    my ($edge, $data_table, $ids_otus) = @_ ;
    #print Dumper($edge);
    #print Dumper($data_table);
    #print Dumper($ids_otus);
    $correl = correl($data_table->{$ids_otus->{$edge->{FROMID}}}, $data_table->{$ids_otus->{$edge->{TOID}}});
  
    my $n=scalar(@{$data_table->{$ids_otus->{$edge->{FROMID}}}}); 
     
    my $t_stat= t_statistic($correl, scalar(@{$data_table->{$ids_otus->{$edge->{FROMID}}}})); 
    #print $edge->{FROMID}. "  ".$edge->{TOID}.">>>".$correl."\n"; 

    #print Statistics::Distributions::tprob($n-2,$t_stat)." ".$t_stat." ". $correl."   ".$edge->{CORREL}."\n"; 
    return Statistics::Distributions::tprob($n-2,$t_stat);
}


sub t_statistic {
   my($correlation, $n ) = @_;
   if( DEBUG==1) {
     print "T-Statistic ". $correlation."\n";
   }

   if($correlation < 0 ) {$correlation = -1*$correlation; }
   if( $correlation > 0.99999) {return 10000;}
   my $stat =  $correlation*sqrt( ($n - 2)/(1 - $correlation*$correlation));
   return $stat ;

}

sub p_value_kendall_statistic{
   my($kendall_correlation, $n) = @_;
   my $z_stat = abs($kendall_correlation)*sqrt((9/2)*($n/(2*$n + 5))*($n-1));
   return Statistics::Distributions::uprob($z_stat);
}

sub p_value_spearman_statistic{
   my($spearman_correlation, $n) = @_;
   if( $spearman_correlation ==1)  {
        return 0.0;
   }
  # my $t_stat = abs($spearman_correlation)*sqrt( ($n-2)/(1 - $spearman_correlation*$spearman_correlation) );
   my $t_stat= t_statistic($spearman_correlation, $n); 
   return Statistics::Distributions::tprob($n-2, $t_stat);
}

sub identify_otu_by_label{
     my($read_vector,$sample_parameters_n_labels, $sample_names, $label, $min_perent) = @_;
     my $identity="";
#     print Dumper($read_vector);
#     print Dumper($sample_parameters_n_labels);
#     print Dumper($sample_names);

     my $total =sum($read_vector);
     my %buckets = map { $sample_parameters_n_labels->{$_}->{$label}=>0 } @{$sample_names};
     my $i=0;


     foreach my $sample (@{$sample_names}) {
       $buckets{$sample_parameters_n_labels->{$sample}->{$label}} += $read_vector->[$i];
       $i++;
     }

    # print Dumper(\%buckets);

     foreach my $bucket ( keys %buckets ) {
         if( $buckets{$bucket}/$total > 0.25 ) {
            $identity = $identity.$bucket.";";
         } 
     }
    # print $identity."\n"; 
     return( $identity);

}

sub get_label_from_limits{
    my $labels =shift;
    my $value =shift;
    foreach my $label ( keys %{$labels} ) {
        if( $labels->{$label}->{'min'} <= $value and  $labels->{$label}->{'max'} > $value ) {
            return( $label);
       }
    }
}
sub create_label {
   my $substrings = shift;
   my $label = shift;
   
   foreach my $string (@{$substrings})  {
      if( $label=~ /$string/ig ) {
          return($label);
      }
   }



   my @fields = split(";",$label);
   @fields = reverse(@fields);
   foreach my $field (@fields) {
      if( !($field =~ /unclassified/ig) ) {
          return($field);
      }
   }
   return($label); 
}

sub read_column_file_data {
   my ($filename, $fs) = @_;
   use strict;
   use warnings;
   unless(open(DATAFILE, $filename))
   {
      print  "Cannot open file $filename \n";
      exit;
   }
   my @filedata = <DATAFILE>;
   
   my %hashtable=(); 
   my $count =0;
   my @headers;
   foreach my $line (@filedata) {
      my @fields = split($fs, trim($line));
      if( scalar(@fields) > 0 && length( $fields[0]) ) {
          if( $count ==0 ) {
               @headers = @fields;
            #   print Dumper(\@headers); 
               $count++;
          }
          if( scalar(@headers) == scalar(@fields) ) {
             %{$hashtable{$fields[0]}}=() ;
             for(my $i=1; $i < scalar(@headers); $i++ ) { 
                $hashtable{$fields[0]}{$headers[$i]} = $fields[$i];
             }
          }
          else {
             print STDERR "Rows with incorrect/inconsistent nuber of fields\n";
         
          }
      }
   }
   close(DATAFILE);
   return %hashtable;
}

##########################################################################################
##########################################################################################
### create_metadata_vector  say oxygen content
sub  create_metadata_vector {
    my ($line, $sample_names, $sample_parameters, $meta_parameters, $metadata_vector) = @_;

    if( !($line =~ /^#OTU/) ) { return; }
    map{  @{$metadata_vector->{$_}} = qw(); } @{$meta_parameters}; 

    foreach my $parameter (@{$meta_parameters}) {
       foreach my $sample_name (@{$sample_names}) {
          push( @{$metadata_vector->{$parameter}}, ($sample_parameters->{$sample_name}->{$parameter}));
       }
    }
}


##########################################################################################
##########################################################################################
### create_metadata label limit has structure
sub create_label_limit_hash {
    my ($METADATA_LABEL_LIMIT_FILE) = @_;
   
    my $temp_hash = {};
    if( !defined($METADATA_LABEL_LIMIT_FILE) || $METADATA_LABEL_LIMIT_FILE eq "") {
       $temp_hash->{'NONE'} = {'min'=>MIN_DOUBLE, 'max'=>MAX_DOUBLE};
       return($temp_hash);
    }
    my %label_seen;
    my  @label_limit_lines = get_file_data($METADATA_LABEL_LIMIT_FILE);
    foreach my $line (@label_limit_lines) {
        if( $line =~ /^#/ ) { next; }
       

        my @fields = split("\t", $line);
        if( scalar(@fields) != 3 ) {
            print "Not 3 fields in some rows in file ".$METADATA_LABEL_LIMIT_FILE."\n";
            exit;
        }
        
        if( exists $label_seen{trim($fields[0])} ) {
            print "Duplicate labels in file ".$METADATA_LABEL_LIMIT_FILE."\n";
            exit;
        }
       $label_seen{trim($fields[0])} =1;
       $temp_hash->{trim($fields[0])} = {'min'=>$fields[1], 'max'=>$fields[2]};
    }
    return($temp_hash);
}
##########################################################################################
##########################################################################################
### create_metadata label map file
sub create_label_maps {
    my ($LABEL_MAP_FILE) = @_;
   
    my $temp_hash = {};
    if( !defined($LABEL_MAP_FILE) || $LABEL_MAP_FILE eq "" || !(-e $LABEL_MAP_FILE))
    {
       print STDERR "Cannot find label file ".$LABEL_MAP_FILE."\n";
       return($temp_hash);
    } 
    my  @label_map_lines = get_file_data($LABEL_MAP_FILE);
    foreach my $line (@label_map_lines) {
        if( $line =~ /^#/ ) { next; }
        $line = trim($line);
       

        my @fields = split("\t", $line);
        if( scalar(@fields) != 2 ) {
            print "Not 2 fields in some rows in label map file ".$LABEL_MAP_FILE."\n";
            exit;
        }
        $temp_hash->{trim($fields[1])} = trim($fields[0]);
    }
if(DEBUG==1) {
    print Dumper($temp_hash);
}
    return($temp_hash);
}

##########################################################################################
##########################################################################################
### get sample names
sub  get_sample_names{
    my ($OTU_TABLE_FILE) = @_;
    my @otu_table_lines = get_file_data($OTU_TABLE_FILE);
    my @sample_names;
    my %sample_name_field=();
    foreach my $line (@otu_table_lines) {
       $line=trim($line);
       @sample_names = split("\t", $line);
       if( $line =~ /^#/  and  scalar(@sample_names) > 2 ) { 
           shift @sample_names;
           splice(@sample_names, scalar(@sample_names)-1, 1);    
           my @temp_sample_names = @sample_names;
           # if it has space trim it
           my $i = 0;
           map{ $sample_names[$i] = trim($_); $i++;} @temp_sample_names;
           last; 
       }
    }
    return(@sample_names);
}

##########################################################################################
##########################################################################################
### create_normalizer vector
sub  create_normalizer_vector {
    my ($lines) = @_;
    my %sample_size = (); 
    my @sample_indices = qw(); 
    my @sample_names ; 
    my @normalizer_vector = qw(); 
    for my  $line ( @{$lines} ) {
         $line=trim($line);
         if( $line =~ /^#/ )  { 
             @sample_names = split("\t", $line);
             if( scalar(@sample_names) < 2) { next; } 
             
             my @tmp_sample_names = @sample_names;
             shift @tmp_sample_names;
             splice(@tmp_sample_names, scalar(@tmp_sample_names)-1, 1);    

             %sample_size =  map { $_=> 0} @tmp_sample_names; 
             my $i = 1;
             map { push(@sample_indices,$i);  $i++;} @tmp_sample_names; 
             %sample_size =  map {$_=> 0} @tmp_sample_names; 
             next;
         }
        # print Dumper(\@sample_names);  exit;
        # print Dumper(\@sample_indices);  exit;
    
         my @values = split("\t", $line);
         foreach my $index (@sample_indices) {
             #print $index."\t". $sample_names[$index]."\n";
             #print $values[$index]."\n";
             $sample_size{$sample_names[$index]} =  $sample_size{$sample_names[$index]} + $values[$index];
         } 
    }
    my @sizes = values %sample_size;
    my $max = max( \@sizes);
    
  #  print Dumper(\%sample_size); exit;
   # print Dumper(\%sample_names); exit;
    foreach my $index (@sample_indices) {
       #print $sample_size{$sample_names[$index]}."    ".$sample_names[$index]."   ".$index."\n";
       if(  $sample_size{$sample_names[$index]} > 0 ) {
          push( @normalizer_vector, ($max/$sample_size{$sample_names[$index]}));
       }
       else {
          print "Zero size on sample ".$sample_names[$index]."\n";
          push( @normalizer_vector, 0);
       }
    }
    return( @normalizer_vector);
}



##########################################################################################
### dot_product   -- dot product of vectors
sub dot_product {
   my ($a,$b)=@_;
   my $nbdata = scalar(@{$a});
   my $sum=0;
   try {
   for ($i=0;$i< $nbdata;$i++){
       $sum += ($a->[$i])*($b->[$i]);
    }
   }
   catch MathException with {
       # print "Multipliction error";
   }
   return($sum);
}

##########################################################################################
### scaling_product   -- scaling product of vectors
sub scaling_product {
   my ($a,$b)=@_;
   #print Dumper($a);
   my $nbdata = scalar(@{$a});
   #print "length a ". scalar(@$a). "    length b ".scalar(@$b). "  \n";
   #print Dumper($a);
   try {
   for ($i=0;$i< $nbdata;$i++){
      #print "a =".$i."   ".$a->[$i]."\n";
      #print "b =".$i."   ".$b->[$i]."\n";
      $a->[$i]= ($a->[$i])*($b->[$i]);
    }
    }
    catch MathException with {
      #print "Error in 457";
   #   exit;
      
    }
   return;
}

##########################################################################################
### SS = sum of squared deviations to the mean

sub SS {
   my ($a,$b)=@_;
   #print Dumper($a);
   my $nbdata = scalar(@{$a});
   my ($i,$sum)=(0,0);
   my $meana = average($a);
   my $meanb = average($b);
   for ($i=0;$i< $nbdata;$i++){
      $sum=$sum+($a->[$i]-$meana)*($b->[$i]-$meanb);
    }
   return $sum;
}

##########################################################################################
### Correlation

sub correl {
   my($X,$Y)=@_;

   my $ssxx = SS($X,$X);
   my $ssyy = SS($Y,$Y);
   my $ssxy = SS($X,$Y);
   my $sign=$ssxy/abs($ssxy);
   my $correl=$sign*sqrt($ssxy*$ssxy/($ssxx*$ssyy));
   return $correl;
}


# Kendall tau test 
sub kendall_correl {
   my($X,$Y)=@_;
   my $NX = scalar(@$X);
   my $concordantCount=0;
   my $discordantCount=0;
   my $X0 = 0;
   my $Y0 = 0;

   for( my $i =0; $i < $NX; $i++) {
     for( my $j =$i+1; $j < $NX; $j++) {
         if( ( ($X->[$i] > $X->[$j]) &&  ($Y->[$i] > $Y->[$j]) ) || ( ($X->[$i] < $X->[$j]) &&  ($Y->[$i] < $Y->[$j]) ) ) { 
               $concordantCount++;
         }

         if( ( ($X->[$i] < $X->[$j]) &&  ($Y->[$i] > $Y->[$j]) ) || ( ($X->[$i] > $X->[$j]) &&  ($Y->[$i] < $Y->[$j]) ) ) { 
              $discordantCount++; 
         }

         if( ( $X->[$i] == $X->[$j]) &&  ($Y->[$i] != $Y->[$j])   ) { 
              $X0++;
         }
         if(  ($X->[$i] != $X->[$j]) &&  ($Y->[$i] == $Y->[$j])  ) { 
              $Y0++;
         }
      }
   }

   my $kendall = ($concordantCount - $discordantCount)/(sqrt($discordantCount + $concordantCount + $X0)*sqrt($discordantCount + $concordantCount + $Y0));

   return $kendall;
}

sub spearman_correl{
     my  ($D, $r1, $r2, $RV) = @_;

     if( !defined( $RV->{$r1}) ) {
        rankCalculation($D, $r1, $RV);
     }

     if( !defined( $RV->{$r2}) ){
        rankCalculation($D, $r2, $RV);
     }
    
     my $correl = correl($RV->{$r1}, $RV->{$r2});
     return $correl;
}
    

sub rankCalculation {
     my ($D, $r, $RV) = @_;
     #print Dumper($D->{$r});

     my @S = sort @{$D->{$r}}; 
     #print Dumper(\@S);

     my $n = scalar(@S);
     my %R = ();
     my %N = ();
     for(my $i= 0; $i< $n; $i++) {
       if( !defined($R{$S[$i]}) ) {
          $R{$S[$i]} = $i + 1;
          $N{$S[$i]} = $i + 1;
       }
       else {
          $R{$S[$i]} = ($R{$S[$i]}*$N{$S[$i]} + ($i+1))/($N{$S[$i]} + 1);
          $N{$S[$i]} = $N{$S[$i]} + 1;
       }
     }

     for(my $i= 0; $i< $n; $i++) {
        ${$RV->{$r}}[$i] = $R{$S[$i]};
     }
     #print Dumper($RV->{$r});

}


sub has_sub_strings($) {
    my $substrings=shift;
    my $target_string=shift;
    foreach my $substring (@{$substrings} ) { 
        if(  $target_string=~/$substring/ig ) { 
           return($substring);
        }   
    }   
    return(0);
}


sub trim_otu_name{
   my $otu_name=shift;
   
   $otu_name=~ s/[a-z]__\;//ig;
   $otu_name=~ s/[a-z]__$//ig;
   $otu_name =~s/[a-z]__Unclassified\;//ig;
   $otu_name =~s/[a-z]__Unclassified\>//ig;
   $otu_name=~ s/\;$//ig;
   return $otu_name;
}


sub tax_level_id($) {
   my $taxassign=shift;
   my $level  = shift;
   my @taxids = split(";",$taxassign);
   $levelid="";
   for($i=0;$i <$level-1; $i++)  {
      $levelid=$levelid.$taxids[$i].";";
   }
   $levelid=$levelid.$taxids[$i];
   return $levelid;
}

sub tax_level($) {
   my $taxassign=shift;
   my @taxids = split(";",$taxassign);
   $level=0;
   foreach my $tax (@taxids) {
      if( $tax=~/[a-z]__[a-zA-Z0-9]/i ) {
         $level++;
      }
      else {
         last;
      }
   }
   return $level;
}


sub newhash($) {
   my $val =shift;
   my %temp_hash = ();
   for(my $i=1; $i <=150; $i++) {
       $temp_hash{$i}= $val;
   }
   return (\%temp_hash);
}
   

sub get_ec_numbers($) {
   my $var = shift;
    my @ecnos = split(" ",$var);
    return( @ecnos);

}



# Perl trim function to remove whitespace from the start and end of the string
 
#open_file  subroutine
#
#


sub newhash_from_array($) {
   my $names =shift;
   my $emptyval =shift;
   my %temp_hash = ();
   foreach my $name ( @$names) {
      if( $emptyval =~ "hash" ) {
         $temp_hash{$name} = {};
      }
      if( $emptyval =~ "array" ) {
         $temp_hash{$name} = \();
      }
      if( $emptyval =~ "double" ) {
         $temp_hash{$name} = 0;
      }
   }

   return (\%temp_hash);
}

sub fillhash_from_array($) {
   my $hashtbl =shift;
   my $names =shift;
   my $values =shift;
   my $i=0;
   foreach my $name ( @$names) {
      $hashtbl->{$name} =  $values->[$i];
      $i++;
   }
}



sub addhash_from_array($) {
   my $hashtbl= shift;
   my $names =shift;
   my $values =shift;
   my %temp_hash = ();
   my $i=0;
   foreach my $name ( @$names) {
      $hashtbl->{$name} =  $hashtbl->{$name} + $values->[$i];
      $i++;
   }
}


# Perl trim function to remove whitespace from the start and end of the string
 sub trim($)
 {
    my $string = shift;
    $string =~ s/^\s+//;
    $string =~ s/\s+$//;
    return $string;
}



sub standard_dev($) {
   my $array=shift;
   my @arraysq = qw();

   if( @${array} == 0 ) { die("Cannot compute standard deviation for an array of size 0"); }
   foreach my $val (@${array} ) {
      push( @arraysq, $val*$val);
   }
   return( sqrt( average(\@arraysq) - average($array)*average($array)) );
}

sub average($) {
   my $array=shift;
   
   if( @${array} == 0 ) { die("Cannot compute average for an array of size 0"); }
   return( sum(${array})/@${array} ); 
}


sub sum($) {
   my $array=shift;
   my $sum = 0;
   foreach my $val ( @${array}) {
      $sum = $sum + $val;
   }
   return($sum);
}

sub max($) {
   my $array=shift;
   my $max = $array->[0];
   foreach my $val ( @${array}) {
      if($val > $max) {
        $max =  $val;
      }
   }
   return($max);
}

sub min($) {
   my $array=shift;
   my $min = $array->[0];
   foreach my $val ( @${array}) {
      if($val < $min) {
        $min =  $val;
      }
   }
   return($min);
}
   

#
sub read_fasta_sequences
{
        my $filename =shift;
        my $line;
        my $first = 0;
        unless(open(DATAFILE, $filename)) 
        { 
	      print  "Cannot open file $filename \n";
          exit;
        }
        my @filedata = <DATAFILE>;
        
        my $seqname;
        my %sequences = ();
        foreach my $line (@filedata) {
              #  print $line;
                if ($line =~ /^>/)
                {
                        $seqname = trim($line);
                        $sequences{$seqname}="";
                        if ($first == 0)
                        {
                                $first = 1;
                        }
                        next;
                }
                else {
                    $sequences{$seqname} = $sequences{$seqname}.trim($line);
                }
                if ($first == 0)
                {
                        die"Not a standard FASTA file. Stop.\n";
                }
        }
        close(DATAFILE);
        return(\%sequences);
}
#
 sub get_file_data {
   my ($filename) = @_;
   use strict;
   use warnings;
   unless(open(DATAFILE, $filename)) 
   { 
	 print  "Cannot open file $filename \n";
     exit;	
   }
   my @filedata = <DATAFILE>;
   close(DATAFILE);

   return @filedata;
 }
 
 
 
#open_file  subroutine
#
#
sub open_file {
	my($filename) = @_;
	my $fh;

        unless(open(DATAFILE, $filename)) 
	{ 
		print  "Cannot open file $filename \n";
        	exit;	
        }
	return $fh;
}

sub binary_string {
      my ($array, $target) = @_;
  
      # $low is first element that is not too low;
      # $high is the first that is too high
      #
      my ( $low, $high ) = ( 0, scalar(@$array) );
   
      # Keep trying as long as there are elements that might work.
      #
      
      while ( $low < $high ) {
         # Try the middle element.
  
         use integer;
         my $cur = ($low+$high)/2;
         #   print $array->[$cur]."    ".$target."\n"; 
         if (quotemeta($array->[$cur]) lt quotemeta($target)) {
            $low  = $cur + 1;                     # too small, try higher
         } 
         else {
            $high = $cur;                         # not too small, try lower
         }
      }
      return $low;
}

sub array_unique($) {
   my $array=shift;
   my %unique_hash = ();
   foreach my $val ( @${array}) {
      $unique_hash{$val} = 1;
   }
   return( keys %unique_hash);
}


1;
__END__
# Below is stub documentation for your module. You'd better edit it!

=head1 NAME

Utilities - Perl extension for for file and string operation

=head1 SYNOPSIS

 use Utilities;

1. $word = trim($word);

2. $std_dev = standard_dev(\@array);

3. $avg = average(\@array) ;

4. $num =  sum(\@array) ;

5. @lines = read_fasta_sequences($filename);

6. @lines =  get_file_data($filename);

7. $file = open_file($filename);

8. $position = binary_string(\@array, $query);

9. $num =  max(\@array) ;

10. @arraynum =  array_unique(\@array) ;

11. $num =  min(\@array) ;

=head1 DESCRIPTION

Stub documentation for Utilities, created by h2xs. It looks like the
author of the extension was negligent enough to leave the stub
unedited.

Day to day utilites for string operation, simple statistics calculations, etc.

=head2 EXPORT

None by default.


=head1 SEE ALSO

Mention other useful documentation such as the documentation of
related modules or operating system documentation (such as man pages
in UNIX), or any relevant external documentation such as RFCs or
standards.

If you have a mailing list set up for your module, mention it here.

If you have a web site set up for your module, mention it here.

=head1 AUTHOR

Kishori Konwar, E<lt>kishori@apple.comE<gt>

=head1 COPYRIGHT AND LICENSE

Copyright (C) 2011 by Kishori Konwar

This library is free software; you can redistribute it and/or modify
it under the same terms as Perl itself, either Perl version 5.10.0 or,
at your option, any later version of Perl 5 you may have available.


=cut

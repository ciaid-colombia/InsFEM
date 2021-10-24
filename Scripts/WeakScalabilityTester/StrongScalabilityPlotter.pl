#!/usr/bin/perl
#Usage: WeakScalabilityTester.pl Casename nproc1 nproc2 nproc3 nproc4

use File::Basename;
use Getopt::Long;
use Digest::MD5 qw(md5 md5_hex md5_base64);
use File::Basename;
use Cwd;
use Recursive qw(fcopy rcopy dircopy fmove rmove dirmove);
use constant SUCESS => 0;        # Code for succesful execution.
use constant ERROR => 1;        # Code for error in execution.
use Term::ANSIColor;
use File::Path qw(mkpath rmtree);
#use Archive::Tar;
use Time::HiRes;
use Scalar::Util qw(looks_like_number);
use Naturally;
use Data::Dumper;
#use Sort::Key::Natural qw( natsort );

my $casename = $ARGV[0];
my $numtests = @ARGV -1;

print "$numtests \n";

   my %table;
   my @arr;

for (my $i = 1;$i <= $numtests; $i++) {
   print "$i \n";
   
   $processors = "Processors";
   push (@{$table{$processors}}, $ARGV[$i]);
   
   $newdir = "$casename.gid.$ARGV[$i]";
   chdir $newdir;
   #chdir "results";
   
   my $filename = "results/$casename"."0.log";
   print "$filename \n";

   
   $olddate = "SUMMARY OF COMPUTING TIMES";
   
   open(IS, $filename);
   $readline = <IS>;
   print "$readline \n";
   

   $check = 0;
   while(($check < 1) && ($readline = <IS>))
   {
       print "$readline";
       if(index($readline,$olddate) != -1){
           $check = 1;
       }
   }
   
   #if ($ARGV[$i] != 1) {
   #   $check = 0;
   #   while(($check < 1) && ($readline = <IS>))
   #   {
   #       print "$readline";
   #       if(index($readline,$olddate) != -1){
   #           $check = 1;
   #       }
   #   }
   #}
   
   #read two additional lines
   $readline = <IS>;
   $readline = <IS>;
   

   
   while ($readline = <IS>) {
      $readline =~ s/^\s+//;
      $readline =~ s/\s+$//;
      if (($readline ne "") && (index($readline,"END OF ANALYSIS") == -1 )) {
         ($concept, $time) = split(":",$readline);
         $time =~ s/^\s+//;
         $concept =~ s/^\s+//;
         $concept =~ s/\s+$//;

         ($time) = split(" ",$time);
         print "The concept $concept and the time $time \n";  
         push (@{$table{$concept}}, $time);
      }
   }
   
   chdir "..";
} 



print "hello \n";



#print Dumper(\%table);

for $family ( keys %table ) {
    print "$family: @{ $table{$family} }\n";
}
print "goodbye \n";

my $file = "gnuplot_file";

my @auxprocessors = @{ $table{"Processors"} };
# Call GNUPLOT for scalability plots
# POSTSCRIPT
for $family ( keys %table ) {
    my @auxarray = @{ $table{$family} };

open (GNUPLOT, "|gnuplot");
print GNUPLOT <<EOPLOT;
set term post color "Courier" 12
set output "$family.ps"
set title '$family'
unset key
set xlabel "Number of Processors" font "Courier,14"
set ylabel "Time" font "Courier,14"
EOPLOT
my $maxproc = 0;
my $maxtime = 0;
for my $i (0 .. $#auxarray)
{
   $auxarray[$i] = 1.0/$auxarray[$i]

   $maxproc = [ $maxproc => $auxprocessors[$i] ] -> [ $maxproc <= $auxprocessors[$i] ] ;
   $maxtime = [ $maxtime => $auxarray[$i] ] -> [ $maxtime <= $auxarray[$i] ] ;
}
$maxtime = $maxtime*1.2;
$maxproc = $maxproc*1.2;
print "maxtime maxproc $maxtime $maxproc \n";
print GNUPLOT <<EOPLOT;
set xrange [0:$maxproc]
set yrange [0:$maxtime]
plot '-' using 1:2 title '$family'
EOPLOT
#plot @{ $table{$concept}};  @{ $table{"Processors"}}
for my $i (0 .. $#auxarray)
{
   print GNUPLOT "$auxprocessors[$i] $auxarray[$i] \n";
   print "$auxprocessors[$i] $auxarray[$i] \n";
  
}
print GNUPLOT <<EOPLOT;
EOPLOT
close(GNUPLOT);

}

system("gs -sDEVICE=pdfwrite -dNOPAUSE -dBATCH -dSAFER -sOutputFile=".$casename."_ScalabilityReport.pdf *.ps");
system("rm *ps");



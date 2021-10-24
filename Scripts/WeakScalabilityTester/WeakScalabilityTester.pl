#!/usr/bin/perl
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
#use Sort::Key::Natural qw( natsort );

my $casename = "$ARGV[0].gid";
my $numtests = @ARGV -1;

print "$numtests \n";

for (my $i = 1;$i <= $numtests; $i++) {
   print "$i \n";
   $newdir = "$casename.$ARGV[$i]";
   print "$newdir \n";
   dircopy($casename,$newdir);
   chdir $newdir;
   
   $filename = "run.sh";
   $tempfile = "b.txt";
   $olddate = "#BSUB -n";
   $newdate = "#BSUB -n $ARGV[$i]";
   

   open(IS, $filename);
   open(OS, ">$tempfile");
   while(<IS>)
   {
       if($_ =~ /^$olddate(.*)$/){
           print OS "$newdate\n";
       } else {
           print OS $_;
       }
   }
   close(IS);
   close(OS);
   unlink($filename);
   rename($tempfile, $filename);
   
   system("bsub < run.sh");
   
   chdir "..";
} 

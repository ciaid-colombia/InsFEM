#!/usr/bin/perl

use File::Find;
use strict; use warnings;

use Cwd;
   
   
   
my $dir = getcwd;
my @content;
find( \&wanted, './');
#do_something_with( @content );




sub wanted {
   if ($File::Find::name =~ /test/ and $File::Find::name =~ /\/data\// and $File::Find::name !~ /\.geo\./) {
      
      if ($File::Find::name =~ /\.dat/ or $File::Find::name =~ /\.dom/) {
         push @content, $File::Find::name;
         print "$File::Find::name\n";
         
         my $filename=$File::Find::name;
         print "$filename\n";
         
         #open(FILE, '<', $filename) or die "File not found $filename";
         #open FILE, '<', "$filename"  or die "File not found $filename";
         
         open FILE, '<', "$dir\/$filename"  or die "File not found $dir\/$filename";
         my @lines = <FILE>;
         close(FILE);

         my @newlines;
         foreach(@lines) {
            $_ =~ s/data\/Fluid\///g;
            $_ =~ s/data\/Solid\///g;            
            $_ =~ s/data\///g;
            push(@newlines,$_);
         }

         open(FILE, ">",  "$dir\/$filename") or die "File not found $dir\/$filename";
         print FILE @newlines;
         close(FILE);
         
      }
   } 
  #print "$File::Find::name\n";
  return;
}



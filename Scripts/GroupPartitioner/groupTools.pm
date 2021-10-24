#!/usr/bin/perl -w
package groupTools;
use Exporter;

our @ISA = ( Exporter );
our @EXPORT_OK = qw (getLine getGroupInfo deletefilesindir trim);

sub getLine {
    my ($fileH, $filename, $word) = @_;
    my (@words,@lineword);

    my $i = 0;
    @dat=<$fileH>;
    foreach (@dat){ 
        $stringaux = $_;
        $stringaux =~ s/^\s+//; 
        @words = (split(/ {1,}/, $stringaux )) ;
        if (exists($words[0])) {
            if ($words[0]=~/$word/i) {

                #push @linenum , $.;
                push @lineword, $_;

                ++$i;
            }     
        }

    }
    undef $words;
    return @lineword;

}
sub  trim { my $s = shift; $s =~ s/^\s+|\s+$//g; return $s }

sub getGroupInfo {

    my $casename = $_[0];
    my $mygroupinfo = "$casename/Group.info";
    my $groupname;
    my $IN;
    open($IN,$mygroupinfo);
    my @lineword = getLine($IN,"$mygroupinfo","GROUPNAME");
    my $i = 0;
    foreach (@lineword){ 
        $lineword[$i] =~ s/^\s+//; 
        @words = (split(/ {1,}/, $lineword[$i])) ;
        $groupname[$i] = trim($words[1]);
        ++$i;
    }
    close($IN);
    undef @lineword;
    return @groupname;

}

sub deletefilesindir
{
   my $dir = $_[0];
    unlink glob($dir.'*.post.res');
   opendir (DIR,$dir);
   @files = readdir (DIR);
   closedir (DIR);
   foreach $file (@files)
   {
       next if ($file eq "."); 
       next if ($file eq "..");
       #unlink $file or die "can't delete $file";
       my $fil = "$dir/$file";
       system ("rm $fil") and die "Couldn't remove $file";
   } 
   return 0;
}


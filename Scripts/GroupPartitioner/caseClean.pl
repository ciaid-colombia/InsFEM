#!/usr/bin/perl -w
use Cwd 'abs_path';
use FindBin;
use lib "$FindBin::Bin/../GroupPartitioner";
use groupTools qw (getLine getGroupInfo deletefilesindir trim);

my $casename   = $ARGV[0];
my @cleanfolders = @ARGV;

my $casepath   = abs_path();

my @groupname;

if (-e "$casepath/Group.info")  { 
    #next if ($casepath eq ".");   # skip the current directory entry
    @groupname = getGroupInfo($casepath);
    #print "groupname: @groupname \n";
}
for (my $j=0; $j<@cleanfolders;$j++){
    if (-d "$casepath/$cleanfolders[$j]") {
        if(@groupname) {
            for (my $i=0; $i<@groupname;$i++){
                my $groupresdir = $casepath.'/'.$cleanfolders[$j].'/'.$groupname[$i];
                #print "should enter $groupname[$i] group in $groupresdir\n";
                if (-d $groupresdir) {
                    deletefilesindir($groupresdir);
                    #print "cleaning $groupname[$i] group results \n";
                }
            }
        }
        else {
            deletefilesindir($cleanfolders[$j]);
            #print "cleaning results \n";
        }
    }
}
if(@groupname) {
    undef @groupname;
}

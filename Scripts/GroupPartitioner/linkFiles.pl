#!/usr/bin/perl -w
use Cwd 'abs_path';
use FindBin;
use lib "$FindBin::Bin/../GroupPartitioner";
use groupTools qw (getLine getGroupInfo deletefilesindir trim);

my $casename   = $ARGV[0];
my $postfolder = $ARGV[1];
my $ext        = $ARGV[2];

my $casepath   = abs_path();
my @groupname;

if (-e "$casepath/Group.info")  { 
    #next if ($casepath eq ".");   # skip the current directory entry
    @groupname = getGroupInfo($casepath);
}
if(@groupname) {
    for (my $i=0; $i<@groupname;$i++){
        my $grouppostdir = './'.$postfolder.'/'.$groupname[$i];
        #print "should enter $groupname[$i] group in $grouppostdir\n";
        if (-d $grouppostdir) {
            #print "entered \n";
            my $fil = $grouppostdir.'/'.$casename.'0'.$ext;
            my $fil_ = './'.$groupname[$i].''.$ext;
            if (-f $fil) {
                #print "linking $groupname[$i] group post \n";
                unlink $fil_;
                symlink ($fil,$fil_);
            }
        }
    }
}
else {
        my $grouppostdir = './'.$postfolder;
        #print "should enter $grouppostdir\n";
        if (-d $grouppostdir) {
            my $fil = "$grouppostdir/$casename.'0'.$ext";
            my $fil_ = "$casename/'.'$ext";
            if (-f $fil) {
                #print "linking post \n";
                unlink $fil_;
                symlink ($fil,$fil_);
            }
        }
}
if(@groupname) {
    undef @groupname;
}


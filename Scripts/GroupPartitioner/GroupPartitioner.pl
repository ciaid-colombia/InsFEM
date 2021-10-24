#!/usr/bin/perl -w
use FileHandle;
use List::Util qw(first);
use FindBin;
use lib "$FindBin::Bin/../GroupPartitioner";
use readGidMods qw(readMods %mod2name %name2mod);
use functionLimits qw(getLine @linenum @lineword);

#-----DESCRIPTION------

#This script divides the domain in the arranged groups through GID.
#It searches through the mesh and copies the group info to another file
#by use of the first element of each group.

#-----USAGE------

# Quit unless we have the correct number of command-line args
$num_args = $#ARGV + 1;
if ($num_args != 2) {
    print "\nUsage: GroupPartitioner.pl folder case_name.gid \n";
    exit;
}

# We have two command line args: probdir problemname

$problem = $ARGV[0];
$name    = $ARGV[1];

#------------------------------------------------------------------------------------
#Directories and files
$localdir = "$problem/data/";
if (-e $localdir) { 
    $datfile     = "$localdir$name.dat" ;   
    $geofile     = "$localdir$name.geo.dat" ;   
    $domdatfile  = "$localdir$name.dom.dat" ;   
    $domfixfile  = "$localdir$name.dom.fix" ;   
    $domsolfile  = "$localdir$name.dom.sol" ;   
    #$rstfile    = "$localdir$name.rst";   
}
else {die "Could'nt find the $localdir folder, did you run the GiD case? \n";}

#dom.dat reading
#open(my $IN,"$domdatfile") or die "Can't open $domdatfile, did you run the GiD case?";
$IN = FileHandle->new("<$domdatfile");
my @groupname;
my @grouptotalel;
my @groupTypeElem;
my @grouptotalnod;

#------Get group info from dot.dat file
getLine($IN,"$domdatfile","GROUPNAME");
my $i = 0;
foreach (@lineword){ 
    $lineword[$i]       = ~ s/^\s+//; 
    @words              = (split(/ {1,}/, $lineword[$i])) ;
    $groupname[$i]      = $words[1];
    #$groupTypeElem[$i] = $words[3];
    $grouptotalel[$i]   = $words[5];
    $grouptotalnod[$i]  = $words[7];
    ++$i;
}
$IN->close;

# Clean first
# modified problem directories and dat files
my $numFiles =@groupname;

#FIX THIS LATER FOR PRETTYNESS
#IF NOT GROUP RUN BUT HAS ROM CRAP
if($numFiles == 0) {

    readMods(1);

    $m2n = \%mod2name;
    $n2m = \%name2mod;

    $localdir = "$problem/data/";
    $groupdir = $localdir;
    $IN = FileHandle->new("<$localdir$name.dat");
    $ROM_OUT = FileHandle->new(">$groupdir/$name.$mod.rom.dat");
    $ROMSOL_OUT = FileHandle->new(">$groupdir/$name.$mod.rom.sol");
    my ($begin,$end,$mygroup,@interp_data);
    while (<$IN>) {
        $stringaux = $_;
        $stringaux =~ s/^\s+//; 
        if (/PROBLEM_DATA/../END_PROBLEM_DATA/) {

            @words = (split('_', $stringaux )) ;
            my $name = $words[0];
            if (defined $name) {
                if($name =~ /PODROM/i) {
                    $mod = ${$n2m} {$words[1]};                   #We find the module nickname from its name
                    $name = $name.'_'.$words[1];
                }
                else {
                    $mod = ${$n2m} {$name};                       #We find the module nickname from its name
                }
                if (defined $mod) {
                    $mygroup = ${$m2g} {$mod};                    #We find the group it belongs to
                    $begin = $name.'_PROBLEM';
                    $end   = 'END_'.$name.'_PROBLEM';
                }     
            }     

            if(defined $begin && defined $end) {                  #If MOD_PROBLEM and END are defined
                if (/$begin/../$end/) {
                    if ($mygroup eq $groupname[$i]) {             #If we are in the correct folder 
                        if (/REDUCED_ORDER_PROBLEM/../END_REDUCED_ORDER_PROBLEM/) {
                            if (/#ROM_SOLVER_OPTIONS/../#END_ROM_SOLVER_OPTIONS/) {
                                print $ROMSOL_OUT "$stringaux";          #So rom.dat can be created
                            }
                            else {print $ROM_OUT "$stringaux";}          #So rom.dat can be created
                        }
                    }
                    if ($_ =~ /$end/i) {
                        undef $end;
                        undef $begin;
                        undef $mygroup;
                    }
                }
            }
        }     
    }
    $IN->close;
    $ROM_OUT->close;
    $ROMSOL_OUT->close;
    if (-z "$groupdir/$name.rom.dat") {system ("rm $groupdir/$name.rom.dat") and die "Couldn't remove emtpy ROM file: $groupdir/$name.rom.dat";}
    if (-z "$groupdir/$name.rom.sol") {system ("rm $groupdir/$name.rom.sol") and die "Couldn't remove emtpy ROMSOL file: $groupdir/$name.rom.sol";}
    exit;
}

my $filname= "Group.info";
if (-e $filname){
    unlink "$problem/$filname";}
my $OUT = FileHandle->new(">$problem/$filname");
for ($i=0; $i<$numFiles;$i++){
    print $OUT "GROUPNAME: $groupname[$i]\n";
}
$OUT->close;
undef $OUT;

for ($i=0; $i<$numFiles;$i++){
    $localdir  = "$problem/data/$groupname[$i]";
    $postdir   = "$problem/post/$groupname[$i]";
    $resultdir = "$problem/results/$groupname[$i]";
    #$rstdir = "$problem/rst/$groupname[$i]";
    if (-e $localdir)  { system ("rm -r $localdir");}
    if (-e $postdir)   { system ("rm -r $postdir");}
    if (-e $resultdir) { system ("rm -r $resultdir");}
    #if (-e $rstdir) { system ("rm -r $rstdir");}
}


# Generic files that dont need modifying
# and Process dom.dat folder
for ($i=0; $i<$numFiles;$i++){
    $localdir = "$problem/data/";
    $IN = FileHandle->new("<$problem/$groupname[$i].FemussGroup");
    @dat=<$IN>;
    foreach (@dat){ 
        if (/ELEMENTS NEWFORMAT RENUMBERED/../END_ELEMENTS NEWFORMAT RENUMBERED/) {
            next if /ELEMENTS NEWFORMAT RENUMBERED/ || /END_ELEMENTS NEWFORMAT RENUMBERED/;
            $stringin = $_;
            $stringaux = $stringin;
            $stringaux =~ s/^\s+//;
            @words = (split(/ {1,}/, $stringaux )) ;
            if (not first{$_ eq $words[1]} @groupTypeElem) {
                push @groupTypeElem,$words[1];
            }
        }
    }
    $IN->close;
    #print "The mesh of $groupname[$i] has @groupTypeElem element types\n";

    $groupdir       = "$problem/data/$groupname[$i]";
    $postgroupdir   = "$problem/post/$groupname[$i]";
    $resultgroupdir = "$problem/results/$groupname[$i]";
    $rstgroupdir = "$problem/rst/$groupname[$i]";
    system ("mkdir $groupdir");
    system ("mkdir $postgroupdir");
    system ("mkdir $resultgroupdir");
    if (-e $rstgroupdir) { } else { 
        system ("mkdir $rstgroupdir");}
    system ("cp $problem/data/$name.dom.sol $groupdir/$name.dom.sol");
    system ("cp $problem/data/$name.dom.fix $groupdir/$name.dom.fix");
    $IN = FileHandle->new("<$localdir$name.dom.dat");
    $OUT = FileHandle->new(">$groupdir/$name.dom.dat");
    @dat=<$IN>;
    foreach (@dat){ 
        $stringin = $_;
        $stringaux = $stringin;
        $stringaux =~ s/^\s+//;
        @words = (split(/ {1,}/, $stringaux )) ;
        if (exists($words[0])) {
            if ($words[0] =~ /ELEMENTS=/i) {
                print $OUT "ELEMENTS=                 $grouptotalel[$i]\n";
            }
            elsif ($words[0] =~ /NODAL_POINTS=/i) {
                print $OUT "NODAL_POINTS=             $grouptotalnod[$i]";
            }
            elsif ($words[0] =~ /NODES_PER_ELEMENT=/i) {
                print $OUT "NODES_PER_ELEMENT=             @groupTypeElem\n";
            }
            elsif ($words[0] =~ /INCLUDE/i) {
                if (index($words[1],".geo.dat") != -1){
                    #print $OUT "INCLUDE  data/$groupname[$i]/$name.geo.dat\n";}
                    print $OUT "INCLUDE  $name.geo.dat\n";}
                elsif (index($words[1],".dom.fix") != -1){
                    #print  $OUT "INCLUDE  data/$groupname[$i]/$name.dom.fix\n";}
                    print  $OUT "INCLUDE  $name.dom.fix\n";}
            }
            else {print $OUT "$stringaux";}
        }
    }
    $IN->close;
    $OUT->close;
    undef @groupTypeElem;

#Move renumbered meshes and boundary files :geo.dat  and dom.fix
#geo.dat 
    $localdir = "$problem/data/$groupname[$i]/";
    system ("mv $problem/$groupname[$i].FemussGroup $localdir/$name.geo.dat" or die "Could'nt move geo.dat to $localdir, did you run the GID case with mesh grouping? ");

    #dom.fix
    if (-e "$problem/$groupname[$i].FemussGroupBoundary") {
        system ("mv $problem/$groupname[$i].FemussGroupBoundary $localdir/$name.dom.fix" or die "Could'nt move dom.fix to $localdir, did you run the GID case with mesh grouping? ");


    }
}

#Read the GID problemtype modules and maps their prefixes
readMods(1);

$m2n = \%mod2name;
$n2m = \%name2mod;

#Check which modules exist in the problem and to which group they 
#belong to
my %mod2group;
my $enter = 0;
$localdir = "$problem/data";
$IN = FileHandle->new("<$datfile");
#getFunctionLimits($IN,"$datfile","PROBLEM_DATA");
my ($mod,$mod2);
while (<$IN>) {
    $stringaux = $_;
    $stringaux =~ s/^\s+//; 
    @words = (split(/ {1,}/, $stringaux )) ;
    if (/PROBLEM_DATA/../END_PROBLEM_DATA/) {
        #next if /"PROBLEM_DATA"/ || /"END_PROBLEM_DATA"/;
        @words = (split(/_/, $stringaux )) ;
        if (defined $words[0]) {
            if (exists ${$n2m} {$words[0]}) {
                $enter = 1;
                my $name = $words[0];
                $mod = ${$n2m} {$name};
                if($name =~ /PODROM/i) {
                    $mod2 = ${$n2m} {$words[1]};
                }
            } 
            if($enter){ 
                @words = (split(/ {1,}/, $stringaux )) ;
                if ($words[0]=~/GROUP/i) {
                    my $group = $words[1];
                    $group =~ s/\R//g;
                    $mod2group{$mod}=$group;
                    if (defined $mod2) {$mod2group{$mod2}=$group;}
                    $enter = 0;
                    undef $mod;
                    undef $mod2;
                }     
            }        
        }
    }
}   
$IN->close;

#Move module files to group folder
$m2g = \%mod2group;
$localdir = "$problem/data/";
@alldirs="$localdir";
foreach $dir (@alldirs) {
    opendir(SDIR,$dir) or die "Couldn't open $localdir";
    @dataentries=readdir(SDIR);
    foreach (@dataentries) {
        if( ! ($_=~/^\./)  ) {
            $suff = $_;
            $filename = $_;
            $suff =~ s/$name//;   #drop suffix
            $suff =~ s/\./ /;   #drop dot
            $suff =~ s/^\s+//; 
            @word = (split(/\./, $suff)) ;
            if (exists ${$m2n} {$word[0]}) { #If file is a module file

                my $modul = $word[0];
                my $group = ${$m2g} {$modul};
                $groupdir = "$localdir$group";

                #Process name.module.dat for boundary conds
                if ($_ =~ "$name.$modul.dat") {
                    $IN = FileHandle->new("<$localdir$name.$modul.dat");
                    $OUT = FileHandle->new(">$groupdir/$name.$modul.dat");
                    @dat=<$IN>;
                    foreach (@dat){  #For each file in the folder
                        $stringin = $_;
                        $stringaux = $stringin;
                        $stringaux =~ s/^\s+//;
                        @words = (split(/ {1,}/, $stringaux )) ;
                        if (/BOUNDARY_CONDITIONS/../END_BOUNDARY_CONDITIONS/) {
                            if (exists($words[0])) {
                                if ($words[0] =~ /INCLUDE/i) {
                                    if (index($words[1],".$modul.fix") != -1){
                                        #print  $OUT "INCLUDE  data/$group/$name.$modul.fix\n";}
                                        print  $OUT "INCLUDE  $name.$modul.fix\n";}
                                }
                                else {print $OUT "$stringaux";}
                            }
                        }
                        else {print $OUT "$stringaux";}
                    }
                    $IN->close;
                    $OUT->close;
                    system ("rm $localdir$name.$modul.dat" or die "Couldn't remove $localdir$name.$modul.dat");
                } 
                else {#If file not a module file we just move it without editing

                    system ("mv $localdir$_ $groupdir/." or die "Could'nt create the file $filename inside $groupdir, did you run the GID case with mesh grouping? ");
                }
            }
        }
    }
    closedir(SDIR);
}

$localdir = "$problem/data/";
for ($i=0; $i<$numFiles;$i++){
    $groupdir = "$problem/data/$groupname[$i]";
    $IN = FileHandle->new("<$localdir$name.dat");
    $OUT = FileHandle->new(">$groupdir/$name.dat");
    $ROM_OUT = FileHandle->new(">$groupdir/$name.$mod.rom.dat");
    $ROMSOL_OUT = FileHandle->new(">$groupdir/$name.$mod.rom.sol");
    my ($begin,$end,$mygroup,@interp_data);
    while (<$IN>) {
        $stringaux = $_;
        $stringaux =~ s/^\s+//; 
        if (/PROBLEM_DATA/../END_PROBLEM_DATA/) {

            @words = (split('_', $stringaux )) ;
            my $name = $words[0];
            if (defined $name) {
                if($name =~ /PODROM/i) {
                    $mod = ${$n2m} {$words[1]};                   #We find the module nickname from its name
                    $name = $name.'_'.$words[1];
                }
                else {
                    $mod = ${$n2m} {$name};                       #We find the module nickname from its name
                }
                if (defined $mod) {
                    $mygroup = ${$m2g} {$mod};                    #We find the group it belongs to
                    $begin = $name.'_PROBLEM';
                    $end   = 'END_'.$name.'_PROBLEM';
                }     
            }     

            if(defined $begin && defined $end) {                  #If MOD_PROBLEM and END are defined
                if (/$begin/../$end/) {
                    if ($mygroup eq $groupname[$i]) {             #If we are in the correct folder 
                        if (/REDUCED_ORDER_PROBLEM/../END_REDUCED_ORDER_PROBLEM/) {
                            if (/#ROM_SOLVER_OPTIONS/../#END_ROM_SOLVER_OPTIONS/) {
                                print $ROMSOL_OUT "$stringaux";          #So rom.dat can be created
                            }
                            else {print $ROM_OUT "$stringaux";}          #So rom.dat can be created
                        }
                        elsif (/INTERPOLATE_DATA/../END_INTERPOLATE_DATA/) {
                            push @interp_data, $stringaux;        #So we can then print at EOF
                        }
                        else {print $OUT "$stringaux";}
                    }
                    if ($_ =~ /$end/i) {
                        undef $end;
                        undef $begin;
                        undef $mygroup;
                    }
                }
            }
            elsif(!defined $mygroup) {
                print $OUT "$stringaux";                          #Print whatever doesnt belong to modules
            }   
        }     
        else {print $OUT "$stringaux";}                           #If not in PROBLEM_DATA we print the rest
    }
    $IN->close;
    $OUT->close;
    if(@interp_data) {
        $OUT = FileHandle->new("$groupdir/$name.dat",O_WRONLY|O_APPEND);
        print $OUT "@interp_data";
        $OUT->close;

    }
    $ROM_OUT->close;
    $ROMSOL_OUT->close;
    if (-z "$groupdir/$name.rom.dat") {system ("rm $groupdir/$name.rom.dat") and die "Couldn't remove emtpy ROM file: $groupdir/$name.rom.dat";}
    if (-z "$groupdir/$name.rom.sol") {system ("rm $groupdir/$name.rom.sol") and die "Couldn't remove emtpy ROMSOL file: $groupdir/$name.rom.sol";}
}

system ("rm $geofile")    and die "Couldn't remove $geofile";
system ("rm $domdatfile") and die "Couldn't remove $domdatfile";
system ("rm $domfixfile") and die "Couldn't remove $domfixfile";
system ("rm $domsolfile") and die "Couldn't remove $domsolfile";
system ("rm $datfile")    and die "Couldn't remove $datfile";


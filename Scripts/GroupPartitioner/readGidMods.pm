#!/usr/bin/perl -w
package readGidMods;
use Exporter;
use experimental 'smartmatch';

our @ISA = ( Exporter );
our @EXPORT_OK = qw (readMods %mod2name %name2mod );

sub readMods {

    my (@savedmods, @modNum);# %mod2name, %name2mod );

    $fem_prb="$ENV{'GID_PRB'}";
    opendir(SDIR,$fem_prb);
    @entries=readdir(SDIR);
    foreach $dir (@entries) {
        @modwords = split m/\-/, $dir;
        my $num = $modwords[0];
        if ("$num" =~ /^[+-]?\d+$/ && $num > 0){
            if ("$num" ~~ @modNum){}
            else {
                push(@modNum,$modwords[0]);
                push(@savedmods,$dir);
            }
        }
    }
    close(SDIR);

    foreach $dir (@savedmods) {

        $localdir = "$fem_prb/$dir";
        open(MDIR,$localdir) or die "Could'nt open $localdir"; 
        @lines = <MDIR>;
        $lin = $lines[2]; 
        $lin =~ s/^\s+//; 
        @words = (split(/ {1,}/, $lin)) ;
        if (defined $words[3] and $words[3] ne '' and defined $words[4] and $words[4] ne ''){
            #print "  prefixes : $words[3], $words[4]\n";
            $mod2name{$words[3]}= $words[4];
            $name2mod{$words[4]}= $words[3];
        }
    }
    close(MDIR);
}

1;

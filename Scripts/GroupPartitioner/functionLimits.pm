#!/usr/bin/perl -w
package functionLimits;
use Exporter;

our @ISA = ( Exporter );
our @EXPORT_OK = qw (getLine @linenum @lineword);

sub getLine {
    my ($fileH, $filename, $word) = @_;
    my @words;

    my $i = 0;
    @dat=<$fileH>;
    foreach (@dat){ 
        $stringaux = $_;
        $stringaux =~ s/^\s+//; 
        @words = (split(/ {1,}/, $stringaux )) ;
        if (exists($words[0])) {
            if ($words[0]=~/$word/i) {

                push @linenum , $.;
                push @lineword, $_;

                ++$i;
            }     
        }

    }

}

1;


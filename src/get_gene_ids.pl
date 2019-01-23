#!/usr/bin/perl
use strict;

# script to read a cDNA file, and pair transcript and gene IDs

my %data = ();


open (inFile, $ARGV[0]) || die "Can't open file $!\n";

thirdTry();

close inFile;


sub thirdTry{
	while (my $line = <inFile>){
		next unless ($line =~ /^\>([\w-.]*)\s+/);     # Don't want sequence lines

		my $transcript = $1;

		if ($line =~ /gene:(\w+)/){
	    print "$transcript\t$1\n";
		} else {
	    die "Something wrong with:\n$line\n\n";
	}
    }
}

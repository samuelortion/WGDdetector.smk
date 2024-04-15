#!/usr/bin/env perl
# ref. phase.mcloutp2orthomcl.format.pl from WGDdetector, Yang et al. 2019
# Format the MCL output file into the OrthoMCL format
use strict;
use warnings;

my ($in,$out)=@ARGV;
# die "perl $0 input_mcl_raw_result output_phased_mcl_result\n" if ! $out;
my $num=0;
open (F,"$in") || die"no such file: $in\n";
open (O,">$out") || die "cannot create file: $out\n";
while (<F>) {
    chomp;
    $num++;
    print O "cluster$num:\t$_\n";
}
close F;
close O;
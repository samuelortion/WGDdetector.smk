#!/usr/bin/env perl
# ref. collect_ks.pl from WGDdetector, Yang et al. 2019
# Write ks.dist

use strict;
use warnings;

my $help = "perl $0 output_dist input_ks.gz.1 input_ks.gz.2 input_ks.gz.3 ...\n";
my $output_dist=shift or die "$help";
my @result = @ARGV;
die "$help" if scalar(@result)<1;

my %ks;
my %id;
for my $result (@result){
    open (F,"zcat \"$result\"|")||die"$!";
    while (<F>) {
        chomp;
        my @a=split(/\s+/,$_);
        $ks{$a[0]}{$a[1]}=$a[2];
        $id{$a[0]}++;
        $id{$a[1]}++;
    }
    close F;
}

my @id=sort keys %id;
open (O,">$output_dist")||die"$!";
print O "\t",join("\t",@id),"\n";
for my $k1 (@id){
    print O "$k1\t";
    for my $k2 (@id){
        my $out=10;
        if (exists $ks{$k1}{$k2}){
            $out=$ks{$k1}{$k2};
        }elsif (exists $ks{$k2}{$k1}){
            $out=$ks{$k2}{$k1};
        }elsif ($k1 eq $k2){
            $out=0;
        }
        print O "$out\t";
    }
    print O "\n";
}
close O;

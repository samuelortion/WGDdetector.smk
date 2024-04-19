#!/usr/bin/env perl
# ref. phase.id.pl from WGDdetector, Yang et al. 2019
#
# Usage example:
# perl phase.id.pl input.cds.fa input.pep.fa output.cds.fa output.pep.fa output.id.table

use strict;
use warnings;
use Bio::SeqIO;

my ($incds, $inpep, $out_filtered_cds, $out_filtered_pep, $out_seqid_gennum_mapping)=@ARGV;
die "perl $0 input_cds input_pep filtered_cds filtered_pep seqid_gennum_mapping \n" if ! $out_seqid_gennum_mapping;

my %cds=&read_fasta("$incds");
my %pep=&read_fasta("$inpep");

open (O1,">$out_filtered_cds");
open (O2,">$out_filtered_pep");
open (O3,">$out_seqid_gennum_mapping");
my $num=0;
for my $k1 (sort keys %cds){
    next if ! exists $pep{$k1};
    $num++;
    print O1 ">gene$num\n$cds{$k1}\n";
    print O2 ">gene$num\n$pep{$k1}\n";
    print O3 "$k1\tgene$num\n";
}
close O1;
close O2;
close O3;

sub read_fasta{
    my ($tmp_in_file)=@_;
    my $tmp_fa=Bio::SeqIO->new(-format=>"fasta",-file=>"$tmp_in_file");
    my %r;
    while (my $tmp_seq=$tmp_fa->next_seq) {
        my $tmp_id=$tmp_seq -> id;
        my $tmp_seq=$tmp_seq -> seq;
        $r{$tmp_id}=$tmp_seq;
    }
    return %r;
}

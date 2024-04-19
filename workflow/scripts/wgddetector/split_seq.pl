#!/usr/bin/env perl
# ref. split_seq.pl from WGDdetector, Yang et al. 2019
# Split the input pep anf cds file into a single fasta file per cluster 

use strict;
use warnings;

use Bio::SeqIO;

my ($cds_file,$pep_file,$cluster_file,$tmpdir) = @ARGV;
die "perl $0 cds_file pep_file cluster_file tmpdir \n" if (! $tmpdir);

`mkdir $tmpdir` if (! -e $tmpdir);
`mkdir $tmpdir/align` if (! -e "$tmpdir/align");
`mkdir $tmpdir/align_large` if (! -e "$tmpdir/align_large");

my $clustermaxnum=50; # TODO: make this variable a parameter
my %list=&read_cluster($cluster_file);
my %cds=&read_fasta($cds_file);
my %pep=&read_fasta($pep_file);

`mkdir -p $tmpdir`;
for my $k1 (sort keys %list) {
    my @gene=sort keys %{$list{$k1}};
    if (scalar(@gene) <= $clustermaxnum){
        open (O1,">$tmpdir/align/$k1.input.cds.file");
        open (O2,">$tmpdir/align/$k1.input.pep.file");
        for my $geneid (@gene){
            if (scalar(@gene) <= $clustermaxnum){
                print O1 ">$geneid\n$cds{$geneid}\n";
                print O2 ">$geneid\n$pep{$geneid}\n";
            }
        }
        close O1;
        close O2;
    }else{
        `mkdir -p $tmpdir/align_large/$k1`;
        open (O1,">$tmpdir/align_large/$k1/input.cds.file");
        open (O2,">$tmpdir/align_large/$k1/input.pep.file");
        for my $geneid (@gene){
            print O1 ">$geneid\n$cds{$geneid}\n";
            print O2 ">$geneid\n$pep{$geneid}\n";
        }
        close O1;
        close O2;
    }
}

sub read_fasta {
    my ($tmp_in_fasta_file)=@_;
    my $tmp_fa=Bio::SeqIO->new(-format=>"fasta",-file=>"$tmp_in_fasta_file");
    my %r_fasta;
    while (my $tmp_seq=$tmp_fa->next_seq) {
        my $tmp_id=$tmp_seq -> id;
        my $tmp_seq=$tmp_seq -> seq;
        $r_fasta{$tmp_id}=$tmp_seq;
    }
    return %r_fasta;
}

sub read_cluster {
    my ($tmp_in_cluster_file)=@_;
    my %r_cluster;
    open (F,"$tmp_in_cluster_file")||die"$!";
    while (<F>) {
        chomp;
        my @a=split(/\s+/,$_);
        $a[0]=~s/\:$//;
        for (my $i=1;$i<@a;$i++){
            $r_cluster{$a[0]}{$a[$i]}++;
        }
    }
    close F;
    return %r_cluster;
}

#!/usr/bin/env perl
# ref. calculate_ks.pl and calculate_ks.single.pl from WGDdetector, Yang et al. 2019
# Compute Ks metrics
# A merge of the original WGDdetector calculate_ks.pl and calculate_ks.single.pl in a single Perl script

# use strict;
use warnings;
use Bio::SeqIO;
use Bio::Align::DNAStatistics;
use Bio::AlignIO;
# use threads; # TODO: Do I really need this?

my ($cds_align, $ks_file) = @ARGV;
die "perl $0 cds_align output_ks_file_gz \n" if (!$ks_file);

# Map a seqid to a DNA sequence
my %seq;
my $cds_aln_obj=Bio::SeqIO->new(-format=>"fasta",-file=>$cds_align);
# Map a seqid to a DNA sequence
while (my $seqobj=$cds_aln_obj->next_seq) {
    my $id=$seqobj->id;
    my $seq=$seqobj->seq;
    $seq{$id}=$seq;
}

my @spid = sort keys %seq; # Sort sequences by sequence ids

# Compute Ks
open(O, "| gzip -c > $ks_file") || die "$!"; # Open a file handle that compresses the output to a .gz archive.

# For each pair of sequence id
for (my $i=0; $i < @spid; $i++) {
    for (my $j=$i+1; $j < @spid; $j++) {
        my $spid1 = $spid[$i];
        my $spid2 = $spid[$j];
        my $spseq1 = $seq{$spid1};
        my $spseq2 = $seq{$spid2};
        my ($tmptwoseqalign, $newseq1, $newseq2) = &two_seq_align($spseq1, $spseq2);
        next if $tmptwoseqalign < 90; # TODO: render this parameter configurable.
        my $Ds = &calculate_ds($spid1, $spid2, $newseq1, $newseq2);
        $Ds = 10 if $Ds !~ /^\d/; # TODO: this seems not appropriate to fall back to 10 whenever the output is not a number.
        print O "$spid1\t$spid2\t$Ds\n";
    }
}


sub two_seq_align {
    my ($tmpseq1,$tmpseq2)=@_;
    my ($ralign,$rseq1,$rseq2)=(0,"","");
    my @tmpseq1=split(//,$tmpseq1);
    my @tmpseq2=split(//,$tmpseq2);
    for (my $i=0;$i<@tmpseq2;$i++){
        if (($tmpseq1[$i] ne '-') && ($tmpseq2[$i] ne '-')){
            $rseq1 .= $tmpseq1[$i];
            $rseq2 .= $tmpseq2[$i];
            $ralign++;
        }
    }
    return ($ralign,$rseq1,$rseq2);
}

sub calculate_ds {
    my ($spid1, $spid2, $newseq1, $newseq2) = @_;
    my $tmpfh;
    open ($tmpfh,"echo \">$spid1\n$newseq1\n>$spid2\n$newseq2\n\"|"); # Open a temporary file handler for a FASTA file with the two sequences only.
    my $cds_aln_obj=Bio::AlignIO->new(-format=>"fasta",-fh=>$tmpfh);
    my $cds_aln=$cds_aln_obj->next_aln();
    my $stats = Bio::Align::DNAStatistics->new(); 
    my $results = $stats->calc_all_KaKs_pairs($cds_aln);
    my $Ds;
    for my $an (@$results){
        $Ds= $an->{D_s};
    }
    return $Ds;
}

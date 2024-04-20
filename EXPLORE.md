# Explore

Infos on the original WGDdetector perl scripts.

## Ks estimate : `02.ks_estimate`

How does this work?

This steps uses GNU parallel.
WGDdetector generates a lot of .sh scripts which are mostly just commandlines which run a subset of the task.

| script | what it does | where it has been generated | 
| ------ | ------------ | --------------------------- |
| run.ks.sh | | split_seq.pl |
| collect.ks.sh | |split_seq.pl |
| run.large.split.merge_cmd.ks.collect.sh | | phase.new_sub_GF_seq.pl|
| run.large.split.merge_cmd.ks.sh | | phase.new_sub_GF_seq.pl |
| run.large.split.sh | | split_seq.pl |

## Command decipher

### `run.ks.sh`


```bash
bin/calculate_ks.pl 4 tmp/02.ks_estimate/align/cluster1.input.cds.file tmp/02.ks_estimate/align/cluster1.input.pep.file /usr/bin/muscle /usr/bin/pal2nal.pl
```

Removing junk software path we have:

```bash
bin/calculate_ks.pl 4 tmp/02.ks_estimate/align/cluster1.input.cds.file tmp/02.ks_estimate/align/cluster1.input.pep.file
```

> What does `calculate_ks.pl` do?

First, It launches the `muscle` CDS alignment 'all against all' on `tmp/02.ks_estimate/align/cluster1.input.pep.file` (resulting in `tmp/02.ks_estimate/align/cluster1.input.pep.file.align` file).
Then it converts the protein alignment into the corresponding nucleotide alignment using PAL2NAL.
Finally, It uses Bioperl module to compute the Ks statistics for each pair of homologous sequence and writes it to the file `tmp/02.ks_estimate/align/cluster1.input.cds.file.align.output.ks.gz`

The relevant part is:

```perl
my %seq;
my $cds_aln_obj=Bio::SeqIO->new(-format=>"fasta",-file=>"$cds.align");
# Map a seqid to a DNA sequence
while (my $seqobj=$cds_aln_obj->next_seq) {
    my $id=$seqobj->id;
    my $seq=$seqobj->seq;
    $seq{$id}=$seq;
}
my @spid=sort keys %seq; # Sort sequences by seq ids
open (O,"| gzip -c >$cds.align.output.ks.gz")||die"$!"; # Open a file, to a compressed gzip archive
for (my $i=0;$i<@spid;$i++){
    for (my $j=$i+1;$j<@spid;$j++){
        my ($spid1,$spid2)=($spid[$i],$spid[$j]);
        my ($spseq1,$spseq2)=($seq{$spid1},$seq{$spid2});
        my ($tmptwoseqalign,$newseq1,$newseq2)=&two_seq_align($spseq1,$spseq2);
        next if $tmptwoseqalign < 90;
        my $Ds=`$basedir/bin/calculate_ks.single.pl \"$spid1\" \"$newseq1\" \"$spid2\" \"$newseq2\"`;
        $Ds=10 if $Ds !~ /^\d+/;
        chomp $Ds;
        print O "$spid1\t$spid2\t$Ds\n";
    }
}
close O;

# 
sub two_seq_align{
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
```

> What does `calculate_ks.single.pl` do?

It creates a sample fasta file with two sequences ids `$spid1`, `$spid2` and their sequences.
Then, it reads the sequence as a fasta file again and extract alignment statistics using Bioperl.

```perl
my $tmpfh;
open ($tmpfh,"echo \">$spid1\n$newseq1\n>$spid2\n$newseq2\n\"|");
my $cds_aln_obj=Bio::AlignIO->new(-format=>"fasta",-fh=>$tmpfh);
my $cds_aln=$cds_aln_obj->next_aln();
my $stats = Bio::Align::DNAStatistics->new(); 
my $results = $stats->calc_all_KaKs_pairs($cds_aln);
my $Ds;
for my $an (@$results){
    $Ds= $an->{D_s};
}
$Ds=10 if $Ds !~ /^\d+/;

print "$Ds\n";
```


### `collect.ks.sh`

```text
bin/collect_ks.pl output/02.ks_estimate cluster1 tmp/02.ks_estimate/align/cluster1.input.cds.file.align.output.ks.gz
```

> What does `collect_ks.pl` do?

```perl
my %ks;
my %id;
for my $result (@result){
    open (F,"zcat $result|")||die"$!";
    while (<F>) {
        chomp;
        my @a=split(/\s+/,$_);
        $ks{$a[0]}{$a[1]}=$a[2];
        $id{$a[0]}++;
        $id{$a[1]}++;
    }
    close F;
}
```

`$result` is valued to the filename `cluster<number>.unput.cds.file.align.output.ks.gz` generated in `calculate_ks.pl`

This file is structured as the following text extract shows:

```
gene6   gene62  1.01430908738512
gene61  gene62  0.020563316867574
```

In this code sample `$a[0]` is a gene id (format `gene<number`), so is `$a[1]`. `$a[2]` is a float, corresponding to the previously computed Ks estimate.

So this perl code sample stores the Ks value for each pair of gene in the .gz file in the `%ks` hashmap.

The `%id` hashmap stores how many times each gene id was found in the gz file, either in the first or the second column.


```perl
my @id=sort keys %id;
open (O,">$outputdir/ks_martix/$cluster_name/ks.dist")||die"$!";
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
        print  O "$out\t";
    }
    print  O "\n";
}
close O;
```

```perl
my @id=sort keys %id;
```
Sort the gene ids by lexicographic order.

It then writes a file `ks.dist` in a cluster subfolder.

For each pair of gene id keys `$k1` and `$k2`, if there is a Ks estimate in `%ks`, in either direction, it extracts the estimate and print it in the distance file so that it looks like the following extract:

```text
	gene15	gene16	gene33	gene34	gene35	gene7	gene71	gene72	gene73	gene76	gene77	gene87
gene15	0	0.123592998642446	2.41061912624257	10	10	10	10	10	10	1.29968922690385	1.14565873648533	1.09752871538553	
gene16	0.123592998642446	0	10	10	10	10	2.38706979103902	2.23714430296334	10	1.13586692077591	1.10324117997687	1.01738170190691	
gene33	2.41061912624257	10	0	0.368905964471772	0.247809729091584	0.697013936392055	1.44293921510591	1.49285542209054	1.80541777398976	10	10	10	
gene34	10	10	0.368905964471772	0	0.21952725262641	0.652598057480441	1.57586306043495	1.63617979203181	1.49014689612897	10	10	2.75521319759289	
gene35	10	10	0.247809729091584	0.21952725262641	0	0.59536268152535	1.79102943115131	1.84832108883335	1.77459233534853	10	10	10	
gene7	10	10	0.697013936392055	0.652598057480441	0.59536268152535	0	1.19209458412741	1.20877143730181	1.38363696955638	10	3.62123530297678	10	
gene71	10	2.38706979103902	1.44293921510591	1.57586306043495	1.79102943115131	1.19209458412741	0	0.0170075315032724	0.260397149988142	10	10	10	
gene72	10	2.23714430296334	1.49285542209054	1.63617979203181	1.84832108883335	1.20877143730181	0.0170075315032724	0	0.258087119049242	10	10	10	
gene73	10	10	1.80541777398976	1.49014689612897	1.77459233534853	1.38363696955638	0.260397149988142	0.258087119049242	0	3.92209638619892	2.87594063660687	3.30087044194534	
gene76	1.29968922690385	1.13586692077591	10	10	10	10	10	10	3.92209638619892	0	0.120499020636339	0.658172427734962	
gene77	1.14565873648533	1.10324117997687	10	10	10	3.62123530297678	10	10	2.87594063660687	0.120499020636339	0	0.785744971120637	
gene87	1.09752871538553	1.01738170190691	10	2.75521319759289	10	10	10	10	3.30087044194534	0.658172427734962	0.785744971120637	0	
```

This file will then be read by R to extract clusters with `hclust` with file  `hclust_ks.pl`.


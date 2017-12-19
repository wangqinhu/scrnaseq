#!/usr/bin/env perl
#===============================================================================
#
# Description: tpm / ercc / filter noise / tsne_mat
#
# Copyright (c) 2017 Northwest A&F University
# Author: Qinhu Wang
# Email: wangqinhu@nwafu.edu.cn
#
#===============================================================================

use strict;
use warnings;
use Data::Dumper;

my $count_file = "data/matrix/count.tsv";
my $length_file = "data/mm10/mm10.genelen.txt";
my $out_norm = "data/matrix/scRNAseq.norm.txt";
my $out_tsne = "data/matrix/scRNAseq.tsne.txt";
my $out_log = "data/log/scRNAseq_log.txt";
my @empty_well = ('W1', 'W384', 'W385', 'W768', 'W769', 'W1152', 'W1153', 'W1536', 'W1537', 'W1920');

my $max_cell = 1920;
my $empty_ercc_threshold = 100;
my $nosie_gene_threshold = 0;

my %Lf = ();
my %Mf = ();
my %rpk = ();
my %tpm = ();

my %gene = ();
my %cell = ();
my %ercc = ();
my %expr = ();

my %clean_gene = ();
my %clean_cell = ();

unlink($out_log);

# load gene length
open (LEN, $length_file) or die $!;
while (<LEN>) {
	chomp;
	next if (/^\s*$/);
	my @w = split /\t/;
	$Lf{$w[0]} = $w[1];
}
close LEN;

# cal rpk and “per million” scaling factor
open (IN, $count_file) or die $!;
while (<IN>) {
	chomp;
	next if (/^\s*$/);
	next if (/^\t/);
	my @w = split ();
	for (my $i = 1; $i <= $max_cell; $i++) {
		my $cell = 'W' . $i;
		my $gene = $w[0];
		$cell{$cell} = 1;
		$gene{$gene} = 1;
		$rpk{$cell}{$gene} = $w[$i] / ($Lf{$gene}/1000);
		$Mf{$cell} += $rpk{$cell}{$gene};
	}
}
close IN;

# cal tpm and ercc_sum
for my $cell (sort by_strnum keys %cell) {
	foreach my $gene (sort keys %gene) {
		$tpm{$cell}{$gene} = $rpk{$cell}{$gene} / ($Mf{$cell}/1000000);
		if ($gene =~ /^ERCC/) {
			$ercc{$cell} += $tpm{$cell}{$gene};
		}
	}
}

# filter empty_ercc and empty cell
for my $cell (sort by_strnum keys %cell) {
	# filter emptu_ercc cell
	next if (is_empty_ercc($cell));
	# filter emptu cell
	next if (is_empty_cell($cell));
	$clean_cell{$cell} = 1;
}

# filter ercc and noise gene
foreach my $gene (sort keys %gene) {
	# filter ercc gene
	next if ($gene =~ /^ERCC/);
	# fitler noise gene
	next if (is_noise($gene));
	$clean_gene{$gene} = 1;
}

# norm by ercc
foreach my $gene (sort keys %clean_gene) {
	for my $cell (sort by_strnum keys %clean_cell) {
		my $norm = $tpm{$cell}{$gene} / $ercc{$cell};
		$norm = int ($norm * 1000 + 0.5) /1000 * 100;
		$expr{$cell}{$gene} = $norm;
	}
}

# output norm matrix
open (NORM, ">$out_norm") or die "Cannot open $out_norm: $!\n";
for my $cell (sort by_strnum keys %clean_cell) {
	print NORM "\t", $cell;
}
print NORM "\n";
foreach my $gene (sort by_strnum keys %clean_gene) {
	print NORM "$gene";
	for my $cell (sort keys %clean_cell) {
		print NORM "\t", $expr{$cell}{$gene};
	}
	print NORM "\n";
}
close NORM;

# output tsne matrix
open (TSNE, ">$out_tsne") or die "Cannot open $out_tsne: $!\n";
for my $gene (sort keys %clean_gene) {
	print TSNE "\t", $gene;
}
print TSNE "\n";
foreach my $cell (sort by_strnum keys %clean_cell) {
	print TSNE "$cell";
	for my $gene (sort keys %clean_gene) {
		print TSNE "\t", $expr{$cell}{$gene};
	}
	print TSNE "\n";
}
close TSNE;


#=======
# subs
#=======

sub by_strnum {
	$a =~ /(\d+)/;
	my $numa = $1;
	$b =~ /(\d+)/;
	my $numb = $1;
	return $numa <=> $numb;
}

sub is_noise {
	my $gene = shift;
	foreach my $cell (@empty_well) {
		if ($tpm{$cell}{$gene} > $nosie_gene_threshold) {
			open (LOG, ">>$out_log") or die $!;
			print LOG $gene, "\t", "$cell\n";
			close LOG;
			return 1;
		}
	}
}

sub is_empty_ercc {
	my $cell = shift;
	if ($ercc{$cell} < $empty_ercc_threshold) {
		open (LOG, ">>$out_log") or die $!;
		print LOG $cell, "\t", "Empty ERCC\n";
		close LOG;
		return 1;
	}
}

sub is_empty_cell {
	my $cell = shift;
	foreach my $empty_well (@empty_well) {
		if ($empty_well eq $cell) {
			return 1;
		}
	}
}

#!/usr/bin/env perl
#===============================================================================
#
# Description: construct the expression matrix for scRNAseq
#
# Copyright (c) 2017 Northwest A&F University
# Author: Qinhu Wang
# Email: wangqinhu@nwafu.edu.cn
#
#===============================================================================

use strict;
use warnings;
use Data::Dumper;

my $dir = $ARGV[0] || "data/count";
my $out = $ARGV[0] || "data/matrix/counts.tsv";

my %count;
my %cell;
my %gene;

opendir (DIR, $dir) or die "Cannot open $dir: $!\n";
foreach my $tiss (readdir DIR) {
	next if ($tiss =~ /^\./);
	opendir (SUBDIR, "$dir/$tiss") or die "Cannot open $dir/$tiss: $!\n";
	foreach my $file (readdir SUBDIR) {
		next if ($file =~ /^\./);
		my $cell = $file;
		$cell =~ s/.count.tsv//;
		$cell{$cell} = 1;
		read_count("$dir/$tiss/$file");
	}
	close SUBDIR;
}
close DIR;

# output matrix
open (OUT, ">$out") or die "Cannot open $out: $!\n";
for my $cell (sort by_strnum keys %cell) {
	print OUT "\t", $cell;
}
print OUT "\n";
foreach my $gene (sort keys %gene) {
	print OUT $gene;
	for my $cell (sort by_strnum keys %cell) {
		if (!exists $count{$gene}{$cell}) {
			$count{$gene}{$cell} = 0;
		}
		print OUT "\t", $count{$gene}{$cell};
	}
	print OUT "\n";
}
close OUT;

sub read_count {
	my $file = shift;
	open (IN, $file) or die "Cannot open $file: $!\n";
	while (<IN>) {
		chomp;
		next if (/gene\tcell\tcount/);
		my @w = split /\t/;
		$gene{$w[0]} = 1;
		$count{$w[0]}{$w[1]} = $w[2];
	}
	close IN;
}

sub by_strnum {
	$a =~ /(\d+)/;
	my $numa = $1;
	$b =~ /(\d+)/;
	my $numb = $1;
	return $numa <=> $numb;
}
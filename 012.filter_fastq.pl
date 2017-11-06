#!/usr/bin/env perl
#===============================================================================
#
# Description: Filter poly-T and short reads, format fastq for umi_tools
#
# Copyright (c) 2017 Northwest A&F University
# Author: Qinhu Wang
# Email: wangqinhu@nwafu.edu.cn
#
#===============================================================================

use strict;
use warnings;
use Data::Dumper;

my $fastq_input = $ARGV[0] || "input.fq";
my $fastq_output = $ARGV[1] || "output.fq";
my $qual_threshold = 20;
my $poly_T_len = 10;
my $seq_len_thresh = 30;

filter_fastq();

sub filter_fastq {
	open (FQI, $fastq_input) or die "Cannot open fastq: $fastq_input, $!\n";

	my @n = split /\//, $fastq_input;
	my $tiss = $n[-2];
	my $cell = $n[-1];
	$cell =~ s/\.fq//;

	my %fastq;
	my ($fq_1, $fq_2, $fq_3, $fq_4, $rmt, $rmt_qual);
	my ($read_count, $short_read_count, $clean_count) = (0, 0, '0000000');
	while ($fq_1 = <FQI>) {
		$fq_2 = <FQI>;
		$fq_3 = <FQI>;
		$fq_4 = <FQI>;
		chomp $fq_2;
		chomp $fq_4;
		$read_count++;
		# find poly T
		if ($fq_2 =~ /(T{$poly_T_len})/) {
			$fq_2 = $`;
			# filter short reads
			if (length $fq_2 < $seq_len_thresh) {
				$short_read_count++;
				next;
			}
			$fq_4 = substr($fq_4, 0, length($fq_2));
		}
		$clean_count++;
		# umi/rmt
		($rmt, $rmt_qual) = split /\s/, $fq_1;
		$rmt =~ s/^\@//;
		$rmt = "$tiss.$clean_count\_$cell\_$rmt";
		$fastq{$rmt}{$fq_2} = $fq_4;
	}
	close FQI;

	print "$cell\t$tiss\t$read_count\t$clean_count\t$short_read_count\n";

	# write fastq
	open (FQO, ">$fastq_output") or die "Cannot open fastq: $fastq_output, $!\n";
	foreach my $rmt (sort keys %fastq) {
		foreach my $read (sort keys $fastq{$rmt}) {
			print FQO "@" . "$rmt\n$read\n+\n", $fastq{$rmt}{$read}, "\n";
		}
	}
	close FQO;
}


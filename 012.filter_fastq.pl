#!/usr/bin/env perl
#===============================================================================
#
# Description: Filter low quanlity fastq
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
my $qual_offset = 33;
my $qual_threshold = 27;
my $poly_T_len = 10;
my $seq_len_thresh = 30;

filter_fastq();

sub filter_fastq {
	open (FQI, $fastq_input) or die "Cannot open fastq: $fastq_input, $!\n";
	my %fastq;
	my ($fq_1, $fq_2, $fq_3, $fq_4, $rmt, $rmt_qual);
	my @rmt_qual;
	my $rmt_flag = 0;
	my ($read_count, $low_rmt_count, $short_read_count, $pcr_dup_count, $clean_count);
	while ($fq_1 = <FQI>) {
		$fq_2 = <FQI>;
		$fq_3 = <FQI>;
		$fq_4 = <FQI>;
		chomp $fq_2;
		chomp $fq_4;
		$read_count++;
		# reset vars
		@rmt_qual = ();
		$rmt_flag = 0;
		# filter low qual RMT
		($rmt, $rmt_qual) = split /\s/, $fq_1;
		@rmt_qual = split //, $rmt_qual;
		foreach my $q (@rmt_qual) {
			my $num = (ord $q) - $qual_offset;
			if ($num < $qual_threshold) {
				$rmt_flag = 1;
			}
		}
		if ($rmt_flag == 1) {
			$low_rmt_count++;
			next;
		}
		# filter polyT
		if ($fq_2 =~ /(T{$poly_T_len})/) {
			$fq_2 = $`;
			if (length $fq_2 < $seq_len_thresh) {
				$short_read_count++;
				next;
			}
			$fq_4 = substr($fq_4, 0, length($fq_2));
		}
		# filter PCR duplication
		if (exists $fastq{$rmt}{$fq_2}) {
			$pcr_dup_count++;
			my ($low1, $sum1) = qual_info($fq_4);
			my ($low2, $sum2) = qual_info($fastq{$rmt}{$fq_2});
			if ($low1 > $low2) {
				$fastq{$rmt}{$fq_2} = $fq_4;
			} else {
				if ($low1 == $low2 && $sum1 > $sum2) {
					$fastq{$rmt}{$fq_2} = $fq_4;
				}
			}
		} else {
			$clean_count++;
			$fastq{$rmt}{$fq_2} = $fq_4;
		}
	}
	close FQI;

	my @n = split /\//, $fastq_input;
	my $cell = $n[-2] . "\t" . $n[-1];
	$cell =~ s/\.fq//;
	print "$cell\t$read_count\t$clean_count\t$low_rmt_count\t$short_read_count\t$pcr_dup_count\n";

	# write fastq
	open (FQO, ">$fastq_output") or die "Cannot open fastq: $fastq_output, $!\n";
	foreach my $rmt (sort keys %fastq) {
		foreach my $read (sort keys $fastq{$rmt}) {
			print FQO "$rmt\n$read\n+\n", $fastq{$rmt}{$read}, "\n";
		}
	}
	close FQO;
}

# quality information
# input: quality string
# output: lowest, sum of quality
sub qual_info {
	my $qual = shift;
	my @w = split //, $qual;
	my $sum = 0;
	my $low = 100;
	my $num = 0;
	for (my $i = 0; $i < @w; $i++) {
		$num = (ord $w[$i]) - $qual_offset;
		$sum += $num;
		if ($num < $low) {
			$low = $num;
		}
	}
	return ($low, $sum);
}

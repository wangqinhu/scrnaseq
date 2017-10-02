#!/usr/bin/env perl
#===============================================================================
#
# Description: Convert the raw fastq of sc-rnaseq to clean fastq
#
# Copyright (c) 2017 Northwest A&F University
# Author: Qinhu Wang
# Email: wangqinhu@nwafu.edu.cn
#
#===============================================================================

use strict;
use warnings;
use Data::Dumper;

my $plate_conf = $ARGV[0] || "data/conf/sc-rnaseq.plate.conf";
my $fastq_conf = $ARGV[1] || "data/conf/sc-rnaseq.fastq.conf";
my $clean_dir = $ARGV[2] || "data/cell";

# experiment design
# position is 0-indexed
my $pool_barcode_start = 3;
my $pool_barcode_len = 4;
my $cell_barcode_start = 0;
my $cell_barcode_len = 7;
my $rmt_start = 7;
my $rmt_len = 8;
my $clean_start = 7;
my $clean_len = 68;

sc_plate2well();

sub sc_plate2well {
	my $well = load_plate_conf($plate_conf);
	my $seqs = load_fastq_conf($fastq_conf);
	foreach my $run (sort keys %{$seqs}) {
		print STDERR "Processing run: $run ...\n";
		my $fq1 = $seqs->{$run}->{'R1'};
		my $fq2 = $seqs->{$run}->{'R2'};
		open (FQ1, $fq1) or die "Cannot open fastq file 1: $fq1, $!\n";
		open (FQ2, $fq2) or die "Cannot open fastq file 2: $fq2, $!\n";
		my ($fq1_1, $fq1_2, $fq1_3, $fq1_4);
		my ($fq2_1, $fq2_2, $fq2_3, $fq2_4);
		while ($fq1_1 = <FQ1>) {
			$fq1_2 = <FQ1>;
			$fq1_3 = <FQ1>;
			$fq1_4 = <FQ1>;
			$fq2_1 = <FQ2>;
			$fq2_2 = <FQ2>;
			$fq2_3 = <FQ2>;
			$fq2_4 = <FQ2>;
			my $pool_barcode = substr($fq1_2, $pool_barcode_start, $pool_barcode_len);
			my $cell_barcode = substr($fq2_2, $cell_barcode_start, $cell_barcode_len);
			my $rmt = substr($fq2_2, $rmt_start, $rmt_len);
			next if (!exists $well->{$run}->{$pool_barcode});
			next if (!exists $well->{$run}->{$pool_barcode}->{$cell_barcode});
			my $sam = $well->{$run}->{$pool_barcode}->{$cell_barcode}->{'sam'};
			my $cell = $well->{$run}->{$pool_barcode}->{$cell_barcode}->{'cell'};
			if (!-d "$clean_dir/$sam") {
				system("mkdir -p $clean_dir/$sam");
			}
			open (WELL, ">>$clean_dir/$sam/$cell.fq") or die "Cannot open $clean_dir/$sam/$well.fq: $!\n";
			my $clean_seq = substr($fq1_2, $clean_start, $clean_len);
			my $clean_qual = substr($fq1_4, $clean_start, $clean_len);
			my $rmt_qual = substr($fq2_4, $rmt_start, $rmt_len);
			print WELL "@" . "$rmt $rmt_qual\n$clean_seq\n+\n$clean_qual\n";
			close WELL;
		}
		close FQ1;
		close FQ2;
	}
	print STDERR "Done\n";
}

sub load_plate_conf {
	my $palte_conf = shift;
	print STDERR "Parsing plate configuration $palte_conf ...\n";
	my %well = ();
	my ($well, $run, $sam, $pool_barcode, $cell_barcode);
	open (PC, $palte_conf) or die "Cannot open plate_conf file: $palte_conf, $!\n";
	while (<PC>) {
		chomp;
		next if (/^\#/);
		next if (/^\s*$/);
		my @w = split /\t/;
		my ($cell, $sam, $run, $pool_barcode, $cell_barcode) = ($w[0], $w[3], $w[4], $w[6], $w[7]);
		$pool_barcode = uc($pool_barcode);
		$cell_barcode = uc($cell_barcode);
		$well{$run}{$pool_barcode}{$cell_barcode}{"sam"} = $sam;
		$well{$run}{$pool_barcode}{$cell_barcode}{"cell"} = $cell;
	}
	close PC;
	return \%well;
}

sub load_fastq_conf {
	my $fastq_conf = shift;
	print STDERR "Parsing fastq configuration $fastq_conf ...\n";
	my %seqs = ();
	open (FC, $fastq_conf) or die "Cannot open fastq_conf file: $fastq_conf, $!\n";
	my ($file, $run);
	while (<FC>) {
		chomp;
		next if (/^\#/);
		next if (/^\s*$/);
		($file, $run) = split /\t/;
		if ($file =~ /R1\.f/) {
			$seqs{$run}{'R1'} = $file;
		} elsif (($file =~ /R2\.f/)) {
			$seqs{$run}{'R2'} = $file;
		} else {
			die "Unknown fastq file: $file in fastq_conf file $fastq_conf\n";
		}
	}
	close FC;
	return \%seqs;
}

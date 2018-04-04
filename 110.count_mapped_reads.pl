#!/usr/bin/env perl
#===============================================================================
#
# Description: 
#
# Copyright (c) 2017 Northwest A&F University
# Author: Qinhu Wang
# Email: wangqinhu@nwafu.edu.cn
#
#===============================================================================

use strict;
use warnings;
use Data::Dumper;

my $dir = shift;
opendir (DIR, $dir) or die $!;
foreach my $file (sort readdir DIR) {
	next if ($file =~ /^\./);
	open (IN, "$dir/$file") or die $!;
	my ($cell, $time) = (undef, 0);
	my $line1 = <IN>;
	my $line2 = <IN>;
	my $line3 = <IN>;
	my $line4 = <IN>;
	my $line5 = <IN>;
	my $line6 = <IN>;
	my $line7 = <IN>;
	if ($line1 =~ /\S\/(\S+)/) {
		$cell = $1;
	}
	if ($line5 =~ /^\s+(\d+)\s/) {
		$time += $1;
	}
	if ($line6 =~/^\s+(\d+)\s/) {
		$time += $1;
	}
	print $cell, "\t", $time, "\n";
	close IN;
}
close DIR;

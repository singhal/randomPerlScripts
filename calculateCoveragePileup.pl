use warnings;
use strict;

my $pileup = shift;

my %coverage; my $c = 'NA';
open(IN, "<$pileup");
while(<IN>) {
	chomp(my $line = $_);
	my @d = split(/\t/,$line);
	if ($d[0] eq $c) {
		$coverage{$d[0]}{'cov'} += length($d[5]);
		$coverage{$d[0]}{'base'}++;
		}
	else {	
		my $cov = $coverage{$c}{'cov'} / $coverage{$c}{'base'} if %coverage;
		print $c, "\t", $cov, "\n";
		%coverage = ();
		$c = $d[0];
		}
	}
	

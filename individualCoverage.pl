use warnings;
use strict;

my %l;
open(IN, "</Users/singhal/Desktop/coverage/library");
while(<IN>) {
	chomp(my $line = $_);
	my @d = split(/\t/,$line);
	push(@{$l{$d[2]}},$d[0]);
	}
close(IN);	
	
foreach my $contact (keys %l) {
	my $out = '/Users/singhal/Desktop/coverage/' . $contact . ".coverage.out";
	open(OUT, ">$out");
	my %cov;
	foreach my $lib (@{$l{$contact}}) {
		my $file = '/Users/singhal/Desktop/coverage/' . $lib . '.coverage.out';
		open(IN, "<$file");
		while(<IN>) {
			chomp(my $line = $_);
			my @d = split(/\t/,$line);
			$cov{$d[0]}{$lib} = $d[1];
			}
		close(IN);	
		}
	foreach my $gene (keys %cov) {
		print OUT $gene, "\t";
		foreach my $lib (sort {$a cmp $b} @{$l{$contact}}) {
			if ($cov{$gene}{$lib}) {
				print OUT $cov{$gene}{$lib}, "\t";
				}
			else {
				print OUT "NA\t";
				}
			}
		print OUT "\n";
		}
	close(OUT);	
	}	
use warnings;
use strict;

my %compare = (
'CarliaN_CarliaS' => {'seq1' => 'Carlia_N', 'seq2' => 'Carlia_S' },
'LamproN_LamproS' => {'seq1' => 'Lampro_N', 'seq2' => 'Lampro_S' },
'LamproC_LamproS' => {'seq1' => 'Lampro_C', 'seq2' => 'Lampro_S' },
'SaproC_SaproS' => {'seq1' => 'Sapro_C', 'seq2' => 'Sapro_S' },
);

my %l;
open(IN, "</Users/singhal/Desktop/coverage/library");
while(<IN>) {
	chomp(my $line = $_);
	my @d = split(/\t/,$line);
	push(@{$l{$d[2]}},$d[0]);
	}
close(IN);	
	

foreach my $compare (keys %compare) {
	my $out = '/Users/singhal/Desktop/coverage/' . $compare . ".coverage.out";
	open(OUT, ">$out");
	my $c1 = $compare{$compare}{'seq1'};
	my $c2 = $compare{$compare}{'seq2'};
	
	my $seq1 = '/Users/singhal/Desktop/coverage/' . $c1 . '_trinity.fa.final.annotated';
	my $seq2 = '/Users/singhal/Desktop/coverage/' . $c2 . '_trinity.fa.final.annotated';
	
	my %seq1; my %seq2;
	open(IN,"<$seq1");
	while(<IN>) {
		chomp(my $line = $_);
		if ($line =~ m/>(\S+).*(ENSAC\S+)/) {
			my $contig = $1; my $gene = $2;
			$seq1{$gene} = $contig;
			}
		}
	close(IN);	
	open(IN,"<$seq2");
	while(<IN>) {
		chomp(my $line = $_);
		if ($line =~ m/>(\S+).*(ENSAC\S+)/) {
			my $contig = $1; my $gene = $2;
			$seq2{$gene} = $contig;
			}
		}
	close(IN);	
	
	my %cov1;
	foreach my $lib (@{$l{$c1}}) {
		my $file = '/Users/singhal/Desktop/coverage/samples/' . $lib . '.coverage.out';
		open(IN, "<$file");
		while(<IN>) {
			chomp(my $line = $_);
			my @d = split(/\t/,$line);
			$cov1{$d[0]}{$lib} = $d[1];
			}
		close(IN);	
		}
	my %cov2;
	foreach my $lib (@{$l{$c2}}) {
		my $file = '/Users/singhal/Desktop/coverage/samples/' . $lib . '.coverage.out';
		open(IN, "<$file");
		while(<IN>) {
			chomp(my $line = $_);
			my @d = split(/\t/,$line);
			$cov2{$d[0]}{$lib} = $d[1];
			}
		close(IN);	
		}	
				
	foreach my $gene (keys %seq1) {
		print OUT $gene, "\t";
		foreach my $lib (sort {$a cmp $b} @{$l{$c1}}) {
			if ($cov1{$seq1{$gene}}{$lib}) {
				print OUT $cov1{$seq1{$gene}}{$lib}, "\t";
				}
			else {
				print OUT "NA\t";
				}
			}
		foreach my $lib (sort {$a cmp $b} @{$l{$c2}}) {
			if ($seq2{$gene}) {	
				if ($cov2{$seq2{$gene}}{$lib}) {
					print OUT $cov2{$seq2{$gene}}{$lib}, "\t";
					}
				else {
					print OUT "NA\t";
					}
				}
			else {
				print OUT "NA\t";
				}
			}	
		print OUT "\n";
		}
	close(OUT);	
	}	
use warnings;
use strict;

my $dir = '/Users/singhal/Desktop/genomics/';
my $covN = $dir . 'selfAlign/Carlia_N.coverage.out';
my $seqN = $dir . 'seqfiles/Carlia_N_trinity.fa.final.annotated';
my $covS = $dir . 'selfAlign/Carlia_S.coverage.out';
my $seqS = $dir . 'seqfiles/Carlia_S_trinity.fa.final.annotated';

my $cN = cov($covN);
my $cS = cov($covS);
my $sN = parseSeq($seqN);
my $sS = parseSeq($seqS);
my %covN = %$cN; my %covS = %$cS; my %seqN = %$sN; my %seqS = %$sS;

foreach my $gene (keys %seqN) {
	if ($seqS{$gene}) {
		my $cov1 = $covN{$seqN{$gene}{'contig'}};
		my $cov2 = $covS{$seqS{$gene}{'contig'}};
		if ($cov1 && $cov2 >= 250) {
			print $gene, "\t", $cov1, "\t", $cov2, "\n";
			}
		}
	}	

sub cov {
	my ($cov) = @_;
	my %cov;
	open(IN, "<$cov");
	foreach(<IN>) {
		chomp(my $line = $_);
		my @d = split(/\t/,$line);
		$cov{'>' . $d[0]} = $d[1];
		}
	close(IN);
	return(\%cov);
	}
	
sub parseSeq {
	my ($seq) = @_;
	
	my %seq;
	
	open(IN, "<$seq");
	while(<IN>) {
		chomp(my $line = $_);
		if ($line =~ m/>/) {
			my @d = split(/\t/,$line);
			chomp(my $seq = <IN>);
			if (scalar(@d) > 1) {
				if ($seq{$d[2]}) {
					$seq{$d[2]} = {'seq' => $seq, 'info' => $d[1], 'contig' => $d[0]} if length($seq) > length($seq{$d[2]}{'seq'});
					}
				else {
					$seq{$d[2]} = {'seq' => $seq, 'info' => $d[1], 'contig' => $d[0]};
					}
				}	
			}
		}	
	
	close(IN);	
	return(\%seq);	
	}	
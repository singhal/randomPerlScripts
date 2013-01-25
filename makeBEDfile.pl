use warnings;
use strict;

my @contacts = (["Lampro_N","Lampro_C"],["Sapro_C","Sapro_S"],["Carlia_N","Carlia_S"], ["Lampro_C", "Lampro_S"]);
my $dir = '/media/DataDrive/sutureGenomics/';
my @cov = qw(1 10 20 50);

foreach my $c (@contacts) {
    foreach my $cov (@cov) {
	my @contact = @{$c};
	my $seqfile = $dir . 'seqfiles/' . $contact[0] . '_trinity.fa.final.annotated';
	my $name = $contact[0] . '_' . $contact[1];
	my $out = $dir . 'snpTesting/' . $name . '.' . $cov . '.bed';

	my $seq = parseSeq($seqfile);
	$seq = coverageCheck($seq, $dir . 'dadi/try1/' . $name . '/' . $contact[0] . '.depth.out',$cov);
	$seq = coverageCheck($seq, $dir . 'dadi/try1/' . $name . '/' . $contact[1] . '.depth.out',$cov);
	$seq = deleteUncovered($seq);
	printBED($seq,$out);
    }
}

	
##########################
# behold the subroutines #
##########################
	
sub printBED {
    my ($var,$out) = @_;
    my %var = %{$var};
    
    open(OUT, ">$out");

    foreach my $c (sort {$a cmp $b} keys %var) {
	my @pos = sort {$a <=> $b} keys %{$var{$c}};
	my $start = $pos[0] - 1;
	my $end = $pos[$#pos] + 2;
	unless ($end == 2) { 
	    print OUT $c, "\t", $start, "\t", $end, "\n";
	}
    }
    close(OUT);
}

sub deleteUncovered {
	my ($var) = @_;
	foreach my $c (keys %{$var}) {
		foreach my $pos (keys %{$var->{$c}}) {
			if ($var->{$c}->{$pos}->{'cov'}) {
				unless ($var->{$c}->{$pos}->{'cov'} == 3) {
					delete($var->{$c}->{$pos});
					}
				}
			else {
				delete($var->{$c}->{$pos});
				}
			}	
		}	
	return($var);	
	}	
	
sub coverageCheck {
	my ($var,$depth,$cov) = @_;
	
	my %var = %{$var};
	open(IN, "<$depth");
	while(<IN>) {
		chomp(my $line = $_);
		my @d = split(/\t/,$line);
		
		if ($var{$d[0]}{$d[1]}) {
			my $covCheck = 0;
			$var{$d[0]}{$d[1]}{'cov'}++;	
			for (my $i = 2; $i < 7; $i++) {
				$covCheck++ if $d[$i] >= $cov;
				}
			unless ($covCheck == 5 ) {
				delete($var{$d[0]}{$d[1]});
				}
			}	
		}
			
	return(\%var);
	}

sub parseSeq {
	my ($s) = @_;
	my %s; 
	open(IN, "<$s");
	while(<IN>) {
		chomp(my $line = $_);
		if ($line =~ m/>(\S+).*ENS.*/) {
			my $c = $1; my @d = split(/\t/,$line);

			chomp(my $seq = <IN>);
			my $lseq = length($seq);
			
			my ($u5,$u3,$ge,$gs);
			
			if ($d[1] =~ m/5u(\d+)/) {
				$u5 = $1;
				}
			if ($d[1] =~ m/gs(\d+)/) {
				$gs = $1;
				}
			if ($d[1] =~ m/ge(\d+)/) {
				$ge = $1;
				}
			if ($d[1] =~ m/3u(\d+)/) {
				$u3 = $1;
				}	
				
			if ($u3) {	
				for (my $i = $u3; $i <= $lseq; $i++) {
					$s{$c}{$i}{'cov'}++;
					}
				}			
			}
		}
	close(IN);	
	return(\%s);	
	}	
	

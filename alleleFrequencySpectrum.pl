use warnings;
use strict;

#need to control for the possibility (wrong though it might be, of a three way hit)

my $dir = '/Users/singhal/Desktop/genomics/snpDiscovery/';
my $pileup1 = $dir . 'CarliaS.mpileup.out';
my $numInd = 5;
my $minInd = 5;
my @cov = (5,10,20,30,40,50,50);
my @af = (2,2,3,3,3,3,5);

for (my $i = 0; $i < scalar(@cov); $i++) {
	my $cov = $cov[$i];
	my $af = $af[$i];
	my $out = $dir . 'alleleFreqSpectrum_cov' . $cov . '_af' . $af . '.out';
	my $p1 = parsePileup($pileup1,$cov,$af,$out);
	}
	
sub parsePileup {
	my ($file,$cov,$af,$out) = @_;
	
	open(OUT, ">$out");
	my %snp;
	open(IN, "<$file");
	while(<IN>) {
		chomp(my $line = $_);
		my @d = split(/\t/,$line);
		my $c = $d[0]; my $pos = $d[1]; my $ref = $d[2];
		#checking to see what the coverage is at each individual
		my $pass = 0; my @reads;
		for (my $i = 0; $i < $numInd; $i++) {
			$pass++ if $d[3 + 3*$i] >= $cov;
			my $read = $d[4 + 3*$i];
			while ($read =~ /([\+|-](\d+))/g) {
				my $n = $2;
				$read =~ s/[\+|\-]\d+[atgcn]{$n}//i;
				}
			push(@reads, $read);
			}
		#have high enough coverage at the number of required individuals	
		if ($pass >= $minInd) {	
			my $aref = parseReads($ref,\@reads,$af);
			#array of arrays - $numInd individuals and two alleles at each
			my @a = @{$aref}; 
			
			my %allele;
			for (my $i = 0; $i < scalar(@a); $i++) {
				$allele{$a[$i][0]}++;
				$allele{$a[$i][1]}++;
				}
			my @alleles = sort {$allele{$a} <=> $allele{$b}} keys %allele;
			my $maf = $allele{$alleles[0]}/($numInd*2);
			print OUT $maf, "\n" if $maf < 1;
			}
		}
	close(OUT);
	close(IN);
	}	
	
sub parseReads {
	my ($ref,$r,$af) = @_;
	my @reads = @{$r}; my @a;
	for (my $i = 0; $i < scalar(@reads); $i++) {
		my %a;
		while ($reads[$i] =~ /([atgcATGC])/g) {
			$a{uc($1)}++;
			}
		foreach my $a (keys %a) {
			if ($a{$a} >= $af) {
				my @ref = $reads[$i] =~ m/([\.|\,])/g;
				if (scalar(@ref) >= $af) {
					#looks to be polymorphic
					$a[$i][0] = $ref; $a[$i][1] = $a;
					}
				else {
					#looks to be fixed
					$a[$i][0] = $a; $a[$i][1] = $a;
					}
				}
			else {
				$a[$i][0] = $ref; $a[$i][1] = $ref;
				}
			}
		unless ($a[$i][0]) {
			$a[$i][0] = $ref; $a[$i][1] = $ref;
			}
		}
	return(\@a);	
	}
	

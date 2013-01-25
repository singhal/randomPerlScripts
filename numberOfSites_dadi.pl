use warnings;
use strict;

#assumes throughout that you are working with 5 individuals per lineage and just 2 lineages
my $dir = '/media/DataDrive/sutureGenomics/lineageSNPs/lamproGillies/';
my $aligndir = $dir . 'align/';
my $final = $dir . 'LamproGillies.seq.out';
my $outgroup = $dir . 'lamprogilliesOutgroup.mpileup.out';
my $depth1 = $dir . 'Lampro_C.depth.out';
my $depth2 = $dir . 'Lampro_S.depth.out';
my $seq = $dir . 'Lampro_C_trinity.fa.final.annotated';
my $cov = '20';
my $qual = '20';
my $numSnp = '2';
my $numWin = '1';
my @seqFiles = </media/DataDrive/sutureGenomics/lineageSNPs/lamproGillies/seqFiles/*annotated>;

########################
# run the subroutines! #
########################

my %var;
#will want to check coverage at that position; at least $cov for every individual
my $var = checkCoverage($depth1,\%var);
print "Checked coverage (pass 1) of all the variants!\n";
$var = checkCoverage($depth2,$var);
print "Checked coverage (pass 2) of all the variants!\n";
$var = checkUTR($seq,$var);
print "Checked position of all the variants!\n";
$var = cleanUp($var);
print "Cleaned up the variants without coverage!\n";
#will need to polarize
$var = polarizePile($outgroup,$var);
print "Polarized all the variants (pass 1)!\n";
#polarizing with blast will recover a decent number of SNPs, but SUPER slow
$var = polarizeBlast($seq,\@seqFiles,$var);
print "Polarized all the variants (pass 2)!\n";
$var = cleanUp2($var);
print "Cleaned up the variants without polarization!\n";

open(FINAL, ">$final");
foreach my $c (keys %{$var}) {
	my $num = scalar(keys %{$var->{$c}});
	print FINAL $c, "\t", $num, "\n";
	}
close(FINAL);	

###########################
# behold the subroutines! #
###########################


sub cleanUp2 {
	my ($var) = @_;
	
	my %var = %{$var};
	
	foreach my $c (keys %var) {
		foreach my $pos (keys %{$var{$c}}) {
			unless($var{$c}{$pos}{'polar'}) {
				delete($var{$c}{$pos});
				}
			}
		}
	
	foreach my $c (keys %var) {
		my @keys = keys %{$var{$c}};
		unless(scalar(@keys) > 0) {
			delete($var{$c});
			}
		}
	return(\%var);	
	}	

sub cleanUp {
	my ($var) = @_;
	
	my %var = %{$var};
	
	foreach my $c (keys %var) {
		foreach my $pos (keys %{$var{$c}}) {
			if($var{$c}{$pos}{'cov'}) {
				unless ($var{$c}{$pos}{'cov'} == 2) {
					delete($var{$c}{$pos});
					}
				}
			}
		}
	
	foreach my $c (keys %var) {
		my @keys = keys %{$var{$c}};
		unless(scalar(@keys) > 0) {
			delete($var{$c});
			}
		}
	return(\%var);	
	}	
	
sub checkCoverage {
	my ($depth,$var) = @_;
	
	my %var = %{$var};
	open(IN, "<$depth");
	while(<IN>) {
		chomp(my $line = $_);
		my @d = split(/\t/,$line);
		
		my $covCheck = 0;	
		for (my $i = 2; $i < 7; $i++) {
			$covCheck++ if $d[$i] >= $cov;
			}
		if ($covCheck == 5 ) {
			$var{$d[0]}{$d[1]}{'cov'}++;
			}	
		}
			
	return(\%var);
	}
	
	
sub checkUTR {
	my ($refseq,$snp) = @_;
		
	my $seq = parseSeq($refseq);
	my %seq = %{$seq};	
	my %snp = %{$snp};
	
	foreach my $c (keys %snp) {
		foreach my $pos (keys %{$snp{$c}}) {
			#is this in cds
			my $loc;
			if ($pos <= $seq{$c}{'ge'} && $pos <= $seq{$c}{'ge'}) {
				$loc = 'cds';
				}
			elsif ($pos > $seq{$c}{'ge'}) {
				if ($seq{$c}{'3u'}) {
					$loc = '3u';
					}
				else {
					$loc = 'undef';
					}
				}	
			elsif ($pos < $seq{$c}{'gs'}) {
				if ($seq{$c}{'5u'}) {
					$loc = '5u';
					}
				else {
					$loc = 'undef';
					}
				}	
			if ($loc =~ m/\du/) {
				$snp{$c}{$pos}{'loc'} = $loc;
				}
			else {
				delete($snp{$c}{$pos});
				}
			}
		}	
		
	return(\%snp);
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
			$s{$c}{'seq'} = $seq;
			$s{$c}{'match'} = $d[2];
			
			if ($d[1] =~ m/5u(\d+)/) {
				$s{$c}{'5u'} = $1;
				}
			if ($d[1] =~ m/gs(\d+)/) {
				$s{$c}{'gs'} = $1;
				}
			if ($d[1] =~ m/ge(\d+)/) {
				$s{$c}{'ge'} = $1;
				}
			if ($d[1] =~ m/3u(\d+)/) {
				$s{$c}{'3u'} = $1;
				}	
			}
		}
	close(IN);	
	return(\%s);	
	}		

sub polarizePile {
	my ($outgroup,$var) = @_;
	
	my %var = %{$var};
	open(IN, "<$outgroup");
	while(<IN>) {
		chomp(my $line = $_);
		my @d = split(/\t/,$line);
		
		if ($var{$d[0]}{$d[1]}) {
			my $snp = $d[4] . $d[7] . $d[10] . $d[13] . $d[16] . $d[19] . $d[22] . $d[25] . $d[28] . $d[31];  
			my %snp;
			
			while ($snp =~ m/(\+\d+|\-\d+)/g) {
				my $match = $1;
				my $num = $1 if $match =~ m/(\d+)/;
				$match = '\\' . $match;	
				$snp =~ s/$match[atgcn]{$num}//i;
				}
			$snp =~ s/\*//g;	

			if (length($snp) > $numSnp) {
				
				while ($snp =~ m/([a|t|c|g])/gi) {
					$snp{uc($1)}++;
					}
			
				if (scalar(keys %snp) >  0 ) {
					my @snps = sort { $snp{$b} <=> $snp{$a} } keys %snp;
					my $winner = $snps[0];
				
					if ($snp{$winner} > $numWin) {
						$var{$d[0]}{$d[1]}{'polar'}++;
						}
					}	
				}	
			}
		}	
	close(IN);	
	return(\%var);	
	}

sub polarizeBlast {
	my ($seq,$seqFiles,$var) = @_;
	my %var = %{$var};
	
	my %seq;
	open(IN, "<$seq");
	while(<IN>) {
		chomp(my $line = $_);
		if ($line =~ m/>(\S+)/) {
			my $id = $1;
			chomp(my $seq = <IN>);
			$seq{$id} = $seq;
			}
		}	
	close(IN);	
	
	my %otherseq;
	my @seqFiles2;
	foreach my $seqFile (@$seqFiles) {
		open(SEQ, "<$seqFile");
		my $seqout = $seqFile . ".blast";
		push(@seqFiles2,$seqout);
		open(SEQOUT, ">$seqout");
		while(<SEQ>) {
			chomp(my $line = $_);
			if ($line =~ m/>(\S+)/) {
				my $id = $1;
				chomp(my $seq = <SEQ>);
				$otherseq{$seqout}{$id} = $seq;
				print SEQOUT ">", $id, "\n$seq\n";
				}
			}	
		close(SEQ); close(SEQOUT);
		my $call = system("formatdb -i $seqout -p F") unless (-f $seqout . ".nin");
		}
		
	#need to do this by variant
	foreach my $c (keys %var) {	
		foreach my $pos (keys %{$var{$c}}) {
			unless ($var{$c}{$pos}{'polar'}) {
				my $aln = $aligndir . $c . ".fa";
				my $aln_out = $aln . ".aln";
				
				unless(-f $aln) {
					open(OUT, ">test.fa");
					print OUT ">" . $c . "\n" . $seq{$c} . "\n";
					close(OUT);
					open(ALN, ">$aln");
					print ALN ">" . $c . "\n" . $seq{$c} . "\n";
					#for each variant, find the best matching contig in the other lineages (watch orientation!)
		
					foreach my $seqFile (@seqFiles2) {
						my @call = `blastall -p blastn -d $seqFile -i test.fa -m 8 -a 4 -b 	1`;
						if (@call) {
							my @d = split(/\t/,$call[0]);
							if ($d[8] < $d[9] && $d[6] < $d[7] || $d[8] > $d[9] && $d[6] > $d[7]) 	{
								print ALN ">" . $d[1] . "\n" . $otherseq{$seqFile}{$d[1]} . "\n";	
								}
							else {
								my $other =  $otherseq{$seqFile}{$d[1]};
								$other = reverse($other);
								$other =~ tr/ATGCatgc/TACGtacg/;
								print ALN ">" . $d[1] . "\n" . $other . "\n";
								}
							}
						}
					close(ALN);
					#align sequences
					my $call = system("muscle -in $aln -out $aln_out");
					}
					
				#parse alignment to find base in other lineages; be careful to watch change in position
				my %aln; my $id; my %pos_ref;
				open(ALN, "<$aln_out");
				while(<ALN>) {
					chomp(my $line = $_);
					if ($line =~ m/>(\S+)/) {
						$id = $1;
						}
					else {
						$aln{$id} .= $line;
						}
					}
				close(ALN);	
				foreach my $id (keys %aln) {
					my @seq = split(//,$aln{$id});
					$aln{$id} = \@seq;
					}
			
				my @ref = @{$aln{$c}}; my $counter = 1;
				for (my $i = 0; $i < scalar(@ref); $i++) {
					unless ($ref[$i] eq '-') {
						$pos_ref{$counter} = $i;
						$counter++;
						}
					}
					
				#allowing for alignment not being found, define the reference 
				my %snp; my $count;
				foreach my $otherc (keys %aln) {
					$snp{$aln{$otherc}[$pos_ref{$pos}]}++ if $aln{$otherc}[$pos_ref{$pos}] =~ m/[ATGC]/i;
					$count++ if $aln{$otherc}[$pos_ref{$pos}] =~ m/[ATGC]/i;
					}
			
				if ($count >  $numSnp ) {
					my @snps = sort { $snp{$b} <=> $snp{$a} } keys %snp;
					my $winner = $snps[0];
					
					if ($snp{$winner} > $numWin) {
						$var{$c}{$pos}{'polar'}++;
					    }
					}
				}	
			}
		}		
	return(\%var);	
	}	

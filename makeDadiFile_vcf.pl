use warnings;
use strict;

#assumes throughout that you are working with 5 individuals per lineage and just 2 lineages

my $dir = '/media/DataDrive/sutureGenomics/lineageSNPs/lamproGillies/';
my $aligndir = $dir . 'align/';
my $final = $dir . 'LamproGillies.dadi.out';
my $vcf = $dir . 'jointLamprogillies.vcf.out';
my $outgroup = $dir . 'lamprogilliesOutgroup.mpileup.out';
my $depth1 = $dir . 'Lampro_C.depth.out';
my $depth2 = $dir . 'Lampro_S.depth.out';
my $seq = $dir . 'Lampro_C_trinity.fa.final.annotated';
my $cov = '20';
my $qual = '20';
my $numSnp = '2';
my $numWin = '1';
my @seqFiles = </media/DataDrive/sutureGenomics/lineageSNPs/lamproGillies/seqFiles/*annotated>;

### need to do some benchmarking, because mpileup / depth does not necessarily put out data for every single base

########################
# run the subroutines! #
########################

mkdir($aligndir) unless(-d $aligndir);
#first need to get all the variant sites read in!
my $var = readVariant($vcf);
print "Read in all the variants!\n";
#will want to make sure it is in the UTR
$var = checkUTR($seq,$var);
print "Checked position of all the variants!\n";
#will want to check coverage at that position; at least $cov for every individual
$var = checkCoverage($depth1,$var);
print "Checked coverage (pass 1) of all the variants!\n";
$var = checkCoverage($depth2,$var);
print "Checked coverage (pass 2) of all the variants!\n";
#because of the nature of the depth file, we need to make sure every position has been checked for coverage twice, and delete if not
$var = deleteUncovered($var);
print "Removed variants without coverage!\n";
#will need to polarize
$var = polarizePile($outgroup,$var);
print "Polarized all the variants (pass 1)!\n";
#polarizing with blast will recover a decent number of SNPs, but SUPER slow
$var = polarizeBlast($seq,\@seqFiles,$var);
print "Polarized all the variants (pass 2)!\n";
#ditto with respect to polarization
$var = deleteUnpolarized($var);
print "Removed variants without polarization!\n";
#will need to put it in the right output for dadi
makeDadi($seq,$var);

###########################
# behold the subroutines! #
###########################

sub makeDadi {
	my ($seq,$var) = @_;
	
	my %seq;
	open(IN, "<$seq");
	while(<IN>) {
		chomp(my $line = $_);
		if ($line =~ m/>(\S+)/) {
			my $id = $1;
			chomp(my $seq = <IN>);
			my @seq = split(//,$seq);
			$seq{$id} = \@seq;
			}
		}	
	close(IN);
	
	open(OUT, ">$final");
	print OUT "Carlia\tOutgroup\tAllele1\tNpop\tSpop\tAllele2\tNpop\tSpop\tGene\tPosition\n";
	my %var = %{$var};
	foreach my $c (keys %var) {
		foreach my $pos (keys %{$var{$c}}) {
			my $bp1 = $seq{$c}[$pos-2];
			my $bp3 = $seq{$c}[$pos];
			my $tri = $bp1 . $var{$c}{$pos}{'ref'} . $bp3;
			my $a1_1 = 10 - $var{$c}{$pos}{'lin1'};
			my $a1_2 = 10 - $var{$c}{$pos}{'lin2'};
			print OUT $tri, "\t", $tri, "\t", $var{$c}{$pos}{'ref'}, "\t", $a1_1, "\t", $a1_2, "\t", $var{$c}{$pos}{'alt'}, "\t",  $var{$c}{$pos}{'lin1'}, "\t",  $var{$c}{$pos}{'lin2'}, "\t", $c, "\t", $pos, "\n";
			}
		}
	close(OUT);	
	}
	
			
sub deleteUncovered {
	my ($var) = @_;
	foreach my $c (keys %{$var}) {
		foreach my $pos (keys %{$var->{$c}}) {
			unless ($var->{$c}->{$pos}->{'lin1'}) {
				$var->{$c}->{$pos}->{'lin1'} = 0;
				}
			unless ($var->{$c}->{$pos}->{'lin2'}) {
				$var->{$c}->{$pos}->{'lin2'} = 0;
				}	
			if ($var->{$c}->{$pos}->{'cov'}) {
				unless ($var->{$c}->{$pos}->{'cov'} == 2) {
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
	
sub deleteUnpolarized {
	my ($var) = @_;
	foreach my $c (keys %{$var}) {
		foreach my $pos (keys %{$var->{$c}}) {
			unless ($var->{$c}->{$pos}->{'polar'}) {
				delete($var->{$c}->{$pos});
				}
			}	
		}	
	return($var);	
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
					    if ($winner eq $var{$c}{$pos}{'ref'}) {
						#don't need to do anything; everything is w.r.t. to derived	already
						$var{$c}{$pos}{'polar'}++;
					    }
					    elsif ($winner eq $var{$c}{$pos}{'alt'}) {
						#hey hey! have it all backwards; need to fix.
						$var{$c}{$pos}{'polar'}++;
						$var{$c}{$pos}{'alt'} = $var{$c}{$pos}{'ref'};
						$var{$c}{$pos}{'ref'} = $winner;
						$var{$c}{$pos}{'lin1'} = 10 - $var{$c}{$pos}{'lin1'};
						$var{$c}{$pos}{'lin2'} = 10 - $var{$c}{$pos}{'lin2'};
						}
					}
					}
				}	
			}
		}		
	return(\%var);	
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
				my $ref = 0;
				$ref++ while ($snp =~ m/[\.|\,]/g);
				$snp{$var{$d[0]}{$d[1]}{'ref'}} = $ref if $ref;
			
				while ($snp =~ m/([a|t|c|g])/gi) {
					$snp{uc($1)}++;
					}
			
				if (scalar(keys %snp) >  0 ) {
					my @snps = sort { $snp{$b} <=> $snp{$a} } keys %snp;
					my $winner = $snps[0];
				
					if ($snp{$winner} > $numWin) {
					if ($winner eq $var{$d[0]}{$d[1]}{'ref'}) {
						#don't need to do anything; everything is w.r.t. to derived already
						$var{$d[0]}{$d[1]}{'polar'}++;
						}
					elsif ($winner eq $var{$d[0]}{$d[1]}{'alt'}) {
						#hey hey! have it all backwards; need to fix.
						$var{$d[0]}{$d[1]}{'polar'}++;
						$var{$d[0]}{$d[1]}{'alt'} = $var{$d[0]}{$d[1]}{'ref'};
						$var{$d[0]}{$d[1]}{'ref'} = $winner;
						$var{$d[0]}{$d[1]}{'lin1'} = 10 - $var{$d[0]}{$d[1]}{'lin1'};
						$var{$d[0]}{$d[1]}{'lin2'} = 10 - $var{$d[0]}{$d[1]}{'lin2'};
						}	
					}
					}	
				}	
			}
		}	
	close(IN);	
	return(\%var);	
	}

sub checkCoverage {
	my ($depth,$var) = @_;
	
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

sub readVariant {
	my ($vcf) = @_;
	
	my %var;
	open(IN, "<$vcf");
	while(<IN>) {
		chomp(my $line = $_);
		if ($line =~ m/contig/) {
			my @d = split(/\t/,$line);
			
			unless ($d[4] =~ m/\,/) {
				$var{$d[0]}{$d[1]}{'ref'} = $d[3];
				$var{$d[0]}{$d[1]}{'qual'} = $d[5];
				$var{$d[0]}{$d[1]}{'alt'}= $d[4];
				for (my $i = 9; $i < 14; $i++) {
					$var{$d[0]}{$d[1]}{'lin1'}++ if $d[$i] =~ m/1\//;
					$var{$d[0]}{$d[1]}{'lin1'}++ if $d[$i] =~ m/\/1/;
					}	
				for (my $i = 14; $i < 19; $i++) {
					$var{$d[0]}{$d[1]}{'lin2'}++ if $d[$i] =~ m/1\//;
					$var{$d[0]}{$d[1]}{'lin2'}++ if $d[$i] =~ m/\/1/;
					}	
				}	
			}
		}	
	close(IN);	
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
	

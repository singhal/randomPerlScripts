#####################################################################
# this script takes a mpileup output files and uses it to find snps,#
#	both defining their type, location, and coding type             #
# assumptions: annotated sequence format is in                      #
#	8annotateTranscriptomes.pl format                               #
# written by sonal.singhal1 [at] gmail.com, 4 March 2012            #
#####################################################################

#need to control for the possibility (wrong though it might be, of a three way hit)

use warnings;
use strict;

my $dir = '/Users/singhal/Desktop/genomics/snpDiscovery/';
my $pileup1 = $dir . 'CarliaS.mpileup.out';
my $pileup2 = $dir . 'CarliaN_toS.mpileup.out';
my $refseq = '/Users/singhal/Desktop/genomics/seqfiles/oldAnnotations/Carlia_S_trinity.fa.final.annotated';
my $numInd = 5;
my $minInd = 4;
#my @cov = (5,10,20,30,40,50,50);
#my @af = (2,2,3,3,3,3,5);
my @cov = (50,50);
my @af = (3,5);
my $tmpdir = $dir . 'tmp/';

#mkdir($tmpdir) unless (-f $tmpdir);
#splitPileup($pileup1,$tmpdir,1);	
#splitPileup($pileup2,$tmpdir,2);
for (my $i = 0; $i < scalar(@cov); $i++) {
	my $cov = $cov[$i];
	my $af = $af[$i];	
	my $out = $dir . 'annotatedSNPs_cov' .  $cov . '_af' . $af . '.out';
	open(OUT, ">$out"); 
	my $seq = parseSeq($refseq);
	my %seq = %$seq;
	foreach my $c (keys %seq) {
		my $tmp1 = $tmpdir . $c . '_1.out';
		my $tmp2 = $tmpdir . $c . '_2.out';
		if (-f $tmp1 && -f $tmp2) {
			my $p1 = parsePileup($tmp1,$cov);
			my $p2 = parsePileup($tmp2,$cov);
		
			my $snp = comparePileup($p1,$p2,$af);
			$snp = annotateSNP($seq,$snp);
	
			my %snp = %$snp;
			foreach my $c (sort {$a cmp $b} keys %snp) {
				foreach my $pos (sort {$a <=> $b} keys %{$snp{$c}}) {
					print OUT $c, "\t", $pos, "\t", $snp{$c}{$pos}{'a1'}, "\t", $snp{$c}{$pos}{'a2'}, "\t", $snp{$c}{$pos}{'type'}, "\t", $snp{$c}{$pos}{'loc'}, "\t", $snp{$c}{$pos}{'coding'}, "\n";
					}
				}
			}	
		}	
	close(OUT);	
	}
	
sub splitPileup {
	my ($pile, $tmpdir,$id) = @_;
	
	my $c = 'NA';
	open(IN, "<$pile");
	while(<IN>) {
		chomp(my $line = $_);
		if ($line =~ m/(contig\S+)/) {
			if ($1 eq $c) {
				print OUT $line, "\n";
				}
			else {
				#new file
				if ($c eq 'NA') {
					$c = $1;
					my $file = $tmpdir . $c . "_" . $id . '.out';
					open(OUT, ">$file");
					print OUT $line, "\n";
					}
				else {
					$c = $1;
					close(OUT);
					my $file = $tmpdir . $c . "_" . $id . '.out';
					open(OUT, ">$file");
					print OUT $line, "\n";
					}
				}
			}
		}
	close(IN);
	}	
	
sub annotateSNP {
	my ($seq,$snp) = @_;
	
	#define location of SNP
	my %snp = %$snp; 
	my %seq = %$seq;
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
			
			my $coding;
			if ($loc eq 'cds') {
				my $start = substr $seq{$c}{'seq'}, 0, $pos - 1;
				my $end = substr $seq{$c}{'seq'}, $pos;
				#need to introduce my mutation1
				my $orf1 = $start . $snp{$c}{$pos}{'a1'} . $end;
				#need to introduce my mutation2
				my $orf2 = $start . $snp{$c}{$pos}{'a2'} . $end;
				
				#this could lead to a SNP change -- let's see!
				my $gs = $seq{$c}{'gs'} - 1;
				my $ln = $seq{$c}{'ge'} - $seq{$c}{'gs'} + 1;
				$orf1 = substr $orf1, $gs, $ln;
				$orf2 = substr $orf2, $gs, $ln;
				my $aa1 = translate($orf1);
				my $aa2 = translate($orf2);
				if ($aa1 eq $aa2) {
					$coding = 'syn';
					}
				else {
					$coding = 'ns';
					}		
				}
			else {
				$coding = 'noncoding';
				}
			$snp{$c}{$pos}{'loc'} = $loc;
			$snp{$c}{$pos}{'coding'} = $coding;
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

sub comparePileup {
	my ($p1,$p2,$af) = @_;
	my %p1 = %$p1; my %p2 = %$p2;
	
	my %snp;
	
	foreach my $c (sort {$a cmp $b} keys %p1) {
		foreach my $pos (sort {$a <=> $b} keys %{$p1{$c}}) {
			if ($p2{$c}{$pos}) {
				foreach my $ref (keys %{$p1{$c}{$pos}}) {		
				 	my $a1 = parseReads($ref,$p1{$c}{$pos}{$ref},$af);
				 	my $a2 = parseReads($ref,$p2{$c}{$pos}{$ref},$af);
				 	
				 	my @a1 = @{$a1}; my @a2 = @{$a2};
				 	my %a1; my %a2;
				 	for (my $i = 0; $i < scalar(@a1); $i++) {
				 		$a1{$a1[$i][0]}++; $a1{$a1[$i][1]}++;	
						}
					for (my $i = 0; $i < scalar(@a2); $i++) {
						$a2{$a2[$i][0]}++; $a2{$a2[$i][1]}++;
						}	
				 		 	
				 	if (scalar(keys %a1) == 1 && scalar(keys %a2) == 1) {
				 		#are they the same base?
				 		my @b1 = keys %a1; my @b2 = keys %a2;
				 		#if not, no SNP
				 		if ($b1[0] ne $b2[0]) {
				 			$snp{$c}{$pos} = {'type' => 'fixed', 'a1'=>$b1[0],'a2'=>$b2[0]};
				 			}
				 		}
				 	elsif (scalar(keys %a1) > 1 && scalar(keys %a2) == 1) {	
				 		my @b1 = sort {$a1{$b} <=> $a1{$a}} keys %a1;
				 		$snp{$c}{$pos} = {'type'=> 'poly1', 'a1'=>$b1[0],'a2'=>$b1[1]};
				 		}
				 	elsif (scalar(keys %a1) == 1 && scalar(keys %a2) > 1) {	
				 		my @b2 = sort {$a2{$b} <=> $a2{$a}} keys %a2;
				 		$snp{$c}{$pos} = {'type'=> 'poly2', 'a1'=>$b2[0],'a2'=>$b2[1]};
				 		}	
				 	else {
				 		my $shared;
				 		foreach (keys %a1) {
				 			$shared++ if $a2{$_};
				 			}
				 		if ($shared > 1) {
				 			my @b1 = sort {$a1{$b} <=> $a1{$a}} keys %a1;
				 			$snp{$c}{$pos} = {'type'=>'shared','a1'=>$b1[0],'a2'=>$b1[1]};
				 			}
				 		}
					}
				}
			}
		}
	return(\%snp);	
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
	
sub parsePileup {
	my ($file,$cov) = @_;
	
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
			$snp{$c}{$pos}{$ref} = \@reads;
			}
		}
	
	close(IN);

	return(\%snp);
	}
	
sub translate {
	my $string = shift;
	$string = uc($string);
	my @codons = $string =~ m/(\S\S\S)/g;
	my %codons = (	'ATG'=>'M','ACG'=>'T','CTG'=>'L','CCG'=>'P','GTG'=>'V','GCG'=>'A','TTG'=>'L','TCG'=>'S',
					'ATA'=>'I','ACA'=>'T','CTA'=>'L','CCA'=>'P','GTA'=>'V','GCA'=>'A','TTA'=>'L','TCA'=>'S',
					'ATC'=>'I','ACC'=>'T','CTC'=>'L','CCC'=>'P','GTC'=>'V','GCC'=>'A','TTC'=>'F','TCC'=>'S',
					'ATT'=>'I','ACT'=>'T','CTT'=>'L','CCT'=>'P','GTT'=>'V','GCT'=>'A','TTT'=>'F','TCT'=>'S',
					'AGG'=>'R','AAG'=>'K','CGG'=>'R','CAG'=>'Q','GGG'=>'G','GAG'=>'E','TGG'=>'W','TAG'=>'*',
					'AGA'=>'R','AAA'=>'K','CGA'=>'R','CAA'=>'Q','GGA'=>'G','GAA'=>'E','TGA'=>'*','TAA'=>'*',
					'AGC'=>'S','AAC'=>'N','CGC'=>'R','CAC'=>'H','GGC'=>'G','GAC'=>'D','TGC'=>'C','TAC'=>'Y',
					'AGT'=>'S','AAT'=>'N','CGT'=>'R','CAT'=>'H','GGT'=>'G','GAT'=>'D','TGT'=>'C','TAT'=>'Y');
	my $translate;
	foreach(@codons) {
		if ($codons{$_}) {
			$translate = $translate . $codons{$_};
			}
		else {
			$translate = $translate . 'X';
			}
		}
	return($translate);
	}	
	
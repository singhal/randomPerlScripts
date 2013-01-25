use warnings;
use strict;

my $seq = '/Users/singhal/Desktop/genomics/seqFiles/Carlia_N_trinity.fa.final.annotated';
my $vcfFile = '/Users/singhal/Desktop/genomics/carliaN.all.vcf';

my $seqInfo = seq($seq);
my $vcf = parseVcf($vcfFile);
defineOrf($seqInfo,$vcf);

sub parseVcf {
	my ($file) = @_;
	
	my %vcf;
	open(IN, "<$file");
	while(<IN>) {
		chomp(my $line = $_);
		if ($line =~ m/^contig/) {	
			my @d = split(/\t/,$line);
			#vcf uses coding of base 1
			my $snp = 'snp';
			$snp = 'indel' if $d[7] =~ m/INDEL/;
			my @snp = split(/,/,$d[4]);
			my $af = $1 if $d[7] =~ m/AF1=([0-9|\.]+)/;
			$vcf{$d[0]}{$snp}{$d[1] - 1} = {'var' => \@snp, 'afreq'=> $af }; 	
			}
		}
	close(IN);
	
	return(\%vcf);
	}	
			
			
sub seq {
	my ($file) = @_;
	my %seq;
	open(IN, "<$file");
	while(<IN>) {
		chomp(my $line = $_);
		if ($line =~ m/>(\S+)/) {
			my $id = $1;
			chomp(my $seq = <IN>);
			my @d = split(/\t/,$line);
			$seq{$id} = {'seq' => $seq, 'info' => $d[1], 'gene' => $d[2]};
			}
		}
	close(IN);
	return(\%seq);
	}

sub defineOrf {
	my ($seq, $vcf) = @_;
	
	my %vcf = %$vcf;
	
	my %seq = %$seq;
	foreach my $c (keys %seq) {
		if ($seq{$c}{'info'}) {
			#$vcf{$d[0]}{$snp}{$d[1] - 1} = {'snp' => \@snp, 'afreq'=> $af }; 
			#do mutations -- start with snps
			my $newseq  = $seq{$c}{'seq'};
			foreach my $pos (keys %{$vcf{$c}{'snp'}}) {
				if ($vcf{$c}{'snp'}{$pos}{'afreq'} > 0.49) {
					my $seq1 = substr $seq{$c}{'seq'}, 0, $pos;
					my $seq2 = substr $seq{$c}{'seq'}, $pos + 1;
					$newseq = $seq1 . $vcf{$c}{'snp'}{$pos}{'var'}[0] . $seq2;
					}
				}	
			
			#then do indels. as you do indels, realize that the frame is shifting
			my $info = $seq{$c}{'info'};
			my $gs = $1 if $info =~ m/gs(\d+)/;
			my $ge = $1 if $info =~ m/ge(\d+)/;
			my $gs0 = $gs - 1;
			my $length = $ge - $gs + 1;
		
			my $s1 = $seq{$c}{'seq'};
			my $orf1 = substr $s1, $gs0, $length;
			my $orf2 = substr $newseq, $gs0, $length;
		
			my $aa1 = translate($orf1);
			my $aa2 = translate($orf2);
		
			if ($aa1 ne $aa2) {
				if ($aa1 =~ m/\*/ && $aa2 !~ m/\*/) {
					print ">", $c, "\n", $aa1, "\n", $aa2, "\n";
					}
				}
			}	
		}
	return(\%seq);	
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
			print "ERROR: ILLEGAL PASS TO CODON TRANSLATION: $_ is not a codon!\n";
			$translate = $translate . 'X';
			}
		}
	return($translate);
	}	
	
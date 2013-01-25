use warnings;
use strict;

my $main_dir = '/Users/singhal/Desktop/eugongylus/finalAssemblies/';
my $carlia_n = $main_dir . 'Carlia_N_trinity.fa.final.annotated';
my $carlia_s = $main_dir . 'Carlia_S_trinity.fa.final.annotated';
my $lampro_n = $main_dir . 'Lampro_N_trinity.fa.final.annotated';
my $lampro_c = $main_dir . 'Lampro_C_trinity.fa.final.annotated';
my $lampro_s = $main_dir . 'Lampro_S_trinity.fa.final.annotated';
my $sapro_c = $main_dir . 'Sapro_C_trinity.fa.final.annotated';
my $sapro_s = $main_dir . 'Sapro_S_trinity.fa.final.annotated';
my $out = $main_dir . 'geneticDivergences.out';
my $np = 8;
my $evalue = 1e-20;

my %compare = (
'CarliaN_CarliaS' => {'seq1' => $carlia_n, 'seq2' => $carlia_s },
'LamproN_LamproC' => {'seq1' => $lampro_n, 'seq2' => $lampro_c },
'LamproC_LamproS' => {'seq1' => $lampro_c, 'seq2' => $lampro_s },
'SaproC_SaproS' => {'seq1' => $sapro_c, 'seq2' => $sapro_s },
'CarliaN_LamproN' => {'seq1' => $carlia_n, 'seq2' => $lampro_n },
'CarliaN_SaproC' => {'seq1' => $carlia_n, 'seq2' => $sapro_c },
'SaproC_LamproN' => {'seq1' => $sapro_c, 'seq2' => $lampro_n }
);

open(FINAL, ">$out");
foreach my $compare (keys %compare) {
	my $seq1 = $compare{$compare}{'seq1'};
	my $seq2 = $compare{$compare}{'seq2'};
	
	my $s1b = parseSeq2($seq1);
	my $s2b = parseSeq2($seq2);
	recipBlastProt($s1b,$s2b,$seq2,$compare,$out);
	}
close(FINAL);
	
sub recipBlastProt {
	my ($s1,$s2,$seq2,$compare, $out) = @_;
	
	my %ti = ('A' => 'G', 'G' => 'A', 'C' => 'T', 'T' => 'C');
	
	my %l;
	my %seq1 = %$s1; my %seq2 = %$s2;
	my $blast = "recipBlastProt.out";
	open(OUT, ">annotatedRef.fa");
	foreach my $c (keys %seq1) {
		if ($seq1{$c}{'info'}) {
			my $gs1 = $1 if $seq1{$c}{'info'} =~ m/gs(\d+)/;
			my $ge1 = $1 if $seq1{$c}{'info'} =~ m/ge(\d+)/;
			my $length = $ge1 - $gs1 + 1;
			my $seq1 = substr $seq1{$c}{'seq'}, $gs1 - 1, $length;
			$seq1{$c}{'cds'} = $seq1;
			print OUT ">$c\n", $seq1, "\n";
			}
		}
	my $newseq = newdb($seq2);	
	my $call1 = system("formatdb -i $newseq -p F");
	my $call2 = system("blastall -p blastn -d $newseq -i annotatedRef.fa -a $np -e $evalue -m 8 -o $blast -b 1");
	close(OUT);
	my $call3 = system("rm $newseq*");		
	my $call = system("rm annotatedRef.fa");	
	
	my %matches;
	open(IN, "<$blast");
	while(<IN>) {
		chomp(my $line = $_);
		my @d = split(/\t/,$line);
		$matches{$d[0]} = $d[1];
		}
	close(IN);
	
	open(OUT, ">>$out");
	
	foreach my $c1 (keys %matches) {
		open(TAR, ">target.fa");
		open(QUER, ">query.fa");
		print TAR ">$c1\n$seq1{$c1}{'cds'}\n";
		print QUER ">$matches{$c1}\n$seq2{$matches{$c1}}{'seq'}\n";
		close(TAR); close(QUER);		
		my @call = `blastn -query query.fa -subject target.fa -outfmt 0`;
		
		if (@call) {
			my $gc; my $ti = 0; my $tv = 0; my $l = -1;
			for (my $i = 0; $i < scalar(@call); $i++) {
				if ($call[$i] =~ m/^Query\s+\d+/) {
					my $line1 = $call[$i];
					my $line2 = $call[$i+2];
					$i = $i + 1;
				
					my @d1 = $line1 =~ m/(\S+)/g;
					my @d2 = $line2 =~ m/(\S+)/g;
					my @b1 = split(//,$d1[2]);
					my @b2 = split(//,$d2[2]);
					
					for (my $j = 0; $j < scalar(@b1); $j++) {
						if (uc($b1[$j]) =~ m/[A|T|G|C]/i) {
							if (uc($b1[$j]) =~ m/[G|C]/i) {
								$gc++;
								}	
							if (uc($b2[$j]) =~ m/[A|T|G|C]/) {	
								$l++;
								unless (uc($b1[$j]) eq uc($b2[$j])) {
									#this is a mutation!
									if ( $ti{uc($b1[$j])} eq uc($b2[$j]) ) {
										$ti++;
										}
									else {	
										$tv++;
										}
									}	
								}
							}
						}
					}		
				}	
			if ($l > 0) {	
				my $p = $ti/$l; my $q = $tv/$l; my $w = $gc/$l;
				my $div = (-2*$w) * ( 1 - $w) *  log(1 - ( $p/ (2 * $w * (1 - $w) ) ) - $q ) - 0.5 * ( 1 - 2 * $w * ( 1 - $w ) ) * log( 1 - 2 * $q );
				print OUT $compare, "\t", $c1, "\t", length($seq1{$c1}{'cds'}), "\t", $l, "\t", $div, "\n";	
				}
			}
		}	
	my $call4 = system("rm query.fa target.fa");
	close(OUT);
	unlink($blast);	
	}
	
sub parseSeq2 {
	my ($seq) = @_;
	
	my %seq;
	
	open(IN, "<$seq");
	while(<IN>) {
		chomp(my $line = $_);
		if ($line =~ m/>(\S+)/) {
			my $c = $1;
			my @d = split(/\t/,$line);
			chomp(my $seq = <IN>);
			if (scalar(@d) > 1) {
				$seq{$c} = {'seq' => $seq, 'info' => $d[1]};
				}	
			else {
				$seq{$c} = {'seq' => $seq};
				}
			}
		}	
	close(IN);	
	return(\%seq);	
	}	
	
sub	newdb {
	my ($s) = @_;
	
	open(IN, "<$s");
	my $newseq = 'newseq.fa';
	open(OUT, ">$newseq");
	while(<IN>) {
		my $line = $_;
		if ($line =~ m/(>\S+)/) {
			print OUT $1, "\n";
			}
		else {
			print OUT $line;
			}
		}	
	close(IN);
	close(OUT);
	return($newseq);
	}
		


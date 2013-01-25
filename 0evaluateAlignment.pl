use warnings;
use strict;

my $dir = '/Users/singhal/Desktop/genomics/alignStats/';
my $timing = $dir . 'alignerTiming.out';
my @contacts = qw(Carlia_N Lampro_C Lampro_N);
my @progs = qw(bowtie bowtie2 bwa novo smalt soap);

my $masterout = $dir . "alignEvaluation.out";

my %master;
foreach my $contact (@contacts) {
	foreach my $prog (@progs) {
		my $bam = $dir . $contact . "." . $prog . ".sorted.bam";
		my $file1 = $dir . $contact . "." . $prog . ".vcf";
		
		print "working on $prog analysis of $contact now!\n";
		
		my $call1 = system("samtools index $bam") unless(-f $bam . ".bai");
		my $pairs = paired($bam);
		my $reads = allReads($bam);
		my $time = timeAlign($timing,$contact,$prog);
		my $vcf1 = parseVcf($file1);
		
		$master{$contact . "_" . $prog}{'pairs'} = $pairs;
		$master{$contact . "_" . $prog}{'reads'} = $reads;
		$master{$contact . "_" . $prog}{'time'} = $time;
		
		foreach my $prog2 (@progs) {
			if ($prog2 ne $prog) {
				print "\tcomparing $prog to $prog2 now!\n";
				my @sort = sort {$a cmp $b} ($prog,$prog2);
				my $id = $sort[0] . '_' . $sort[1];
				unless ($master{$contact . "_" . $prog}{$id . 'uniqSnp'}) {
					my $bam2 = $dir . $contact . "." . $prog2 . ".sorted.bam";
					my $file2 = $dir . $contact . "." . $prog2 . ".vcf";
					my $vcf2 = parseVcf($file2);
					my $afCorr = afCompare($vcf1,$vcf2);
					my $uniqSnp = snpCompare($file1,$file2);
					my $uniqIndel = indelCompare($file1,$file2);
					my $covCorr = covCompare($bam,$bam2);					
					$master{$contact . "_" . $prog}{$id . 'uniqSnp'} = $uniqSnp;
					$master{$contact . "_" . $prog}{$id . 'uniqIndel'} = $uniqIndel;
					$master{$contact . "_" . $prog}{$id . 'afCorr'} = $afCorr;
					$master{$contact . "_" . $prog}{$id . 'covCorr'} = $covCorr;
					}
				}
			}
		}
	}
	
my %header;
foreach my $k1 (keys %master) {
	foreach my $k2 (keys %{$master{$k1}}) {
		$header{$k2}++;
		}
	}

open(MASTER, ">$masterout");
print MASTER "id\t";
foreach my $d (sort {$a cmp $b} keys %header) {	
	print MASTER $d, "\t";
	}
print MASTER "\n";	

foreach my $id (keys %master) {
	print MASTER $id, "\t";
	foreach my $d (sort {$a cmp $b} keys %header) {	
		if ($master{$id}{$d}) {
			print MASTER $master{$id}{$d}, "\t";
			}
		else {
			print MASTER "NA\t";
			}
		}
	print MASTER "\n";
	}
close(MASTER);	
	
#num of paired reads mapped
sub paired {
	my ($bam) = @_;
	my @call = `cat $bam | samtools view -f 0x0002 - | wc`;
	my $reads = $1 if $call[0] =~ m/^\s+(\d+)/;
	$reads = 'NA' unless($reads);
	return($reads);
	}

#num of reads total mapped
sub allReads {
	my ($bam) = @_;
	my @call = `samtools idxstats $bam`;
	my $totReads;
	foreach my $line (@call) {
		my @d = split(/\t/,$line);
		$totReads += $d[2];
		}
	return($totReads);
	}
	
#time it takes to run
sub timeAlign {
	my ($file,$contact,$prog) = @_;
	
	my %time;
	open(IN, "<$file");
	while(<IN>) {
		chomp(my $line = $_);
		my @d = split(/\t/,$line);
		$time{$d[0]}{$d[1]} = $d[2];
		}
	close(IN);	
	
	return($time{$prog}{$contact});
	}

#comparison of found snps
sub snpCompare {
	my ($vcf1,$vcf2) = @_;
	
	my $total;
	
	my $call1 = system("grep -v 'INDEL' $vcf1 | grep 'contig' | cut -f 1,2 > out1");
	my $call2 = system("grep -v 'INDEL' $vcf2 | grep 'contig' | cut -f 1,2 > out2");

	my @call1 = `wc out1`;
	my @call2 = `wc out2`;
	my @call3 = `cat out1 out2 | sort | uniq | wc`;
	my $tot1 = $1 if $call1[0] =~ m/^\s+(\d+)/;
	my $tot2 = $1 if $call2[0] =~ m/^\s+(\d+)/;
	my $uniq = $1 if $call3[0] =~ m/^\s+(\d+)/;
	my $tot = $tot1 + $tot2;
	my $unique = sprintf("%.3f", $uniq/$tot);
	my $call = system("rm out1 out2");
	
	return($unique);
	}

#compare indels
sub indelCompare {
	my ($vcf1,$vcf2) = @_;
	
	my $total;
	
	my $call1 = system("grep 'INDEL' $vcf1 | grep 'contig' | cut -f 1,2 > out1");
	my $call2 = system("grep 'INDEL' $vcf2 | grep 'contig' | cut -f 1,2 > out2");

	my @call1 = `wc out1`;
	my @call2 = `wc out2`;
	my @call3 = `cat out1 out2 | sort | uniq | wc`;
	my $tot1 = $1 if $call1[0] =~ m/^\s+(\d+)/;
	my $tot2 = $1 if $call2[0] =~ m/^\s+(\d+)/;
	my $uniq = $1 if $call3[0] =~ m/^\s+(\d+)/;
	my $tot = $tot1 + $tot2;
	my $unique;
	if ($tot) {
		$unique = sprintf("%.3f", $uniq/$tot);
		}
	else {
		$unique = 'NA';
		}
	my $call = system("rm out1 out2");	
		
	return($unique);
	}

#comparison of found snps allele frequencies
sub afCompare {
	my ($vcf1,$vcf2) = @_;

	my $tracker = 0;
	my $x = [];
	my %vcf1 = %$vcf1; my %vcf2 = %$vcf2;
	foreach my $c (keys %vcf1)	{
		foreach my $pos (keys %{$vcf1{$c}{'snp'}}) {
			if ($vcf2{$c}{'snp'}{$pos}) {
				$x -> [$tracker][1] = $vcf1{$c}{'snp'}{$pos}{'afreq'};
				$x -> [$tracker][2] = $vcf2{$c}{'snp'}{$pos}{'afreq'};
				$tracker++;
				}
			}
		}
		
	my $corr = correlation($x);
	return($corr);
	}
		
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

#comparison of coverage
sub covCompare {
	my ($bam1,$bam2) = @_;
	
	unless (-f $bam2 . ".bai") {
		my $call1 = system("samtools index $bam2");
		}
		
	my $cov1 = coverageDepth($bam1);	
	my $cov2 = coverageDepth($bam2);
	my %cov1 = %$cov1; my %cov2 = %$cov2;
	
	my $tracker = 0;
	my $x = [];
	
	foreach my $c1 (keys %cov1) {
		if ($cov2{$c1}) {
			if ($tracker < 10000) {
				$x->[$tracker][1] = $cov1{$c1}; 
				$x->[$tracker][2] = $cov2{$c1};
				}
			$tracker++;
			}
		}
			
	my $correlation = correlation($x);
	return($correlation);
	}
		
sub coverageDepth {
	my ($bam) = @_;
	
	my %final;
	
	my $out = $bam . "depth.out";
	my $call = system("samtools depth $bam > $out") unless (-f $out);
	
	my %coverage; my $c = 'NA';
	open(OUT,"<$out");
	while(<OUT>) {
		chomp(my $line = $_);
		my @d = split(/\t/,$line);
		if ($d[0] eq $c) {
			$coverage{$d[0]}{'cov'} += $d[2];
			$coverage{$d[0]}{'base'}++;
			}
		else {	
			if (%coverage) {
				my $cov = $coverage{$c}{'cov'} / $coverage{$c}{'base'}; 
				$final{$c} = $cov;
				%coverage = ();
				}
			$c = $d[0];
			}
		}
	#unlink($out);	
	return(\%final);	
	}
	
###################################################################################
# Pearson correlation code from                                                   #  
# http://davetang.org/muse/2010/11/29/calculating-pearson-correlation-using-perl/ #	
###################################################################################	

sub mean {
	my ($x)=@_;
	my $num = scalar(@{$x}) - 1;
	my $sum_x = '0';
	my $sum_y = '0';
	for (my $i = 0; $i < scalar(@{$x}); ++$i){
		$sum_x += $x->[$i][1];
		$sum_y += $x->[$i][2];
   		}
   	my $mu_x = $sum_x / $num;
   	my $mu_y = $sum_y / $num;
   	return($mu_x,$mu_y);
	}

sub ss {
	my ($x,$mean_x,$mean_y,$one,$two)=@_;
	my $sum = '0';
	for (my $i=0;$i<scalar(@{$x});++$i){
		$sum += ($x->[$i][$one]-$mean_x)*($x->[$i][$two]-$mean_y);
		}
	return $sum;
	}
 
sub correlation {
	my ($x) = @_;
	my ($mean_x,$mean_y) = mean($x);
	my $ssxx=ss($x,$mean_x,$mean_y,1,1);
	my $ssyy=ss($x,$mean_x,$mean_y,2,2);
	my $ssxy=ss($x,$mean_x,$mean_y,1,2);
	my $xcorrel;
	if ($ssxx == 0 || $ssyy == 0 || $ssxy == 0) {
		$xcorrel = 'NA';
		}
	else {	
		my $correl=correl($ssxx,$ssyy,$ssxy);
		$xcorrel=sprintf("%.4f",$correl);
		}
	return($xcorrel);
	}
 
sub correl {
	my($ssxx,$ssyy,$ssxy)=@_;
	my $sign=$ssxy/abs($ssxy);
	my $correl=$sign*sqrt($ssxy*$ssxy/($ssxx*$ssyy));
	return $correl;
	}		
		

		
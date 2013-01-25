####################################################################
# this script takes a sam output file and uses it to find snps,    #
#	both defining their type and location on the read to evaluate  #
#	error profiling                                                #
# assumptions: SAM output includes optional MD tag                 #
# written by sonal.singhal1 [at] gmail.com, 1 March 2012           #
####################################################################

use warnings;
use strict;

#hash of hashes of hashes: with reference first, then snp, then pos
my %type;

my $file = shift;

open(IN, "<$file");
while(<IN>) {
	chomp(my $line1 = $_);
	my @d1 = split(/\t/,$line1);
	if ($d1[0] =~ m/_F/) {
		chomp(my $line2 = <IN>);
		my @d2 = split(/\t/,$line2);
		#there are matches both here
		if ($d1[3] && $d2[3]) {
			if ($d1[3] < $d2[3]) {
				forwardSNP($line1,\@d1);
				reverseSNP($line2,\@d2);	
				}
			else {
				forwardSNP($line2,\@d2);
				reverseSNP($line1,\@d1);
				}
			}
		}
	}
close(IN);	
foreach my $ref (sort {$a cmp $b} keys %type) {
	foreach my $snp (sort {$a cmp $b} keys %{$type{$ref}}) {
		unless ($ref =~ m/N/ || $snp =~ m/N/) {
			foreach my $pos (sort {$a <=> $b} keys %{$type{$ref}{$snp}}) {
				print $ref, "_", $snp, "\t", $pos, "\t", $type{$ref}{$snp}{$pos}, "\n";
				}
			}	
		}
	}
			
sub forwardSNP {
	my ($line, $d) = @_;
	my @d = @{$d};

	if ($line =~ m/MD:Z:(\S+)/) {
		my $cigar = $d[5];
		my $seq = $d[9];
		my $error = $1;
		
		my $adj = 0;
		if ($cigar =~ m/^(\d+)S/) {
			$adj = $1;
			}

		#there is a deletion here
		if ($cigar =~ m/D/ && $cigar !~ m/I/) {
			while ($error =~ /(?=(\d+[ATCGN])\d+)/g) {		
				my $errorpos = $+[0] + length $1;
				my $suberr = substr $error, 0, $errorpos;
				my @dist = $suberr =~ m/(\d+)/g;
				my @snp = $suberr =~ m/(?=\d+([ATCGN])\d+)/g;
				my $pos;
				foreach (@dist) {
					$pos += $_;
					}
				$pos += scalar(@snp)  + 1 + $adj;
				my $relpos = sprintf("%.2f", $pos/length($seq));	
					
				my @seq = split(//,$seq);
				my $ref = $1 if $suberr =~ m/([ATCGN])$/;
				$type{$ref}{$seq[$pos - 1]}{$relpos}++;
				}
			}
		#there is an insertion & maybe deletion too....	
		elsif ($cigar =~ m/I/) {
			while ($error =~ /(?=(\d+[ATCGN])\d+)/g) {	
				my $errorpos = $+[0] + length $1;
				my $suberr = substr $error, 0, $errorpos;
				my @dist = $suberr =~ m/(\d+)/g;
				my @snp = $suberr =~ m/(?=\d+([ATCGN])\d+)/g;
				my $pos;
				foreach (@dist) {
					$pos += $_;
					}
				$pos += scalar(@snp)  + 1;
				
				#now need to correct the position for the insertion...
				my %shift; my $loc = 1; my $totInsert = 0;
				while ($cigar =~ /(?=(\d+)([M|I]))/g) {	
					my $type = $2; my $num = $1;
					if ($type eq 'M') {
						while ($num > 0) {	
							$shift{$loc} = $loc + $totInsert;
							$num--; $loc++;
							}
						}
					else {
						$totInsert += $num;
						}
					}
				$pos = $shift{$pos} + $adj;
				
				my $relpos = sprintf("%.2f", $pos/length($seq));	
					
				my @seq = split(//,$seq);
				my $ref = $1 if $suberr =~ m/([ATCGN])$/;
				$type{$ref}{$seq[$pos - 1]}{$relpos}++;
				}
			}
		#straight up alignment
		else {
			#this is a case where there is actually a sequence mismatch
			unless (length($seq) =~ m/$error/) {
				my @dist = ($error =~ m/(\d+)/g);
				my @snp = ($error =~ m/([ATCGN])/g);
				
				#this position is with start referenced to 1
				for (my $i = 0; $i < scalar(@snp); $i++) {
					my $j = $i; my $pos = 0;
					while ($j >=0) {
						$pos = $pos + $dist[$j] + 1;
						$j--;
						}
					$pos = $pos + $adj;
					my $relpos = sprintf("%.2f", $pos/length($seq));	
					
					my @seq = split(//,$seq);
					$type{$snp[$i]}{$seq[$pos - 1]}{$relpos}++;
					}
				}
			}
		}
	}

sub reverseSNP {
	my ($line, $d) = @_;
	my @d = @{$d};

	if ($line =~ m/MD:Z:(\S+)/) {
		my $cigar = $d[5];
		my $seq = $d[9];
		my $error = $1;

		my $adj = 0;
		if ($cigar =~ m/^(\d+)S/) {
			$adj = $1;
			}

		#there is a deletion here
		if ($cigar =~ m/D/ && $cigar !~ m/I/) {
			while ($error =~ /(?=(\d+[ATCGN])\d+)/g) {		
				my $errorpos = $+[0] + length $1;
				my $suberr = substr $error, 0, $errorpos;
				my @dist = $suberr =~ m/(\d+)/g;
				my @snp = $suberr =~ m/(?=\d+([ATCGN])\d+)/g;
				my $pos;
				foreach (@dist) {
					$pos += $_;
					}
				$pos += scalar(@snp)  + 1 + $adj;
				my $relpos = sprintf("%.2f", $pos/length($seq));	
				$relpos = 1 - $relpos;				
	
				my @seq = split(//,$seq);
				my $ref = $1 if $suberr =~ m/([ATCGN])$/;
				$type{$ref}{$seq[$pos - 1]}{$relpos}++;
				}
			}
		#there is an insertion & maybe deletion too....	
		elsif ($cigar =~ m/I/) {
			while ($error =~ /(?=(\d+[ATCGN])\d+)/g) {	
				my $errorpos = $+[0] + length $1;
				my $suberr = substr $error, 0, $errorpos;
				my @dist = $suberr =~ m/(\d+)/g;
				my @snp = $suberr =~ m/(?=\d+([ATCGN])\d+)/g;
				my $pos;
				foreach (@dist) {
					$pos += $_;
					}
				$pos += scalar(@snp)  + 1 ;
				
				#now need to correct the position for the insertion...
				my %shift; my $loc = 1; my $totInsert = 0;
				while ($cigar =~ /(?=(\d+)([M|I]))/g) {	
					my $type = $2; my $num = $1;
					if ($type eq 'M') {
						while ($num > 0) {	
							$shift{$loc} = $loc + $totInsert;
							$num--; $loc++;
							}
						}
					else {
						$totInsert += $num;
						}
					}
				$pos = $shift{$pos} + $adj;
				
				my $relpos = sprintf("%.2f", $pos/length($seq));
				$relpos = 1 - $relpos;		
					
				my @seq = split(//,$seq);
				my $ref = $1 if $suberr =~ m/([ATCGN])$/;
				$type{$ref}{$seq[$pos - 1]}{$relpos}++;
				}
			}
		#straight up alignment
		else {
			#this is a case where there is actually a sequence mismatch
			unless (length($seq) =~ m/$error/) {
				my @dist = ($error =~ m/(\d+)/g);
				my @snp = ($error =~ m/([ATCGN])/g);
				
				#this position is with start referenced to 1
				for (my $i = 0; $i < scalar(@snp); $i++) {
					my $j = $i; my $pos = 0;
					while ($j >=0) {
						$pos = $pos + $dist[$j] + 1;
						$j--;
						}
					$pos = $pos + $adj;
					my $relpos = sprintf("%.2f", $pos/length($seq));
					$relpos = 1 - $relpos;		
					
					my @seq = split(//,$seq);
					$type{$snp[$i]}{$seq[$pos - 1]}{$relpos}++;
					}
				}
			}
		}
	}


	
	


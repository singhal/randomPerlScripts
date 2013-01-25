use warnings;
use strict;

my $sdir = '/media/DataDrive/sutureGenomics/seqfiles/';
my $dir = '/media/DataDrive/sutureGenomics/snpTesting/';
my @contact = ('Lampro_C_Lampro_S');

foreach my $contact (@contact) {
	my @pileup = <$dir*$contact*pile*out>;
	my $seq = $sdir . $1 . '_annotated.fa' if $contact =~ m/^([a-z]+_[a-z])/i;
	my $other = $dir . $contact . '.seq.out';
	my $out = $dir . $contact . '.fa';
	open(OUT, ">$out");
	my $newdir = $dir . $contact . 'TMP/';
	mkdir($newdir) unless(-d $newdir);
	my $bed = $dir . $contact . '.bed';

	my %bed;
	open(IN, "<$bed");
	while(<IN>) {
	    chomp(my $line = $_);
	    my @d = split(/\t/,$line);
	    for (my $i = $d[1]; $i < $d[2]; $i++) {
		$bed{$d[0]}{$i + 1}++;
	    }
	}
	close(IN);
			
	open(IN, "<$other");
	my $c = 'NA';
	while(<IN>) {
		chomp(my $line = $_);
		if ($line =~ m/(contig\S+)/) {
			if ($1 eq $c) {
				print TMPOUT $line, "\n";
				}
			else {
				#new file
				if ($c eq 'NA') {
					$c = $1;
					my $file = $newdir . $c . '.anc.out';
					open(TMPOUT, ">$file");
					print TMPOUT $line, "\n";
					}
				else {
					$c = $1; 
					close(TMPOUT);
					my $file = $newdir . $c . '.anc.out';
					open(TMPOUT, ">$file");
					print TMPOUT $line, "\n";
					}
				}
			}
		}		
	close(IN);

	my %names; my %c;
	foreach my $pile (@pileup) {	
		open(IN, "<$pile");
		my $name = $1 if $pile =~ m/mpileup(\d)/i;		
		$names{$name}++;
		
		my $junk = <IN>;
		my $c = 'NA';
		while(<IN>) {
			chomp(my $line = $_);
			if ($line =~ m/(contig\S+)/) {
				if ($1 eq $c) {
					print TMPOUT $line, "\n";
					}
				else {
					#new file
					if ($c eq 'NA') {
						$c = $1; $c{$c}++;
						my $file = $newdir . $c . "_" . $name . '.out';
						open(TMPOUT, ">$file");
						print TMPOUT $line, "\n";
						}
					else {
						$c = $1; $c{$c}++;
						close(TMPOUT);
						my $file = $newdir . $c . "_" . $name . '.out';
						open(TMPOUT, ">$file");
						print TMPOUT $line, "\n";
						}
					}
				}
			}
		close(IN);
		}		

	my %master;
	open(IN, "<$seq");
	while(<IN>) {
	    chomp(my $line = $_);
	    if ($line =~ m/>(\S+)/) {
		my $c = $1;
		chomp(my $seq = <IN>);
		my @seq = split(//,$seq);
		$master{$c} = \@seq;
	    }
	}
	close(IN);

	foreach my $c (sort {$a cmp $b} keys %c) {
	    my %seq; my %anc; my %var;
		
		my $ancout = $newdir . $c . '.anc.out';
		if (-f $ancout) {
		open(IN, "<$ancout");
		while(<IN>) {
		    chomp(my $line = $_);
		    my @d = split(/\t/,$line);
			for (my $i = 2; $i <= 6; $i++) {
				$anc{$d[0]}{$d[1]}{$d[$i]}++ unless $d[$i] =~ m/N/;
				}
			}
		}
		
		foreach my $name (keys %names) {
			my $in = $newdir . $c . '_' . $name . '.out';
			if (-f $in) {
				open(IN, "<$in");	

				while(<IN>) {
					chomp(my $line = $_);
					my @d = split(/\t/,$line);
	    				
					my $snp = '.';
					for (my $i = 4; $i < scalar(@d); $i = $i + 3) {
					    $snp .= $d[$i];
						}
							
					while ($snp =~ m/(\+\d+|\-\d+)/g) {
						my $match = $1;
						my $num = $1 if $match =~ m/(\d+)/;
						$match = '\\' . $match; 
						$snp =~ s/$match[atgcn]{$num}//i;
						}
                             
					my $ref = 0; $ref++ while $snp =~ m/[\.|\,]/g;
					my %snp; my $alt = 0;
					while ($snp =~ m/([a|t|c|g])/gi) {
					$snp{uc($1)}++;
					$alt++;
					}
           	
					my $totcov = $alt + $ref;	
					my @snps = sort {$snp{$b} <=> $snp{$a}} keys %snp;	            	    
        
					my $af = 0;
	    		    $af = sprintf("%.3f",$snp{$snps[0]}/$totcov) if @snps;
	                
	     		   if ($totcov > 2) {
	       			 	if ($af > 0.5) {
    	    				$seq{$d[1]}{$name} = $snps[0];
					$var{$d[1]}++;
    	    				}
    	    			else {
    	    				$seq{$d[1]}{$name} = $d[2];
    	    				}
    	    			}
				}
			}
		}
        
		print OUT ">", $c, "\n";
		for (my $i = 0; $i < scalar(@{$master{$c}}); $i++) {
		    my $pos = $i + 1;
		    if ($var{$pos}) {
			if ($seq{$pos}{'1'} && $seq{$pos}{'2'}) {
				if ($seq{$pos}{'1'} eq $seq{$pos}{'2'}) {
				    my @snps = sort {$anc{$c}{$pos}{$b} <=> $anc{$c}{$pos}{$a}} (keys %{$anc{$c}{$pos}});
				    if (@snps) {
					if ($anc{$c}{$pos}{$snps[0]} > 2) {
					    print OUT $snps[0];
						}
					else {
					    print OUT $seq{$pos}{'1'};
						}
			    	}
			    else {
					print OUT $seq{$pos}{'1'};
			    }
			}
			else {
			    my @snps = sort {$anc{$c}{$pos}{$b} <=> $anc{$c}{$pos}{$a}} (keys %{$anc{$c}{$pos}});
			    if (@snps) {
					if ($anc{$c}{$pos}{$snps[0]} > 2) {
					    if ($seq{$pos}{'1'} eq $snps[0] || $seq{$pos}{'2'} eq $snps[0]) {
							print OUT $snps[0];
				    		}
				   		 else {
							if ($master{$c}[$i] eq $seq{$pos}{'1'} || $master{$c}[$i] eq $seq{$pos}{'2'}) {
					  		  #non ideal
					  		  print OUT $master{$c}[$i];
					  		  }
							else {
					  	 	 	#non ideal
					    		print OUT $master{$c}[$i];
								}
						 }
					}
				else {
				    #non ideal
				    print OUT $master{$c}[$i];
					}
			    }
			else {
				#non ideal
				print OUT $master{$c}[$i];
			    }
			}
			}
		#cannot confidently call seq here based on sequence
		else {
			#non ideal
			print OUT $master{$c}[$i];
		}
		    }
		    else {
			print OUT $master{$c}[$i];
		    }
		}
	print OUT "\n";
	}
 	close(OUT);	
 	my $call = system("rm -r $newdir");
}
	

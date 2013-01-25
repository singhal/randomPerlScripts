use warnings; 
use strict;

#defines the reassembled
my %a;
my @assemblies = </Users/singhal/Desktop/assemblies/*>;
foreach my $a (@assemblies) {
	my $c = $1 if $a =~ m/(contig\d+)/;
	open(IN, "<$a");
	while(<IN>) {
		if ($_ =~ m/>/)  {
			chomp(my $seq = <IN>);
			if ($a{$c}) {
				$a{$c} = $seq if length($seq) > length($a{$c});
				}
			else {
				$a{$c} = $seq;
				}
			}
		}
	close(IN);	
	}	
	
	
my %s;
my $s = '/Users/singhal/Desktop/genomics/seqfiles/Carlia_N_trinity.fa.final.annotated';
open(IN, "<$s");
while(<IN>) {
	chomp(my $line = $_);
	if ($line =~ m/>(\S+)/) {
		my $c = $1;
	#	my $info = $2;
		chomp(my $seq = <IN>);
		$s{$c}{'seq'} = $seq;
	##	$s{$c}{'info'} = $info;
		}
	}	
close(IN);

foreach my $c (keys %a) {
	my $diff = length($a{$c}) - length($s{$c}{'seq'});
	print $diff, "\n";
	}
use strict;
use warnings;

my $l;
my $data = shift;
open(IN, "<$data");
while(<IN>) {
	chomp(my $line = $_);
	if ($line =~ m/>/) {
		}
	else {
		$l += length($line);
		}
	}
close(IN);
print $l, "\n";
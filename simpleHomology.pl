use warnings;
use strict;

my @files = <*>;

my %id;
	
foreach my $file (@files) {
	open(IN, "<$file");
	while(<IN>) {
		my $line = $_;
		if ($line =~ m/(ENSAC\S+)/) {
			$id{$1}++;
			}
		}
	close(IN);
	}
		
foreach (keys %id) {
	print $_, "\n" if $id{$_} >= scalar(@files);
	}

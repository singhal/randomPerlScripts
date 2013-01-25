use warnings;
use strict;

my $file = shift;
open(IN, "<$file");

while(<IN>) {
    chomp(my $line = $_);
    if ($line =~ m/^\@(HS\S+)/) {
	my $id = $1;
	my $seq = <IN>;
	print ">", "$1", "\n", $seq;
    }
}
close(IN);

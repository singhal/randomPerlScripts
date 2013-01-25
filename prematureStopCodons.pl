use warnings;
use strict;

#check for premature stop codons
#my @assemblies = </Users/singhal/Desktop/genomics/seqfiles/*>;
my $s = shift;

#foreach my $s (@assemblies) {
my %seq;
open(IN, "<$s");
while(<IN>) {
    chomp(my $line = $_);
    if ($line =~ m/>(\S+)/) {
my $id = $1;
chomp(my $seq = <IN>);
my @d = split(/\t/,$line);
if ($d[1]) {
    $seq{$id} = {'seq' => $seq, 'info' => $d[1], 'gene' => $d[2]};
}
    }
}
close(IN);

extendOrf($s, \%seq);
#}

sub extendOrf {
    my ($name,$seq) = @_;

    my $count;
    my %seq = %$seq;
    foreach my $c (keys %seq) {
	if ($seq{$c}{'info'}) {
	    my $info = $seq{$c}{'info'};
	    my $gs = $1 if $info =~ m/gs(\d+)/;
	    my $ge = $1 if $info =~ m/ge(\d+)/;
	    my $gs0 = $gs - 1;
	    my $length = $ge - $gs + 1;
	    
	    my $s = $seq{$c}{'seq'};
	    
	    my $def_orf = substr $s, $gs0, $length;
	    my $def_aa = translate($def_orf);
	    
	    if ($def_aa) {
		my $new_aa;
		if ($def_aa =~ m/\*/) {
		    $new_aa = $1 if $def_aa =~ m/^([A-Z]+)\*/;
		}
		else {
		    $new_aa = $def_aa;
		}
		
		
		#if (length($def_aa) - length($new_aa) > 10 && (length($def_aa) - length($new_aa))/length($def_aa) > 0.1) {
		
		if (length($def_aa) ne length($new_aa)) {
		    if (length($def_aa) - length($new_aa) > 1) {
			print $c, "\t", length($new_aa)/length($def_aa), "\n";
			#print $def_aa, "\n";
		    }
		    #$count++;
		}
	    }
	}
    }
#my $per = $count/scalar(keys %seq);
#print $name, "\t", $count, "\t", $per, "\n";
}


sub translate {
    my $string = shift;
    $string = uc($string);
    my @codons = $string =~ m/(\S\S\S)/g;
    my %codons = ('ATG'=>'M','ACG'=>'T','CTG'=>'L','CCG'=>'P','GTG'=>'V','GCG'=>'A','TTG'=>'L','TCG'=>'S',
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
#	    print "ERROR: ILLEGAL PASS TO CODON TRANSLATION: $_ is not a codon!\n";
	    $translate = $translate . 'X';
	}
    }
    return($translate);
}

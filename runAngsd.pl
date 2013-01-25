use warnings;
use strict;

my $dir = '/media/DataDrive/sutureGenomics/';
my $results_dir =  $dir . 'snpTesting/';
my $results_dir2 = $results_dir . 'dadi/';
my $library = $dir . 'library';

my @contacts = (["Lampro_N","Lampro_C"],["Sapro_C","Sapro_S"],["Carlia_N","Carlia_S"],["Lampro_C","Lampro_S"]);
mkdir($results_dir2) unless(-d $results_dir2);

#defines the library
my %lib;
open(IN, "<$library");
while (<IN>) {
	chomp(my $line = $_);
	my @d = split(/\t/, $line);
	my $lib = $d[2] . '/' . $d[0] . '/';
	$lib{$d[2]}{$d[0]} = {'dir' => $lib, 'index' => $d[1]}; 
    } 
close(IN);


foreach my $contact (@contacts) {
	my @contact = @{$contact};
	my $name = $contact[0] . '_' . $contact[1];
	my $bed = $results_dir . $name . '.bed';
	my $seq = $results_dir . $name . '.fa';
	my $bamdir = $results_dir . $name . '/';

	my @posterior;

	my $all;
	foreach my $lineage (@contact) {
		
		my $bamfiles;
		foreach my $id (sort {$a cmp $b} keys %{$lib{$lineage}}) {
			$bamfiles .= $bamdir . $id . '.sorted.bam' . "\t";
			}
		$all .= $bamfiles;

		my $out = $results_dir2 . $lineage . '.' . $name;
#		my $call1 = system("samtools mpileup -d 20000 -A -E -I -l $bed -uf $seq $bamfiles | bcftools view -bvcg -> raw.bcf");
#		my $call2 = system("bcftools view raw.bcf > $out" . ".vcf");

#		my $call1 = system("/home/singhal/programs/angsd0.03/angsd.g++ -anc $seq -realSFS 1 -outfiles $out mpileup -f $seq -l $bed -d 20000 -A -I -g $bamfiles > NULL");
#		my $sfsout = $out . '.sfs';
#		my $call2 = system("/home/singhal/programs/angsd0.03/misc/optimSFS.gcc -nChr 10 -binput $sfsout");
		my $posterior = $out . '_posterior_probabilities.txt';
#		my $sfsml = $out . '.sfs.ml';
#		my $call3 = system("/home/singhal/programs/angsd0.03/misc/sfstools.g++ -nChr 10 -priorFile $sfsml -sfsFile $sfsout -dumpBinary 0 > $posterior");
		
		push(@posterior,$posterior);
		}

#	my $out = $results_dir2 . $name . '.joint.vcf';
#	my $call2 = system("samtools mpileup -d 20000 -A -E -I -l $bed -uf $seq $all | bcftools view -bvcg - > var_out.raw.bcf");
#	my $call3 = system("bcftools view var_out.raw.bcf > $out");
		
	my $sfs2d = $results_dir2 . $name . '.2DSFS';	
	my $call = system("perl /media/DataDrive/sutureGenomics/pipeline/random/2dsfs.pl -k 5 $posterior[0] $posterior[1] > $sfs2d");	
#	my $call2 = system("rm NULL");
	}	
	
			
	

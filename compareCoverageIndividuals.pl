use warnings;
use strict;

#home directory with all my files; make sure to end with a backslash
my $dir = '/media/DataDrive/sutureGenomics/';
#this file has all the info about my library, and creates a directory of all my files
my $lib = $dir . 'library';
my $np = 2;

#defines the library
open(IN, "<$lib");
my %lib;
while (<IN>) {
	chomp(my $line = $_);
	my @d = split(/\t/, $line);
	my $lib = $d[2] . '/' . $d[0] . '/';
	$lib{$d[2]}{$d[0]} = {'dir' => $lib, 'index' => $d[1]};	
	}
close(IN);

###########################
# run the subroutines     #
###########################	

foreach my $contact (keys %lib) {
    my $oldseq = $dir . 'oldAnnotations/' . $contact . '_trinity.fa.final.annotated';
    my $seq = $dir . 'oldAnnotations/' . $contact . '_annotated.fa';
    open(IN, "<$oldseq");
    open(OUT, ">$seq");
    while(<IN>) {
	my $line = $_;
	if ($line =~ m/>(\S+).*ENS.*/) {
	    my $c = $1;
	    my $fa = <IN>;
	    print OUT ">", $c, "\n", $fa;
	}
    }
    close(IN); close(OUT);
  #  my $call1 = system("bowtie2-build -q $seq $seq");
    foreach my $lib (sort {$a cmp $b} keys %{$lib{$contact}}) {
		my $subdir = $dir . $lib{$contact}{$lib}{'dir'};
	    print "doing $lib now!\n";
	    my $file1gz = $subdir . $lib . '_1p_final.fastq.gz';
	    my $file2gz = $subdir . $lib . '_2p_final.fastq.gz';
	    my $fileugz = $subdir . $lib . '_u_final.fastq.gz';

	    #my $file1 = unzip($file1gz);
	    #my $file2 = unzip($file2gz);
	    #my $fileu = unzip($fileugz);
		#runMapping($seq,$file1,$file2,$fileu,$subdir,$dir) unless (-f $subdir . $lib . '.sorted.bam') ;

		print "Coverage for $lib now!\n";
		runCoverage($seq,$subdir,$lib);

	 #   my @files = ($file1,$file2,$fileu);
	 #   zip(\@files);
		}
    }


sub runCoverage {
	my ($seq,$subdir,$lib) = @_;
	my $call = system("samtools faidx $seq");
	my $out = $subdir . $lib . "mpileup.out";
	my $call1 = system("samtools mpileup -A -f $seq $subdir" . $lib . ".sorted.bam > $out");
	}

sub runMapping {
	my ($seq,$file1,$file2,$fileu,$subdir,$dir) = @_;
	my $call1 = system("bowtie2 -x $seq -1 $file1 -2 $file2 -S bowtie1.sam -5 5 -3 5 --sensitive -k 10 -X 300 -p $np");
	my $call2 = system("bowtie2 -x $seq $fileu -S bowtie2.sam -5 5 -3 5 --sensitive -k 10  -p $np");
	my $call3 = system("samtools view -bS bowtie1.sam > bowtie1.bam");
	my $call4 = system("samtools view -bS bowtie2.sam > bowtie2.bam");
	my $call5 = system("samtools merge bowtie.bam bowtie1.bam bowtie2.bam");
	my $call6 = system("samtools sort bowtie.bam $subdir" . $lib . ".sorted");
	my $call7 = system("rm bowtie1.sam bowtie2.sam bowtie.bam bowtie1.bam bowtie2.bam");
	}

sub zip {
    my ($file) = @_;
    foreach my $_ (@$file) {
	my $call = system("gzip -1 $_");
    }
}

sub unzip {
    my ($file) = @_;
    my $unzip = $1 if $file =~ m/(\S+)\.gz/;
    my $call = system("gunzip $file");
    return($unzip);
}


		       

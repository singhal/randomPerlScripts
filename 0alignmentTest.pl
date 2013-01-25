use warnings;
use strict;

my @assemblies = </media/DataDrive/sutureGenomics/*/assemblies/*trinity*final>;
my $rootdir = '/media/DataDrive/sutureGenomics/';
my $np = 4;
my $time_out = $rootdir . "alignerTiming.out"; 
open(OUT, ">>$time_out");

foreach my $assembly (@assemblies) {
	my $contact = $1 if $assembly =~ m/$rootdir([a-z|_]+)/i;
	if ($contact =~ m/Lampro_N/) {
	my $outdir =  $rootdir . $contact . '/';
	
	my $reads1gz = $outdir . $contact . '_1.fastq.gz';
	my $reads2gz = $outdir . $contact . '_2.fastq.gz';
	my $readsugz = $outdir . $contact . '_u.fastq.gz';

	my $reads1 = unzipFiles($reads1gz);
	my $reads2 = unzipFiles($reads2gz);
	my $readsu = unzipFiles($readsugz);
	
	my $files = runAligners($assembly, $np, $reads1, $reads2,$readsu,$contact);
	
	rezipFiles($reads1);
	rezipFiles($reads2);
	rezipFiles($readsu);
	}
	}
	
sub rezipFiles {
	my ($file) = @_;
	
	my $call = system("gzip -1 $file");
	}

sub runAligners {
	my ($assembly, $np, $forward,$reverse,$unpaired,$contact) = @_;
	
	my @files;

	#BWA
	unless (-f $rootdir . $contact . ".bwa.bam") {
		my $call1 = system("bwa index $assembly");
		my $start1 = time;
		my $call2 = system("bwa aln -n 0.05 -t $np $assembly $forward > $forward" . ".sai");
		my $call3 = system("bwa aln -n 0.05 -t $np $assembly $reverse > $reverse" . ".sai");
		my $call4 = system("bwa aln -n 0.05 -t $np $assembly $unpaired > $unpaired" . ".sai");
		my $call5 = system("bwa sampe -n 10 $assembly $forward.sai $reverse.sai $forward $reverse > bwa1.sam");
		my $call6 = system("bwa samse -n 10 $assembly $unpaired.sai $unpaired > bwa2.sam");
		my $end1 = int((time - $start1)/60);
		print OUT "bwa\t$contact\t$end1\n";
		my $call7 = system("samtools view -bS bwa1.sam > bwa1.bam");
		my $call8 = system("samtools view -bS bwa2.sam > bwa2.bam");
		my $call9 = system("samtools merge $rootdir" . "$contact" . ".bwa.bam bwa1.bam bwa2.bam");
		push(@files, $rootdir . $contact . ".bwa.bam");
		my $call10 = system("rm *sam bwa1.bam bwa2.bam $forward*sai $reverse*sai");
		}
	
	
	#novoalign
	unless (-f $rootdir . $contact . ".novo.bam") {
		my $call11 = system("novoindex $assembly.novo $assembly");
		my $start2 = time;
		my $call12 = system("novoalign -d $assembly.novo -F ILM1.8 -o SAM -f $forward $reverse > novo1.sam");
		my $call13 = system("novoalign -d $assembly.novo -F ILM1.8 -o SAM -f $unpaired > novo2.sam");
		my $end2 = int((time - $start2)/60);
		print OUT "novo\t$contact\t$end2\n";
		my $call14 = system("samtools view -bS novo1.sam > novo1.bam");
		my $call15 = system("samtools view -bS novo2.sam > novo2.bam");
		my $call16 = system("samtools merge $rootdir" . "$contact" . ".novo.bam novo1.bam novo2.bam");
		push(@files, $rootdir . $contact . ".novo.bam");
		my $call17 = system("rm *sam novo1.bam novo2.bam");
		}

	#smalt
	unless (-f $rootdir . $contact . ".smalt.bam") {
		my $call18 = system("smalt index $assembly.smalt $assembly");
		my $start3 = time;
		my $call19 = system("smalt map -f sam -l pe -n $np -o smalt1.sam $assembly.smalt $forward $reverse");
		my $call20 = system("smalt map -f sam -n $np -o smalt2.sam $assembly.smalt $unpaired");
		my $end3 = int((time - $start3)/60);
		print OUT "smalt\t$contact\t$end3\n";
		
		my $calla = system("samtools faidx $assembly");
		my $call21 = system("samtools view -bS -t $assembly.fai smalt1.sam > smalt1.bam");
		my $call22 = system("samtools view -bS -t $assembly.fai smalt2.sam > smalt2.bam");
		my $call23 = system("samtools merge $rootdir" . "$contact" . ".smalt.bam smalt1.bam smalt2.bam");
		push(@files, $rootdir . $contact . ".smalt.bam");
		my $call24 = system("rm *sam smalt1.bam smalt2.bam");
		}

	#bowtie1
	unless (-f $rootdir . $contact . ".bowtie.bam") {
		my $call33 = system("bowtie-build $assembly $assembly.bowtie1");
		my $start5 = time;
		my $call34 = system("bowtie $assembly.bowtie1 -1 $forward -2 $reverse -5 5 -3 5 -X 300 -a --best -e 200 -S bowtie1.sam");
		my $call35 = system("bowtie $assembly.bowtie1 $unpaired -5 5 -3 5 -a --best -e 200 -S bowtie2.sam");
		my $end5 = int((time - $start5)/60);
		print OUT "bowtie\t$contact\t$end5\n";
		my $call36 = system("samtools view -bS bowtie1.sam > bowtie1.bam");
		my $call37 = system("samtools view -bS bowtie2.sam > bowtie2.bam");
		my $call38 = system("samtools merge $rootdir" . "$contact" . ".bowtie.bam bowtie1.bam bowtie2.bam");
		push(@files, $rootdir . $contact . ".bowtie.bam");
		my $call39 = system("rm *sam bowtie1.bam bowtie2.bam");
		}

	#bowtie2
	unless (-f $rootdir . $contact . ".bowtie2.bam") {
		my $call40 = system("bowtie2-build $assembly $assembly.bowtie2");
		my $start6 = time;
		my $call41 = system("bowtie2 -x $assembly.bowtie2 -1 $forward -2 $reverse -S bow2tie1.sam -5 5 -3 5 --sensitive -k 10 -X 300 -p $np --no-mixed");
		my $call42 = system("bowtie2 -x $assembly.bowtie2 $unpaired -S bow2tie2.sam -5 5 -3 5 --sensitive -k 10 -p $np");
		my $end6 = int((time - $start6)/60);
		print OUT "bowtie2\t$contact\t$end6\n";
		my $call43 = system("samtools view -bS bow2tie1.sam > bow2tie1.bam");
		my $call44 = system("samtools view -bS bow2tie2.sam > bow2tie2.bam");
		my $call45 = system("samtools merge $rootdir" . "$contact" . ".bowtie2.bam bow2tie1.bam bow2tie2.bam");
		push(@files, $rootdir . $contact . ".bowtie2.bam");
		my $call46 = system("rm *sam bow2tie1.bam bow2tie2.bam");
		}

	#stampy -- cannot run on that machine
	#unless (-f $rootdir . $contact . ".stampy.bam") {
	#	my $call25 = system("/home/singhal/programs/stampy/stampy.py --species=lizard --assembly=trinity -G $assembly $assembly");
	#	my $start4 = time;
	#	my $call26 = system("/home/singhal/programs/stampy/stampy.py -g $assembly -H $assembly");
	#	my $call27 = system("/home/singhal/programs/stampy/stampy.py -g $assembly -h $assembly -o stampy1.sam --substitutionrate=0.05 --insertsize=110 -M $forward $reverse");
	#	my $call28 = system("/home/singhal/programs/stampy/stampy.py -g $assembly -h $assembly -o stampy2.sam --substitutionrate=0.05 -M $unpaired");
	#	my $end4 = int((time - $start4)/60);
	#	print OUT "stampy\t$contact\t$end4\n";
	#	my $call29 = system("samtools view -bS stampy1.sam > stampy1.bam");
	#	my $call30 = system("samtools view -bS stampy2.sam > stampy2.bam");
	#	my $call31 = system("samtools merge $rootdir" . "$contact" . ".stampy.bam stampy1.bam stampy2.bam");
	#	push(@files, $rootdir . $contact . ".stampy.bam");
	#	my $call32 = system("rm *sam stampy1.bam stampy2.bam");
	#	}

	#soap2 -- cannot run on that machine
	#unless (-f $rootdir . $contact . ".soap.bam") {
	#	my $call47 = system("2bwt-builder $assembly");
	#	my $start7 = time;
	#	my $call48 = system("soap -a $forward -b $reverse -D $assembly.index -o soap1 -2 out.sam -m 0 -x 500 -r 2 -v 4 -p $np");
	#        my $call49 = system("soap -a $unpaired -D $assembly.index -o soap2 -r 2 -v 4 -p $np");
	#	my $callb = system("soap2sam.pl soap1 > soap1.sam");
	#	my $callc = system("soap2sam.pl soap2 > soap2.sam");
	#	my $end7 = int((time - $start7)/60);
	#	print OUT "soap2\t$contact\t$end7\n";
	#	my $calla = system("samtools faidx $assembly");
	#	my $call50 = system("samtools view -bS -t $assembly.fai soap1.sam > soap1.bam");
	#	my $call51 = system("samtools view -bS -t $assembly.fai soap2.sam > soap2.bam");
	#	my $call52 = system("samtools merge $rootdir" . "$contact" . ".soap.bam soap1.bam soap2.bam");
	#	push(@files, $rootdir . $contact . ".soap.bam");
	#	my $call53 = system("rm *sam soap1.bam soap2.bam soap1 soap2");
	#	}
	
	return(\@files);
	}

sub unzipFiles {
	my ($file) = @_;
	my $call = system("gunzip $file");
	my $file2 = $1 if $file =~ m/(.*)\.gz/;
	return($file2);
	}

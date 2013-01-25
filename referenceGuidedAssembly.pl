use warnings;
use strict;
use Devel::Size qw(total_size);
use FileHandle;

#reference guided assembly
my $dir = '/Users/singhal/Desktop/genomics/';
my $file_dir = $dir . 'contigFiles/';
my $seq = $dir . 'seqFiles/Carlia_N_trinity.fa.final.annotated';
my $bam = $dir . 'carliaN.sorted.bam';
my $reads1 = $dir . 'Carlia_N_1.fastq';
my $reads2 = $dir . 'Carlia_N_2.fastq';
my $readsu = $dir . 'Carlia_N_u.fastq';

makeContigFiles();
#doAssembly();

sub doAssembly {
	my @files = <$file_dir*u.fastq>;
	foreach my $u (@files) {
		my $p1 = $1 . '1.fastq' if $u =~ m/(\S+_)u/;
		my $p2 = $1 . '2.fastq' if $u =~ m/(\S+_)u/;
	
		my $c = $1 if $u =~ m/(contig\d+)/;
		
		my $abyss = system("abyss-pe n=5 s=130 e=0 E=0 c=2 in=\'$p1 $p2\' se=$u name=$c" . ".assembly k=51");
		}
	}

sub makeContigFiles {
	my %seq;
	open(IN, "<$seq");
	while(<IN>) {
		if ($_ =~ m/>(\S+).*gs/) {
			$seq{$1}++ unless(-f $file_dir . $1 . '_u.fastq');
			}
		}
	close(IN);
		
	my %reads;	
	foreach my $c (keys %seq) {	
		#this hash is getting too big! time to dump it	
		if (total_size(\%reads) > 3e7) {
			print "dumping!\n";
		
			#take care of paired reads
			open(IN1, "<$reads1");
			open(IN2, "<$reads2");
			while(<IN1>) {
				chomp(my $read_id1 = $_);
				chomp(my $seq1 = <IN1>);
				chomp(my $qual_id1 = <IN1>);
				chomp(my $qual1 = <IN1>);
			
				chomp(my $read_id2 = <IN2>);
				chomp(my $seq2 = <IN2>);
				chomp(my $qual_id2 = <IN2>);
				chomp(my $qual2 = <IN2>);
			
				my $id = $1 if $read_id1 =~ m/\@(\S+)_[F|R]/;
				if ($reads{'p'}{$id}) {
					foreach my $contig (keys %{$reads{'p'}{$id}}) {
						my $file1 = $file_dir . $contig . '_1.fastq';
						my $file2 = $file_dir . $contig . '_2.fastq';
						open(OUT1, ">>$file1");
						open(OUT2, ">>$file2");
						print OUT1 $read_id1 . "\n" . $seq1 . "\n" . $qual_id1 . "\n" . $qual1 . "\n";
						print OUT2 $read_id2 . "\n" . $seq2 . "\n" . $qual_id2 . "\n" . $qual2 . "\n";
						close(OUT1); close(OUT2);
						}
					}	
				}
			close(IN1); close(IN2);
		
			#do unpaired reads
			open(IN, "<$readsu");
			while(<IN>) {
				chomp(my $read_id = $_);
				chomp(my $seq = <IN>);
				chomp(my $qual_id = <IN>);
				chomp(my $qual = <IN>);
			
				my $id = $1 if $read_id =~ m/\@(\S+)/;
				if ($reads{'up'}{$id}) {
					foreach my $contig (keys %{$reads{'up'}{$id}}) {
						my $file = $file_dir . $contig . '_u.fastq';
						open(OUT, ">>$file");
						print OUT $read_id . "\n" . $seq . "\n" . $qual_id . "\n" . $qual . "\n";
						close(OUT);
						}
					}
				}
			close(IN);	
			
			%reads = ();
			}
		else {	
			my @call = `samtools view $bam $c`;
			foreach my $line (@call) {
				my @d = split(/\t/,$line);
				#paired end!
				if ($d[0] =~ m/(\S+)_[F|R]/) {
					$reads{'p'}{$1}{$c}++;
					}
				else {
					$reads{'up'}{$d[0]}{$c}++;
					}
				}	
			}
		}
	}	
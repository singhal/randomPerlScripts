use Bio::Tools::Run::RemoteBlast;
use Bio::SeqIO;
use strict;
use warnings;

#http://doc.bioperl.org/releases/bioperl-current/bioperl-live/Bio/Tools/Run/RemoteBlast.html#CODE1

my $prog = 'blastn';
my $db   = 'nr';
my $e_val= '1e-10';
my $result_dir = '/media/DataDrive/sutureGenomics/annoStats/orphans/';
my $tracker = 1;
my $numBlast = 70000;
my $out = $result_dir . "annotation.out";
my $seq = '/media/DataDrive/sutureGenomics/Carlia_N/assemblies/Carlia_N_trinity.fa.final';
my $blast = '/media/DataDrive/sutureGenomics/Carlia_N/assemblies/Carlia_N_trinity.fa.final.blast.out';

open(OUT, ">$out");

mkdir($result_dir) unless (-d $result_dir);
my $file = $result_dir . 'leftovers.fa';
open(IN, "<$seq"); my %seq;
#open(SEQ, ">$file"); 
my $c;
while(<IN>) {
    chomp(my $line = $_);
    if ($line =~ m/>(\S+)/) {
	$c = $1;
    }
    else {
	$seq{$c} .= $line;
    }
}
close(IN);
open(IN, "<$blast");
while(<IN>) {
    my $c = $1 if $_ =~ m/^(\S+)/;
    delete($seq{$c}) if $seq{$c};
}
close(IN);
#foreach my $c (keys %seq) {
#    print SEQ ">", $c, "\n", $seq{$c}, "\n";
#}
#close(SEQ);

my @params = ( '-prog' => $prog,
         '-data' => $db,
         '-expect' => $e_val,
         '-readmethod' => 'SearchIO' );
my $factory = Bio::Tools::Run::RemoteBlast->new(@params);

#$v is just to turn on and off the messages
my $v = 1;
my $str = Bio::SeqIO->new(-file=>$file , -format => 'fasta' );

while (my $input = $str->next_seq()){
#Blast a sequence against a database:

    if ($tracker < $numBlast && $tracker > 54000) {
	my $r = $factory->submit_blast($input);

	print STDERR "waiting..." if( $v > 0 );
	while ( my @rids = $factory->each_rid ) {
	    foreach my $rid ( @rids ) {
		my $rc = $factory->retrieve_blast($rid);
		if( !ref($rc) ) {
		    if( $rc < 0 ) {
			$factory->remove_rid($rid);
		    }
		    print STDERR "." if ( $v > 0 );
		    sleep 5;
		} else {
		    my $result = $rc->next_result();
		    #save the output
		    my $filename = $result_dir . $tracker."\.out";
		    $factory->save_output($filename);
		    $factory->remove_rid($rid);
		    print OUT "\nQuery Name: ", $result->query_name(), "\n";
		    while ( my $hit = $result->next_hit ) {
			next unless ( $v > 0);
			print OUT "\thit name is ", $hit->name, "\n";
			while( my $hsp = $hit->next_hsp ) {
			    print OUT "\t\tscore is ", $hsp->score, "\n";
			}
		    }
		}
	    }
	}
    }
    $tracker++;
}

close(OUT);

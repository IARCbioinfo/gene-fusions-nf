
###############################################################################
# Author: Alex Di Genova 
# Laboratory: ERABLE/INRIA
# Copyright (c)
# year: 2017
###############################################################################
use Data::Dumper;
use Getopt::Std;
use Bio::SeqIO;
use strict;

sub usage {
   print "$0 usage : -a <fasta transcripts>  -b <saples prefix:T000>  -c <number of samples,def:10>\n";
   print "Error in use\n";
   exit 1;
}

my %opts = ();
getopts( "a:b:c:", \%opts );
if ( !defined $opts{a} or !defined $opts{b} ) {
   usage;
}

#control the number of samples selected
my $n_samples=11;
if(defined $opts{c}){
	$n_samples=$opts{c};
	$n_samples++;
}
#read the transcript sequences
my $file=$opts{a};
my $in=Bio::SeqIO->new(-format=>'Fasta',-file=>$file) or die "cannot open seq file $file\n";
my $hash=();
while(my $seq=$in->next_seq()){
     $hash->{$seq->display_id()}=$seq->seq();
}
#print Dumper($hash);
my @k=sort keys %{$hash};
open(LOG,">".$opts{b}.".samples.log") or die "cannot open file $opts{b}.samples.log\n";

#we simulate the reads using wgssim 
for(my $s=1; $s<$n_samples; $s++){
my $sample=$opts{b}.$s;

my $used=();
my $chimeric=();
my $max_sim=int(rand()*10+1);
print LOG join(" ","LOG","SAMPLE=$sample","Number of fusion transcripts = $max_sim")."\n";
print join(" ","LOG","SAMPLE=$sample","Number of fusion transcripts = $max_sim")."\n";
#we create 5 gene fusion
for(my $i=0; $i<$max_sim; $i++){
	my $rand1=0;
	my $rand2=0;
	do{
		$rand1=int(rand()*$#k)+1;
		$rand2=int(rand()*$#k)+1;
	}while($rand1 == $rand2 or length($hash->{$k[$rand1]}) <500 or length($hash->{$k[$rand2]})< 500);

	#print join(" ",$rand1,$rand2)."\n";
	my $start1=1;
	my $stop1=int(length($hash->{$k[$rand1]})/2);
	my $start2=int(length($hash->{$k[$rand2]})/2);
	my $stop2=int(length($hash->{$k[$rand2]}));
	my $p_1=substr($hash->{$k[$rand1]},$start1,$stop1);
	my $p_2=substr($hash->{$k[$rand2]},$start2,abs($stop2-$start2));
	my $c_name=join("__",$k[$rand1],$k[$rand2]);
	$used->{$k[$rand1]}=1;
	$used->{$k[$rand2]}=1;
	#print join(" ",$c_name,length($p_1),length($p_2),length($p_1.$p_2),$p_1.$p_2)."\n";
	#print join("\n",">".$c_name,$p_1.$p_2)."\n";
	$chimeric->{$c_name}=$p_1.$p_2;
}

foreach my $t(@k){
	next if (defined $used->{$t});
	#print join("\n",">".$t,$hash->{$t})."\n";
	$chimeric->{$t}=$hash->{$t};
}

foreach my $t(keys %{$chimeric}){
	my $n_reads=rand()*5000;
	$n_reads=50 if($n_reads < 50);
	print LOG join(" ","LOG","SAMPLE=$sample","T=$t","N=".int($n_reads),"C=".int($n_reads*200/length($chimeric->{$t})))."\n";
	open(FILE,">ref.tmp") or die "cannot open file ref.tmp\n";
	print FILE join("\n",">".$t,$chimeric->{$t})."\n";
	#my $cmd="./wgsim -d 200 -1 100 -2 100 "
	my $cmd=" ./wgsim/wgsim  -d 200 -N $n_reads -1 100 -2 100 ref.tmp r1.tmp.fastq r2.tmp.fastq >wgsim.log 2>wgsim.err";
	system($cmd);
	system("cat r1.tmp.fastq >> $sample.R1.fastq");
	system("cat r2.tmp.fastq >> $sample.R2.fastq");
	close(FILE);
  }
  #we remove the temporary files
  system("rm -f ref.tmp r1.tmp.fastq r2.tmp.fastq wgsim.log wgsim.err");		
  #we compress the simulated sample, overrride if file exist
  system("gzip -f $sample.R1.fastq");
  system("gzip -f $sample.R2.fastq");
	
}

close(LOG);


use strict;
use constant COIMAP => " files/cois_mapping.txt"; # CSV that maps COI name to taxonomic group
use constant COILEN => " files/db.len"; # map each COI to it's sequence length
use constant PARSESCRPIT => " scripts/parseBlastResults.pl";
use constant isCOIREADS => 1;
use constant DATA_SUFFIX => ".final";
use constant SUBMIT_JOB => "bsub -q sorek"; 
open (D,"dirs"); #TXT file with the directory names (one per sample)
my @dirs = <D>;
close(D);
chomp @dirs;
open (OUT,">blastUnmappedReads".DATA_SUFFIX);
open(S,">sortBlastResults".DATA_SUFFIX.".sh");
foreach my $dir (@dirs){
	my $combainedBlast = "$dir/$dir" ."blastn" . DATA_SUFFIX;
	if(isCOIREADS) {
		print S SUBMIT_JOB . " sort $dir/*coi.reads.blastn.COI" .DATA_SUFFIX. " -k1,1 -k3,3nr -o $combainedBlast;\n";
	}
	else {
		print S SUBMIT_JOB . " sort $dir/*.0N.fasta.blastn.COI -k 1,1 -k3,3nr -o $combainedBlast;\n";
	}
	my $outFile = "$dir/$dir.pipe.blast.parsed" .DATA_SUFFIX;
	print SUBMIT_JOB .PARSESCRPIT . COIMAP . COILEN . " $combainedBlast $outFile;\n";
	#Extract unmapped reads;
	my $unmappedFile = "$outFile.unmapped";
	my $unmappedReads = "$unmappedFile.reads.name";
	my $unmappedSequence = "$unmappedFile.reads";
	print OUT "cut -f 1 $unmappedFile > $unmappedReads;\n";
	if (isCOIREADS){
		print OUT "cat $dir/*coi.reads |grep -A 1 -f $unmappedReads |grep -v '\\-\\-' > $unmappedSequence;\n";
	}
	else {
	print OUT "cat $dir/*.fasta |grep -A 1 -f $unmappedReads |grep -v '\\-\\-' > $unmappedSequence;\n";
	}
}
close(S);
close(OUT);

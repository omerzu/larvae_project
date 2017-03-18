use strict;
use constant DBPATH => "db.fa ";
use constant BLAST_ARGS => "-F F -e 1e-10 -m 8";
use constant isCOIreads => 1; # Use previously identified COI reads instead of working on the entire raw reads
use constant DATA_SUFFIX => ".final";
use constant SUBMIT_JOB => "bsub -q sorek blastall -p blastn -d "; # Replace this line tor work with distributed system.
my $path = "fa.files"; #line separated file with all paths the all .fa files 
if (isCOIreads) {$path = "coi.reads";} # replace the original path file with the file the points to COI reads only.
open (F,$path);
my @fs = <F>;
chomp @fs;
foreach my $fa (@fs){
	my $resFile = $fa.".blastn.COI".DATA_SUFFIX;
	print SUBMIT_JOB .DBPATH .BLAST_ARGS." -i $fa -o $resFile;\n";
}
  

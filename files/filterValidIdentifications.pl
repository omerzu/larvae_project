#Read the outfile with the coverage and print the fish with sufficient coverage
use strict;
use constant MIN_COV => 50;
my ($in) = @ARGV;
open (IN, $in) || die "Error while opening $in\n";
my $sample;
#print "Sample\tOrder\Family\tGenus\Species\tReads\tCOV\tMismatchs\tCOI\n";
print "Sample\tOrder\tFamily\tGenus\tSpecies\tReads\tAmbRead\tCoverage\tMismatches\tCOI\tMutations\n";
while (my $l = <IN>){
	if ($l =~ m/Sample: (\d+[A,B]?),/) {
		$sample = $1;
	}
	if ($l =~/,,/) { #line with larvae report
		my @data = split(/,/,$l);
		my($sp,$reads,$amb,$cov,$mismatch,$coi,$misDetails) = ($data[2],$data[3],$data[5],$data[8],$data[9],$data[12],$data[13]);
		my ($o,$f,$g,@as) = split(/_/,$sp);
		my $s = join("_",@as);
		if ($cov >= MIN_COV) {
			#print "$sample\t$o\t$f\t$g\t$s\t$reads\t$cov\t$mismatch\t$coi\n";
			print "$sample\t$o\t$f\t$g\t$s\t$reads\t$amb\t$cov\t$mismatch\t$coi\t$misDetails";
		}
	}
}
close (IN);

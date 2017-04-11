#The script merge the two main results file: the file that summerize the identification of all samples (createSummeryOfSamples) and the reports from the pipeline
use strict;
my ($pipeResults, $mainOutput) = @ARGV;
open (R,$pipeResults) || die "Error while opening $pipeResults\n";
my %sampleAndFishTocoverageAndMismath;
my $line = <R>; #ignore header
while ($line = <R>) {
	chomp $line;
	my ($sample,$fish,$coi,$coveredbps,$coilength,$coverage,$mismatch,$bpsCovered,$mm) = split (/\s+/,$line);
	$sampleAndFishTocoverageAndMismath{lc("$sample\_$fish")} = ",$coverage,$mismatch,$bpsCovered,$coilength,$coi,$mm";
}
close(R);
#read the results file and identify the fish
open (RES,$mainOutput) ||die "Error while opening $mainOutput\n";
my $sample = "";
while (my $line = <RES>) {
	chomp $line;
	if ($line =~ m/Sample: (\d+)/) {
		$sample = $1;
	}
	if ($line =~ m/,,([\w,\.]*?),/) {
	#if ($line =~ m/,,(.*?),/) {
		my $sampleFish = lc("$sample\_$1");
		print $line  . $sampleAndFishTocoverageAndMismath{$sampleFish}."\n";
	}
	else {
		print "$line\n";
	}
}
close(RES);

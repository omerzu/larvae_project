#The script reads the mapping results and print the following data: % unmapped reads (mean and total)
# counts the total #reads mapped to each fish and the #of samples it appears
#Also print for each fish the list of samples that it appears in
use String::Util 'trim';
use List::Util qw(max);
use List::MoreUtils qw(uniq);
use List::Util qw(sum);
use Math::Round qw(:all);
use constant false => 0;
use constant true  => 1;
use constant DEBUG_AMB  => 1;
sub readManualData;
sub printBlastFishResults; 
sub printManualResults;
sub readRedSeaAndGulfData; 
sub intersect;
use constant DATA_SUFFIX => "";
my ($resultsPath,$manualResults,$samplesPath,$sampleToTotalFishFile,$GulfPath,$RedSeaPath) = @ARGV;
@ARGV == 6 or die "Usage:<path to file with the results>\t<path to maual(Naama) results>\t<path to file with the samples name>\t<sample number to total fish in tube>\t<path to file with the fish in the Gulf>\t<path to file with the fish in the red sea>\n";
my (%fishToReads,%fishToSamples,%missedFish,%manualMissed,$totalReads,$totalUnmappedReads,@unmappedFractions);
#Parse manual results
my ($SampleToManualFamilyRef,$Sample_FishToCountRef) = readManualData($manualResults,$samplesPath);
my %SampleToManualFamily = %{$SampleToManualFamilyRef};
my %Sample_FishToCount = %{$Sample_FishToCountRef};
my ($allMappedRead,$allUnmappedReads) = (0,0);
my ($fishInGulfRef,$fishInRedRef) = readRedSeaAndGulfData($GulfPath,$RedSeaPath);
my (%hallMappedRead,%hallUnmappedReads,%hSampleToMappedReads);
my %spToSamples;
#Get the files paths
open(F, $resultsPath) or die "Error while opening $resultsPath\n";
my @files = <F>;
chomp @files;
close(F);
my %sampleToTotalFishCount;
open(S,$sampleToTotalFishFile) or die "Error while opening $sampleToTotalFishFile\n";
while (my $l = <S>){chomp$l; my ($sample,$count)=split(/\s+/,$l); $sampleToTotalFishCount{$sample} = $count;}
close(S);
#Iterate over the samples
foreach my $file (@files) {
	my %htotalReads;
	my @path = split(/\//,$file);
	my $samplePath = pop(@path);
	my @splitedPath = split(/\./,$samplePath);
	my $sampleID = $splitedPath[0];
	open(SAMPLE,$file) or die "Erorr while opening $file\n";
	my $line = <SAMPLE>; #Ignore header
	my %fishToReadsInSample;
	my %fishToPIDsInSample;
	my %fish_ambiguousToReadsInSample;
	my %ambiguousReadsToFishInSample;
	my %fishFamilyToFishs; #maps Family to lower taxonomic units (genus,sp.)
	my $totalReads = 0;
	while (my $line = <SAMPLE>) { #working on specific sample
		chomp $line;
		if (!(substr($line,0,1) eq "#")){ #ignore the ambiguous reads, will handle them later 
		my ($Read,$COI,$pid,$Order,$Family,$Genus,$Species)	= split (/\s+/,$line);
		if (substr($line,0,1) eq "*"){ #Same read mapped to more than one COI of the fish, count it only once
			next if exists ($hallMappedRead{$Read});
		}
		$hallMappedRead{$Read} = 1;
		$htotalReads{$Read} = 1;
		my @sampleMappedReads; if(exists ($hSampleToMappedReads{$sampleID})){@sampleMappedReads = @{$hSampleToMappedReads{$sampleID}}};
		push(@sampleMappedReads,$Read);
		$hSampleToMappedReads{$sampleID} = \@sampleMappedReads; 
		my $fishFamily = "$Order" ."_$Family";
		my @fishTax = ($Order,$Family,$Genus,$Species);
		my $fishID = join("_",@fishTax);
		#Store the fish in the right family and update the reads that mapped to it
		my @ids;
		if (exists $fishFamilyToFishs{$fishFamily}){@ids = @{$fishFamilyToFishs{$fishFamily}}};
		push (@ids, $fishID);
		@ids = uniq(@ids);
		$fishFamilyToFishs{$fishFamily} = \@ids;
		my $fishReads = 1;
		my @pids;		
		if (exists $fishToReadsInSample{$fishID}) {
			$fishReads = $fishReads + $fishToReadsInSample{$fishID};
			@pids = @{$fishToPIDsInSample{$fishID}};
		}
		push(@pids,$pid);
		$fishToPIDsInSample{$fishID} = \@pids;
		$fishToReadsInSample{$fishID} = $fishReads;
		my $fishReads = 1;
		if (exists($fishToReads{$fishID})) {$fishReads = 1 + $fishToReads{$fishID};}
		$fishToReads{$fishID} = $fishReads;
		my @fishSamples;
		if (exists($fishToSamples{$fishID})) {@fishSamples = @{$fishToSamples{$fishID}};}
		push(@fishSamples,$sampleID);
		@fishSamples = uniq(@fishSamples);
		$fishToSamples{$fishID} = \@fishSamples;
		}
		else { #Ambiguous Hit
			my ($Read,$COI,$pid,$Order,$Family,$Genus,$Species)	= split (/\s+/,$line);
			$htotalReads{$Read} = 1;
			$hallMappedRead{$Read} = 1;
			my $fishFamily = "$Order" ."_$Family";
			my @fishTax = ($Order,$Family,$Genus,$Species);
			my $fishID = join("_",@fishTax);
			my @fishReads;
			if (exists $fish_ambiguousToReadsInSample{$fishID}) {
				@fishReads = @{$fish_ambiguousToReadsInSample{$fishID}};
			}
			push(@fishReads, $Read); 
			@fishReads = uniq(@fishReads);
			#print "$fishID\t". join(",",@fishReads);
			$fish_ambiguousToReadsInSample{$fishID} = \@fishReads;
			my @fishes;
			if (exists $ambiguousReadsToFishInSample{$Read}) {
				@fishes = @{$ambiguousReadsToFishInSample{$Read}};
			}
			push(@fishes,$fishID);
			@fishes = uniq(@fishes);
			#print "$Read\t".join(",",@fishes)."\n";
			$ambiguousReadsToFishInSample{$Read} = \@fishes;
		}
	}
	close(SAMPLE);
	#Convert the PID is sample to the mean
	foreach my $fish (keys (%fishToPIDsInSample)) {
		my @pids = @{$fishToPIDsInSample{$fish}};
		$fishToPIDsInSample{$fish} = nearest(.1,sum(@pids) / scalar (@pids));
	}
	my $unmappedPath = "$file.unmapped";
	open(UNM, $unmappedPath) or die "Error while openinng $unmappedPath\n";
	my @unmappedReads = <UNM>;
	close(UNM);
	my $unmappedReadsCount = scalar (@unmappedReads) - 1; #-1 is due to the header
	$totalReads = scalar (keys %htotalReads) +  scalar (@unmappedReads) - 1; #add the unmapped reads to the reads count of the current file
	$allUnmappedReads = $allUnmappedReads + $unmappedReadsCount;
	$allMappedRead = scalar (keys %hallMappedRead);
	my $sampleMapped = scalar (uniq(@{$hSampleToMappedReads{$sampleID}}));
	open(RC,"$sampleID/$sampleID.all.reads.cnt") or die "Error while opening $sampleID/$sampleID.all.reads.cnt\n";
	my $allReads = <RC>;chomp $allReads;
	close(RC);
	open(my $OUT,">$file.out" . DATA_SUFFIX) or die "Error while opening $file.out\n";
	my $sid = $sampleID;
	$sid =~ s/\D//g;
	for my $fish (keys %fishToReadsInSample) {
		my @samples; if (exists $spToSamples{$fish}){@samples = @{$spToSamples{$fish}};}
		push(@samples,$sid);
		$spToSamples{$fish} = \@samples;
	}
	#Print the mapping information
	#Add the unmapped reads
	my $mappedFraction = nearest(.01,(($totalReads-$unmappedReadsCount)/$totalReads *100));
	my $identifedSp = scalar(keys %fishToReadsInSample);
	if (exists $fishToReadsInSample{"Perciformes_Cichlidae_Tilapia_zillii"}) {$identifedSp = $identifedSp-1;}
	print $OUT "Sample: $sid,Larvae in Sample: $sampleToTotalFishCount{$sid},Number of species identified: $identifedSp\n";
	print $OUT "All Reads:,$allReads,COI Reads:,$totalReads,Percent COI from total:,".nearest(.001,100*$totalReads/$allReads).",mapped reads:,$sampleMapped,unmapped reads:,$unmappedReadsCount,Percent mapped:,$mappedFraction\n";
	print $OUT "Morphology_Counts,Larvae_Family,Larvae_id,#MappedReads,MeanPID,#AmbiguousReads,isGulf,isRedSea,%COIcovered,#mismatch,TotalBpsCovered,COIlengeh\n";
	#find the common families and print them first
	
	my @manualFamilies = (); if (exists ($SampleToManualFamily{$sid})){@manualFamilies = @{$SampleToManualFamily{$sid}};}
	my @blastFamilies = keys(%fishFamilyToFishs);
	my @commonFamilies = sort (intersect(@manualFamilies,@blastFamilies));
	#print "sampleID = $sampleID\t". join (",",@blastFamilies) . "\n";
	my @blastUniq = sort(array_minus(@blastFamilies,@manualFamilies));
	#Store the fish that were missed manually
	foreach my $family(@blastUniq) {
		my @samples;
		if (exists($manualMissed{$family})) {@samples = @{$manualMissed{$family}}};
		push(@samples,$sampleID);
		$manualMissed{$family} = \@samples;
	}
	my @manUniq = sort(array_minus(@manualFamilies,@blastFamilies));
	foreach my $manfamily (@manUniq){my @samples; if (exists($missedFish{$manfamily})){@samples = @{$missedFish{$manfamily}}};push(@samples,$sampleID);@samples = uniq(@samples); $missedFish{$manfamily} = \@samples;}
	#printResults(\@manualFamilies,\%fishToReadsInSample,\%fish_ambiguousToReadsInSample,\%fishFamilyToFishs,\%fishToPIDsInSample,\%ambiguousReadsToFishInSample,$sampleID,$OUT,$fishInGulfRef,$fishInRedRef);
	printResults(\@manualFamilies,\%fishToReadsInSample,\%fish_ambiguousToReadsInSample,\%fishFamilyToFishs,\%fishToPIDsInSample,\%ambiguousReadsToFishInSample,$sid,$OUT,$fishInGulfRef,$fishInRedRef);
	close($OUT);
	#Print the list of samples for each fish;
}
#for my $fish (keys %spToSamples) {
#		print "$fish\t" . join (",", @{$spToSamples{$fish}}) . "\n";
#	}
#@in -path to csv file with the manual identification:New Name,Number of larvae,Order,Family
#@out - hash sample->@fish_family
sub readManualData {
	my ($path,$samplesPath)  = @_;
	open(S,$samplesPath) or die "Error while opening $samplesPath\n";
	my @samples = <S>;
	close (S);
	my %samplesHash; #samplesHash ->make sure that reads data corresponds to sequenced samples
	chomp @samples;
	foreach my $s (@samples){my $sid = $s;$sid =~ s/\D//g;$samplesHash{$sid} = 1}; #handle non-numeric samples name
	open (F,$path) or die "Error while opening $path\n";
	my (%SampleToFamily,%Sample_FishToObserved); #%SampleToFamily: maps sample ID to the fish family in it. #Sample_FishToCount ->Counts for each fis(Order_family) in the samples how many fishes Observed
	my $line = <F>; #ignore header
	while ($line = <F>) {
		chomp $line;
		my ($sample,$count,$Order,$Family) = split(/\t/,$line);
		next if (!(exists ($samplesHash{$sample})));
		$Order = trim($Order);
		$Family = trim($Family);
		my $fish = "$Order"."_$Family";
		if (my $Family eq "UNK") {$fish = "$Order"}; #Use Order instead of Order_UNK
		my @fishInSample;
		if (exists $SampleToFamily{$sample}) {@fishInSample = @{$SampleToFamily{$sample}}};
		push(@fishInSample,$fish);
		@fishInSample = uniq (@fishInSample);
		$SampleToFamily{$sample} = \@fishInSample;
		if (exists $Sample_FishToObserved{"$sample"."_$fish"}) {$count = $count + $Sample_FishToObserved{"$sample"."_$fish"};}#Sometime the same family appears in different lines-> combine counts
		$Sample_FishToObserved{"$sample"."_$fish"} = $count;
	}
	close(F);
	return(\%SampleToFamily,\%Sample_FishToObserved);
}
#calculate and print for each fish the data:fishID\t%mappedReads\ttotalReads\tmeanPID

sub printResults {
	my ($manualFamiliesRef,$fishToReadsInSampleRef,$fish_ambiguousToReadsInSampleRef,$fishFamilyToFishsRef,$fishToPIDsInSampleRef,$ambiguousReadsToFishInSampleRef,$sample,$OUT,$fishInGulfRef,$fishInRedRef) = @_;
	my @manualFamilies = @{$manualFamiliesRef};
	my %fish_ambiguousToReadsInSample = %{$fish_ambiguousToReadsInSampleRef};
	my %ambiguousReadsToFishInSample = %{$ambiguousReadsToFishInSampleRef};
	my %fishToReadsInSample = %{$fishToReadsInSampleRef};
	my %fishFamilyToFishs = %{$fishFamilyToFishsRef};
	my %fishToPIDsInSample = %{$fishToPIDsInSampleRef};
	my %fishInGulf = %{$fishInGulfRef};
	my %fishInRed = %{$fishInRedRef};
	my %order_unk;
	my %printedFamily;
	foreach my $family(@manualFamilies) {
		my ($o,$f) = split(/_/,$family);
		if ($f eq "UNK") {
			$order_unk{$o}=1;
			next;
		}
		my $sampleFam = "$sample" ."_$family";
		my $manualCount = $Sample_FishToCount{"$sampleFam"};
		print $OUT "$manualCount,$family\n";
		$printedFamily{$family} = 1;
		#Print the relevant blast fish
		if (exists $fishFamilyToFishs{$family}) {
			foreach my $fish (@{$fishFamilyToFishs{$family}}) {
				if (!(exists $fish_ambiguousToReadsInSample{$fish})) {my @t = ();$fish_ambiguousToReadsInSample{$fish} = \@t;};
				my @arrfish = split("_",$fish);
				my $fishName = join("_",($arrfish[2],$arrfish[3]));
				my $isGulf = 0; $isGulf = 1 if  exists ($fishInGulf{lc $fishName});
				my $isRed = 0;$isRed = 1 if exists ($fishInRed{lc $fishName});
				print $OUT ",,$fish,$fishToReadsInSample{$fish},$fishToPIDsInSample{$fish},".scalar(@{$fish_ambiguousToReadsInSample{$fish}}).",$isGulf,$isRed\n";
				#Remove the ambigous reads;
				foreach my $read (@{$fish_ambiguousToReadsInSample{$fish}}){
					if (DEBUG_AMB) {
					delete $ambiguousReadsToFishInSample{$read};
					}
				}
			}
		}
		else {
			print $OUT ",,NA,NA\n";
		}
	}
	#find blast identified fish that can be linked with order_unk
	$printedFamily{"Perciformes_Cichlidae"} = 1; #make sure not to print the cichlid data
	foreach my $family(keys(%fishFamilyToFishs)) {
		next if exists $printedFamily{$family};
		my ($order,$fam) = split(/_/,$family);
		if (exists $order_unk{$order}) {
			$order_unk{$order} = 0;
			my $sampleFam = "$sample" ."_$order\_UNK";
			my $manualCount = $Sample_FishToCount{"$sampleFam"};
			print $OUT "$manualCount,$family\n";
			$printedFamily{$family} = 1;
			#Print the relevant blast fish
			foreach my $fish (@{$fishFamilyToFishs{$family}}) {
			if (!(exists $fish_ambiguousToReadsInSample{$fish})) {my @t = ();$fish_ambiguousToReadsInSample{$fish} = \@t;};
				my @arrfish = split("_",$fish);
				my $fishName = join("_",($arrfish[2],$arrfish[3]));
				my $isGulf = 0; $isGulf = 1 if  exists ($fishInGulf{lc($fishName)});
				my $isRed = 0;$isRed = 1 if exists ($fishInRed{lc($fishName)});
				print $OUT ",,$fish,$fishToReadsInSample{$fish},$fishToPIDsInSample{$fish},".scalar(@{$fish_ambiguousToReadsInSample{$fish}}).",$isGulf,$isRed\n";
				foreach my $read (@{$fish_ambiguousToReadsInSample{$fish}}){
					if (DEBUG_AMB) {
					delete $ambiguousReadsToFishInSample{$read};
					}
				}
			}
		}
	}
	#Print the Order_unk that were not found
	foreach my $order (keys %order_unk) {
		next if $order eq "UNK"; 
		if ($order_unk{$order} == 1) {
			my $sampleFam = "$sample\_" .$order ."_UNK";
			my $manualCount = $Sample_FishToCount{"$sampleFam"};
			print $OUT "$manualCount,$order\n";
			print $OUT ",,NA,NA\n";
		}
	}
	#Print all other families that were not identified manually
	foreach my $family(keys(%fishFamilyToFishs)) {
		next if exists $printedFamily{$family};
		print $OUT "0,$family\n";
		foreach my $fish (@{$fishFamilyToFishs{$family}}) {
			if (!(exists $fish_ambiguousToReadsInSample{$fish})) {my @t = ();$fish_ambiguousToReadsInSample{$fish} = \@t;};
			my @arrfish = split("_",$fish);
			my $fishName = join("_",($arrfish[2],$arrfish[3]));
			my $isGulf = 0; $isGulf = 1 if  exists ($fishInGulf{lc($fishName)});
			my $isRed = 0;$isRed = 1 if exists ($fishInRed{lc($fishName)});
			print $OUT ",,$fish,$fishToReadsInSample{$fish},$fishToPIDsInSample{$fish},".scalar(@{$fish_ambiguousToReadsInSample{$fish}}).",$isGulf,$isRed\n";
			foreach my $read (@{$fish_ambiguousToReadsInSample{$fish}}){
					if (DEBUG_AMB) {
					delete $ambiguousReadsToFishInSample{$read};
					}
			}
		}
	}
	#print ambiguous reads
	my @ambFish;
	foreach my $read (keys %ambiguousReadsToFishInSample) {
		#print "$read\n";
		my @tmpFish = @{$ambiguousReadsToFishInSample{$read}};
		#print "$read\t@tmpFish\n";
		push(@ambFish,@tmpFish);
		@ambFish = uniq(@ambFish);
	}
	foreach my $fish (@ambFish) {
			if (!(exists $fish_ambiguousToReadsInSample{$fish})) {my @t = ();$fish_ambiguousToReadsInSample{$fish} = \@t;};
			my @arrfish = split("_",$fish);
			my $fishName = join("_",($arrfish[2],$arrfish[3]));
			my $isGulf = 0; $isGulf = 1 if  exists ($fishInGulf{lc($fishName)});
			my $isRed = 0;$isRed = 1 if exists ($fishInRed{lc($fishName)});
			print $OUT "AMB,AMB,$fish,$fishToReadsInSample{$fish},$fishToPIDsInSample{$fish},".scalar(@{$fish_ambiguousToReadsInSample{$fish}}).",$isGulf,$isRed\n";
			foreach my $read (@{$fish_ambiguousToReadsInSample{$fish}}){
					if (DEBUG_AMB) {
					delete $ambiguousReadsToFishInSample{$read};
					}
			}
	}
}
# Parse text files that were used as initial ref for being RS / Gulf sp. The final data was curated manually later in the project
sub readRedSeaAndGulfData {
	my ($GulfPath,$RedSeaPath)  = @_;
	my (%fishInGulf,%fishInRed);
	open(G,"$GulfPath") || die "Error while opening $GulfPath\n";
	while (my $genus_species =<G>) {
		chomp $genus_species;
		$fishInGulf{lc($genus_species)} = 1;
	}
	close(G);
	open(R,$RedSeaPath) || die "Eror while opening $RedSeaPath\n";
	while (my $fish = <R>) {
		chomp $fish;
		my ($fishNumber,$genus_species) = split(/\t/,$fish);
		$fishInRed{lc($genus_species)} = 1;
	}
	close(R);
	return (\%fishInGulf,\%fishInRed);
}
# Array.utils import is broken, adding required methods
sub intersect(\@\@) {
	my %e = map { $_ => undef } @{$_[0]};
	return grep { exists( $e{$_} ) } @{$_[1]};
}
sub array_minus(\@\@) {
	my %e = map{ $_ => undef } @{$_[1]};
	return grep( ! exists( $e{$_} ), @{$_[0]} ); 
}

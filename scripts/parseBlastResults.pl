use strict;
use String::Util 'trim';
use List::MoreUtils qw(uniq);
use constant false => 0;
use constant true  => 1;
use constant LENGTH_TH  => 90;
use constant PID_TH  => 98;
use constant SHORT_PID_TH  => 98;
use constant MIN_LEN_TH  => 50;
use constant START_COI  => 10;
use constant END_COI => 10;
use constant SHORT_HIT => -1;
use constant WEAK_HIT => 0;
use constant STRONG_HIT => 1;
use constant EDGE_HIT => 2;
use constant DIFF_LEN_FACTOR => 1.02;
sub printMappedReads;
sub isStrongHit;
sub upadteReadsHash;
sub updateTotalReads;
sub updateEdgeReads;
my ($COImappingFile,$coiLengthFile,$blastFile,$outputFile) = @ARGV;
@ARGV == 4 or die "Usage:<path to the COI mapping file>\t<path to the COI length file>\t<path to blast results>\t<output path>\n";
my (%hCOI_Order,%hCOI_Family,%hCOI_Genus,%hCOI_Species,%hCOI_taxa,%hCOI_Len,%hReadToBestHits,%hTotalReads,%hShortReads,%hReadToLen,%hReadToEdgeLen,%hEdgeReads,%hEdgeReadsToBestHits);

#read COI mapping file
open(COIMAP,$COImappingFile) or die "Error while opening $COImappingFile\n";
my $line = <COIMAP>; #Ignore header;
while ($line = <COIMAP>) {
	chomp $line;
	my ($COI,$order,$family,$genus,$species) = split (/\t/,$line);
	#make sure that there are no white space issues;
	$COI = trim($COI);
	$order = trim($order);
	$family = trim($family);
	$genus = trim($genus);
	$species = trim($species);
	#Create mapping
	$hCOI_Order{$COI} = $order;
	$hCOI_Family{$COI} = $family;
	$hCOI_Genus{$COI} = $genus;
	$hCOI_Species{$COI} = lc($species);
	$hCOI_taxa{$COI} = $order.$family.$genus.$species;
}
close (COIMAP);
open(COILEN,$coiLengthFile) or die "Error while opening $coiLengthFile\n";
while ($line = <COILEN>) {
	chomp $line;
	my ($COI,$len) = split (/\t/,$line);
	#make sure that there are no white space issues;
	$COI = trim($COI);
	$len = trim($len);
	$hCOI_Len{$COI} = $len;
}
close (COILEN);
open(BLAST,$blastFile) or die "Error while opening $blastFile\n";
my ($read,$pid,$len,$start,$end,$isLongRead) = ("-1","-1","-1","-1","-1",false);
#Read BLAST results and update the hashes
while(my $hit = <BLAST>) {
	chomp $hit;
	my @arrLine = split (/\t/,$hit);
	my ($cur_read,$cur_COI,$cur_pid,$cur_len,$cur_start,$cur_end) = ($arrLine[0],$arrLine[1],$arrLine[2],$arrLine[3],$arrLine[8],$arrLine[9]);
	next if ($cur_len < MIN_LEN_TH); #Ignore very short hits.
	#First, Add the read to the total reads hash.
	%hTotalReads = %{updateTotalReads($cur_read, $cur_len,$hit,\%hTotalReads)};
	my $strongHitCode = isStrongHit($cur_pid,$cur_len,$cur_start,$cur_end,$cur_COI,\%hCOI_Len,$cur_read);
	if ($strongHitCode == SHORT_HIT) {$hShortReads{$cur_read} = 1;next;} #the is unmapped due to short hit, we will not consider it as unmapped afterwards
	if ($strongHitCode == EDGE_HIT) {
		$hEdgeReads{$cur_read} = 1;
		my @refs = upadteReadsHash($hit,$cur_read,$cur_pid,$cur_len,$cur_start,$cur_end,\%hEdgeReadsToBestHits,\%hReadToEdgeLen,\%hCOI_taxa,$cur_COI);
		%hEdgeReadsToBestHits = %{$refs[0]};
		%hReadToEdgeLen = %{$refs[1]};
		next;
	} 
	next if($strongHitCode == WEAK_HIT);#the is unmapped due to short hit, we will not consider it as unmappped afterwards
	my @refs = upadteReadsHash($hit,$cur_read,$cur_pid,$cur_len,$cur_start,$cur_end,\%hReadToBestHits,\%hReadToLen,\%hCOI_taxa,$cur_COI);
	%hReadToBestHits = %{$refs[0]};
	%hReadToLen = %{$refs[1]};
}

my %uniqlyMappedCOIs;
my %coisWithLongHit;
#Find the cois with unique reads:
foreach my $key_read (keys(%hReadToBestHits)){
	my @hits = @{$hReadToBestHits{$key_read}};
	if (scalar (@hits) == 1) {
		my @arrLine = split (/\t/,$hits[0]);
		my ($cur_read,$cur_COI,$cur_pid,$cur_len,$cur_start,$cur_end) = ($arrLine[0],$arrLine[1],$arrLine[2],$arrLine[3],$arrLine[8],$arrLine[9]);
		if ($cur_len > LENGTH_TH) { #not short hit..
			$uniqlyMappedCOIs{$cur_COI} = 1;
		}
	}
	else{
		my @taxa;
		my @tmpCOIs;
		for my $h (@hits) {
			my @arrLine = split (/\t/,$h);
			my ($cur_read,$cur_COI,$cur_pid,$cur_len,$cur_start,$cur_end) = ($arrLine[0],$arrLine[1],$arrLine[2],$arrLine[3],$arrLine[8],$arrLine[9]);
			if ($cur_len > LENGTH_TH) { #not short hit..
				my $taxon = $hCOI_taxa{$cur_COI};
				push(@taxa,$taxon);
				push(@tmpCOIs,$cur_COI);
			}
		}
		if (scalar (uniq(@taxa)) == 1 ) {
			for my $c (@tmpCOIs) {
				$uniqlyMappedCOIs{$c} = 1;
			}
		}
	}
}
#Iterate again over the best hits and remove ambiguity
foreach my $key_read (keys(%hReadToBestHits)){
	my @hits = @{$hReadToBestHits{$key_read}};
	if (scalar (@hits) > 1) {
		my @hitsToUniqueCOI;
		foreach my $hit (@hits) {
			my @arrHit = split(/\t/,$hit);
			my $COI = $arrHit[1];
			#Found unique COI, remove the other reads
			if (exists $uniqlyMappedCOIs{$COI}) {
				push (@hitsToUniqueCOI,$hit);
			}
		}
		#Keep only hits to COIs from a coi that have support from unique reads.
		if (scalar @hitsToUniqueCOI > 0) {
			$hReadToBestHits{$key_read} = \@hitsToUniqueCOI;
		}
	}
}
#print results - full hits
my %taxonToReport; #Mark the reported taxon, in order to report their edge reads later.
open(OUT, ">$outputFile") or die "Error while opening $outputFile\n";
print OUT "Read\tCOI\tpid\tOrder\tFamily\tGenus\tSpecies\n";
foreach my $key_read (keys(%hReadToBestHits)){
	my @hits = @{$hReadToBestHits{$key_read}};
	my (@cois,@order,@family,@genus,@species,@ids,@sp_ids);
	my($COI,$pid);
	foreach my $hit (@hits) {
		my @arrHit = split(/\t/,$hit);
		($COI,$pid,$len) = ($arrHit[1],$arrHit[2],$arrHit[3]);
		my $t = lc($hCOI_taxa{$COI});
		$taxonToReport{$t} = 1;
		if ($len <= LENGTH_TH) { #print it only if it the coi is supported by long hit
			next if (!(exists $coisWithLongHit{$COI}));
		}
		push (@cois,$COI);
		push(@order,$hCOI_Order{$COI});
		push(@family,$hCOI_Family{$COI});
		push(@genus,$hCOI_Genus{$COI});
		push(@species,$hCOI_Species{$COI});
		push(@ids,"$COI\t$pid\t$hCOI_Order{$COI}\t$hCOI_Family{$COI}\t$hCOI_Genus{$COI}\t$hCOI_Species{$COI}");
		push (@sp_ids, "$hCOI_Order{$COI}\t$hCOI_Family{$COI}\t$hCOI_Genus{$COI}\t$hCOI_Species{$COI}"); 
		@sp_ids = uniq(@sp_ids);
		@ids = uniq(@ids);
	}
	if (scalar (@ids) == 1) {	#unique mapping
		print OUT "$key_read\t$ids[0]\n";
	}
	#ambiguous mappings print all options and add # to mark it
	# when the ambiguity is in the coi only and not in the mapping to species print it with *
	else{
		if (scalar @sp_ids == 1){
			foreach my $id (@ids){
				print OUT "*$key_read\t$id\n";
			}
		}
		else {
			foreach my $id (@ids){
				print OUT "#$key_read\t$id\n";
			}
		}
	}
}
#Print the edge reads which are supported by full hits
foreach my $key_read (keys(%hEdgeReadsToBestHits)){
	next if exists $hReadToBestHits{$key_read}; #This edge read was found also full read
	my @hits = @{$hEdgeReadsToBestHits{$key_read}};
	my (@cois,@order,@family,@genus,@species,@ids,@sp_ids);
	my($COI,$pid);
	foreach my $hit (@hits) {
		my @arrHit = split(/\t/,$hit);
		($COI,$pid,$len) = ($arrHit[1],$arrHit[2],$arrHit[3]);
		my $t = lc($hCOI_taxa{$COI});
		next if (!(exists($taxonToReport{$t})));
		push (@cois,$COI);
		push(@order,$hCOI_Order{$COI});
		push(@family,$hCOI_Family{$COI});
		push(@genus,$hCOI_Genus{$COI});
		push(@species,$hCOI_Species{$COI});
		push(@ids,"$COI\t$pid\t$hCOI_Order{$COI}\t$hCOI_Family{$COI}\t$hCOI_Genus{$COI}\t$hCOI_Species{$COI}");
		push (@sp_ids, "$hCOI_Order{$COI}\t$hCOI_Family{$COI}\t$hCOI_Genus{$COI}\t$hCOI_Species{$COI}"); 
		@sp_ids = uniq(@sp_ids);
		@ids = uniq(@ids);
	}
	if (scalar (@ids) == 1) {	#unique mapping
		print OUT "$key_read\t$ids[0]\n";
	}
	#ambiguous mappings print all options and add # to mark it
	# when the ambiguity is in the coi only and not in the mapping to species print it with *
	else{
		if (scalar @sp_ids == 1){
			foreach my $id (@ids){
				print OUT "*$key_read\t$id\n";
			}
		}
		else {
			foreach my $id (@ids){
				print OUT "#$key_read\t$id\n";
			}
		}
	}
}
close(OUT);
open(OUTUNM, ">$outputFile.unmapped") or die "Error while opening $outputFile\n";
print OUTUNM "Read\tCOI\tpid\tOrder\tFamily\tGenus\tSpecies\n";
#print the reads with unknown COI
foreach my $read (keys %hTotalReads) {
	next if exists $hReadToBestHits{$read};
	next if exists $hShortReads{$read};
	next if exists $hEdgeReads{$read};
	my @arrHit = split(/\t/,$hTotalReads{$read});
	my($COI,$pid) = ($arrHit[1],$arrHit[2]);
	print OUTUNM "$read\t$COI\t$pid\tUNK\tUNK\tUNK\tUNK\n";
}
close(OUTUNM);
#print "Total mapped: " . keys (%hReadToBestHits) . "Total reads = " .keys (%hTotalReads) . "\n";
sub isStrongHit{
	my ($cur_pid,$cur_len,$cur_start,$cur_end,$cur_coi,$hCOIlenRef,$cur_read) = @_;
	my %hCOIlen = %{$hCOIlenRef};
	my $coiLenTH = $hCOIlen{$cur_coi} - END_COI;	
	if ($cur_pid < PID_TH) {
		return WEAK_HIT;
	}
	if ($cur_len > LENGTH_TH) {
		return STRONG_HIT;
	}
	if ($cur_start < START_COI || $cur_end < START_COI || $cur_start > $coiLenTH || $cur_end > $coiLenTH) { #good short hit
		if ($cur_pid > SHORT_PID_TH) {
			return EDGE_HIT;
		}
		else {
			return SHORT_HIT;
		}
	}
	return WEAK_HIT;
}
sub updateTotalReads{
	my ($cur_read, $cur_len,$hit,$hReadToBestHitsRef) = @_;
	my %hReadToBestHits = %{$hReadToBestHitsRef};
	if (exists ($hReadToBestHits{$cur_read})){
		return(\%hReadToBestHits);
	}	
	$hReadToBestHits{$cur_read} = $hit;
	return(\%hReadToBestHits);
}
sub upadteReadsHash{
	my ($cur_hit,$cur_read,$cur_pid,$cur_len,$cur_start,$cur_end,$hReadToBestHitsRef,$hReadToLenRef,$hCOI_taxaRef,$cur_COI) = @_;
	my %hReadToBestHits = %{$hReadToBestHitsRef};
	my %hReadToLen = %{$hReadToLenRef};
	my %hCOI_taxa = %{$hCOI_taxaRef};
	if (!exists $hReadToBestHits{$cur_read}) { #first time this read is mapped
		my @hits = ();
		push (@hits, $cur_hit);
		$hReadToBestHits{$cur_read} = \@hits;
		my @arrLine = split (/\t/,$cur_hit);
		my ($best_COI,$best_pid,$best_len,$best_start,$best_end) = ($arrLine[1],$arrLine[2],$arrLine[3],$arrLine[8],$arrLine[9]);
		$hReadToLen{$cur_read} = $best_len;
	}
	else {
		my $best_hit = (@{$hReadToBestHits{$cur_read}})[0];
		my @arrLine = split (/\t/,$best_hit);
		my ($best_COI,$best_pid,$best_len,$best_start,$best_end) = ($arrLine[1],$arrLine[2],$arrLine[3],$arrLine[8],$arrLine[9]);
		$best_len = $hReadToLen{$cur_read};
		if (($cur_len > LENGTH_TH) && ($best_len < LENGTH_TH)) { #full hit instead of short hit, replace;
			my @hits = ();
			push (@hits, $cur_hit);
			$hReadToBestHits{$cur_read} = \@hits;
			$hReadToLen{$cur_read} = $cur_len;
		}
		next if (($cur_len < LENGTH_TH) && ($best_len > LENGTH_TH)); #current read is short and the previous was long -> ignore short read
		#Both reads are either short or long, update if the pids are equal, break equality if there if the alignment length factor if large (SHORT_LEN_FACTOR)
		if ((($cur_len < LENGTH_TH) && ($best_len < LENGTH_TH)) ||(($cur_len > LENGTH_TH) && ($best_len > LENGTH_TH))){
		next if $best_pid > $cur_pid; #better PID was observed
			if ($cur_len > $best_len) {
				if ($cur_len > ($best_len * DIFF_LEN_FACTOR)) {
					my @hits = ();
					push (@hits, $cur_hit);
					$hReadToBestHits{$cur_read} = \@hits;
					$hReadToLen{$cur_read} = $cur_len;
				}
				else {
					my @hits = @{$hReadToBestHits{$cur_read}};
					push(@hits,$cur_hit);
					$hReadToBestHits{$cur_read} = \@hits;
				}
			}
			else {
				next if $best_len > ($cur_len * DIFF_LEN_FACTOR);
				my @hits = @{$hReadToBestHits{$cur_read}};
				push(@hits,$cur_hit);
				$hReadToBestHits{$cur_read} = \@hits;
			}
		}
	}	
	return(\%hReadToBestHits,\%hReadToLen);
}

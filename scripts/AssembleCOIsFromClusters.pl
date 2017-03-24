#!/usr/bin/perl -w
my $clusters = "clusters_of_reads_mapping_with_100id_at_least_80bp";
my $reads = "all_unmapped_reads.fasta";
my $dir = "denovo_recreate/";
my $dbPath = "COI_db.fa";
open (READS, $reads) || die "cannot open reads file";
while (<READS>){
    if (/^>(\S+)\s+(\d):N/){
	my ($name,$mate) = ($1,$2);
	my $seq = <READS>;
	
	$mate = "R".$mate;

	$data{$name}{$mate} = $seq;

    }
}

open (CL,$clusters) || die "cannot open cluster file";

$clust_num=0;
while (<CL>){
    next if /clustering/;
    $clust_num++;
    chomp;
    my @reads = split;
    
    # taking only clusters with 10 reads or more
    last if (scalar (@reads) <10);
    
    # making a directory for each cluster
    my $dir_name = $dir."cluster_".$clust_num;
    unless (mkdir $dir_name){
	warn "cannot make directory $dir_name";
    }

    my %used_reads;
    my $r1 = $dir_name."/R1";
    my $r3 = $dir_name."/R3";
    open (R1,">$r1") || die "cannot open $r1 for writting";
    open (R3,">$r3") || die "cannot open $r3 for writting";

    foreach my $read (@reads){
	next if (exists $used_reads{$read});
	$used_reads{$read} = 1;
	
	print R1 ">",$read,"\n";
	print R1 $data{$read}{R1};
	print R3 ">",$read,"\n";
	print R3 $data{$read}{R3};
    }
    close R1;
    close R3;
    
    # assembling 
    my $command1 = "velveth $dir_name 25 -fasta -shortPaired -separate $r1 $r3";
    my $command2 = "velvetg $dir_name -min_contig_lgth 100";

    system $command1;
    system $command2;

    # finding which part of each contig maps to COI
    my $contigs_file = $dir_name."/contigs.fa";
    my $blast_out = $contigs_file.".round1_blast_COI";

    my $command3 = "blastall -i $contigs_file -d $dbPath -F F -b 1 -v 1 -e 0.0001 -m 8 -p blastn -o  $blast_out";

    system $command3;
    
    open (CONTIGS,$contigs_file) || die "cannot open contigs file $contigs";
    while (<CONTIGS>){
	if (/>(\S+)/){
	    if ($flag){
		$contigs{$con_name} = $con_seq;
		undef $con_seq;
	    }
	    $con_name = $1;
	    $flag = 1;
	}
	else {
	    chomp;
	    $con_seq .= $_;
	}
    }
    $contigs{$con_name} = $con_seq;

    close CONTIGS;

    # trimming contigs to fit just the COI
	    
    my $new_contigs_file = $contigs_file."._round1_trimmed_for_COI";
    open (NEW_CONTIGS,">$new_contigs_file") || die "cannot open new contigs file $new_contigs_file";
    
    open (BLAST, $blast_out) || die "cannot open file $blast_out";
    while (<BLAST>){
	my @line = split;
	my ($curr_contig,$fr,$to) = ($line[0],$line[6],$line[7]);
	
	my $new_contig_seq = substr ($contigs{$curr_contig},$fr-1,$to-$fr+1);
	print NEW_CONTIGS ">",$curr_contig,"_$fr","_$to COI_fr $line[8] COI_to $line[9]\n",$new_contig_seq,"\n";	
    }
    close BLAST;
    close NEW_CONTIGS;

    # looking for all the reads that map to the trimmed contigs

    my $command4 = "formatdb -i $new_contigs_file -p F";
    system $command4;

    my $second_blast_out = $new_contigs_file.".blast_all_unmapped_reads";
    my $command5 = "blastall -i $reads -d $new_contigs_file -p blastn -e 0.00001 -F F -m 8 -o $second_blast_out";
    system $command5;
    
    my %used_reads2;
    my $r1_second = $dir_name."/R1_second_round";
    my $r3_second = $dir_name."/R3_second_round";
    open (R1_2,">$r1_second") || die "cannot open $r1_second for writting";
    open (R3_2,">$r3_second") || die "cannot open $r3_second for writting";

    open (SECOND_BLAST,$second_blast_out) || die "cannot open second blast out file $second_blast_out";
    while (<SECOND_BLAST>){
	my @line = split;
	if ($line[2] eq "100.00" && $line[3]>=80){
	    my $read = $line[0];
	    next if (exists $used_reads2{$read});
	    $used_reads2{$read} = 1;
	    print R1_2 ">",$read,"\n";
	    print R1_2 $data{$read}{R1};
	    print R3_2 ">",$read,"\n";
	    print R3_2 $data{$read}{R3};	    
	}
    }

    close R1_2;    
    close R3_2;

    my $command6 = "velveth $dir_name 25 -fasta -shortPaired -separate $r1_second $r3_second";
    my $command7 = "velvetg $dir_name -min_contig_lgth 100";

    system $command6;
    system $command7;

#    ====== second round of trimming ====

    # finding which part of each contig maps to COI
    my $contigs_file = $dir_name."/contigs.fa";
    $blast_out = $contigs_file.".round2_blast_COI";

    my $command3 = "blastall -i $contigs_file -d $dbPath -F F -b 1 -v 1 -e 0.0001 -m 8 -p blastn -o  $blast_out";

    system $command3;
    
    open (CONTIGS,$contigs_file) || die "cannot open contigs file $contigs";
    while (<CONTIGS>){
	if (/>(\S+)/){
	    if ($flag2){
		$contigs2{$con_name} = $con_seq;
		undef $con_seq;
	    }
	    $con_name = $1;
	    $flag2 = 1;
	}
	else {
	    chomp;
	    $con_seq .= $_;
	}
    }
    $contigs2{$con_name} = $con_seq;
    close CONTIGS;

    # trimming contigs to fit just the COI
	    
    my $new_contigs_file2 = $contigs_file."._round2_trimmed_for_COI";
    open (NEW_CONTIGS,">$new_contigs_file2") || die "cannot open new contigs file $new_contigs_file2";
    
    open (BLAST, $blast_out) || die "cannot open file $blast_out";
    while (<BLAST>){
	my @line = split;
	my ($curr_contig,$fr,$to) = ($line[0],$line[6],$line[7]);
	
	my $new_contig_seq = substr ($contigs2{$curr_contig},$fr-1,$to-$fr+1);
	print NEW_CONTIGS ">",$curr_contig,"_$fr","_$to COI_fr $line[8] COI_to $line[9]\n",$new_contig_seq,"\n";	
    }
    close BLAST;
    close NEW_CONTIGS;



 #  exit (0);

}
    
    
    
    
    
    

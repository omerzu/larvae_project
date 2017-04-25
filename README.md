# larvae_project
This is the source code for analyzing larvae sequncing data.
The suggested structrue is to have a directory for each sequnced sample.

## preparations
1. Build the local COI db (formatdb -i db.fa -p F)
2. Create a text file with a line for each directory that contain a sequenced sample (see dirs file for example)  

## Running the pipeline
1. Create the blast calls:
 `perl scripts/createBlastCalls.pl > callBlast.sh` 
2. Parse Blast results:
 `perl scripts/createParseBlastCalls.pl >ParseBlastRes.sh` and execute the newly created files:
   *  `bash sortBlastResults.sh`
   *  `bash ParseBlastRes.sh`
3. Building initial sample summay:
   * `ls */*.pipe.blast.parsed.final > pipe.blast.parsed.files`
   * count total number of reads (using grep -c for example) in each sample and save in to all.reads.cnt file
   * `perl scripts/createInitialReport.pl pipe.blast.parsed.files files/initial_morph_data.txt dirs files/sampleNameToLarvaeNumber.txt files/GOA.txt  files/RS.txt` 
   *  **NOTE:** This script was initially used to summarize the parsed data for easy comparison with the initial morphological data.
    The script also added initial annotation of Red Sea and Gulf species (using RS.txt and GOA.txt respectively), the species list were manually curated later so those files may contain partial and inaccurate data.
 4. combine the files from the previous step: `for f in $(cat dirs);do cat $f/$f*.out >>pipe.initial.report.csv;done;`
 5. Calculate coverage for each identified COI: `scripts/simple_coverage_calculator.py --out coverage_data.txt`
6. Add coverage data to the report:
	 `perl scripts/addCoverageData.pl coverage_data.txt pipe.initial.report.csv >pipe.report.csv`
7. Filter out identification with less than 50% COI coverage:
	`perl scripts/filterValidIdentifications.pl pipe.report.csv > All_fish_identified_in_all_samples_with_COI_coverage_of_at_least_50_percent.csv`

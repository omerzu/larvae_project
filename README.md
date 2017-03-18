# larvae_project
This is the source code for analyzing larvae sequncing data.
The suggested structrue is to have a directory for each sequnced sample.

## Running the pipeline
1. Build the local COI db (formatdb -i db.fa -p F)
2. Create the blast calls:
 `perl scripts/createBlastCalls.pl > callBlast.sh` 
3. Parse Blast results:
 `perl scripts/createParseBlastCalls.pl >ParseBlastRes.sh` and execute the newly created files:
   *  `bash sortBlastResults.sh`
   *  `bash ParseBlastRes.sh`

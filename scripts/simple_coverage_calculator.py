"""
This Code replaces usage of another pipeline that was also detected SNPs between the referrnce COI to sequncing data.
For backward-compatibility it outputs few NA columns
"""
from Bio.SeqIO import parse
from optparse import OptionParser
def Main():
    parser = OptionParser()
    parser.add_option("--db_file",
                      default="db.fa",
                      help="COI DB to calculate coverage")
    parser.add_option("--blast",
                      default="coi.reads.blastn.COI.final",
                      help="suffix of blast res file")
    parser.add_option("--hits",
                      default="pipe.blast.parsed.final",
                      help="Suffix of the mapping file between the reads to the COI from previous steps")
    parser.add_option("--samples_file",
                      default="dirs",
                      help="A file with the names of all the direcetories")
    parser.add_option("-o","--out",
                      help="Output file")
    options, _ = parser.parse_args()
    with open(options.samples_file)  as f:
        samples = f.readlines()
    with open(options.db_file) as f:
        cois = list(parse(f, "fasta"))

    out_file = open(options.out, 'w')
    out_file.write("Sample\tTaxonomy\tCOI\tCoveredBPs\tCOI_Length\tcoverage\tNA\tNA\n")
    for sample in samples:
        sample = sample.rstrip("\n")
        coiToCoverdBases, coiTax = SampleCoverage("{}/{}.{}".format(sample,sample, options.hits), "{}/{}.{}".format(sample,sample, options.blast), cois)
        for coi, coverage in coiToCoverdBases.iteritems():
            covered_bps = len([bp for bp in coverage if bp])
            if covered_bps:
                out_file.write("{}\t{}\t{}\t{}\t{}\t{}\tNA\tNA\n".format(sample,coiTax[coi], coi, covered_bps, len(coverage), round(float(covered_bps)/len(coverage)*100,2)))

def SampleCoverage(hits_file, blast_file, cois):
    coiToCoverdBases = {}
    coiToTaxonomy = {}
    for coi in cois:
        coiToCoverdBases[coi.id] = [False] * len(coi.seq)
    with open(hits_file) as f:
        hits = f.readlines()
    readsToCOI = {}
    # Parse the hits file, the first raw is header
    for hit in hits[1:]:
        read, COI, _, order, family, genus, species = hit.split("\t")
        # hits start with '*' are reads that mapped to difference COI of the same species so it OK to count it
        read = read.replace("*", "", 1)
        readsToCOI[read] = COI
        taxon = "{}_{}_{}_{}".format(order, family, genus, species)
        coiToTaxonomy[COI] = taxon.rstrip("\n")
    with open(blast_file) as f:
        blast_lines = f.readlines()
    for blast_line in blast_lines:
        query, subj, pid, aln_len, mismathces, gaps, qs, qe, ss, se, _, _ = blast_line.split("\t")
        ss = int(ss)
        se = int(se)
        matching_COI = readsToCOI.get(query)
        if not matching_COI or matching_COI != subj:
            continue
        covered_indices = xrange(ss - 1, se) if ss < se else xrange(se - 1, ss, 1)
        for i in covered_indices:
            coiToCoverdBases[matching_COI][i] = True
    return coiToCoverdBases, coiToTaxonomy

if __name__ == "__main__":
    Main()

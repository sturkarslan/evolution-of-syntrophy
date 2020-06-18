# Collect stats from all bam files and write them into a file
import glob, sys, os, string, os.path
from os import path

picardPath = "/Users/sturkars/Applications/picard.jar"
epds = ['results_09_2']
organisms = ['dvh']
for epd in epds:
    for organism in organisms:

        if organism == 'dvh':
            resultsDir = '/Volumes/omics4tb/sturkarslan/dvh-coculture-rnaseq/dvh-single-cells/%s/%s'%(epd,organism)
            referenceFasta = '/Volumes/omics4tb/sturkarslan/dvh-coculture-rnaseq/dvh-single-cells/reference/Desulfovibrio_vulgaris_str_hildenborough.GCA_000195755.1.30.dna.genome.fasta'
        elif organism == 'mmp':
            resultsDir = '/Volumes/omics4tb/sturkarslan/singleCell-UA3-mmp/%s/%s'%(epd,organism)
            referenceFasta = '/Volumes/omics4tb/sturkarslan/dvh-coculture-rnaseq/dvh-single-cells/reference/Methanococcus_maripaludis_s2.GCA_000011585.1.30.dna.genome.fasta'

        outfile = '%s/collected_alignment_summary_metrics.txt'%resultsDir
        print('Processing: EPD: %s, Organism: %s Outfile: %s\n' %(epd,organism,outfile))

        bamFiles = glob.glob('%s/*/alignment_results/*_marked.bam' %(resultsDir))

        # Open output file
        g = open(outfile, 'w')

        # create and write header
        header = 'Organism\tSample\tCATEGORY\tTOTAL_READS\tPF_READS\tPCT_PF_READS\tPF_NOISE_READS\tPF_READS_ALIGNED\tPCT_PF_READS_ALIGNED\tPF_ALIGNED_BASES\tPF_HQ_ALIGNED_READS\tPF_HQ_ALIGNED_BASES\tPF_HQ_ALIGNED_Q20_BASES\tPF_HQ_MEDIAN_MISMATCHES\tPF_MISMATCH_RATE\tPF_HQ_ERROR_RATE\tPF_INDEL_RATE\tMEAN_READ_LENGTH\tREADS_ALIGNED_IN_PAIRS\tPCT_READS_ALIGNED_IN_PAIRS\tPF_READS_IMPROPER_PAIRS\tPCT_PF_READS_IMPROPER_PAIRS\tBAD_CYCLES\tSTRAND_BALANCE\tPCT_CHIMERAS\tPCT_ADAPTER\tSAMPLE\tLIBRARY\tREAD_GROUP\n'
        g.write(header)

        # process each BAM file with Picard CollectAlignmentSummaryMetric and write to output and then read the output
        for file in bamFiles:
            #print('File is: %s' %file)
            sampleName = file.split('_marked.bam')[-2].split('/')[-1].split('-')[0]
            picardSummaryOutput = file.split('_marked.bam')[0] + '_picardSummaryMetrics.txt'
            print('sample Name is: %s' %(sampleName))
            #print('PicardOutput File is: %s' %(picardSummaryOutput))
        
            # run picard    
            picard_cmd = 'java -jar %s CollectAlignmentSummaryMetrics R=%s I=%s O=%s' %(picardPath, referenceFasta, file, picardSummaryOutput)
            #print('Running Picard with command:\n %s' %picard_cmd)
            os.system(picard_cmd)    
            
            # open picard output to read the stats if exists
            try:   
                f = open(picardSummaryOutput, 'r')
                # Skip first lines to get the pair data
                line = f.readlines()[9]
                #print(line)
                vector = line.split('\t')
                #print(vector)
                vectorJoined = '%s\t%s\t'%(organism,sampleName) + '\t'.join(vector)
                #print(vectorJoined)
            except:
                vectorJoined = '\t'.join('NA')

            #g.write('%s\t'%organism)
            #g.write('%s\t'%sampleName)    
            g.write(vectorJoined)
            #g.write("\n")
            print    
            #sys.exit()

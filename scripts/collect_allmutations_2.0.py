
############# Serdar Turkarslan / Institute for Systems Biology ###############
# Last update: 06/05/2020
###############################################################################
# Collect all mutations in EPDs, clonal isolate, 1000-gen, clonal isolates
# and ancestors into a single attributes, collated and matrix files
###############################################################################
from __future__ import division
import glob, sys
import datetime
#import pandas as pd

now = datetime.datetime.now()
runDate = '%s%s%s' %(now.month,now.day,now.year)

runDate = "04172019"
outfile = '/proj/omics4tb/sturkarslan/dvh-mutation-verifications/Dvh_mutations_allsamples_attributes_' + runDate + '.txt'
outfile2 = '/proj/omics4tb/sturkarslan/dvh-mutation-verifications/Dvh_mutations_allsamples_collated_' + runDate + '.txt'
organism = 'dvh'
# paths for all sample folders
paths = [#"/proj/omics4tb/sturkarslan/EPD/evolved_lines/after_300g/results/dvh/*/",
         "/proj/omics4tb/sturkarslan/EPD/EPD_seq/results/%s/*/" %organism,
         "/proj/omics4tb/sturkarslan/clonal-isolates/results/%s/*/" %organism,
         "/proj/omics4tb/sturkarslan/syntrophy_raw_sequences/1000-gen/results/%s/*/" %organism,
         "/proj/omics4tb/sturkarslan/syntrophy_raw_sequences/Ancestors/results/%s/*/" %organism,
         "/proj/omics4tb/sturkarslan/syntrophy_raw_sequences/Early-gen/results/%s/*/" %organism]

# create a list of all folders from these paths
folders = []
for path in paths:
    folder = glob.glob(path)
    folders.append(folder)

## Remove unwanted folders from the list
exception = ["AK_43", "AK_44", "AK_47", "AK_48", "AK_49"]
exceptionfull = ["/proj/omics4tb/sturkarslan/clonal-isolates/results/%s/"%organism + i + "/" for i in exception ]
# create a flat list out of nested lists
folderlistfull = [item for sublist in folders for item in sublist]
folderlist = iter([j for j in folderlistfull if j not in exceptionfull])
# file to get list of all variants in all single cells
# count files
print("Procesing counts files...\n")
#print(folderlistfull)
#sys.exit()
# function to calculate means for frequencies
def mean(numbers):
    return round(float(sum(numbers)) / max(len(numbers), 1), 2)

h = open(outfile, 'w')
headernames = ["variant_id","experiment","sample","source","position","initial_nt_1st","initial_nt","changed_nt","mutation","effect","type","codon_change","aa_change","gene_id","region","freq","predictor","freq_predictor","read_number","\n"]
headerline = "\t".join(headernames)
h.write(headerline)
samplenames = []
#countfiles = []
for folder in folderlist:
    # remove EPD HA2, HR2 and UA3 lines that are actuall 1000-gen samples
    folderExceptions = ['HA2','HR2','UA3','CI_36','CI_37','old_results_label_swap']
    excludedFolders = []
    for exception in folderExceptions:
        print(exception)
        if exception == ('CI_36' or 'CI_37'):
            excludedFolders.append("/proj/omics4tb/sturkarslan/clonal-isolates/results/%s/%s/" %(organism, exception))
        elif exception == ('HA2' or 'HR2' or'UA3'):
            excludedFolders.append("/proj/omics4tb/sturkarslan/EPD/EPD_seq/results/%s/%s/" %(organism,exception))
        else:
             excludedFolders.append("/proj/omics4tb/sturkarslan/Early-gen/results/%s/%s/" %(organism,exception))
    print
    print('exceludedFolders:%s' %excludedFolders)
    print
    sys.exit()
    if folder in excludedFolders:
        print('Found exception folder: %s' %folder)
        next(folderlist)
        continue
    #folder = '/proj/omics4tb/sturkarslan/syntrophy_raw_sequences/Early-gen/results/dvh/HA2_15/'
    folder_type = folder.split('/')[-5]
    print
    print('folder:%s' %folder)
    print('folder_type:%s' %folder_type)
    #print('%scombined_variants/*_merged_variants_final.txt'%folder)

    if folder_type == 'Early-gen':
        variantfile = glob.glob( ('%scombined_variants/*_merged_variants_final.txt'%folder))
        samplename = variantfile[0].split("/")[-3]
        experiment = variantfile[0].split("/")[-6]
    else:
        variantfile = glob.glob( ('%s*.consensus.variants.FINAL.txt'%folder))
        samplename = variantfile[0].split("/")[-2]
        experiment = variantfile[0].split("/")[-5]
    #sys.exit()
    #variantfile = glob.glob((folder + 'combined_variants/' + "*.merged_variants_final.txt") or (folder + "*.consensus.variants.FINAL.txt"))
    #variantfile = glob.glob( ('%scombined_variants/*_merged_variants_final.txt'%folder))
    print('variantfile:%s'%variantfile)
    print

    if experiment == "after_300g":
        samplenameCut = samplename.split("_")[0]
        samplenameClean = "EG_" + samplenameCut
    elif experiment == "EPD_seq" and samplename == "00_S7":
        samplenameClean = "AN_" + samplename + "_Stahl"
        experiment = "Ancestors"
    elif experiment == "EPD_seq":
        samplenameClean = "EP_" + samplename
    elif experiment == "clonal-isolates":
        samplenameClean = "CI_" + samplename
    elif experiment == "1000-gen":
        samplenameClean = "TG_" + samplename
    elif experiment == "Ancestors":
        samplenameClean = "AN_" + samplename
    elif experiment == "Early-gen":
        print(samplename)
        samplenameCut = samplename.split("_")[1]
        samplenameClean = "EG_" + samplename
    else:
        samplenameClean = samplename

    print('sampleName:%s, experiment:%s' %(samplenameClean, experiment))
    #sys.exit()

    f = open(variantfile[0], 'r')
    # loop through each line
    for line in f:
        fields = line.split("\t")
        chromosome = fields[0]
        chromosome = chromosome.replace("Chromosome", "Chr")
        coordinate = fields[1]
        alternative = fields[2]

        sysname = fields[10]
        if sysname == "":
            locus = "IG"
        else:
            locus = sysname.replace("DVU_", "DVU")

        programs = fields[13]

        variantname = "%s-%s-%s" %(chromosome,locus,coordinate)
        programcount = len(programs.split(":"))
        #print(programs, programcount)

        #remove first reference base position column
        #del fields[2]
        #print(fields)
        #line2 = str(("\t").join(fields))

        if programcount >= 2:
            line2write = variantname + "\t" + experiment + "\t" + samplenameClean + "\t" + line
            h.write(line2write)
            samplenames.append(samplenameClean)
h.close()
sys.exit()


## create a matrix of mutation frequencies per sample
matrixfile = '/proj/omics4tb/sturkarslan/dvh-mutation-verifications/Mmp_mutations_allsamples_matrix_' + runDate + '.txt'
f = open(outfile, 'r')
g = open(matrixfile, 'w')
next(f)
variants = []
samples = []
frequencies = []
variantsDict = {}
for line in f:
    fields = line.split("\t")
    variantid = fields[0]
    samplename = fields[2]

    frequencies = fields[17]
    frequencies = frequencies.replace("%", "")
    frequencies = list(frequencies.split(":"))
    print(frequencies)
    #convert to float
    frequencies = [float(i) for i in frequencies]
    meanFreq = mean(frequencies)
    variantsDict[(variantid, samplename)] = meanFreq
    variants.append(variantid)
    samples.append(samplename)
# get list of keys to check if sample/variant exists
keylist = list(variantsDict.keys())

# write sample names as column names
headers = list(sorted(set(samples), reverse=True))
headers2write = "Variantid" + "\t" + "\t".join(headers) + "\n"
g.write(headers2write)

for variant in list(sorted(set(variants))):
    print(variant)
    frequencylist = []
    frequencylist.append(variant)
    for sample in list(sorted(set(samples), reverse=True)):
        # check to see if key exst for given variantid sample name combination
        if (variant, sample) in keylist:
            frequency = variantsDict[(variant, sample)]
            frequencylist.append(frequency)
        else:
            print("not in keylist")
            frequencylist.append("")
    frequencylist.append("\n")
    print(frequencylist)
    # convert to string for writing
    frequencylist = [str(i) for i in frequencylist]
    line2write = "\t".join(frequencylist)
    #print(line2write)
    g.write(line2write)
g.close()
f.close()


print("################ Part 3 #################")

t = open(outfile2, 'w')

for variant in list(sorted(set(variants))):
    s = open(outfile, 'r')
    samplelist = []
    linelist = []
    experimentlist = []
    print("Working with variant: %s"  %variant)
    for line in s:
        fields = line.split("\t")
        print("This is fields: %s" %fields)
        variantid = fields[0]
        samplename = fields[2]
        experiment = fields[1]
        print("this is variant: %s" %variantid)
        if variant == variantid:
            samplelist.append(samplename)
            print(samplelist)
            experimentlist.append(experiment)
            linelist.append(line)
    #uniqueline = "\t".join(linelist[0])
    uniquesamples = ":".join(samplelist)
    line2write = variant + "\t" + uniquesamples + "\t" +  "\n"#uniqueline + "\n"
    t.write(line2write)
t.close()

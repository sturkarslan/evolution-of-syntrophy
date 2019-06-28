
############# Serdar Turkarslan / Institute for Systems Biology ###############
# Last update: 05/13/2019
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








#
## Open variant file and loop through each variant to create a variantList
#print("Procesing variants file...\n")
#f = open(variantFile, 'r')
#variantList = []
#for line in f:
#    fields = line.split("\t")
#    one = fields[0]
#    two = fields[1]
#    three = fields[3]
#    locus = fields[6]
#    names = "%s-%s-%s" %(one, two,locus)
#    variantID = "%s-%s-%s_%s" %(one, two, three,locus)
#    variantList.append(variantID)
#
#
#v = open(mutationmatrix, 'w') ## open matrix output only file for writing
#h = open(outfile, 'w') ##open output file for matrix and cel and mutation names
#headerlist = []
#headerlist.append("Variant")
#
#for cfile in countfiles:
#    if cfile[0].split("/")[6] == "dvh":
#        headername = cfile[0].split("/")[7].split("_allsamples_bamreadcount_parsed.txt")[0] + "-" + cfile[0].split("/")[4].split("_allsamples_bamreadcount_parsed.txt")[0]
#
#    if cfile[0].split("/")[7] == "dvh":
#        headername = cfile[0].split("/")[8].split("_allsamples_bamreadcount_parsed.txt")[0].split("-")[0]
#
#    if cfile[0].split("/")[6] == "after_300g":
#        headername = cfile[0].split("/")[9].split("_allsamples_bamreadcount_parsed.txt")[0].split("_")[0]
#
#    #print(headername)
#    headerlist.append(headername)
#headerlist.append("\n")
#header = "\t".join(headerlist)
#h.write(header)
#
### write variants into mutation names file
#t = open(mutationnames, 'w')
#for mutation in variantList:
#    mutationShort = mutation.replace("Chromosome", "")
#    mutationShort = mutationShort.replace("pDV", "p")
#    mutationShort = mutationShort.replace("DVU_", "DVU")
#    mutationShort = mutationShort.split("_")[1] + "." + mutationShort.split("-")[1] + mutationShort.split("-")[0]
#    t.write(mutationShort)
#    t.write("\n")
##twrite = "\n".join(variantList)
#t.close()
#
### write cell names into cell names fileName
#s = open(cellnames, 'w')
#swrite = "\t".join(headerlist)
#s.write(swrite)
#s.close()
#
#n = 1
## loop through each variant and search each count file to collect mutations status
#for variant2 in variantList:
#    totalvariants = len(variantList)
#    print("--------------------%s/%s--------------------" %(n,totalvariants))
#    print("Processing %s..." %(variant2))
#    print
#    statusList = []
#    matrixList = []
#    variant = variant2.split("_")[0]
#    variant2 = variant2.replace("\n", "")
#    statusList.append(variant2)
#
#    # loop through each single cell folder to read bmreadcount files and parse
#    m = 1
#    totalfolders = len(folderlist)
#    for folder in folderlist:
#            print("----------%s/%s----------" %(m,totalfolders))
#            print("Processing %s..." %(folder))
#            print
#            # get count file
#            bamFile = glob.glob(folder + "*_marked.bam")[0]
#            consensusVariant = glob.glob(folder + "*consensus.variants.FINAL.txt")[0]
#            sample = bamFile.split("_marked.bam")[0]
#            countFile = sample + "_allsamples_bamreadcount.txt"
#            resultFile = countFile.split("_allsamples_bamreadcount.txt")[0] + "_allsamples_bamreadcount_parsed.txt"
#
#            # create a list of variants from consesnus variant calls for each folder
#            consensusList = []
#            k = open(consensusVariant, 'r')
#            for line in k:
#                row = line.split("\t")
#                chromosome1 = row[0]
#                if chromosome1 == "NC_002937":
#                    chromosome = "Chromosome"
#                elif chromosome1 == "NC_005863":
#                    chromosome = "pDV"
#                else:
#                    chromosome = chromosome1
#                coord = row[1]
#                reference = row[2]
#                #alt = row[15]
#                joint = "%s-%s-%s" %(chromosome, coord, reference)
#                consensusList.append(joint)
#
#
#            # create a dictionary of mutations as keys and status as values
#            countDict = {}
#            g = open(resultFile, 'r')
#            for line in g:
#                row1 = line.split("\t")
#                chromosome1 = row1[0]
#                if chromosome1 == "NC_002937":
#                    chromosome = "Chromosome"
#                elif chromosome1 == "NC_005863":
#                    chromosome = "pDV"
#                else:
#                    chromosome = chromosome1
#                coord = row1[1]
#                reference = row1[2]
#                alt = row1[10]
#                status = row1[11]
#                ins = row1[7]
#                joint2 = "%s-%s-%s" %(chromosome, coord, reference)
#                rcount = str(reference.split(":")[1])
#                joint2 = joint2.split(":")[0]
#                summary = reference + "|" + alt + "|" + ins + "|" + rcount + "|" + status
#                countDict[joint2] = summary
#            keylist = list(countDict.keys())
#            #print(keylist)
#
#            # check variants if they are in consensus, if not check in bam counts if not assign 3 (NA)
#            print(variant)
#            ## mutations is called by callers
#            if variant in consensusList:
#                status = countDict[variant] + "|1"
#                print("%s in consensusList" %variant)
#
#            # mutation is not
#            elif variant in keylist:
#                print("%s in keyList" %variant)
#                verify = int(countDict[variant].split("|")[3])
#                if verify >= 5:
#                    print("%s is higher than 5" %verify)
#                    status = countDict[variant] + "|0"
#                else:
#                    status = countDict[variant] + "|3"
#                    print("%s is LOWER than 5" %verify)
#
#            # mutation can not be called
#            else:
#                status = "3"
#
#            ## Do mutation presence/absence calls agree?
#            if status == "3":
#                status = "3"
#            else:
#                verify1 = status.split("|")[4]
#                verify2 = status.split("|")[5]
#            if verify1 == verify2:
#                status = status + "|OK"
#            else:
#                status = status + "|CHECK"
#
#            m = m + 1
#            statusList.append(status)
#            matrixList.append(status)
#    #sys.exit()
#    n = n + 1
#    #print(statusList)
#    #statusList.append("\n")
#    #print(statusList[0:3])
#    line2write = "\t".join(statusList)
#    line2write = line2write + "\n"
#    #print(line2write)
#    h.write(line2write)
#    ## for matrix and names
#    vwrite = "\t".join(matrixList) + "\n" ## for matrix ony
#
#    v.write(vwrite)
#h.close()
#v.close()

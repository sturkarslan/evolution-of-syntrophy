############# Serdar Turkarslan / Institute for Systems Biology ###############
# Last update: 06/05/2020
###############################################################################
# DNASeq variant calling pipeline for read quality control (trim_galore), 
# read alignments to reference (bwa) and variant calling by GATK, Samtools and Varscan.
###############################################################################
import glob
import sys
import os
import re
#from tqdm import tqdm

############# Programs #############
fastqc = "/users/sturkars/FastQC/fastqc" # path to fastqc executable
javaPath = "/usr/bin/java"
piccardPath = "/users/sturkars/picard-tools-1.139/picard.jar"
gatkPath = "/users/sturkars/gatk/GenomeAnalysisTK.jar"
varscanPath = "/users/sturkars/VarScan.v2.3.9.jar"
snpeff_path = "/users/sturkars/snpEffv43/snpEff/"
############# Globals ##############
organism = "mmp"
# data and results directories
run_dir = "/proj/omics4tb/sturkarslan/syntrophy_raw_sequences/Early-gen"
data_dir = "%s/data" %run_dir
data_trimmed_dir = "%s/trimmed" %data_dir
genome_dir = "%s/reference" %run_dir
fastqc_dir = "%s/data/fastqc_results" %run_dir
genomeGff = '%s/%s.GCA_000195755.1.30.gtf' %(genome_dir, organism)
known_sites = '%s-variants-compiled_sorted.vcf' %(organism)

######### Annotation databases #######
# snpEff databases
if organism == "mmp":
    snpeff_db = "Methanococcus_maripaludis_s2"
    genome_fasta = '%s/Methanococcus_maripaludis_s2.GCA_000011585.1.30.dna.genome.fasta'  %(genome_dir)
if organism == "dvh":
    snpeff_db = "Desulfovibrio_vulgaris_str_hildenborough"
    genome_fasta = '%s/Desulfovibrio_vulgaris_str_hildenborough.GCA_000195755.1.30.dna.genome.fasta' %(genome_dir)


############# Functions ##############
def get_data():
    arguments = sys.argv
    print('Arguments: %s' %arguments)
    data_folders = glob.glob('%s/%s' %(data_dir,arguments[1]))
    data_folders = [element for element in data_folders if element not in ('%s,%s')%(data_trimmed_dir,fastqc_dir)]
    print ('data_folders: %s,%s' %(len(data_folders),data_folders))
    return data_folders

def create_dirs(samtools_results,gatk_results,varscan_results,data_trimmed_dir,fastqc_dir,alignment_results,combined_variants):
    dirs = [samtools_results,gatk_results,varscan_results,data_trimmed_dir,fastqc_dir,alignment_results,combined_variants]
    for dir in dirs:
        # create results folder
        print(dir)
        if not os.path.exists('%s' %(dir)):
            os.makedirs('%s' %(dir))
            print ('\033[31m %s directory doesn NOT exists. I am creating it. \033[0m' %(dir))
        else:
            print ('\033[31m %s directory exists. Not creating. \033[0m' %(dir))


####################### Trimgalore for quality and trimming ###############################
def trimgalore(first_pair_file,second_pair_file,folder_name,sample_id,file_ext):
    #print("1stpair:%s, 2ndpair:%s, folder_name:%s, sample_name:%s")%(first_pair_file,second_pair_file,folder_name,sample_name)
    print
    print("\033[34m Running TrimGalore \033[0m")
    # create sample spepcific trimmed directory
    if not os.path.exists('%s/%s' %(data_trimmed_dir,folder_name)):
        os.makedirs('%s/%s' %(data_trimmed_dir,folder_name))
    # create sample spepcific fastqcdirectory
    if not os.path.exists('%s/%s' %(fastqc_dir,folder_name)):
        os.makedirs('%s/%s' %(fastqc_dir,folder_name))
    # run Command
    cmd = 'trim_galore --fastqc_args "--outdir %s/%s/" --paired --output_dir %s/%s/ %s %s' %(fastqc_dir,folder_name,data_trimmed_dir,folder_name,first_pair_file, second_pair_file)
    print
    print ('++++++ Trimgalore Command:', cmd)
    print
    #os.system(cmd)


####################### BWA for alignment ###############################
def runBWA(alignment_results,file_ext,first_file_name,second_file_name,lane,folder_name,sample_id,RGId, RGSm, RGLb, RGPu,files_2_delete):
    print
    print ("\033[34m Running BWA alignment... \033[0m")
    # create sample spepcific results directory
    if not os.path.exists('%s' %(alignment_results)):
        os.makedirs('%s' %(alignment_results))
    # define result files
    if file_ext == "gz":
        first_pair_trimmed = '%s/%s/%s_val_1.fq.gz'%(data_trimmed_dir,folder_name,first_file_name)
        second_pair_trimmed = '%s/%s/%s_val_2.fq.gz'%(data_trimmed_dir,folder_name,second_file_name)
    else:
        first_pair_trimmed = '%s/%s/%s_val_1.fq'%(data_trimmed_dir,folder_name,first_file_name)
        second_pair_trimmed = '%s/%s/%s_val_2.fq'%(data_trimmed_dir,folder_name,second_file_name)
    print ('Trimmed Files:\n 1st:%s \n 2nd:%s' %(first_pair_trimmed,second_pair_trimmed))
    print

    # modify read group information
    read_group = "@RG\\tID:%s\\tPL:ILLUMINA\\tSM:%s\\tLB:%s\\tPU:%s" %(RGId, RGSm, RGLb, RGPu)

    base_file_name = '%s/%s'%(alignment_results,sample_id)
    print('Base Filename: %s' %(base_file_name))

    # bwa run command
    cmd = "bwa mem -t 4 -R '%s' %s %s %s > %s.sam" %(read_group, genome_fasta, first_pair_trimmed, second_pair_trimmed,base_file_name)
    print( "++++++ Run BWA Command", cmd)
    os.system(cmd)

    files_2_delete.append('%s.sam'%(base_file_name))
    return base_file_name


####################### samtools fixmate, sort and index to cleanup read pair info and flags ###############################
def run_samtools_fixmate(base_file_name,sample_id,files_2_delete):
    print
    print( "\033[34m Running SAMtools fixmate... \033[0m")
    cmd1 = 'samtools fixmate -O bam %s.sam %s_fixmate.bam' %(base_file_name,base_file_name)
    cmd2 = 'samtools sort  -@ 8 -O bam -o %s_sorted.bam -T /proj/omics4tb/sturkarslan/tmp/%s_temp %s_fixmate.bam' %(base_file_name, sample_id, base_file_name)
    cmd3 = 'samtools index %s_sorted.bam' %(base_file_name)

    print ("++++++ Samtools Fixmate Command: ", cmd1)
    os.system(cmd1)
    print ("++++++ Samtools Sort Command: ", cmd2)
    os.system(cmd2)
    print ("++++++ Samtools Index Command: ", cmd3)
    os.system(cmd3)
    # add  temp files to list to delete
    temp_files = ['%s_sorted.bam'%base_file_name, '%s_sorted.bam.bai'%base_file_name, '%s_fixmate.bam'%base_file_name]
    for temp_file in temp_files:
        files_2_delete.append(temp_file)


####################### GATK 1st Pass ###############################
def runGATK(base_file_name,files_2_delete):
    print
    print ("\033[34m Running GATK Realigner.. \033[0m")
    #run Target Interval Creater command
    cmd1 = '%s -Xmx128m -jar %s -T RealignerTargetCreator -R %s -I %s_sorted.bam -o %s.intervals' %(javaPath, gatkPath, genome_fasta, base_file_name,base_file_name)
    #run Indel Realigner command
    cmd2 = '%s -Xmx4G -jar %s -T IndelRealigner -R %s -I %s_sorted.bam -targetIntervals %s.intervals -o %s_realigned.bam' %(javaPath, gatkPath, genome_fasta, base_file_name, base_file_name, base_file_name)
    # index bam file
    cmd3 = 'samtools index %s_realigned.bam' %(base_file_name)
    # Detect covariates
    cmd4 = '%s -Xmx4G -jar %s -T BaseRecalibrator -R %s -I %s_realigned.bam -knownSites %s-variants-compiled_sorted.vcf -o %s_recal.table' %(javaPath, gatkPath, genome_fasta, base_file_name, organism, base_file_name)
    # Adjust quality scores
    cmd5 = '%s -Xmx4G -jar %s -T PrintReads -R %s -I %s_realigned.bam -BQSR %s_recal.table -o %s_recal.bam' %(javaPath, gatkPath, genome_fasta, base_file_name, base_file_name, base_file_name)

    print("++++++ Command GATK Interval Creater: ", cmd1)
    print
    os.system(cmd1)
    print
    print( "++++++ Command GATK Realigner: ", cmd2)
    print
    os.system(cmd2)
    print
    print( "++++++ Command GATK index BAM: ", cmd3)
    print
    os.system(cmd3)
    print
    print("++++++ Command GATK BaseRecalibrator: ", cmd4)
    print
    os.system(cmd4)
    print
    print( "++++++ Command GATK PrintReads: ", cmd5)
    print
    os.system(cmd5)
    print
    temp_files = ['%s.intervals'%base_file_name,'%s_realigned.bam'%base_file_name,'%s_realigned.bai'%base_file_name,'%s_recal.table'%base_file_name,'%s_recal.bam'%base_file_name,'%s_recal.bai'%base_file_name]
    for temp_file in temp_files:
        files_2_delete.append(temp_file)


####################### Mark duplicates with GATK ###############################
def runMarkDuplicates(alignment_results,exp_name, base_file_name):
    print
    print("\033[34m Running Mark Duplicates.. \033[0m")
    print( "Exp Name:" + exp_name)
    print
    # collect list of recalibrated bam files
    recal_bams = glob.glob('%s*_recal.bam' %(base_file_name))
    print('recal_bams:%s' %(recal_bams))
    # creaate command line parameter for each file
    print( 'Input Recalibrated BAM Files...')
    bamList = []
    for i in recal_bams:
        inputAdd = 'INPUT=%s' %(i)
        bamList.append(inputAdd)
        print(i)
    bam_list_joined = " ".join(bamList)
    print
    print( 'Output BAM File...')
    marked_bam_name = '%s/%s_marked.bam' %(alignment_results,exp_name)
    metrics_file = '%s/%s.metrics' %(alignment_results,exp_name)
    alignment_files_path = '%s/%s'%(alignment_results,exp_name)
    print (marked_bam_name)
    print('bam_list_joined:%s'%(bam_list_joined))

    # MArk Duplicates
    cmd1 = '%s -Xmx4G -jar %s MarkDuplicates %s VALIDATION_STRINGENCY=LENIENT METRICS_FILE=%s OUTPUT=%s' %(javaPath, piccardPath, bam_list_joined, metrics_file, marked_bam_name)
    print
    print( "++++++ Mark Duplicated Command:... ", cmd1)
    os.system(cmd1)

    # index bam file
    cmd2 = 'samtools index %s' %(marked_bam_name)
    print
    print( "++++++ Index bamfile Command:... ", cmd2)
    os.system(cmd2)
    return alignment_files_path


####################### Samtools Variant Calling ###############################
def samtools_variants(samtools_results,alignment_files_path,exp_name,folder_name): # With samtools
    print
    print( "\033[34m Running SAMtools Variant Calling.. \033[0m")
    # create samtools results specific results directory
    samtools_files_path = '%s/%s'%(samtools_results,exp_name)

    # Produce BCF file with all locations in the genome
    cmd1 = 'bcftools  mpileup --threads 8 -Ou -f %s %s_marked.bam | bcftools call -vmO z --ploidy 1 -o %s_samtools.vcf.gz' %(genome_fasta,alignment_files_path, samtools_files_path)
    # Prepare vcf file for querying
    cmd2 = 'tabix -p vcf %s_samtools.vcf.gz' %(samtools_files_path)
    #Filtering
    percentageString = "%"
    cmd3 = "bcftools filter -O z -o %s_samtools_final.vcf.gz -s LOWQUAL -i '%sQUAL>30' %s_samtools.vcf.gz" %(samtools_files_path, percentageString, samtools_files_path)
    print( "++++++ Variant Calling mpileup: ", cmd1)
    os.system(cmd1)
    print
    print( "++++++ Variant Calling tabix: ", cmd2)
    os.system(cmd2)
    print
    print( "++++++ Variant Calling filtering: ", cmd3)
    os.system(cmd3)
    print
    #files_2_delete.append('%s_sorted.bam, %s_sorted.bam.bai, %s_fixmate.bam'%(base_file_name,base_file_name,base_file_name))
    return samtools_files_path


####################### Varscan Variant Calling ###############################
def varscan_variants(alignment_files_path,varscan_results,exp_name,folder_name,files_2_delete): # with varscan
    print
    print( "\033[34m Running Varscan.. \033[0m")
    # create varscan results specific results directory
    varscan_files_path = '%s/%s'%(varscan_results,exp_name)

    # varscan mpileup
    cmd1 = 'samtools mpileup --input-fmt-option nthreads=8 -B -f %s -o %s.pileup %s_marked.bam' %(genome_fasta, varscan_files_path, alignment_files_path)
    # varscan for snps
    cmd2 = '%s -Xmx128m -jar %s mpileup2snp %s.pileup --output-vcf 1 --min-coverage 8 --min-reads2 2 --min-avg-qual 30 --strand-filter 0 > %s_varscan_snps_final.vcf' %(javaPath, varscanPath, varscan_files_path, varscan_files_path)
    # varscan for indels
    cmd3 = '%s -Xmx128m -jar %s mpileup2indel %s.pileup --output-vcf 1 --min-coverage 8 --min-reads2 2 --min-avg-qual 30 --strand-filter 0 > %s_varscan_inds_final.vcf' %(javaPath, varscanPath, varscan_files_path, varscan_files_path)
    print( "++++++ samtools Mpileup: ", cmd1)
    print
    os.system(cmd1)
    print( "++++++ Varscan for SNPs: ", cmd2)
    print
    os.system(cmd2)
    print
    print( "++++++ Varscan for INDELS: ", cmd3)
    os.system(cmd3)

    files_2_delete.append('%s.pileup' %(varscan_files_path))

    return varscan_files_path



####################### GATK Variant Calling ###############################
def gatk_variants(alignment_files_path,gatk_results,exp_name,folder_name): #with GATK HaploTypeCaller
    print
    print( "\033[34m Running GATK Haplotype Variant Caller.. \033[0m")
    # create varscan results specific results directory
    gatk_files_path = '%s/%s'%(gatk_results,exp_name)

    # haplotype command
    cmd1 = '%s -Xmx4G -jar %s -T HaplotypeCaller -R %s -I %s_marked.bam --genotyping_mode DISCOVERY -ploidy 1 -stand_emit_conf 30 -stand_call_conf 30 -o %s_gatk_raw.vcf' %(javaPath, gatkPath, genome_fasta, alignment_files_path, gatk_files_path)
    # Select snp variants
    cmd2 = '%s -Xmx128m -jar %s -T SelectVariants -R %s -V %s_gatk_raw.vcf -selectType SNP -o %s_gatk_snps.vcf' %(javaPath, gatkPath, genome_fasta, gatk_files_path,gatk_files_path)
    # Apply filters to SNPs
    cmd3 = "%s -Xmx128m -jar %s -T VariantFiltration -R %s -V %s_gatk_snps.vcf --filterExpression 'QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0' --filterName 'my_snp_filter' -o %s_gatk_snps_filtered.vcf" %(javaPath, gatkPath, genome_fasta, gatk_files_path, gatk_files_path)
    # Select indel variants
    cmd4 = '%s -Xmx128m -jar %s -T SelectVariants -R %s -V %s_gatk_raw.vcf -selectType INDEL -o %s_gatk_inds.vcf' %(javaPath, gatkPath, genome_fasta, gatk_files_path,gatk_files_path)
    # Apply filters to indels
    cmd5 = "%s -Xmx128m -jar %s -T VariantFiltration -R %s -V %s_gatk_inds.vcf --filterExpression 'QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0' --filterName 'my_indel_filter' -o %s_gatk_inds_filtered.vcf" %(javaPath, gatkPath, genome_fasta, gatk_files_path, gatk_files_path)
    # Merge vcf files
    cmd6 = "%s -Xmx128m -jar %s -T CombineVariants -R %s --variant %s_gatk_snps_filtered.vcf --variant %s_gatk_inds_filtered.vcf -o %s_gatk_final.vcf  -genotypeMergeOptions UNSORTED" %(javaPath, gatkPath, genome_fasta, gatk_files_path, gatk_files_path, gatk_files_path)

    print( "++++++ GATK HaplotypeCaller Comnand: ", cmd1)
    print
    os.system(cmd1)
    print( "++++++ Select SNP Variants: ", cmd2)
    print
    os.system(cmd2)
    print( "++++++ Applying filters for SNPs: ", cmd3)
    print
    os.system(cmd3)
    print( "++++++ Select IndelVariants: ", cmd4)
    print
    os.system(cmd4)
    print( "++++++ Applying filters for Indelss: ", cmd5)
    print
    os.system(cmd5)
    print( "++++++ Merging vcf files: ", cmd6)
    print
    os.system(cmd6)

    return gatk_files_path


####################### Run SNPEff annotations ###############################
def run_snpeff(samtools_files_path,varscan_files_path,gatk_files_path,combined_variants,exp_name):
    print
    print( "\033[34m Running SNPEff Annotations.. \033[0m")

    vcf_files = ['%s_samtools_final.vcf.gz'%(samtools_files_path), '%s_varscan_inds_final.vcf'%(varscan_files_path), '%s_varscan_snps_final.vcf'%(varscan_files_path), '%s_gatk_final.vcf'%(gatk_files_path)]
    #print('vcf_files:%s' %(vcf_files))
    varscan_files = ['%s_varscan_inds_final.vcf'%(varscan_files_path), '%s_varscan_snps_final.vcf'%(varscan_files_path)]
    #print('varscan_files:%s' %(varscan_files))

    # create output file for combined variants output
    combined_variants_output = '%s/%s_combined_variants.txt' %(combined_variants,exp_name)
    print('combined_variants_output:%s' %(combined_variants_output))

    # open the final output file for writing before combining
    run_snpeff.t = open(combined_variants_output, 'w')

    for vcf_file in vcf_files:
        print(vcf_file + '\n')
        # creat file names for output
        snpeff_vcf = re.split('final.', vcf_file)[0] + 'snpeff.' + re.split('final.', vcf_file)[1]
        snpeff_filtered_vcf = re.split('final.', vcf_file)[0] + 'snpeff_filtered.' + re.split('final.', vcf_file)[1]
        snpeff_stats = re.split('final.', vcf_file)[0] + 'snpeff_stats.txt'
        snpeff_final = re.split('final.', vcf_file)[0] + 'snpeff_final.txt'

        # snpeff formateff commands
        cmd1 ='%s -Xmx2g -jar %ssnpEff.jar -ud 0 -classic -csvStats %s -geneId -lof -v -formatEff -o gatk %s %s > %s'  %(javaPath, snpeff_path, snpeff_stats, snpeff_db, vcf_file, snpeff_vcf)
        # snpeff filtering command
        cmd2 = 'cat %s | %s -Xmx128m -jar %sSnpSift.jar filter "(FILTER = \'PASS\') & (EFF[*].CODING != \'NON_CODING\')" > %s' %(snpeff_vcf, javaPath, snpeff_path, snpeff_filtered_vcf)
        # create final one line variant file
        cmd3 ='cat %s | perl %sscripts/vcfEffOnePerLine.pl | %s -Xmx128m -jar %sSnpSift.jar extractFields - CHROM POS REF ALT AF AC DP MQ "(FILTER = \'PASS\')" "EFF[*].EFFECT" "EFF[*].IMPACT" "EFF[*].FUNCLASS" "EFF[*].CODON" "EFF[*].AA" "EFF[*].AA_LEN" "EFF[*].GENE" "EFF[*].CODING" "EFF[*].RANK" "EFF[*].DISTANCE"> %s'%(snpeff_filtered_vcf, snpeff_path, javaPath, snpeff_path, snpeff_final)

        print( "++++++ Running SNPEff Formateff command: ", cmd1)
        print
        os.system(cmd1)
        print( "++++++ Running SNPEff Filtering: ", cmd2)
        print
        os.system(cmd2)
        print( "++++++ Running SNPEff Oneline final formatter: ", cmd3)
        print
        os.system(cmd3)

        # Run function to combine vcf files into a single file from 3 callers
        combine_variants(snpeff_filtered_vcf,snpeff_final,combined_variants)
    run_snpeff.t.close()

####################### Combine variants from snpeff outputs ###############################
def combine_variants(snpeff_filtered_vcf,snpeff_final,combined_variants):
    print
    print( "\033[34m Running Combine variants.. \033[0m")

    vcf_file_program = snpeff_filtered_vcf.split('/')[-2]
    print('snpeff_filtered_vcf:%s' %snpeff_filtered_vcf)

    ## open vcf file for processing and skip comment lines
    infile = open(snpeff_filtered_vcf, 'r')
    if os.stat(snpeff_filtered_vcf).st_size > 0:

        for num, line in enumerate(infile, 1):
            if "#CHROM" in line:
                N = num
                print( 'N ='+ str(N))

        f = open(snpeff_filtered_vcf, 'r')
        for i in range(N):
            try:
                f.readline()
            except NameError:
                var_exists = False
                print( "N is not defined")
            else:
                var_exists = True
                #print( "var_exists")
                pass
                #next(f)
    else:
        print( "Empty File")

    ## we want to collect freq column for varsab vcf file
    frequencies = []
    #frequencies.append('') # skip the headerline
    for line in f:
        vector = line.split("\t")
        #print(vector)
        ## if vcf file is varscan file
        if vcf_file_program == "varscan_results":
            ## collect vcf field for frequency
            freq = vector[9].split(':')[6]

            #print('++++++++Freq: %s+++++++++' %(freq))
            chm = line.split('\t')[0]
            my_fields = vector[7].split(';')
            #adp = myFields.split("ADP=")[1]
            ## append freq to list of freq
            frequencies.append(freq)

            ## if vcf file is samtools file
        elif vcf_file_program == "samtools_results":
            #collect list of all fields from vcf file
            my_fields = vector[7].split(';')
            #print(my_fields)
            # grab DP4 field
            dp4 = next(filter(lambda x:'DP4' in x, my_fields))
            #print dp4
            dp4_fields = dp4.split("DP4=")[1]
            fr = float(dp4_fields.split(",")[0])
            rr = float(dp4_fields.split(",")[1])
            fa = float(dp4_fields.split(",")[2])
            ra = float(dp4_fields.split(",")[3])
            freq = (fa + ra) / (fr + rr + fa + ra)
            freq = round(freq * 100, 2)
            freq = str(freq) + "%"

            frequencies.append(freq)

            ## if vcf file is gatk file
        elif vcf_file_program == "gatk_results":
            ## collect vcf field
            my_fields = vector[9].split(':')
            #print myFields
            if my_fields[0] != ".":
                ad = float(my_fields[1].split(',')[1])
                dp = float(my_fields[2])
                freq = ad / dp
                freq = round(freq * 100, 2)
                freq = str(freq) + "%"
            ## append freq to list of freq
                frequencies.append(freq)
            else:
                frequencies.append('')
        else:
            frequencies.append('')

    print( 'Length of Frequencies= '+ str(len(frequencies)))
    #print(frequencies)
    #sys.exit()
    # filename for the output from converting vcf to oneliner with frequency and program added
    outfile_w_freq = combined_variants + '/' + snpeff_filtered_vcf.split('/')[-1].split('_snpeff_filtered')[0] + '_outfile_w_freq.txt'

    s = open(outfile_w_freq, 'w')
    g = open(snpeff_final, 'r')
    next(g)
    myList = []
    index = 0
    for line in g:
        # change alternative chromosome names
        chromosome = line.split('\t')[0]
        if chromosome == "NC_002937":
            print( "Found alternative chromosome name: %s" %chromosome)
            chm = "Chromosome"
        elif chromosome == "NC_005863":
            print( "Found alternative chromosome name: %s" %chromosome)
            chm = "pDV"
        else:
            chm = line.split('\t')[0]

        pos = line.split('\t')[1]
        refs = list(line.split('\t')[2])[0] # grab only first character in reference sequence
        ref = line.split('\t')[2]
        alt = line.split('\t')[3]
        eff = line.split('\t')[9]
        imp = line.split('\t')[10]
        fnc = line.split('\t')[11]
        cdn = line.split('\t')[12]
        aac = line.split('\t')[13]
        loc = line.split('\t')[15]
        cod = line.split('\t')[16]
        dep = line.split('\t')[6]
        fre = frequencies[index]
        pro = vcf_file_program.split('_')[0]

        lineToWrite = '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' %(chm,pos,refs,ref,alt,eff,imp,fnc,cdn,aac,loc,cod,fre,pro,dep)
        s.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' %(chm,pos,refs, ref,alt,eff,imp,fnc,cdn,aac,loc,cod,fre,pro,dep))

        myList.append(lineToWrite)
        index = index + 1
    # write into final output file
    for element in myList:
        run_snpeff.t.write('%s' %element)

####################### Collate variants from three programs ###############################
def salomon(aggregatedList):
    #this function takes a list of lists of variants and find the consensus list

    # f.1. finding the unique positions
    uniqueLocations=[]
    for variants in aggregatedList:
        for variant in variants:
            uniqueLocation=variant[:3]
            if uniqueLocation not in uniqueLocations:
                uniqueLocations.append(uniqueLocation)

    # f.2. building the full consensus list
    consensus_list=[]
    for uniqueLocation in uniqueLocations:
        callers=[]
        body=[]
        freqs = []
        dps = []
        freq=''
        dp =''
        freqFloats = []

        for variants in aggregatedList:
            for variant in variants:
                if uniqueLocation == variant[:3]:
                    body=variant[:-2]
                    callers.append(variant[-2])
                    if variant[-2] == 'varscan':
                        freq=variant[-3]
                        freqs.append(freq)
                        freqFloat = freq.split('%')[0]
                        freqFloats.append(freqFloat)
                        #print variant[:-2]
                        dp = "NA"
                        dps.append(dp)

                        if freq == '':
                            print( 'WARNING varscan did not provide frequency value for variant')
                            stringVariant='\t'.join(variant)
                            #print stringVariant

                    if variant[-2] == 'samtools':
                        freq=variant[-3]
                        freqs.append(freq)
                        freqFloat = freq.split('%')[0]
                        freqFloats.append(freqFloat)
                        dp =variant[-1]
                        dps.append(dp)
                        if freq == '':
                            print( 'WARNING samtools did not provide frequency value for variant')
                            stringVariant='\t'.join(variant)
                            #print stringVariant

                    if variant[-2] == 'gatk':
                        dp = variant[-1]
                        dps.append(dp)
                        freq=variant[-3]
                        freqs.append(freq)
                        freqFloat = freq.split('%')[0]
                        freqFloats.append(freqFloat)

        # incorporating what we found
        callersString=':'.join(callers)
        freqsString=':'.join(freqs)
        dpsString=':'.join(dps)
        freqFloatString = max(freqFloats)

        # Get the samtools frequency as final frequency if available otherwise pick varscan freq.
        print(callers)
        print(freqs)
        if 'samtools' in callers and len(callers) > 1:
            final_freq = freqs[callers.index('samtools')]
        elif 'varscan' in callers and 'samtools' not in callers:
            final_freq = freqs[callers.index('varscan')]
        else:
            final_freq = freqs[0]

        body.append(callersString)
        body.append(freqsString)
        body.append(dpsString)
        body.append(freqFloatString)
        body.append(final_freq)
        consensus_list.append(body)

    return consensus_list


####################### Variant Retriever ###############################
def variantRetriever(combined_output_file):

    # this function retrieves the variants for each caller
    varscan_list = []
    gatk_list = []
    samtools_list = []
    print('combined_output')
    print(os.stat(combined_output_file).st_size)

    # Start reading each line and append to appropriate program list
    with open(combined_output_file) as f:
        for line in f:
            vector = line.split('\t')
            vector[-1]=vector[-1].replace('\n','')
            program=vector[-2]
            #print('program is:%s' %(program))
            if program == 'varscan':
                varscan_list.append(vector)
            if program == 'gatk':
                gatk_list.append(vector)
            if program == 'samtools':
                samtools_list.append(vector)
        f.close()

    return varscan_list,gatk_list,samtools_list


####################### Collate variants ###############################
def collate_variants(combined_output_file,merged_variants_file):
    print
    print( "\033[34m Running Collate variants.. \033[0m")

    # 2. recovering list of variants
    varscan_list,gatk_list,samtools_list = variantRetriever(combined_output_file)
    print( 'detected variants',len(varscan_list),len(gatk_list),len(samtools_list))

    # 3. finding consensus list of variants
    print( 'merging...')
    consensus_list = salomon([varscan_list,gatk_list,samtools_list])
    print( 'final set',len(consensus_list))

    # 4. writing a consensus list
    print( 'writing file...')
    #print(merged_variants_file)

    g = open(merged_variants_file, 'w')
    for element in consensus_list:
        line2write='\t'.join(element)
        line2write=line2write+'\n'
        g.write(line2write)
    g.close()
    #return collate_variants.combined_output_file


####################### Delete temporary files ###############################
def delete_temp_files(files_2_delete):
    print
    print( "\033[34m Deleting Temporry files.. \033[0m")
    for file in files_2_delete:
        cmd = 'rm %s' %(file)
        print(cmd)
        os.system(cmd)


####################### Running the Pipeline ###############################
def run_pipeline():
    folder_count = 1
    data_folders = get_data()
    #sys.exit()
    files_2_delete = [] # create a list of files to delete later and keep adding to the list
    #data_folders = tqdm(data_folders)
    # Loop through each data folder
    for data_folder in data_folders:
        #data_folders.set_description('Processing %s\n' %data_folder)
        folder_name = data_folder.split('/')[-1]
        print
        print
        print( '\033[33mProcessing Folder: %s of %s (%s)\033[0m' %(folder_count, len(data_folders), folder_name))

        # get the list of first file names in paired end sequences
        first_pair_files = glob.glob('%s/*_1.fastq*' %(data_folder))
        second_pair_files = glob.glob('%s/*_2.fastq*' %(data_folder))
        #print(first_pair_files, second_pair_files)

        # Program specific results directories
        samtools_results = "%s/results/%s/%s/samtools_results" %(run_dir,organism,folder_name)
        gatk_results = "%s/results/%s/%s/gatk_results" %(run_dir,organism,folder_name)
        varscan_results = "%s/results/%s/%s/varscan_results" %(run_dir,organism,folder_name)
        alignment_results = "%s/results/%s/%s/alignment_results" %(run_dir,organism,folder_name)
        combined_variants = "%s/results/%s/%s/combined_variants" %(run_dir,organism,folder_name)
        # final results files
        combined_output_file = '%s/%s_combined_variants.txt' %(combined_variants,folder_name)
        #print('combined_output_file:%s' %(combined_output_file))
        merged_variants_file = '%s/%s_merged_variants_final.txt' %(combined_variants,folder_name)
        #print('combined_output_file:%s' %(combined_output_file))

        # Loop through each file and create filenames
        file_count = 1
        for first_pair_file in first_pair_files:
            first_file_name_full = first_pair_file.split('/')[-1]
            second_pair_file = first_pair_file.replace('_1.', '_2.')
            second_file_name_full = second_pair_file.split('/')[-1]
            file_ext = first_pair_file.split('.')[-1]

            print( '\033[32m Processing File: %s of %s (%s)\033[0m' %(file_count, len(first_pair_files), first_file_name_full ))

            first_file_name = re.split('.fastq|.fastq.gz',first_file_name_full)[0]
            second_file_name = re.split('.fastq|.fastq.gz',second_file_name_full)[0]
            print('first_file_name:%s, second_file_name:%s' %(first_file_name,second_file_name))

            # Collect Sample attributes
            exp_name = folder_name
            print("exp_name: %s" %(exp_name))
            lane = "L001" #fileName.split("/")[-1].split("_")[2]
            print("Lane: %s" %(lane))
            #sample_id = file_name.split('_R1')[0]
            sample_id = re.split('.fastq|.fastq.gz', first_file_name)[0]
            #sample_id = sample_id.replace('_R1', '')
            print("sample_id: %s"  %(sample_id))

            # create readgroup info
            RGId = sample_id
            RGSm = exp_name
            RGLb = exp_name
            RGPu = lane
            print( "RG: ID: %s SM: %s LB: %s PU: %s" %(RGId, RGSm, RGLb, RGPu))
            #sys.exit()

            # 00. Get directories
            create_dirs(samtools_results,gatk_results,varscan_results,data_trimmed_dir,fastqc_dir,alignment_results,combined_variants)
            # 01. Run TrimGalore
            trimgalore(first_pair_file,second_pair_file,folder_name,sample_id,file_ext)
            # 02. Run bwa alignment to produce SAM file.
            base_file_name = runBWA(alignment_results,file_ext,first_file_name,second_file_name,lane,folder_name,sample_id,RGId, RGSm, RGLb, RGPu,files_2_delete)
            # 03. Run samtools fixmate
            run_samtools_fixmate(base_file_name,sample_id,files_2_delete)
            #sys.exit()
            # 0.4 Run GATK 1st PASS
            runGATK(base_file_name,files_2_delete)
            file_count = file_count + 1

        # 05. Run Mark duplicates
        alignment_files_path = runMarkDuplicates(alignment_results,exp_name,base_file_name)
        # 06. Run Samtools variant calling
        samtools_files_path = samtools_variants(samtools_results,alignment_files_path,exp_name,folder_name)
        # 07. Run varscan variant calling
        varscan_files_path = varscan_variants(alignment_files_path,varscan_results,exp_name,folder_name,files_2_delete)
        # 08. Run GATK variant Calling
        gatk_files_path = gatk_variants(alignment_files_path,gatk_results,exp_name,folder_name)
        # 09. Run SNPEff annotations
        vcf_file = run_snpeff(samtools_files_path,varscan_files_path,gatk_files_path,combined_variants,exp_name)
        # 10. Collate variants into single file from 3 callers and unify them
        collate_variants(combined_output_file,merged_variants_file)
        # 11. Delete temporary files_2_delete
        #delete_temp_files(files_2_delete)

        folder_count = folder_count + 1
        #sys.exit()

run_pipeline()

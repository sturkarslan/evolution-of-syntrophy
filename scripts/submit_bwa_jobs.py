import glob
import sys
import os
import string

# data and results directories
run_dir = "/proj/omics4tb/sturkarslan/syntrophy_raw_sequences/Early-gen"
data_dir = "%s/data" %run_dir
data_trimmed_dir = "%s/trimmed" %data_dir
fastqc_dir = "%s/data/fastqc_results" %run_dir
jobscripts_dir = "%s/scripts/jobscripts_mmp" %run_dir
jobscripts_logs = "%s/logs" %jobscripts_dir
# create sample spepcific results directory
if not os.path.exists('%s' %(jobscripts_dir)):
    os.makedirs('%s' %(jobscripts_dir))

if not os.path.exists('%s' %(jobscripts_logs)):
    os.makedirs('%s' %(jobscripts_logs))


############# Functions ##############
data_folders = glob.glob('%s/*' %(data_dir))
data_folders = [element for element in data_folders if '_15' in element or '_76' in element]

data_folders = [element for element in data_folders if element not in ('%s,%s,%s/HA2_15,%s/HA2_76,%s/UR1_76,%s/UR1_45, %s/old')%(data_trimmed_dir,fastqc_dir,data_dir,data_dir,data_dir,data_dir,data_dir)]    

print(data_folders)
#sys.exit()

folderCount = 1
for data_folder in data_folders:
    folder_name = data_folder.split('/')[-1]
    jobscript = '%s/%s_bwa_pipeline.csh' %(jobscripts_dir,folder_name)
    # write to job file
    with open(jobscript,'w') as g:
      g.write('#!/bin/bash\n\n')
      g.write('#$ -N %s\n'%(folder_name))
      g.write('#$ -o %s/%s_log.txt\n' %(jobscripts_logs,folder_name))
      g.write('#$ -e %s/%s_log.txt\n' %(jobscripts_logs,folder_name))
      g.write('#$ -pe smp 40\n')
      g.write('#$ -S /bin/bash\n\n')

      g.write('#Sample: %s\n' %(data_folder))

      # changing terminal to bash
      g.write('bash\n\n')
      g.write('source /users/sturkars/.bashrc \n\n')

      # change directory
      g.write('cd %s/scripts\n\n' %(run_dir))

      job_cmd = 'python bwa_pipeline_2.0.py %s' %(folder_name)
      print(job_cmd)
      # spades command
      g.write('%s' %(job_cmd))

    g.close()
    folderCount = folderCount + 1

    # submit each job with qsub
    cmd = 'qsub %s' %jobscript
    print(cmd)
    print
    os.system(cmd)
    #sys.exit()

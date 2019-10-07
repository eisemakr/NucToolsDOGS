#!/usr/bin/python

#===============================================================================
# This scripts starts the Nuctools pipeline 
# PRE: The MySQL Database must contain the 'samplesheet'-table
#===============================================================================

import mysql.connector
import os
from optparse import OptionParser
import subprocess


#===============================================================================
# Parsing arguments
#===============================================================================
parser = OptionParser(usage="./Nuctoolsdev.py [-f] [-u] [-a] [-p] [-d] [-j] [-t]", version="FunGen Nuctools 2019")
parser.set_defaults(verbose=True)
parser.add_option('-f', '--projectid', dest='projectid', help='project id', metavar='string')
parser.add_option('-u', '--url', dest='url', help='url of sql-server  e.g. 10.11.22.10"', metavar='url')
parser.add_option('-a', '--acc', dest='acc', help='sql server', metavar='string')
parser.add_option('-p', '--pw', dest='pw', help='password', metavar='string')
parser.add_option('-d', '--db', dest='db', help='sql database', metavar='string')
parser.add_option('-j', '--project', dest='project', help='folder project folder where to save data', metavar='string')
parser.add_option('-t', '--threads', type='int', dest='ncores', default=20, help='number of compute threads for NGM. [default: 20]', metavar='string' )

(options, args) = parser.parse_args()
projects_dir = str(options.project)
projectid = str(options.projectid)
ncores = str(options.ncores)

h = projects_dir.replace("/mnt/fungen/", "/fungen-srv/") + "/Project_" + projectid + "/logs"
#old version!####
#logs_path = projects_dir + "/Project_" + projectid + "/logs"
#print projectid
#logs_path.replace("/mnt/fungen/","/fungen-srv/")
#################
#===============================================================================
# Step 1: Checking project director, if it does not exist -> exit
# We need to replace the /mnt/fungen to the path of the filer
#===============================================================================
if not os.path.exists(projects_dir + "/Project_" + projectid):
    print "\n#Can not find project folder " + projects_dir + "/Sample_" + projectid
    exit()

#===============================================================================
# Step 2: Connecting to MySQL-Database
#===============================================================================

conn = mysql.connector.Connect(host=options.url, user=options.acc, password=options.pw, database=options.db)
c = conn.cursor()


#===============================================================================
# Step 3: Check if tophat output is already there
#===============================================================================


command = "SELECT DISTINCT sample_sheet.sampleProject,sample_sheet.sampleID,sample_sheet.sampleRef,sample_sheet.restInsertSize FROM sample_sheet WHERE sample_sheet.sampleProject = '" + str(projectid) + "' ORDER BY 'sample_sheet.sampleID' ASC"
c.execute("" + str(command) + "")

for row in c:
	out = str(row).replace("(u'", "").replace("', u'", ",").replace("')", "")
	project_name = out.split(",")[0]
	sample_name = out.split(",")[1]
	sample_path = str(projects_dir + "/Project_" + project_name).replace(" ", "") + "/Sample_" + sample_name
	reference_name = out.split(",")[2].split("+")
	rest_insertsize = out.split(",")[3]
	
	for i in reference_name:
		conn2 = mysql.connector.Connect(host=options.url, user=options.acc, password=options.pw, database=options.db)
		c2 = conn2.cursor()
		temp = "SELECT genome_reference.path, genome_reference.extension FROM genome_reference WHERE genome_reference.sampleRef LIKE '" + str(i) + "'"
		c2.execute("" + str(temp) + "")
		genome_path = ""
		for genome in c2:
			res = str(genome).replace("(u'", "").replace("', u'", ",").replace("')", "")
			genome_path = res.split(",")[0] + res.split(",")[1][1:]
		c2.close()
		conn2.close()
	
		#if(os.path.exists(sample_path + "/ngm_" + i)):
		#    continue
	
		# TODO get insert length of paired end reads from db ...also only map RNA-Seq reads
		# R1 and R2 exist
		extension = ""
		qualExt =""
		ext=False
		if((os.path.exists(sample_path + "/" + sample_name + "_R1.fastq") and
			os.path.exists(sample_path + "/" + sample_name + "_R2.fastq")) or
			(os.path.exists(sample_path + "/" + sample_name + "_R1.fastq.gz") and
			os.path.exists(sample_path + "/" + sample_name + "_R2.fastq.gz"))):
			if(os.path.exists(sample_path + "/" + sample_name + "_R1.fastq.gz") and os.path.exists(sample_path + "/" + sample_name + "_R2.fastq.gz")):
				extension =  ".gz"
				if(os.path.exists(sample_path + "/" + sample_name + "_R1_AQ20.fastq.gz") and os.path.exists(sample_path + "/" + sample_name + "_R2_AQ20.fastq.gz")):
					qualExt="_AQ20"
			if(os.path.exists(sample_path + "/ngm_" + i)):
				print "\n#pid=$(echo \"mkdir " + sample_path + "/ngm_" + i + "\" | qsub -N mkdir_" + sample_name + " -l nodes=1:ppn=1 -l walltime=120:00:00 " + " -e " + logs_path + " -o " + logs_path + ")"
				print "#pid=$(echo \"ngm -t " + ncores + " --no-progress -b -o " + sample_path + "/ngm_" + i +"/accepted_hits.bam -r " + genome_path+ " -1 " + sample_path + "/" + sample_name + "_R1"+qualExt+".fastq"+extension+" -2 " + sample_path + "/" + sample_name + "_R2"+qualExt+".fastq"+extension +" 2> " + sample_path + "/ngm_" + i +"/ngm_mapping_reports.txt \" | qsub -N ngm_" + sample_name + "  -l nodes=1:ppn=" + ncores + "  -l walltime=120:00:00 -W depend=afterok:$pid -e " + logs_path + " -o " + logs_path + ")"
				print "#pid=$(echo \"/mnt/fungen/db01/software/Python-2.7.5/bin/python2.7 /mnt/fungen/db01/software/pipeline/report_R/bam_stat.py -i " + sample_path + "/ngm_" + i +"/accepted_hits.bam > " + sample_path + "/ngm_" + i + "/mapping_stat.txt 2>&1\" | qsub -N Stat_" + sample_name + " -l nodes=1:ppn=1 -l walltime=4:00:00 -W depend=afterok:$pid -e " + logs_path + " -o " + logs_path + ")"
				print "#pid=$(echo \"Rscript /mnt/fungen/db01/software/pipeline/report_R/mapping_report_ngm.R " + sample_name + " " + sample_path + "/ngm_" + i + " " + sample_path + "/ngm_" + i + "/mapping_stat.txt \" | qsub -N mappingReport_" + sample_name + " -l nodes=1:ppn=1 -l walltime=4:00:00 -W depend=afterok:$pid -e " + logs_path + " -o " + logs_path + ")"
				print "#pid=$(echo \"samtools view -b -f 4 " + sample_path + "/ngm_" + i +"/accepted_hits.bam > "+sample_path + "/ngm_" + i +"/"+ sample_name + "_unmapped_ngm.bam \" | qsub -N extractUnmapped" + sample_name + " -l nodes=1:ppn=1 -l walltime=4:00:00 -W depend=afterok:$pid -e " + logs_path + " -o " + logs_path + ")"
				#NEW make index for bam
				print "#pid=$(echo \"samtools sort " + sample_path + "/ngm_" + i + "/accepted_hits.bam "+ sample_path + "/ngm_" + i + "/accepted_hits_sorted \" | qsub -N indexBam" + sample_name + " -l nodes=1:ppn=1 -l walltime=4:00:00 -W depend=afterok:$pid -e " + logs_path + " -o " + logs_path + ")"
				print "#pid=$(echo \"samtools index " + sample_path + "/ngm_" + i + "/accepted_hits_sorted.bam \" | qsub -N indexBam" + sample_name + " -l nodes=1:ppn=1 -l walltime=4:00:00 -W depend=afterok:$pid -e " + logs_path + " -o " + logs_path + ")"
				print "#pid=$(echo \"bamtools convert -in "+sample_path + "/ngm_" + i +"/"+ sample_name + "_unmapped_ngm.bam -format fastq > "+sample_path + "/ngm_" + i +"/"+ sample_name + "_unmapped_ngm.fq \" | qsub -N bam2fq" + sample_name + " -l nodes=1:ppn=1 -l walltime=4:00:00 -W depend=afterok:$pid -e " + logs_path + " -o " + logs_path + ")"
				print "#pid=$(echo \"pigz -p 8 "+sample_path + "/ngm_" + i +"/"+ sample_name + "_unmapped_ngm.fq \" | qsub -N gzip_fq" + sample_name + " -l nodes=1:ppn=8 -l walltime=4:00:00 -W depend=afterok:$pid -e " + logs_path + " -o " + logs_path + ")"
				print "#pid=$(echo \"samtools sort -n " +sample_path + "/ngm_" + i +"/"+ sample_name + "_unmapped_ngm.bam " +sample_path + "/ngm_" + i +"/"+ sample_name + "_unmapped_ngm_sorted \" | qsub -N sortUnmappedBAM" + sample_name + " -l nodes=1:ppn=1 -l walltime=4:00:00 -W depend=afterok:$pid -e " + logs_path + " -o " + logs_path + ")"
				print "#pid=$(echo \"bamToFastq -i " +sample_path + "/ngm_" + i +"/"+ sample_name + "_unmapped_ngm_sorted.bam -fq "+sample_path + "/ngm_" + i +"/"+ sample_name + "_unmapped_ngm_R1.fq -fq2 "+sample_path + "/ngm_" + i +"/"+ sample_name + "_unmapped_ngm_R2.fq\" | qsub -N rmUnmappedBAM" + sample_name + " -l nodes=1:ppn=1 -l walltime=4:00:00 -W depend=afterok:$pid -e " + logs_path + " -o " + logs_path + ")"
				print "#pid=$(echo \"pigz -p 8 "+sample_path + "/ngm_" + i +"/"+ sample_name + "_unmapped_ngm_R1.fq \" | qsub -N gzip_fqUM1" + sample_name + " -l nodes=1:ppn=8 -l walltime=4:00:00 -W depend=afterok:$pid -e " + logs_path + " -o " + logs_path + ")"
				print "#pid=$(echo \"pigz -p 8 "+sample_path + "/ngm_" + i +"/"+ sample_name + "_unmapped_ngm_R2.fq \" | qsub -N gzip_fqUM2" + sample_name + " -l nodes=1:ppn=8 -l walltime=4:00:00 -W depend=afterok:$pid -e " + logs_path + " -o " + logs_path + ")"
				print "#pid=$(echo \"rm " +sample_path + "/ngm_" + i +"/"+ sample_name + "_unmapped_ngm.bam \" | qsub -N rmUnmappedBAM" + sample_name + " -l nodes=1:ppn=1 -l walltime=4:00:00 -W depend=afterok:$pid -e " + logs_path + " -o " + logs_path + ")"
				print "#pid=$(echo \"rm " + sample_path + "/ngm_" + i + "/accepted_hits.bam \" | qsub -N rmUnsortedBAM" + sample_name + " -l nodes=1:ppn=1 -l walltime=4:00:00 -W depend=afterok:$pid -e " + logs_path + " -o " + logs_path + ")"
			else:
				print "\npid=$(echo \"mkdir " + sample_path + "/ngm_" + i + "\" | qsub -N mkdir_" + sample_name + " -l nodes=1:ppn=1 -l walltime=120:00:00 " + " -e " + logs_path + " -o " + logs_path + ")"
				print "pid=$(echo \"ngm -t " + ncores + " --no-progress -b -o " + sample_path + "/ngm_" + i +"/accepted_hits.bam -r " + genome_path+ " -1 " + sample_path + "/" + sample_name + "_R1"+qualExt+".fastq"+extension+" -2 " + sample_path + "/" + sample_name + "_R2"+qualExt+".fastq"+extension +" 2> " + sample_path + "/ngm_" + i +"/ngm_mapping_reports.txt \" | qsub -N ngm_" + sample_name + "  -l nodes=1:ppn=" + ncores + "  -l walltime=120:00:00 -W depend=afterok:$pid -e " + logs_path + " -o " + logs_path + ")"
				print "pid=$(echo \"/mnt/fungen/db01/software/Python-2.7.5/bin/python2.7 /mnt/fungen/db01/software/pipeline/report_R/bam_stat.py -i " + sample_path + "/ngm_" + i +"/accepted_hits.bam > " + sample_path + "/ngm_" + i + "/mapping_stat.txt 2>&1\" | qsub -N Stat_" + sample_name + " -l nodes=1:ppn=1 -l walltime=4:00:00 -W depend=afterok:$pid -e " + logs_path + " -o " + logs_path + ")"
				print "pid=$(echo \"Rscript /mnt/fungen/db01/software/pipeline/report_R/mapping_report_ngm.R " + sample_name + " " + sample_path + "/ngm_" + i + " " + sample_path + "/ngm_" + i + "/mapping_stat.txt \" | qsub -N mappingReport_" + sample_name + " -l nodes=1:ppn=1 -l walltime=4:00:00 -W depend=afterok:$pid -e " + logs_path + " -o " + logs_path + ")"
				print "pid=$(echo \"samtools view -b -f 4 " + sample_path + "/ngm_" + i +"/accepted_hits.bam > "+sample_path + "/ngm_" + i +"/"+ sample_name + "_unmapped_ngm.bam \" | qsub -N extractUnmapped" + sample_name + " -l nodes=1:ppn=1 -l walltime=4:00:00 -W depend=afterok:$pid -e " + logs_path + " -o " + logs_path + ")"
				#NEW make index for bam
				print "pid=$(echo \"samtools sort " + sample_path + "/ngm_" + i + "/accepted_hits.bam "+ sample_path + "/ngm_" + i + "/accepted_hits_sorted \" | qsub -N indexBam" + sample_name + " -l nodes=1:ppn=1 -l walltime=4:00:00 -W depend=afterok:$pid -e " + logs_path + " -o " + logs_path + ")"
				print "pid=$(echo \"samtools index " + sample_path + "/ngm_" + i + "/accepted_hits_sorted.bam \" | qsub -N indexBam" + sample_name + " -l nodes=1:ppn=1 -l walltime=4:00:00 -W depend=afterok:$pid -e " + logs_path + " -o " + logs_path + ")"
				print "pid=$(echo \"bamtools convert -in "+sample_path + "/ngm_" + i +"/"+ sample_name + "_unmapped_ngm.bam -format fastq > "+sample_path + "/ngm_" + i +"/"+ sample_name + "_unmapped_ngm.fq \" | qsub -N bam2fq" + sample_name + " -l nodes=1:ppn=1 -l walltime=4:00:00 -W depend=afterok:$pid -e " + logs_path + " -o " + logs_path + ")"
				print "pid=$(echo \"pigz -p 8 "+sample_path + "/ngm_" + i +"/"+ sample_name + "_unmapped_ngm.fq \" | qsub -N gzip_fq" + sample_name + " -l nodes=1:ppn=8 -l walltime=4:00:00 -W depend=afterok:$pid -e " + logs_path + " -o " + logs_path + ")"
				print "pid=$(echo \"samtools sort -n " +sample_path + "/ngm_" + i +"/"+ sample_name + "_unmapped_ngm.bam " +sample_path + "/ngm_" + i +"/"+ sample_name + "_unmapped_ngm_sorted \" | qsub -N sortUnmappedBAM" + sample_name + " -l nodes=1:ppn=1 -l walltime=4:00:00 -W depend=afterok:$pid -e " + logs_path + " -o " + logs_path + ")"
				print "pid=$(echo \"bamToFastq -i " +sample_path + "/ngm_" + i +"/"+ sample_name + "_unmapped_ngm_sorted.bam -fq "+sample_path + "/ngm_" + i +"/"+ sample_name + "_unmapped_ngm_R1.fq -fq2 "+sample_path + "/ngm_" + i +"/"+ sample_name + "_unmapped_ngm_R2.fq\" | qsub -N rmUnmappedBAM" + sample_name + " -l nodes=1:ppn=1 -l walltime=4:00:00 -W depend=afterok:$pid -e " + logs_path + " -o " + logs_path + ")"
				print "pid=$(echo \"pigz -p 8 "+sample_path + "/ngm_" + i +"/"+ sample_name + "_unmapped_ngm_R1.fq \" | qsub -N gzip_fqUM1" + sample_name + " -l nodes=1:ppn=8 -l walltime=4:00:00 -W depend=afterok:$pid -e " + logs_path + " -o " + logs_path + ")"
				print "pid=$(echo \"pigz -p 8 "+sample_path + "/ngm_" + i +"/"+ sample_name + "_unmapped_ngm_R2.fq \" | qsub -N gzip_fqUM2" + sample_name + " -l nodes=1:ppn=8 -l walltime=4:00:00 -W depend=afterok:$pid -e " + logs_path + " -o " + logs_path + ")"
				print "pid=$(echo \"rm " +sample_path + "/ngm_" + i +"/"+ sample_name + "_unmapped_ngm.bam \" | qsub -N rmUnmappedBAM" + sample_name + " -l nodes=1:ppn=1 -l walltime=4:00:00 -W depend=afterok:$pid -e " + logs_path + " -o " + logs_path + ")"
				print "pid=$(echo \"rm " + sample_path + "/ngm_" + i + "/accepted_hits.bam \" | qsub -N rmUnsortedBAM" + sample_name + " -l nodes=1:ppn=1 -l walltime=4:00:00 -W depend=afterok:$pid -e " + logs_path + " -o " + logs_path + ")"
			# Only single end R1 exists
		elif(os.path.exists(sample_path + "/" + sample_name + "_R1.fastq") and
			not os.path.exists(sample_path + "/" + sample_name + "_R2.fastq") or
			(os.path.exists(sample_path + "/" + sample_name + "_R1.fastq.gz") and
			not os.path.exists(sample_path + "/" + sample_name + "_R2.fastq.gz"))):
				if((os.path.exists(sample_path + "/" + sample_name + "_R1.fastq.gz") and not os.path.exists(sample_path + "/" + sample_name + "_R2.fastq.gz"))):
					extension = ".gz"
					if(os.path.exists(sample_path + "/" + sample_name + "_R1_AQ20.fastq.gz") and not os.path.exists(sample_path + "/" + sample_name + "_R2_AQ20.fastq.gz")):
						qualExt="_AQ20"
				if(os.path.exists(sample_path + "/ngm_" + i)):
					print "\n#pid=$(echo \"mkdir " + sample_path + "/ngm_" + i + "\" | qsub -N mkdir_" + sample_name + " -l nodes=1:ppn=1 -l walltime=120:00:00 " + " -e " + logs_path + " -o " + logs_path + ")"
					print "#pid=$(echo \"ngm -t " + ncores + " --no-progress -b -o " + sample_path + "/ngm_" + i +"/accepted_hits.bam -r " + genome_path + " -q " + sample_path + "/" + sample_name + "_R1"+qualExt+".fastq"+extension +" 2> " + sample_path + "/ngm_" + i +"/ngm_mapping_reports.txt \" | qsub -N ngm_" + sample_name + " -l nodes=1:ppn=" + ncores + "  -l walltime=120:00:00 -W depend=afterok:$pid -e " + logs_path + " -o " + logs_path + ")"
					print "#pid=$(echo \"/mnt/fungen/db01/software/Python-2.7.5/bin/python2.7 /mnt/fungen/db01/software/pipeline/report_R/bam_stat.py -i " + sample_path + "/ngm_" + i +"/accepted_hits.bam > " + sample_path + "/ngm_" + i + "/mapping_stat.txt 2>&1\" | qsub -N Stat_" + sample_name + " -l nodes=1:ppn=1 -l walltime=4:00:00 -W depend=afterok:$pid -e " + logs_path + " -o " + logs_path + ")"
					print "#pid=$(echo \"Rscript /mnt/fungen/db01/software/pipeline/report_R/mapping_report_ngm.R " + sample_name + " " + sample_path + "/ngm_" + i + " " + sample_path + "/ngm_" + i + "/mapping_stat.txt \" | qsub -N mappingReport_" + sample_name + " -l nodes=1:ppn=1 -l walltime=4:00:00 -W depend=afterok:$pid -e " + logs_path + " -o " + logs_path + ")"
					print "#pid=$(echo \"samtools view -b -f 4 " + sample_path + "/ngm_" + i +"/accepted_hits.bam > "+sample_path + "/ngm_" + i +"/"+ sample_name + "_unmapped_ngm.bam \" | qsub -N extractUnmapped" + sample_name + " -l nodes=1:ppn=1 -l walltime=4:00:00 -W depend=afterok:$pid -e " + logs_path + " -o " + logs_path + ")"
					#NEW make index for bam
					print "#pid=$(echo \"samtools sort " + sample_path + "/ngm_" + i + "/accepted_hits.bam "+ sample_path + "/ngm_" + i + "/accepted_hits_sorted \" | qsub -N indexBam" + sample_name + " -l nodes=1:ppn=1 -l walltime=4:00:00 -W depend=afterok:$pid -e " + logs_path + " -o " + logs_path + ")"
					print "#pid=$(echo \"samtools index " + sample_path + "/ngm_" + i + "/accepted_hits_sorted.bam \" | qsub -N indexBam" + sample_name + " -l nodes=1:ppn=1 -l walltime=4:00:00 -W depend=afterok:$pid -e " + logs_path + " -o " + logs_path + ")"
					print "#pid=$(echo \"bamtools convert -in "+sample_path + "/ngm_" + i +"/"+ sample_name + "_unmapped_ngm.bam -format fastq > "+sample_path + "/ngm_" + i +"/"+ sample_name + "_unmapped_ngm.fq \" | qsub -N bam2fq" + sample_name + " -l nodes=1:ppn=1 -l walltime=4:00:00 -W depend=afterok:$pid -e " + logs_path + " -o " + logs_path + ")"
					print "#pid=$(echo \"pigz -p 8  "+sample_path + "/ngm_" + i +"/"+ sample_name + "_unmapped_ngm.fq \" | qsub -N gzip_fq" + sample_name + " -l nodes=1:ppn=8 -l walltime=4:00:00 -W depend=afterok:$pid -e " + logs_path + " -o " + logs_path + ")"
					print "#pid=$(echo \"rm " +sample_path + "/ngm_" + i +"/"+ sample_name + "_unmapped_ngm.bam \" | qsub -N rmUnmappedBAM" + sample_name + " -l nodes=1:ppn=1 -l walltime=4:00:00 -W depend=afterok:$pid -e " + logs_path + " -o " + logs_path + ")"
					print "#pid=$(echo \"rm " + sample_path + "/ngm_" + i + "/accepted_hits.bam \" | qsub -N rmUnsortedBAM" + sample_name + " -l nodes=1:ppn=1 -l walltime=4:00:00 -W depend=afterok:$pid -e " + logs_path + " -o " + logs_path + ")"
				else:
					print "\npid=$(echo \"mkdir " + sample_path + "/ngm_" + i + "\" | qsub -N mkdir_" + sample_name + " -l nodes=1:ppn=1 -l walltime=120:00:00 " + " -e " + logs_path + " -o " + logs_path + ")"
					print "pid=$(echo \"ngm -t " + ncores + " --no-progress -b -o " + sample_path + "/ngm_" + i +"/accepted_hits.bam -r " + genome_path + " -q " + sample_path + "/" + sample_name + "_R1"+qualExt+".fastq"+extension +" 2> " + sample_path + "/ngm_" + i +"/ngm_mapping_reports.txt \" | qsub -N ngm_" + sample_name + " -l nodes=1:ppn=" + ncores + " -l walltime=120:00:00 -W depend=afterok:$pid -e " + logs_path + " -o " + logs_path + ")"
					print "pid=$(echo \"/mnt/fungen/db01/software/Python-2.7.5/bin/python2.7 /mnt/fungen/db01/software/pipeline/report_R/bam_stat.py -i " + sample_path + "/ngm_" + i +"/accepted_hits.bam > " + sample_path + "/ngm_" + i + "/mapping_stat.txt 2>&1\" | qsub -N Stat_" + sample_name + " -l nodes=1:ppn=1 -l walltime=4:00:00 -W depend=afterok:$pid -e " + logs_path + " -o " + logs_path + ")"
					print "pid=$(echo \"Rscript /mnt/fungen/db01/software/pipeline/report_R/mapping_report_ngm.R " + sample_name + " " + sample_path + "/ngm_" + i + " " + sample_path + "/ngm_" + i + "/mapping_stat.txt \" | qsub -N mappingReport_" + sample_name + " -l nodes=1:ppn=1 -l walltime=4:00:00 -W depend=afterok:$pid -e " + logs_path + " -o " + logs_path + ")"
					print "pid=$(echo \"samtools view -b -f 4 " + sample_path + "/ngm_" + i +"/accepted_hits.bam > "+sample_path + "/ngm_" + i +"/"+ sample_name + "_unmapped_ngm.bam \" | qsub -N extractUnmapped" + sample_name + " -l nodes=1:ppn=1 -l walltime=4:00:00 -W depend=afterok:$pid -e " + logs_path + " -o " + logs_path + ")"
					#NEW make index for bam
					print "pid=$(echo \"samtools sort " + sample_path + "/ngm_" + i + "/accepted_hits.bam "+ sample_path + "/ngm_" + i + "/accepted_hits_sorted \" | qsub -N indexBam" + sample_name + " -l nodes=1:ppn=1 -l walltime=4:00:00 -W depend=afterok:$pid -e " + logs_path + " -o " + logs_path + ")"
					print "pid=$(echo \"samtools index " + sample_path + "/ngm_" + i + "/accepted_hits_sorted.bam \" | qsub -N indexBam" + sample_name + " -l nodes=1:ppn=1 -l walltime=4:00:00 -W depend=afterok:$pid -e " + logs_path + " -o " + logs_path + ")"
					print "pid=$(echo \"bamtools convert -in "+sample_path + "/ngm_" + i +"/"+ sample_name + "_unmapped_ngm.bam -format fastq > "+sample_path + "/ngm_" + i +"/"+ sample_name + "_unmapped_ngm.fq \" | qsub -N bam2fq" + sample_name + " -l nodes=1:ppn=1 -l walltime=4:00:00 -W depend=afterok:$pid -e " + logs_path + " -o " + logs_path + ")"
					print "pid=$(echo \"pigz -p 8 "+sample_path + "/ngm_" + i +"/"+ sample_name + "_unmapped_ngm.fq \" | qsub -N gzip_fq" + sample_name + " -l nodes=1:ppn=8 -l walltime=4:00:00 -W depend=afterok:$pid -e " + logs_path + " -o " + logs_path + ")"
					print "pid=$(echo \"rm " +sample_path + "/ngm_" + i +"/"+ sample_name + "_unmapped_ngm.bam \" | qsub -N rmUnmappedBAM" + sample_name + " -l nodes=1:ppn=1 -l walltime=4:00:00 -W depend=afterok:$pid -e " + logs_path + " -o " + logs_path + ")"
					print "pid=$(echo \"rm " + sample_path + "/ngm_" + i + "/accepted_hits.bam \" | qsub -N rmUnsortedBAM" + sample_name + " -l nodes=1:ppn=1 -l walltime=4:00:00 -W depend=afterok:$pid -e " + logs_path + " -o " + logs_path + ")"
		
		#Insert Mapping Stats
		print "pid=$(echo \"python /mnt/fungen/db01/software/pipeline/report_R/insert_bam_stat_db.py --force -f "+projectid+" -u "+options.url+" -a "+options.acc+" -p "+options.pw+" -d "+options.db+"\" | qsub -N ins_stat" + projectid + " -l nodes=1:ppn=1 -W depend=afterok:$pid -e "+logs_path+" -o "+logs_path + ")"
		#Generate Mapping report (SQL)
		projects_path = str(projects_dir + "/Project_" + project_name).replace(" ", "")
		print "pid=$(echo \"Rscript /mnt/fungen/db01/software/pipeline/report_R/Mapping_report.SQL.R "+projectid+" "+projects_path+" 50 \" | qsub -N SQL.report." + projectid + " -l nodes=1:ppn=1 -W depend=afterok:$pid -e "+logs_path+" -o "+logs_path + ")"
		#Generate Mapping report from individual mapping statistics files
		print "pid=$(echo \"find "+projects_path+" -name \"mapping_stat.txt\" -print -exec cat \'{}\' \\; | grep -e \"mapping_stat.txt\" -e \\\"Total records\\\" -e \\\"Unmapped reads\\\" -e \\\"non-unique\\\" -e \\\"unique\\\" > "+projects_path+"/all_stat.txt \" | qsub -N mapping_report." + projectid + " -l nodes=1:ppn=1 -W depend=afterok:$pid -e "+logs_path+" -o "+logs_path + ")"
		print "pid=$(echo \"perl -w /mnt/fungen/db01/software/perl_scripts/convert_report_list2tab.pl --input="+projects_path+"/all_stat.txt --output="+projects_path+"/"+projectid+".mapping_stat.tab.txt \" | qsub -N ConvertMappingReport." + projectid + " -l nodes=1:ppn=1 -W depend=afterok:$pid -e "+logs_path+" -o "+logs_path + ")"

c.close()
conn.close()
exit()

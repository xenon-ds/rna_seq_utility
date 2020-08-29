#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 15 22:08:51 2019

@author: shuba
"""
import subprocess
import os

def createDirectory(directory):
    try:
    # Create target Directory
        os.mkdir(directory)
        print("Directory " , directory ,  " Created ") 
    except FileExistsError:
        print("Directory " , directory ,  " already exists")
    
 
def runsubprocess(commands):
    print("start runsubprocess")
    print("executing : \n"+commands)
    p = subprocess.Popen(commands, stdout=subprocess.PIPE, shell=True)
    (output, err) = p.communicate()
    # This makes the wait possible
    p_status = p.wait()
    print("end runsubprocess")
    rv=p.returncode
    while(rv is None):
        {
                #do nothing
        }



hisatPath=""

refFasta=input("Enter the path of the reference fasta : ")
#
indexPath=input("Enter the path (with name ) of the index you want : ")
nThreads=input("enter the number of threads (<20) : ")

cmd="hisat2-build -p "+nThreads+" "+refFasta+" "+indexPath
#rv=runsubprocess(cmd)

fqDir=input("enter the directory containing the fastq files : ")
files = []
extn=input("enter the file extn fastq or fq.gz : ")
nthreads=input("how many threads ( simultaneous files to be processed in fastqc : ")
fastqcOutDir=input("enter the directory path where fastqc result will be stored : ")
createDirectory(fastqcOutDir)
# Create directory

   
# r=root, d=directories, f = files
for r, d, f in os.walk(fqDir):
    for file in f:
        if ((extn in file) and (".txt" not in file)):
            print(file)
            files.append(os.path.join(r, file))

files.sort()
print(files)
print("will process the files mentioned below : ")
for f in files:
    print(f)

cmd="fastqc -o "+fastqcOutDir+" "+"-t "+str(nthreads)+" "
allFastqFiles=""
for f in files:
    allFastqFiles=allFastqFiles+f+" "
    
    
for f in files:
    cmd=cmd+f+" "

print(cmd)
#cmd="fastqc -h"
runsubprocess(cmd)

#############################
#trimming reads
#############################

trimGaloreOutDir=input("enter the trim galore output directory : ")
createDirectory(trimGaloreOutDir)
nCores=input("enter the number of cores you want to use : ")

#*********
#for single end reads 
#*********
cmd="trim_galore --dont_gzip --cores "+str(nCores)+" --output_dir "+trimGaloreOutDir+" "+allFastqFiles

#*********
#for paired end reads 
#*********
cmd="trim_galore --dont_gzip --cores "+str(nCores)+" --output_dir "+trimGaloreOutDir+" --paired --retain_unpaired "+allFastqFiles

print(cmd)
runsubprocess(cmd)


#############################
#aligning reads to the index
#############################
samOutputDir=input("enter the sam output directory : ")
createDirectory(samOutputDir)
cmd="hisat2 --dta-cufflinks -q -p "+nThreads+" -x "+indexPath+" -U "
trimmedFastqDir=trimGaloreOutDir
fastqFiles = []
samFiles=[]
extn="trimmed.fq"
fileNames=[]
for r, d, f in os.walk(trimmedFastqDir):
    for file in f:
        if extn in file:
            fileNames.append(file)
            fastqFiles.append(os.path.join(r, file))
            temp=file.split('.')
            tempname="sam_"+temp[0]+".sam"
            samFiles.append(os.path.join(samOutputDir, tempname))

totalFiles = len(files)

for i in range(totalFiles):
    tempcmd=cmd+fastqFiles[i]+" -S "+samFiles[i]
    print(tempcmd)
    runsubprocess(tempcmd)
    
#############################
#sorting the sam files 
#############################

cmd = 'samtools sort -O sam --threads 25 -o '
extn='.sam'
samFileNames=[]
for r, d, f in os.walk(samOutputDir):
    for file in f:
        print(file)
        if extn in file:
            samFileNames.append(file)
            samFiles.append(os.path.join(samOutputDir, file))
totalSamFiles=len(samFileNames)
sortedSamFilesOutputDir=input("enter the sorted sam files output directory : ")
for i in range(totalSamFiles):
    sortedSamName="sorted_"+samFileNames[i]
    tempcmd=cmd+os.path.join(sortedSamFilesOutputDir,sortedSamName)+" "+samFiles[i]
    print(tempcmd)
    runsubprocess(tempcmd)

##########################
#running cufflinks 
##########################
    
nThreads=input("enter the number of threads that you want for cufflinks : ")
cufflinksOutDir=input("enter the cufflinks output directory : ")
refGFF=input("Enter the path of the reference gff path : ")
cmd = "cufflinks -p "+nThreads+" -o "
extn="sorted_"
sortedSamFilesDir=sortedSamFilesOutputDir
for r, d, f in os.walk(sortedSamFilesDir):
    for file in f:
        if extn in file:
            print(file)
            tempOutDir=file.split(".")[0]
            print(os.path.join(cufflinksOutDir,tempOutDir))
            createDirectory(os.path.join(cufflinksOutDir,tempOutDir))
            tempcmd=cmd+os.path.join(cufflinksOutDir,tempOutDir)+" "+os.path.join(sortedSamFilesDir,file)
            print(tempcmd)
            runsubprocess(tempcmd)
            

#########################################
#cuffmerge to merge gtf files generated
#########################################
cufflinksOutDir=input("enter the cufflinks output directory : ")
extn="transcripts.gtf"
nThreads=input("enter the number of threads for cuffmerge : ")
refGenomeFasta=input("enter the reference genome fasta path : ")
refGTF=input("enter the reference genome gtf path : ")
cuffmergeOutDir=input("enter the cuffmerge ouput dir : ")
allTranscriptsGTF=[]
for r, d, f in os.walk(cufflinksOutDir):
    for file in f:
        if extn in file:
            print(os.path.join(r,file))
            allTranscriptsGTF.append(os.path.join(r,file))

assembliesFileName="assemblies.txt"
assembliesFile=open(os.path.join(cufflinksOutDir,assembliesFileName),"w")
for filePath in allTranscriptsGTF:
    assembliesFile.write(filePath+"\n")
assembliesFile.close() 
cmd="cuffmerge -g "+refGTF+" -s "+refGenomeFasta+" -p "+nThreads+" -o "+cuffmergeOutDir+" "+os.path.join(cufflinksOutDir,assembliesFileName)
print(cmd)
runsubprocess(cmd)


##########################
#htseq-count
##########################

cmd = "htseq-count --help"
runsubprocess(cmd)


############################
#STAR aligner alignment
############################
starIndexDir="/home/cidr/Documents/work/test_analysis/homoSapien/forSTAR "
nThreads=25
readFilesDir="/media/cidr/d7416dce-cdf6-43ed-9df7-978f8a9438d8/databackup_2019_07_11/work/asma_whole_blood_rnaseq/asma_durbar_analysis_raw_reads_26_jun_2020/trimmed_reads"
cmd="STAR --genomeLoad LoadAndExit --genomeDir "+starIndexDir
print(cmd)
runsubprocess(cmd)

readFiles=[]
prefix=" --outFileNamePrefix "
cmd = "STAR --genomeLoad LoadAndKeep --genomeDir "+starIndexDir+" --runThreadN "+str(nThreads)+" --readFilesIn "
starOutDir="/media/cidr/d7416dce-cdf6-43ed-9df7-978f8a9438d8/databackup_2019_07_11/work/asma_whole_blood_rnaseq/asma_durbar_analysis_raw_reads_26_jun_2020/star_results"
extn="SCR"
i=0

for r, d, f in os.walk(readFilesDir):
    for file in f:
        print(file)
        tempPrefix=os.path.join(starOutDir,"sam_"+file.split(".")[0]+"_")
        tempcmd=cmd+os.path.join(r,file)+prefix+tempPrefix
        print(tempcmd)
        runsubprocess(tempcmd)

######################################################################################
#run the code below if there is a selection based on extension while aligning with STAR. look at the extn variable
######################################################################################
for r, d, f in os.walk(readFilesDir):
    for file in f:
        if extn in file:
            #if(i==0):
              #tempcmd=cmd+" --genomeLoad LoadAndKeep"
            print(file)
            tempPrefix=os.path.join(starOutDir,"sam_"+file.split(".")[0]+"_")
            tempcmd=cmd+os.path.join(r,file)+prefix+tempPrefix
            print(tempcmd)
            #runsubprocess(tempcmd)


###################################################
#for 2 hrs
###################################################
readFilesDir="/home/cidr/Documents/work/SO_7352_8152/trimmed/igra_negative/2HR/"

readFiles=[]
prefix=" --outFileNamePrefix "
cmd = "STAR --genomeLoad LoadAndKeep --genomeDir "+starIndexDir+" --runThreadN "+str(nThreads)+" --readFilesIn "
starOutDir="/media/cidr/d7416dce-cdf6-43ed-9df7-978f8a9438d8/databackup_2019_07_11/work/trimmedSamFiles/staroutput"
extn="SCR"
i=0


for r, d, f in os.walk(readFilesDir):
    for file in f:
        if extn in file:
            #if(i==0):
              #tempcmd=cmd+" --genomeLoad LoadAndKeep"
            print(file)
            tempPrefix=os.path.join(starOutDir,"sam_"+file.split(".")[0]+"_")
            tempcmd=cmd+os.path.join(r,file)+prefix+tempPrefix
            print(tempcmd)
            runsubprocess(tempcmd)

readFilesDir="/home/cidr/Documents/work/SO_7352_8152/trimmed/igra_negative/2HR/"


###################################################
#for 24 hrs
###################################################
readFilesDir="/home/cidr/Documents/work/SO_7352_8152/trimmed/igra_negative/24HR/"

readFiles=[]
prefix=" --outFileNamePrefix "
cmd = "STAR --genomeLoad LoadAndKeep --genomeDir "+starIndexDir+" --runThreadN "+str(nThreads)+" --readFilesIn "
starOutDir="/media/cidr/d7416dce-cdf6-43ed-9df7-978f8a9438d8/databackup_2019_07_11/work/trimmedSamFiles/staroutput"
extn="SCR"
i=0


for r, d, f in os.walk(readFilesDir):
    for file in f:
        if extn in file:
            #if(i==0):
              #tempcmd=cmd+" --genomeLoad LoadAndKeep"
            print(file)
            tempPrefix=os.path.join(starOutDir,"sam_"+file.split(".")[0]+"_")
            tempcmd=cmd+os.path.join(r,file)+prefix+tempPrefix
            print(tempcmd)
            runsubprocess(tempcmd)

############################
#STAR aligner alignment for B3 Samples
############################
starIndexDir="/home/cidr/Documents/work/test_analysis/homoSapien/forSTAR "
nThreads=25
readFilesDir="/home/cidr/Documents/work/SO_7352_8152/trimmed/igra_negative/0HR/"
#cmd="STAR --genomeLoad LoadAndExit --genomeDir "+starIndexDir
#print(cmd)
#runsubprocess(cmd)

readFiles=[]
prefix=" --outFileNamePrefix "
cmd = "STAR --genomeLoad LoadAndKeep --genomeDir "+starIndexDir+" --runThreadN "+str(nThreads)+" --readFilesIn "
starOutDir="/media/cidr/d7416dce-cdf6-43ed-9df7-978f8a9438d8/databackup_2019_07_11/work/trimmedSamFiles/staroutput"
extn="B3"
i=0


for r, d, f in os.walk(readFilesDir):
    for file in f:
        if extn in file:
            #if(i==0):
              #tempcmd=cmd+" --genomeLoad LoadAndKeep"
            print(file)
            tempPrefix=os.path.join(starOutDir,"sam_"+file.split(".")[0]+"_")
            tempcmd=cmd+os.path.join(r,file)+prefix+tempPrefix
            print(tempcmd)
            runsubprocess(tempcmd)


###################################################
#for 2 hrs
###################################################
readFilesDir="/home/cidr/Documents/work/SO_7352_8152/trimmed/igra_negative/2HR/"

readFiles=[]
prefix=" --outFileNamePrefix "
cmd = "STAR --genomeLoad LoadAndKeep --genomeDir "+starIndexDir+" --runThreadN "+str(nThreads)+" --readFilesIn "
starOutDir="/media/cidr/d7416dce-cdf6-43ed-9df7-978f8a9438d8/databackup_2019_07_11/work/trimmedSamFiles/staroutput/B3/2HR"
extn="B3"
i=0


for r, d, f in os.walk(readFilesDir):
    for file in f:
        if extn in file:
            #if(i==0):
              #tempcmd=cmd+" --genomeLoad LoadAndKeep"
            print(file)
            tempPrefix=os.path.join(starOutDir,"sam_"+file.split(".")[0]+"_")
            tempcmd=cmd+os.path.join(r,file)+prefix+tempPrefix
            print(tempcmd)
            runsubprocess(tempcmd)

###################################################
#for 24 hrs
###################################################
readFilesDir="/home/cidr/Documents/work/SO_7352_8152/trimmed/igra_negative/24HR/"

readFiles=[]
prefix=" --outFileNamePrefix "
cmd = "STAR --genomeLoad LoadAndKeep --genomeDir "+starIndexDir+" --runThreadN "+str(nThreads)+" --readFilesIn "
starOutDir="/media/cidr/d7416dce-cdf6-43ed-9df7-978f8a9438d8/databackup_2019_07_11/work/trimmedSamFiles/staroutput"
extn="B3"
i=0


for r, d, f in os.walk(readFilesDir):
    for file in f:
        if extn in file:
            #if(i==0):
              #tempcmd=cmd+" --genomeLoad LoadAndKeep"
            print(file)
            tempPrefix=os.path.join(starOutDir,"sam_"+file.split(".")[0]+"_")
            tempcmd=cmd+os.path.join(r,file)+prefix+tempPrefix
            print(tempcmd)
            runsubprocess(tempcmd)
            
cmd="STAR --genomeLoad Remove --genomeDir "+starIndexDir
runsubprocess()



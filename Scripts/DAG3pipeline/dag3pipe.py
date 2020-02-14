# ===================================
# MG data analysis pipeline
# by R.Gacesa, UMCG (2019)
# ===================================
#
# includes:
# -> raw data qc (fastQC -> kneaddata -> fastQC)
# -> profiling (metaphlan, humann)
# -> targeted analysis:
#       - virulence factors (VFDB)
#       - resistance genes (CARD (WIP), ResFinder DB)
# -> in development:
#       - strain profiling (strainphlan)
#
# -> job management:
#       - data finding
# ==================================

# helper functions
# ==================================
# buffered counter for read number
def bufcount(filename):
    f = open(filename)
    lines = 0
    buf_size = 1024 * 1024
    read_f = f.read
    buf = read_f(buf_size)
    while buf:
        lines += buf.count('\n')
        buf = read_f(buf_size)
    return lines


# imports
# =================================
import copy
import argparse
import os
import csv
import shutil
import glob
import ConfigParser
import subprocess

# ===================================
#              MAIN
# ===================================

# process CL args
print ' ---------------------------------'
print '    DAG3 pipeline (V2) invoked'
print ' ---------------------------------'
parser = argparse.ArgumentParser()
parser.add_argument("cmd",
   help=' command for pipeline, as follows:\n' +
    'finddata    : recursively finds fastq files (or fq.gz) in input folder (--tgt) as root, writes _datalist.tsv in current folder; REQUIRED for writejobs\n' +
    'checkjobs  : checks which input files in _datalist.tsv are done by checking (--tgt) for results\n' +
    'writejobs   : prepares jobs for _datalist.tsv\n' +
    'sortresults : \n' +
    'saveresults : \n' +
    'cleantmp    : \n' +
    'cleanall    : cleans all temporary files\n' +
    'deljobs     : removes all data associated with joblist file [requires datafile input]\n')
parser.add_argument('tgt',help='where [def= current folder (".") ]',default='.',nargs='?')
parser.add_argument('--dataFile',help='which data file to use', default='')
parser.add_argument('--nFiles',help='how many data file entries use', default='-1')
parser.add_argument('--jobsAllSteps',help='if invoked, writes all steps when making jobs, by default writes only needed steps and skips jobs for which results already exist')
parser.add_argument('--outPath',help='where to put results',default='.')
parser.add_argument('--inPath',help='where are results',default='.')
parser.add_argument('--jobConfig',help='master config for making jobs',default='dag3pipe.cfg')

args = parser.parse_args()
#if args.cmd not in set(['finddata','writejobs','sortresults','saveresults','cleantmp','cleanall']):
#    print (' -> ERROR: command <',args.cmd,'> not recognised: run with --help for details')
#    exit(1)
cmd = args.cmd
tgt = args.tgt
dfile = args.dataFile

print ' -> running command <',cmd,'> on target <',tgt,'> and datafile <',dfile,'>'

# ======================================
# Tmp file cleaner:
# > goes through <tgt> and cleans stuff
# -> does not clean tmp files of jobs in progress
# -> otherwise:
# --> cleans humann temp
# --> cleans clean_reads
# --> cleans filtered reads
# ======================================
if args.cmd == 'cleantmp':
    print ' > Executing "cleantmp" command: cleaning temporary files in ',tgt
    print ' --> looking for jobs in progress'
    jlist = subprocess.check_output('sacct --format="JobName%30" -n --allusers', shell=True)
    jhesh = set()
    for j in jlist.split('\n'):
        j = j.strip()
        if not j == '': jhesh.add(j)
    for l in os.listdir(tgt):
        if os.path.isdir(l):
            runna = False
            for j in jhesh:
                if l.strip() in j:
                    print ' W: ',l,'is running (',j,')'
                    runna = True
                    break
            if not runna:
                print ' > ',l,'not running, cleaning'
                if os.path.isdir(l+'/clean_reads'): os.system('rm -r '+l+'/clean_reads')
                if os.path.isdir(l+'/filtering_data'): os.system('rm -r '+l+'/filtering_data')
                if os.path.isdir(l+'/humann2'):
                    for hl in os.listdir(l+'/humann2'):
                        if 'temp' in hl:
                            os.system('rm -r '+l+'/humann2/'+hl)

# =======================================
# Job deleter (deljobs):
# =======================================
if args.cmd == 'deljobs':
    print ' > Executing "deljobs" command'
    if tgt == '.':
        print ' ERROR: deljobs must be executed on datafile (usually _done_hu.tsv)'
        exit(1)
    if not os.path.isfile(tgt):
        print ' ERROR: input file',tgt,'does not exist!'
        exit(1)
    print ' --> preparing delete list (',tgt,')'
    toKill = []
    toKillRel = []
    with open(tgt) as iF:
        rdr = csv.reader(iF,delimiter='\t',quotechar='"')
        for l in rdr:
            toKill.append(l[0])
            toKillRel.append(l[0])
    print 'This command will delete ',len(toKill),'files listed in ',tgt,'!'
    print ' and related .sh, .err and .out files in . '
    cnf = raw_input (' Confirm by entering Y: ')
    if cnf == 'Y':
        for k in toKill:
            print ' -> deleting',k
            if os.path.isfile(k):
                os.remove(k)
            if os.path.isfile(k+'_1.fq.gz'):
                os.remove(k+'_1.fq.gz')
            if os.path.isfile(k+'_2.fq.gz'):
                os.remove(k+'_2.fq.gz')
        for k in toKillRel:
            if os.path.isdir(k):
                shutil.rmtree(k)
            for fl in glob.glob(k+'*.sh'):
                os.remove(fl)
            for fl in glob.glob(k+'*.out'):
                os.remove(fl)
            for fl in glob.glob(k+'*.err'):
                os.remove(fl)
# =====================================================
# Data finder:
# - does recursive search from target location (max depth = 5),
# looks for fastq/fq, .fq.gz / fastq.gz, / .bam; also check if samples are paired or not
# - paired samples are assumed to end with _1 & _2 of same prefix
# - writes output in csv form with location, type (fastq/zip/bam), paired(T/F)
# - default output is "_datalist.csv"
# =====================================================
if args.cmd == 'finddata':
    print ' > Executing "finddata" command'
    print ' --> looking for fastq and/or bam files in',tgt
    dataFnd = [] # each entry will be [name,absolute path & name, type, paired ]
    for root, dirs, files in os.walk(tgt):
        for fl in files:
            aPathFl = os.path.abspath(root)+'/'+fl
            # zipped fastq
            if fl.endswith('.fq.gz') or fl.endswith('.fastq.gz'):
                if fl.endswith('.fq.gz'): fl = fl[:-6]
                elif fl.endswith('.fastq.gz'): fl=fl[:-9]
                dataFnd.append( [fl,aPathFl,"fastq.gz",None] )
            # fastq
            elif fl.endswith('.fastq') or fl.endswith('.fq'):
                if fl.endswith('.fq'): fl = fl[:-3]
                elif fl.endswith('.fastq'): fl=fl[:-6]
                dataFnd.append( [fl,aPathFl,"fastq",None] )
            # bam
            elif fl.endswith('.bam'):
                fl = fl[:-4]
                dataFnd.append( [fl,aPathFl,"bam",None] )
    print '   --> found',len(dataFnd),'files'
    if len(dataFnd) > 0:
       print ' --> identifying paired samples'
       dataFnd = sorted(dataFnd)
       dic = {}
       for l in dataFnd:
           lN = l[0]
           if lN.endswith('_1') or lN.endswith('_2'):
              try:
                 dic[lN[:-2]] += 1
              except:
                 dic[lN[:-2]] = 1
       for l in dataFnd:
           lN = l[0]
           if (lN.endswith('_1') or lN.endswith('_2')) and dic[lN[:-2]] == 2:
              l[3] = True
              l[0] = l[0][:-2]
           else:
              l[3] = False

    if len(dataFnd) > 0:
       print ' --> saving out to _datalist.tsv'
       with open('_datalist.tsv','w') as oF:
           writ = csv.writer(oF,delimiter='\t',quotechar='"',quoting=csv.QUOTE_MINIMAL)
           for l in dataFnd: writ.writerow(l)

# =========================================================
# Find done:
# - takes target (either folder or datafile) and compares to target folder(s) to identify
# a) which samples have been processed [def output: _done_mp.tsv, _done_hu.tsv]
# b) which samples have not been processed [def output: _ndone_mp.tsv,ndone_hu.tsv']
# NYI: b) which step was done for which sample
# ========================================================
if args.cmd == 'checkjobs':
    print ' > Executing "checkjobs" command'
    inPath = args.inPath
    print ' --> identifying which samples in <',inPath,'> have results in <',tgt,'>'
    smpls = set()
    smpls_all = {}
    dFileIn = False
    if os.path.isfile(inPath):
        with open(inPath) as iF:
            dFileIn = True
            rdr = csv.reader(iF,delimiter='\t',quotechar='"')
            for l in rdr: smpls.add(l[0]); smpls_all[l[0]] = l
    elif os.path.isdir(inPath):
        for l in os.listdir(inPath):
            if os.path.isdir(l):
                smpls.add(l)

    if len(smpls) == 0:
        print ' ERROR: no samples found!'
        exit(1)
    print '   --> found',len(smpls),'samples in',inPath
    print ' --> scanning ',tgt
    #dataFnd = [] # each entry will be [name,absolute path & name, type, paired ]
    humannDone = set()
    metaPhlanDone = set()
    humannND = set()
    metaPhlanND = set()
    for root, dirs, files in os.walk(tgt):
        for fl in files:
            aPathFl = os.path.abspath(root)+'/'+fl
            # zipped fastq
            if fl.endswith('pathabundance.tsv'):
                if fl.replace('_kneaddata_merged_pathabundance.tsv','') in smpls:
                    humannDone.add(fl.replace('_kneaddata_merged_pathabundance.tsv',''))
            if fl.endswith('_metaphlan.txt'):
                if fl.replace('_metaphlan.txt','') in smpls:
                    metaPhlanDone.add(fl.replace('_metaphlan.txt',''))
    humannND = smpls - humannDone
    metaPhlanND = smpls - metaPhlanDone

    print '  -> Humann done for ',len(humannDone),'samples [saving to _done_hu.tsv]'
    print '     Humann not done for',len(humannND),'samples [saving to _ndone_hu.tsv]'
    print '  -> Metaphlan done for ',len(metaPhlanDone),'samples [saving to _done_mp.tsv]'
    print '     Metaphlan not done for',len(metaPhlanND),'samples [saving to _ndone_mp.tsv]'
    if not dFileIn:
        with open('_done_hu.tsv','w') as oF:
            for s in humannDone: oF.write(s+'\n')
        with open('_ndone_hu.tsv','w') as oF:
            for s in humannND:  oF.write(s+'\n')
        with open('_done_mp.tsv','w') as oF:
            for s in metaPhlanDone:  oF.write(s+'\n')
        with open('_ndone_mp.tsv','w') as oF:
            for s in metaPhlanND:  oF.write(s+'\n')
    else:
        with open('_done_hu.tsv','w') as oF:
            writ = csv.writer(oF,delimiter='\t',quotechar='"',quoting=csv.QUOTE_MINIMAL)
            for s in humannDone: writ.writerow(smpls_all[s])
        with open('_ndone_hu.tsv','w') as oF:
            writ = csv.writer(oF,delimiter='\t',quotechar='"',quoting=csv.QUOTE_MINIMAL)
            for s in humannND: writ.writerow(smpls_all[s])
        with open('_done_mp.tsv','w') as oF:
            writ = csv.writer(oF,delimiter='\t',quotechar='"',quoting=csv.QUOTE_MINIMAL)
            for s in metaPhlanDone: writ.writerow(smpls_all[s])
        with open('_ndone_mp.tsv','w') as oF:
            writ = csv.writer(oF,delimiter='\t',quotechar='"',quoting=csv.QUOTE_MINIMAL)
            for s in metaPhlanND: writ.writerow(smpls_all[s])

# ====================================================
# copydata:
# - takes input file <example: _datalist.tsv> and copies data files to target
# - requires extra cl parameter (--dataFile)
# ==========================================
if args.cmd == 'copydata':
    if dfile == '':
        print 'ERROR: requires --dataFile argument'
        exit(1)
    print 'copying ',dfile,'to',tgt
    smpls = []
    if os.path.isfile(dfile):
        with open(dfile) as iF:
            rdr = csv.reader(iF,delimiter='\t',quotechar='"')
            for l in rdr: smpls.append(l[1])
    if len(smpls) == 0:
        print ' ERROR: ',dfile,' does not exist, check --dataFile argument!'
        exit(1)

    if args.nFiles == '-1' or int(args.nFiles) > len(smpls): nToProc = len(smpls)
    else: nToProc = int(args.nFiles)

    print '  --> loaded',len(smpls),'entries; copying',nToProc,'files to',tgt
    c = 0
    if not os.path.isdir(tgt):
        print 'ERROR: ',tgt,'folder does not exist!'
        exit(1)
    for s in smpls:
        c += 1
        if c > nToProc: break
        #print s
        if os.path.isfile(tgt+'/'+os.path.basename(s)):
            print 'WARNING:',s,'is already in',tgt,'skipping!'
            continue
        if os.path.isfile(s):
            print 'copying',s,'->',tgt
            shutil.copy2(s, tgt)
        else:
            print 'WARNING:',s,'does not exist!'

# ==================================================
# Result mover:
# - looks for completed jobs in <tgt> [def = .]
# - moves all completed samples to --outPath
# ===================================================
if args.cmd == 'movedone':
    outPath = args.outPath
    if outPath == '.':
        print 'ERROR: incorrect parameters: call with <tgt> --outPath; <tgt> = where to find results, --outPath = where to put results'
        exit(1)
    if not os.path.isdir(outPath):
        print 'ERROR: folder',outPath,'does not exist!'
        exit(1)
    print ' > Executing "movedone" command'
    print ' --> searching for samples in',tgt
    sFnd = []
    sDone = []
    for l in os.listdir(tgt):
        if os.path.isdir(l) and os.path.isdir(l+'/humann2') and os.path.isdir(l+'/metaphlan'):
             sFnd.append(l)
             humannDone = False
             metaDone = False
             for ll in os.listdir(l+'/humann2'):
                 if '_pathabundance.tsv' in ll: humannDone = True
             for ll in os.listdir(l+'/metaphlan'):
                 if '_metaphlan.txt' in ll: metaDone = True
             if metaDone and humannDone:
                 sDone.append(l)
    print '   --> found',len(sFnd),'samples; ',len(sDone),'have humann and metaphlan done'
    print ' --> moving done samples to',outPath
    for s in sDone:
        print 'moving',tgt+'/'+s,'->',outPath+'/'+s
        shutil.move(tgt+'/'+s,outPath+'/'+s)


# ============================================================
# Result sorter:
# - looks for pieces of completed jobs in <tgt> [def = .]
# - puts results into appropriate subfolders under --outPath
# example: all metaphlan results are put into outPath/metaphlan
# =============================================================
if args.cmd == 'sortresults':
    outPath = args.outPath
    outFolder = outPath
    inFolder = tgt
    print ' > Executing "sortresults" command'
    if outPath == '.':
        print 'ERROR: incorrect parameters: call with <tgt> --outPath; <tgt> = where to find results, --outPath = where to put results'
        exit(1)
    if not os.path.isdir(outPath):
        print 'ERROR: folder',outPath,'does not exist!'
        exit(1)
    print ' --> output folders, making new ones as necessary'
    if not os.path.isdir(outFolder+'/humann2'):
        os.mkdir(outFolder+'/humann2')
    if not os.path.isdir(outFolder+'/humann2/gene_families'):
        os.mkdir(outFolder+'/humann2/gene_families')
    if not os.path.isdir(outFolder+'/humann2/path_abundances'):
        os.mkdir(outFolder+'/humann2/path_abundances')
    if not os.path.isdir(outFolder+'/humann2/path_coverage'):
        os.mkdir(outFolder+'/humann2/path_coverage')
    if not os.path.isdir(outFolder+'/qc_preclean'):
        os.mkdir(outFolder+'/qc_preclean')
    if not os.path.isdir(outFolder+'/qc_postclean'):
        os.mkdir(outFolder+'/qc_postclean')
    if not os.path.isdir(outFolder+'/metaphlan'):
        os.mkdir(outFolder+'/metaphlan')
    if not os.path.isdir(outFolder+'/logs'):
        os.mkdir(outFolder+'/logs')
    print ' --> Preparing to copy files'
    c = 0
    for l in os.listdir(inFolder):
        c+=1
        print '   --> sorting ',l
        # sort metaphlan results
        for f in glob.glob(inFolder+'/'+l+'/metaphlan/*'):
            shutil.copy2(f, outFolder+'/metaphlan')
        # sort qc_preclean
        for f in glob.glob(inFolder+'/'+l+'/qc_preclean/*'):
            shutil.copy2(f, outFolder+'/qc_preclean')
        # sort qc_postclean
        for f in glob.glob(inFolder+'/'+l+'/qc_postclean/*'):
            shutil.copy2(f, outFolder+'/qc_postclean')
        # sort humann gene families
        for f in glob.glob(inFolder+'/'+l+'/humann2/*_genefamilies.tsv'):
            shutil.copy2(f, outFolder+'/humann2/gene_families')
        # sort humann gene pathway abundances
        for f in glob.glob(inFolder+'/'+l+'/humann2/*_pathabundance.tsv'):
            shutil.copy2(f, outFolder+'/humann2/path_abundances')
        # sort humann gene pathway coverage
        for f in glob.glob(inFolder+'/'+l+'/humann2/*_pathcoverage.tsv'):
            shutil.copy2(f, outFolder+'/humann2/path_coverage')
        # logs
        for f in glob.glob(inFolder+'/'+l+'/*.log'):
            shutil.copy2(f, outFolder+'/logs')
    print 'ALL DONE, sorted through ',c,'results!'

# =====================================================================
# Job writer:
# - writes a job for each input file in --inPath which might have output files in --outPath
#    -> if --inPath points to file, then it parses through this as if it was _datafile.tsv
#    -> if it points to folder, then it makes jobs for fastq/fq OR fastq.gz/fq.gz OR .bam (pairs) in folder
# TODO: - job files go to --jobs
# TODO: - results are put into --results
# TODO: - by default only creates steps that are needed unless overridden
# TODO: - input fastq/bam/fastq.gz files are NOT copied unless --jobCopyFiles is 1 (note: sbatch cannot copy from prm,
# so they HAVE TO BE copied locally for it to work)

# NOTES:
#  IN-DEV: job writing
#  fastq PE
#  fastq.gz PE
#  TODO: bam
#  TODO: check for which pieces to write and which results exist!
# ======================================================================
if args.cmd == 'writejobs':
    print ' > Executing "writejobs" command'
# check for config
    print '  --> loading config file',args.jobConfig
    if not os.path.isfile(args.jobConfig):
        print " ERROR: config file not found! make sure --jobConfig points to config file"
        exit(1)
    else:
        cfg = ConfigParser.RawConfigParser()
        cfg.read(args.jobConfig)

# here taar be dragon (samples actually):
    dataFnd = [] # each entry will be [name,absolute path & name, type, paired ]

# check for input: input is data file (generated by writejobs)
    if os.path.isfile(tgt):
        print ' --> searching for samples in input file ',tgt
        with open(tgt) as iF:
            rdr = csv.reader(iF,delimiter='\t',quotechar='"')
            for row in rdr: dataFnd.append(row)

# input is folder instead: make jobs for all samples in folder
    elif os.path.isdir(tgt):
        print ' --> searching for samples in input folder ',tgt
        for fl in os.listdir(tgt):
            aPathFl = os.path.abspath(fl)
            #print aPathFl
            # zipped fastq
            if fl.endswith('.fq.gz'):
                if fl.endswith('.fq.gz'): fl = fl[:-6]
                dataFnd.append( [fl,aPathFl,"fq.gz",None] )
            elif fl.endswith('.fastq.gz'):
                if fl.endswith('.fastq.gz'): fl=fl[:-9]
                dataFnd.append( [fl,aPathFl,"fastq.gz",None] )
            # fastq
            elif fl.endswith('.fastq'):
                fl=fl[:-6]
                dataFnd.append( [fl,aPathFl,"fastq",None] )
            elif fl.endswith('.fq'):
                fl = fl[:-3]
                dataFnd.append( [fl,aPathFl,"fq",None] )
            # bam
            elif fl.endswith('.bam'):
                fl = fl[:-4]
                dataFnd.append( [fl,aPathFl,"bam",None] )

        if len(dataFnd) > 0:
            #print ' --> identifying paired samples'
            dataFnd = sorted(dataFnd)
            dic = {}
            for l in dataFnd:
                lN = l[0]
                if lN.endswith('_1') or lN.endswith('_2'):
                    try: dic[lN[:-2]] += 1
                    except: dic[lN[:-2]] = 1
            for l in dataFnd:
                lN = l[0]
                if (lN.endswith('_1') or lN.endswith('_2')) and dic[lN[:-2]] == 2:
                    l[3] = True
                    l[0] = l[0][:-2]
                else:
                    l[3] = False
    else:
        print ' ERROR: something is wrong with input; check <tgt> cl parameter'
        exit()

    # sample list is ready, bit of cleaning
    dataUn = {} # job -> paired or not
    dataFmt = {} # job -> format
    nrPairs = 0.0
    for s in dataFnd:
        dataUn[s[0]] = False
        dataFmt[s[0]] = s[2]
        if s[3]:
            nrPairs+=0.5
            dataUn[s[0]] = True

    print '   --> found',len(dataFnd),'files [',len(dataUn),'unique samples / ',nrPairs,'pairs]'
    #for f in dataFnd:
    #   print f
    #print dataUn

    # iterate over unique samples and write jobs:
    for smpl,paired in dataUn.items():
        fmt = dataFmt[smpl]
        print '  --> Writing job for',smpl,'; paired:',paired,'; format',fmt
        # if sample is paired, extract pairs
        if paired or fmt == 'bam':
            if dataFmt[smpl] == 'bam':
                bamFile = smpl+'.bam'
                pair1 = smpl+'_1.fastq'
                pair2 = smpl+'_2.fastq'
            elif dataFmt[smpl] == 'fq.gz':
                pair1 = smpl+'_1.fq.gz'
                pair2 = smpl+'_2.fq.gz'
            elif dataFmt[smpl] == 'fastq.gz':
                pair1 = smpl+'_1.fastq.gz'
                pair2 = smpl+'_2.fastq.gz'
            elif dataFmt[smpl] == 'fastq':
                pair1 = smpl+'_1.fastq'
                pair2 = smpl+'_2.fastq'
            elif dataFmt[smpl] == 'fq':
                pair1 = smpl+'_1.fq'
                pair2 = smpl+'_2.fq'
            pair1o = str(pair1)
            pair2o = str(pair2)
            print smpl,pair1,pair2
            # write paired fastq jobs
            # part 1: this is kneaddata job, should always be included
            with open(smpl+'_p1.sh','w') as oF:
                oF.write('#!/bin/bash\n')
                oF.write('#SBATCH --job-name=kn_'+smpl+'\n')
                oF.write('#SBATCH --error=__kn_'+smpl+'.err\n')
                oF.write('#SBATCH --output=__kn_'+smpl+'.out\n')
                oF.write('#SBATCH --mem='+cfg.get('KNEAD','memory')+'\n')
                oF.write('#SBATCH --time='+cfg.get('KNEAD','time')+'\n')
                oF.write('#SBATCH --cpus-per-task='+cfg.get('KNEAD','cpus')+'\n')
                oF.write('#SBATCH --open-mode=truncate\n')
                oF.write('# --- LOAD MODULES --- \n')
                for m in cfg.get('KNEAD','modules').split(','):
                    oF.write('module load '+m+'\n')
                if fmt == 'bam':
                    for m in cfg.get('KNEAD','modulesBAM').split(','):
                        oF.write('module load '+m+'\n')
                oF.write('# --- MAKE FOLDERS ---- \n')
                oF.write('mkdir '+smpl+'\n')
                if cfg.get('PIPELINE','doQC') == '1':
                    oF.write('mkdir ./'+smpl+'/qc_preclean \n')
                    oF.write('mkdir ./'+smpl+'/qc_postclean \n')
                oF.write('mkdir ./'+smpl+'/filtering_data \n')
                oF.write('mkdir ./'+smpl+'/clean_reads \n')
                if cfg.get('PIPELINE','doHumann') == '1':
                    oF.write('mkdir ./'+smpl+'/humann2 \n')
                if cfg.get('PIPELINE','doMeta') == '1' and not cfg.get('PIPELINE','doMetaAndStrainPhlan') == '1':
                    oF.write('mkdir ./'+smpl+'/metaphlan \n')
                if cfg.get('PIPELINE','doMetaAndStrainPhlan') == '1':
                    oF.write('mkdir ./'+smpl+'/metaphlan \n')
                    oF.write('mkdir ./'+smpl+'/strainphlan \n')
                if cfg.get('PIPELINE','doVFDB') == '1':
                    oF.write('mkdir ./'+smpl+'/VFDB \n')
                if cfg.get('PIPELINE','doResfinder') == '1':
                    oF.write('mkdir ./'+smpl+'/ResFinder \n')
                if cfg.get('PIPELINE','doVFDB_SB') == '1':
                    oF.write('mkdir ./'+smpl+'/VFDB_SB \n')
                if cfg.get('PIPELINE','doResfinderSB') == '1':
                    oF.write('mkdir ./'+smpl+'/ResFinderSB \n')
                if cfg.get('PIPELINE','doCARD_SB') == '1':
                    oF.write('mkdir ./'+smpl+'/CARD_SB \n')
                if cfg.get('PIPELINE','doGrowth') == '1':
                    oF.write('mkdir ./'+smpl+'/PTR \n')
                # BAM to FASTQ
                if fmt == 'bam':
                    #print bamFile
                    oF.write('#--- BAM to FastQ conversion ---- \n')
                    oF.write('java -jar ${EBROOTPICARD}/picard.jar SamToFastq I='+bamFile+' F='+smpl+'/filtering_data/'+pair1+' F2='+smpl+'/filtering_data/'+pair2+'\n')
                # PRE-QC
                if cfg.get('PIPELINE','doQC') == '1':
                    oF.write('# --- SUBMIT FASTQC JOB (PRE-knead) ---- \n')
                    oF.write('echo "Submitting FastQC (pre-kneaddata)" \n')
                    oF.write('sbatch '+smpl+'_q1.sh \n')

                if cfg.get('PIPELINE','doKnead') == '1':
                    # do trim galore for nextera
                    if cfg.get('KNEAD','dotrimgalore') == '1':
                       oF.write('# --- LOAD MODULES --- \n')
                       for m in cfg.get('KNEAD','modulesTrimGalore').split(','):
                          oF.write('module load '+m+'\n')
                       oF.write('# --- RUN trim-galore / remove nextera --- \n')
                       oF.write('echo "Running Trim-Galore" \n')                       
                       if not fmt == 'bam':
                          oF.write(cfg.get('KNEAD','trimgalore')+' --nextera --output_dir '+smpl+'/filtering_data/'+ ' --dont_gzip --paired '+pair1+' '+pair2+'\n')
                          pair1 = pair1.replace('.fq.gz','.fq')
                          pair2 = pair2.replace('.fq.gz','.fq')
                          oF.write('mv '+smpl+'/filtering_data/'+pair1.replace('_1.','_1_val_1.')+' '+smpl+'/filtering_data/'+pair1+'\n')
                          oF.write('mv '+smpl+'/filtering_data/'+pair2.replace('_2.','_2_val_2.')+' '+smpl+'/filtering_data/'+pair2+ '\n')
                    # converted bam
                       if fmt == 'bam':
                          oF.write(cfg.get('KNEAD','trimgalore')+' --nextera --paired '+smpl+'/filtering_data/'+pair1+' '+smpl+'/filtering_data/'+pair2+' -o '+smpl+'/filtering_data/'+'\n')
                          pair1 = pair1.replace('.fastq','.fq')
                          pair2 = pair2.replace('.fastq','.fq')
                          oF.write('mv '+ smpl+'/filtering_data/'+pair1.replace('_1.','_1_val_1.') +' '+ smpl+'/filtering_data/'+pair1+'\n')
                          oF.write('mv '+ smpl+'/filtering_data/'+pair2.replace('_2.','_2_val_2.') +' '+ smpl+'/filtering_data/'+pair2+'\n')

                    # done with nextera, do knead
                    # ===============================================
                    oF.write('# --- RUN KNEADDATA ---- \n')  # PAIRED END KNEAD
                    oF.write('# --- LOAD MODULES --- \n')
                    # module load
                    for m in cfg.get('KNEAD','modules').split(','):
                        oF.write('module load '+m+'\n')
                    if fmt == 'bam':
                        for m in cfg.get('KNEAD','modulesBAM').split(','):
                            oF.write('module load '+m+'\n')
                    # run knead
                    oF.write('echo "Running Kneaddata" \n')
                    # paired-end fastqs
                    if not fmt == 'bam' and not cfg.get('KNEAD','dotrimgalore') == '1':                        
                        oF.write('kneaddata --input '+pair1+' --input '+pair2+' --threads '+cfg.get('KNEAD','threads')+' --processes '+cfg.get('KNEAD','threads')+' --output '+smpl+'/filtering_data/ --log '+smpl+'/'+smpl+'_kneaddata.log -db '+cfg.get('KNEAD','db')+'\n')
                    # converted bam
                    else: # fmt == 'bam' and not cfg.get('KNEAD','dotrimgalore') == '1':
                        oF.write('kneaddata --input '+smpl+'/filtering_data/'+pair1+' --input '+smpl+'/filtering_data/'+pair2+' --threads '+cfg.get('KNEAD','threads')+' --processes '+cfg.get('KNEAD','threads')+' --output '+smpl+'/filtering_data/ --log '+smpl+'/'+smpl+'_kneaddata.log -db '+cfg.get('KNEAD','db')+'\n')
                    # merge & move stuff, remove contaminants and unpaired reads
                    oF.write('#   -->  clean kneaddata results: \n')
                    if cfg.get('KNEAD','mergeResults') == '1':
                        oF.write('cat '+smpl+'/filtering_data/'+smpl+'_1_kneaddata_paired_1.fastq > '+smpl+'/filtering_data/'+smpl+'_kneaddata_merged.fastq\n')
                        oF.write('cat '+smpl+'/filtering_data/'+smpl+'_1_kneaddata_paired_2.fastq >> '+smpl+'/filtering_data/'+smpl+'_kneaddata_merged.fastq\n')
                    oF.write('mv '+smpl+'/filtering_data/*kneaddata_paired_1.fastq '+smpl+'/clean_reads\n')
                    oF.write('mv '+smpl+'/filtering_data/*kneaddata_paired_2.fastq '+smpl+'/clean_reads\n')
                    # cleanup - save human if required
                    if cfg.get('KNEAD','keephuman') == '1':
                        oF.write('mkdir '+smpl+'/human_reads \n')
                        oF.write('mv '+smpl+'/filtering_data/*Homo_sapiens*.fastq '+smpl+'/human_reads \n')
                    # cleanup - merge data if required
                    if cfg.get('KNEAD','mergeResults') == '1':
                        oF.write('mv '+smpl+'/filtering_data/*kneaddata_merged.fastq '+smpl+'/clean_reads\n')
                    
                    # cleanup - remove junk
                    if cfg.get('KNEAD','cleantmp') == '1':
                        oF.write('rm -r '+smpl+'/filtering_data\n')
                # POST-QC
                if cfg.get('PIPELINE','doQC') == '1':
                    oF.write('# --- SUBMIT FASTQC JOB (POST-knead) ---- \n')
                    oF.write('echo "Submitting FastQC job (post-knead)" \n')
                    oF.write('sbatch '+smpl+'_q2.sh \n')
                # SUBMIT Antibiotic resistance genes finder (if required)
                if cfg.get('PIPELINE','doResfinder') == '1':
                    oF.write('# --- SUBMIT Resfinder job ---- \n')
                    oF.write('echo "Submitting Resistome job (post-knead)" \n')
                    oF.write('sbatch '+smpl+'_a.sh \n')
                # SUBMIT Antibiotic resistance genes finder (shortBRED) (if required)
                if cfg.get('PIPELINE','doResfinderSB') == '1':
                    oF.write('# --- SUBMIT ResfinderSB job ---- \n')
                    oF.write('echo "Submitting Resistome SB job (post-knead)" \n')
                    oF.write('sbatch '+smpl+'_as.sh \n')
                # SUBMIT Virulence factor finder (if required)
                if cfg.get('PIPELINE','doVFDB') == '1':
                    oF.write('# --- SUBMIT Virulence finder job ---- \n')
                    oF.write('echo "Submitting Virulence finder job" \n')
                    oF.write('sbatch '+smpl+'_v.sh \n')
                # SUBMIT Virulence factor finder (shortBRED) (if required)
                if cfg.get('PIPELINE','doVFDB_SB') == '1':
                    oF.write('# --- SUBMIT Virulence finder SB job ---- \n')
                    oF.write('echo "Submitting Virulence finder SB job" \n')
                    oF.write('sbatch '+smpl+'_vs.sh \n')
                # SUBMIT CARD antibiotic resistance (shortBRED) (if required)
                if cfg.get('PIPELINE','doCARD_SB') == '1':
                    oF.write('# --- SUBMIT CARD resistome SB job ---- \n')
                    oF.write('echo "Submitting CARD SB job" \n')
                    oF.write('sbatch '+smpl+'_ac.sh \n')
                # SUBMIT METAPHLAN (if required)
                if cfg.get('PIPELINE','doMeta') == '1':
                    oF.write('# --- SUBMIT METAPHLAN JOB ---- \n')
                    oF.write('echo "Submitting Metaphlan job" \n')
                    oF.write('sbatch '+smpl+'_p2.sh \n')
                # SUBMIT PTR (if required)
                if cfg.get('PIPELINE','doGrowth') == '1':
                    oF.write('# --- SUBMIT PTR/Growth rates JOB ---- \n')
                    oF.write('echo "Submitting PTR/Growth rates job" \n')
                    oF.write('sbatch '+smpl+'_g.sh \n')

                # END OF PART 1

            # WRITE QC 1
            if cfg.get('PIPELINE','doQC') == '1':
                with open(smpl+'_q1.sh','w') as oF:
                    oF.write('#!/bin/bash\n')
                    oF.write('#SBATCH --job-name=q1_'+smpl+'\n')
                    oF.write('#SBATCH --error=__q1_'+smpl+'.err\n')
                    oF.write('#SBATCH --output=__q1_'+smpl+'.out\n')
                    oF.write('#SBATCH --mem='+cfg.get('QC','memory')+'\n')
                    oF.write('#SBATCH --time='+cfg.get('QC','time')+'\n')
                    oF.write('#SBATCH --cpus-per-task='+cfg.get('QC','cpus')+'\n')
                    oF.write('#SBATCH --open-mode=truncate\n')
                    oF.write('# --- LOAD MODULES --- \n')
                    for m in cfg.get('QC','modules').split(','):
                        oF.write('module load '+m+'\n')
                    oF.write('# --- RUN QC --- \n')
                    oF.write('echo "Running FastQC on unclean reads" \n')
                    if fmt == 'bam':
                        oF.write('fastqc -t'+cfg.get('QC','threads')+ '-q -o '+smpl+'/qc_preclean '+ smpl+'/filtering_data/'+smpl+'_1.fastq'+'\n')
                        oF.write('fastqc -t'+cfg.get('QC','threads')+ '-q -o '+smpl+'/qc_preclean '+ smpl+'/filtering_data/'+smpl+'_2.fastq'+'\n')
                    elif not fmt == 'bam':
                        oF.write('fastqc -t '+cfg.get('QC','threads')+' -q -o '+smpl+'/qc_preclean '+pair1o+'\n')
                        oF.write('fastqc -t '+cfg.get('QC','threads')+' -q -o '+smpl+'/qc_preclean '+pair2o+'\n')
            # WRITE QC 2
                with open(smpl+'_q2.sh','w') as oF:
                    oF.write('#!/bin/bash\n')
                    oF.write('#SBATCH --job-name=q2_'+smpl+'\n')
                    oF.write('#SBATCH --error=__q2_'+smpl+'.err\n')
                    oF.write('#SBATCH --output=__q2_'+smpl+'.out\n')
                    oF.write('#SBATCH --mem='+cfg.get('QC','memory')+'\n')
                    oF.write('#SBATCH --time='+cfg.get('QC','time')+'\n')
                    oF.write('#SBATCH --cpus-per-task='+cfg.get('QC','cpus')+'\n')
                    oF.write('#SBATCH --open-mode=truncate\n')
                    oF.write('# --- LOAD MODULES --- \n')
                    for m in cfg.get('QC','modules').split(','):
                        oF.write('module load '+m+'\n')
                    oF.write('# --- RUN QC --- \n')
                    oF.write('echo "Running FastQC on cleaned reads" \n')
                    oF.write('fastqc -t '+cfg.get('QC','threads')+' -q -o '+smpl+'/qc_postclean '+smpl+'/clean_reads/'+smpl+'_kneaddata_merged.fastq'+'\n')
            # WRITE PART 2 (METAPHLAN and NOT STRAINPHLAN)
            if cfg.get('PIPELINE','doMeta') == '1' and not cfg.get('PIPELINE','doMetaAndStrainPhlan') == 1:
                with open(smpl+'_p2.sh','w') as oF:
                    oF.write('#!/bin/bash\n')
                    oF.write('#SBATCH --job-name=m_'+smpl+'\n')
                    oF.write('#SBATCH --error=__m_'+smpl+'.err\n')
                    oF.write('#SBATCH --output=__m_'+smpl+'.out\n')
                    oF.write('#SBATCH --mem='+cfg.get('METAPHLAN','memory')+'\n')
                    oF.write('#SBATCH --time='+cfg.get('METAPHLAN','time')+'\n')
                    oF.write('#SBATCH --cpus-per-task='+cfg.get('METAPHLAN','cpus')+'\n')
                    oF.write('#SBATCH --open-mode=truncate\n')
                    oF.write('# --- LOAD MODULES --- \n')
                    for m in cfg.get('METAPHLAN','modules').split(','):
                        oF.write('module load '+m+'\n')
                    # metaphlan command
                    oF.write('echo "Running Metaphlan" \n')
                    # set DB (if needed)
                    mdb = ''
                    if not cfg.get('METAPHLAN','metaphlanDB') == '':
                        mdb = ' --bowtie2db '+cfg.get('METAPHLAN','metaphlanDB')
                    oF.write('# --- RUN METAHPLAN --- \n')
                    oF.write(cfg.get('METAPHLAN','metaphlan')+' '+smpl+'/clean_reads/'+smpl+'_kneaddata_merged.fastq'+' --input_type multifastq --nproc '+cfg.get('METAPHLAN','threads') + mdb + ' -o '+smpl+'/metaphlan/'+smpl+'_metaphlan.txt'+' --tmp_dir '+smpl+'/metaphlan_tmp' + ' 2>&1 | tee '+smpl+'/'+smpl+'_metaphlan.log'+'\n')
                    # SUBMIT HUMANN2 (if required)
                    if cfg.get('PIPELINE','doHumann') == '1':
                        oF.write('# --- SUBMIT HUMANN2 JOB ---- \n')
                        oF.write('echo "Submitting humann2 job" \n')
                        oF.write('sbatch '+smpl+'_p3.sh \n')
            # END OF PART 2

            # WRITE PART 2/b (METAPHLAN & STRAINPHLAN) (overrides METAPHLAN alone)
            if cfg.get('PIPELINE','doMetaAndStrainPhlan') == '1':
                # METAPHLAN
                with open(smpl+'_p2.sh','w') as oF:
                    oF.write('#!/bin/bash\n')
                    oF.write('#SBATCH --job-name=ms_'+smpl+'\n')
                    oF.write('#SBATCH --error=__ms_'+smpl+'.err\n')
                    oF.write('#SBATCH --output=__ms_'+smpl+'.out\n')
                    oF.write('#SBATCH --mem='+cfg.get('METAPHLAN','memory')+'\n')
                    oF.write('#SBATCH --time='+cfg.get('METAPHLAN','time')+'\n')
                    oF.write('#SBATCH --cpus-per-task='+cfg.get('METAPHLAN','cpus')+'\n')
                    oF.write('#SBATCH --open-mode=truncate\n')
                    oF.write('# --- LOAD MODULES --- \n')
                    for m in cfg.get('METAPHLAN','modules').split(','):
                        oF.write('module load '+m+'\n')
                    oF.write('echo "Running Metaphlan" \n')
                    # run metaphlan
                    oF.write('# --- RUN METAPHLAN --- \n')
                    # set DB (if needed)
                    mdb = ''
                    if not cfg.get('METAPHLAN','metaphlanDB') == '':
                        mdb = ' --bowtie2db '+cfg.get('METAPHLAN','metaphlanDB')
                    oF.write(cfg.get('METAPHLAN','metaphlan')+' '+smpl+'/clean_reads/'+smpl+'_kneaddata_merged.fastq'+' --input_type multifastq --nproc '+cfg.get('METAPHLAN','threads') + mdb + ' -o '+smpl+'/metaphlan/'+smpl+'_metaphlan.txt'+' --tmp_dir '+smpl+'/metaphlan_tmp' + ' --bowtie2out ' +smpl+'/metaphlan/'+smpl+'_metaphlan_bowtie2.txt' + ' --samout ' +smpl+'/metaphlan/'+smpl+'_metaphlan.sam.bz2'+' 2>&1 | tee '+smpl+'/'+smpl+'_metaphlan.log'+'\n')
                    # SUBMIT STRAINPHLAN
                    oF.write('# --- SUBMIT STRAINPHLAN JOB ---- \n')
                    oF.write('echo "Submitting strainphlan job" \n')
                    oF.write('sbatch '+smpl+'_p2sp.sh \n')
                    # SUBMIT HUMANN2 (if required - note: it can go in parallel with strainphlan)
                    if cfg.get('PIPELINE','doHumann') == '1':
                        oF.write('# --- SUBMIT HUMANN2 JOB ---- \n')
                        oF.write('echo "Submitting humann2 job" \n')
                        oF.write('sbatch '+smpl+'_p3.sh \n')
            # WRITE PART 2/sp (Strainphlan) (if required):
                # STRAINPHLAN
                with open(smpl+'_p2sp.sh','w') as oF:
                    oF.write('#!/bin/bash\n')
                    oF.write('#SBATCH --job-name=sp_'+smpl+'\n')
                    oF.write('#SBATCH --error=sp_'+smpl+'.err\n')
                    oF.write('#SBATCH --output=sp_'+smpl+'.out\n')
                    oF.write('#SBATCH --mem='+cfg.get('STRAINPHLAN','memory')+'\n')
                    oF.write('#SBATCH --time='+cfg.get('STRAINPHLAN','time')+'\n')
                    oF.write('#SBATCH --cpus-per-task='+cfg.get('STRAINPHLAN','cpus')+'\n')
                    oF.write('#SBATCH --open-mode=truncate\n')
                    oF.write('# --- LOAD MODULES --- \n')
                    for m in cfg.get('STRAINPHLAN','modules').split(','):
                        oF.write('module load '+m+'\n')
                    oF.write('echo "Running Strainphlan marker generation" \n')
                    oF.write('# --- RUN STRAINPHLAN --- \n')
                    # add bcftools and samtools if needed
                    bcf = ''
                    if not cfg.get('STRAINPHLAN','bcftools') == '':
                        bcf = ' --bcftools_exe '+cfg.get('STRAINPHLAN','bcftools')

                    stools = ''
                    if not cfg.get('STRAINPHLAN','samtools') == '':
                        stools = ' --samtools_exe '+cfg.get('STRAINPHLAN','samtools')

                    oF.write(cfg.get('STRAINPHLAN','samples2markers')+' --input_type sam'+bcf+stools+' --ifn_samples '+smpl+'/metaphlan/'+smpl+'_metaphlan.sam.bz2'+ ' --output_dir '+smpl+'/strainphlan'+' 2>&1 | tee '+smpl+'/'+smpl+'_metaphlan.log'+'\n')


            # WRITE PART 3 (HUMANN) (if required):
            if cfg.get('PIPELINE','doHumann') == '1':
                with open(smpl+'_p3.sh','w') as oF:
                    oF.write('#!/bin/bash\n')
                    oF.write('#SBATCH --job-name=h_'+smpl+'\n')
                    oF.write('#SBATCH --error=__h_'+smpl+'.err\n')
                    oF.write('#SBATCH --output=__h_'+smpl+'.out\n')
                    oF.write('#SBATCH --mem='+cfg.get('HUMANN','memory')+'\n')
                    oF.write('#SBATCH --time='+cfg.get('HUMANN','time')+'\n')
                    oF.write('#SBATCH --cpus-per-task='+cfg.get('HUMANN','cpus')+'\n')
                    oF.write('#SBATCH --open-mode=truncate\n')
                    oF.write('# --- LOAD MODULES --- \n')
                    for m in cfg.get('HUMANN','modules').split(','):
                        oF.write('module load '+m+'\n')
                    oF.write('echo "Running Humann2" \n')
                    oF.write('# --- RUN HUMANN2 --- \n')
                    oF.write(cfg.get('HUMANN','humann')+' --input '+smpl+'/clean_reads/'+smpl+'_kneaddata_merged.fastq --output '+smpl+'/humann2/ --taxonomic-profile '+smpl+'/metaphlan/'+smpl+'_metaphlan.txt --threads '+cfg.get('HUMANN','threads')+' --o-log '+smpl+'/'+smpl+'_humann2.log --remove-temp-output'+'\n')
                    oF.write('# --- CLEANUP --- \n')
                    oF.write('echo "Cleaning redundant data"\n')
                    oF.write('rm -r '+smpl+'/clean_reads\n')
                    #oF.write('rm -r '+smpl+'/filtering_data\n')
                    oF.write('echo " --> ALL DONE !!! <-- "\n')
            # END OF PART 3

            # WRITE PART A (antibiotics resistome) (if required):
            if cfg.get('PIPELINE','doResfinder') == '1':
                with open(smpl+'_a.sh','w') as oF:
                    oF.write('#!/bin/bash\n')
                    oF.write('#SBATCH --job-name=r_'+smpl+'\n')
                    oF.write('#SBATCH --error=__r_'+smpl+'.err\n')
                    oF.write('#SBATCH --output=__r_'+smpl+'.out\n')
                    oF.write('#SBATCH --mem='+cfg.get('RESFINDER','memory')+'\n')
                    oF.write('#SBATCH --time='+cfg.get('RESFINDER','time')+'\n')
                    oF.write('#SBATCH --cpus-per-task='+cfg.get('RESFINDER','cpus')+'\n')
                    oF.write('#SBATCH --open-mode=truncate\n')
                    oF.write('# --- LOAD MODULES --- \n')
                    for m in cfg.get('RESFINDER','modules').split(','):
                        oF.write('module load '+m+'\n')
                    oF.write('echo "Running ResFinder analysis"\n')
                    oF.write('# --- RUN BOWTIE vs ResFinder DB --- \n')
                    oF.write('bowtie2 --threads '+cfg.get('RESFINDER','threads')+' '+cfg.get('RESFINDER','bowtie_mode')+' -x '+cfg.get('RESFINDER','bowtie_index')+' --no-discordant --no-mixed '+'-1 '+smpl+'/clean_reads/'+smpl+'_1_kneaddata_paired_1.fastq'+' -2 '+smpl+'/clean_reads/'+smpl+'_1_kneaddata_paired_2.fastq'+' | samtools view -SF4 > '+smpl+'/ResFinder/alignment.sam'+'\n')
                    oF.write('# --- PARSE ResFinder results\n')
                    # count reads
                    oF.write("NLINES=$(grep -c '' "+smpl+'/clean_reads/'+smpl+'_1_kneaddata_paired_1.fastq'+")\n")
                    oF.write('NREADS=$(( NLINES / 4 ))\n')
                    oF.write('python '+cfg.get('RESFINDER','parser') +' --sam '+smpl+'/ResFinder/alignment.sam'+ ' --annot '+cfg.get('RESFINDER','annotation')+' --fasta '+cfg.get('RESFINDER','fasta')+' --readN $NREADS --out '+smpl+'/ResFinder/'+smpl+'_resfinder_out'+'\n')
                    # kill sam file
                    oF.write('# --- CLEANUP\n')
                    oF.write('rm '+smpl+'/ResFinder/alignment.sam\n')
            # ---- END OF ANTIBIOTICS RESISTOME (part A) ----

            # WRITE PART AS (antibiotics resistome shortBRED) (if required):
            if cfg.get('PIPELINE','doResfinderSB') == '1':
                with open(smpl+'_as.sh','w') as oF:
                    oF.write('#!/bin/bash\n')
                    oF.write('#SBATCH --job-name=rs_'+smpl+'\n')
                    oF.write('#SBATCH --error=__rs_'+smpl+'.err\n')
                    oF.write('#SBATCH --output=__rs_'+smpl+'.out\n')
                    oF.write('#SBATCH --mem='+cfg.get('RESFINDER_SB','memory')+'\n')
                    oF.write('#SBATCH --time='+cfg.get('RESFINDER_SB','time')+'\n')
                    oF.write('#SBATCH --cpus-per-task='+cfg.get('RESFINDER_SB','cpus')+'\n')
                    oF.write('#SBATCH --open-mode=truncate\n')
                    oF.write('# --- LOAD MODULES --- \n')
                    for m in cfg.get('RESFINDER_SB','modules').split(','):
                        oF.write('module load '+m+'\n')
                    oF.write('echo "Running ResFinderSB analysis"\n')
                    oF.write('# --- RUN shortBRED vs ResFinder DB Markers --- \n')
                    oF.write(cfg.get('RESFINDER_SB','shortbred')+ ' --markers '+cfg.get('RESFINDER_SB','markers')+' --wgs '+smpl+'/clean_reads/'+smpl+'_kneaddata_merged.fastq'+' --threads '+cfg.get('RESFINDER_SB','threads')+' --results '+smpl+'/ResFinderSB/'+smpl+'_resfinder_sb.txt'+' --tmp '+smpl+'/ResFinderSB/tmp'+' --usearch '+cfg.get('RESFINDER_SB','usearch')+'\n')
                    # count reads
                    oF.write("NLINES=$(grep -c '' "+smpl+'/clean_reads/'+smpl+'_1_kneaddata_paired_1.fastq'+")\n")
                    oF.write('NREADS=$(( NLINES / 4 ))\n')
                    # run parser
                    oF.write('# --- PARSE shortBRED results\n')
                    oF.write('python '+cfg.get('RESFINDER_SB','parser') +' --inFile '+smpl+'/ResFinderSB/'+smpl+'_resfinder_sb.txt'+' --annot '+cfg.get('RESFINDER_SB','annotation')+' --out '+smpl+'/ResFinderSB/'+smpl+'_resfindersb_out'+' --fasta '+cfg.get('RESFINDER_SB','fasta')+' --readN $NREADS'+'\n')
                    # CLEANUP
                    oF.write('rm -r '+smpl+'/ResFinderSB/tmp'+'\n')
            # ---- END OF ANTIBIOTICS RESISTOME/SB (part AS) ----


            # WRITE PART AC (CARD shortBRED) (if required):
            if cfg.get('PIPELINE','doCARD_SB') == '1':
                with open(smpl+'_ac.sh','w') as oF:
                    oF.write('#!/bin/bash\n')
                    oF.write('#SBATCH --job-name=cs_'+smpl+'\n')
                    oF.write('#SBATCH --error=__cs_'+smpl+'.err\n')
                    oF.write('#SBATCH --output=__cs_'+smpl+'.out\n')
                    oF.write('#SBATCH --mem='+cfg.get('CARD_SB','memory')+'\n')
                    oF.write('#SBATCH --time='+cfg.get('CARD_SB','time')+'\n')
                    oF.write('#SBATCH --cpus-per-task='+cfg.get('CARD_SB','cpus')+'\n')
                    oF.write('#SBATCH --open-mode=truncate\n')
                    oF.write('# --- LOAD MODULES --- \n')
                    for m in cfg.get('CARD_SB','modules').split(','):
                        oF.write('module load '+m+'\n')
                    oF.write('echo "Running CARD/SB analysis"\n')
                    oF.write('# --- RUN shortBRED vs CARD DB Markers --- \n')
                    oF.write(cfg.get('CARD_SB','shortbred')+ ' --markers '+cfg.get('CARD_SB','markers')+' --wgs '+smpl+'/clean_reads/'+smpl+'_kneaddata_merged.fastq'+' --threads '+cfg.get('CARD_SB','threads')+' --results '+smpl+'/CARD_SB/'+smpl+'_CARD_sb.txt'+' --tmp '+smpl+'/CARD_SB/tmp'+' --usearch '+cfg.get('CARD_SB','usearch')+'\n')
                    # CLEANUP
                    oF.write('rm -r '+smpl+'/CARD_SB/tmp'+'\n')


            # ---- END OF ANTIBIOTICS RESISTOME/SB (part AS) ----

            # WRITE PART V (virulence factors) (if required):
            if cfg.get('PIPELINE','doVFDB') == '1':
                with open(smpl+'_v.sh','w') as oF:
                    oF.write('#!/bin/bash\n')
                    oF.write('#SBATCH --job-name=v_'+smpl+'\n')
                    oF.write('#SBATCH --error=__v_'+smpl+'.err\n')
                    oF.write('#SBATCH --output=__v_'+smpl+'.out\n')
                    oF.write('#SBATCH --mem='+cfg.get('VFDB','memory')+'\n')
                    oF.write('#SBATCH --time='+cfg.get('VFDB','time')+'\n')
                    oF.write('#SBATCH --cpus-per-task='+cfg.get('VFDB','cpus')+'\n')
                    oF.write('#SBATCH --open-mode=truncate\n')
                    oF.write('# --- LOAD MODULES --- \n')
                    for m in cfg.get('VFDB','modules').split(','):
                        oF.write('module load '+m+'\n')
                    oF.write('# --- RUN BOWTIE vs VFDB --- \n')
                    oF.write('echo "Running VFDB analysis"\n')
                    oF.write('bowtie2 --threads '+cfg.get('VFDB','threads')+' '+cfg.get('VFDB','bowtie_mode')+' -x '+cfg.get('VFDB','bowtie_index')+' --no-discordant --no-mixed '+'-1 '+smpl+'/clean_reads/'+smpl+'_1_kneaddata_paired_1.fastq'+' -2 '+smpl+'/clean_reads/'+smpl+'_1_kneaddata_paired_2.fastq'+' | samtools view -SF4 > '+smpl+'/VFDB/alignment.sam'+'\n')
                    oF.write('# --- PARSE VFDB results\n')
                    # count reads
                    oF.write("NLINES=$(grep -c '' "+smpl+'/clean_reads/'+smpl+'_1_kneaddata_paired_1.fastq'+")\n")
                    oF.write('NREADS=$(( NLINES / 4 ))\n')
                    oF.write('python '+cfg.get('VFDB','parser') +' --sam '+smpl+'/VFDB/alignment.sam'+ ' --annot '+cfg.get('VFDB','annotation')+' --fasta '+cfg.get('VFDB','fasta')+' --readN $NREADS --out '+smpl+'/VFDB/'+smpl+'_VFDB_out'+'\n')
                    # kill sam file
                    oF.write('# --- CLEANUP\n')
                    oF.write('rm '+smpl+'/VFDB/alignment.sam\n')
            # ---- END OF PART V ----

            # WRITE PART VS (VFDB shortBRED) (if required):
            if cfg.get('PIPELINE','doVFDB_SB') == '1':
                with open(smpl+'_vs.sh','w') as oF:
                    oF.write('#!/bin/bash\n')
                    oF.write('#SBATCH --job-name=vs_'+smpl+'\n')
                    oF.write('#SBATCH --error=__vs_'+smpl+'.err\n')
                    oF.write('#SBATCH --output=__vs_'+smpl+'.out\n')
                    oF.write('#SBATCH --mem='+cfg.get('VFDB_SB','memory')+'\n')
                    oF.write('#SBATCH --time='+cfg.get('VFDB_SB','time')+'\n')
                    oF.write('#SBATCH --cpus-per-task='+cfg.get('VFDB_SB','cpus')+'\n')
                    oF.write('#SBATCH --open-mode=truncate\n')
                    oF.write('# --- LOAD MODULES --- \n')
                    for m in cfg.get('VFDB_SB','modules').split(','):
                        oF.write('module load '+m+'\n')
                    oF.write('echo "Running VFDB_SB analysis"\n')
                    oF.write('# --- RUN shortBRED vs VFDB DB Markers --- \n')
                    oF.write(cfg.get('VFDB_SB','shortbred')+ ' --markers '+cfg.get('VFDB_SB','markers')+' --wgs '+smpl+'/clean_reads/'+smpl+'_kneaddata_merged.fastq'+' --threads '+cfg.get('VFDB_SB','threads')+' --results '+smpl+'/VFDB_SB/'+smpl+'_vfdb_sb.txt'+' --tmp '+smpl+'/VFDB_SB/tmp'+' --usearch '+cfg.get('VFDB_SB','usearch')+'\n')
                    # count reads
                    oF.write("NLINES=$(grep -c '' "+smpl+'/clean_reads/'+smpl+'_1_kneaddata_paired_1.fastq'+")\n")
                    oF.write('NREADS=$(( NLINES / 4 ))\n')
                    # run parser
                    oF.write('# --- PARSE shortBRED results\n')
                    oF.write('python '+cfg.get('VFDB_SB','parser') +' --inFile '+smpl+'/VFDB_SB/'+smpl+'_vfdb_sb.txt'+' --annot '+cfg.get('VFDB_SB','annotation')+' --out '+smpl+'/VFDB_SB/'+smpl+'_vfdb_sb_out'+' --fasta '+cfg.get('VFDB_SB','fasta')+' --readN $NREADS'+'\n')
                    # CLEANUP
                    oF.write('rm -r '+smpl+'/VFDB_SB/tmp'+'\n')
            # ---- END OF VFDB/SB (part VS) ----

            # WRITE PART G (growth factors) (if required):
            if cfg.get('PIPELINE','doGrowth') == '1':
                with open(smpl+'_g.sh','w') as oF:
                    oF.write('#!/bin/bash\n')
                    oF.write('#SBATCH --job-name=g_'+smpl+'\n')
                    oF.write('#SBATCH --error=__g_'+smpl+'.err\n')
                    oF.write('#SBATCH --output=__g_'+smpl+'.out\n')
                    oF.write('#SBATCH --mem='+cfg.get('GROWTH','memory')+'\n')
                    oF.write('#SBATCH --time='+cfg.get('GROWTH','time')+'\n')
                    oF.write('#SBATCH --cpus-per-task='+cfg.get('GROWTH','cpus')+'\n')
                    oF.write('#SBATCH --open-mode=truncate\n')
                    oF.write('# --- LOAD MODULES --- \n')
                    for m in cfg.get('GROWTH','modules').split(','):
                        oF.write('module load '+m+'\n')
                    oF.write('# --- SET PATHS --- \n')
                    for p in cfg.get('GROWTH','paths').split(','):
                        oF.write('export PATH=$PATH:'+p+'\n')
                    oF.write('# --- RUN PTRC --- \n')
                    oF.write('echo Running PTRC CA\n')
                    oF.write('python '+cfg.get('GROWTH','ptrc')+' -pe'+ ' -i1 '+smpl+'/clean_reads/*_kneaddata_paired_1.fastq'+' -i2 '+smpl+'/clean_reads/*_kneaddata_paired_2.fastq'+' -db_path_name '+cfg.get('GROWTH','dbpath')+' -outfol '+smpl+'/PTR'+' -m '+smpl+'/PTR/mapping.map'+' CA'+'\n')
                    oF.write('echo PTRC MAPPING COMPLETE!\n')
                    # remove mapping file - they are not necessary for postprocessing
                    oF.write('echo rm '+smpl+'/PTR/mapping.map\n')
                    oF.write('echo note: all mapped files have to be bundled and processed with PTR after all samples are mapped\n')

        # if sample is unpaired, find it
        else:
            unpaired = False
            for d in dataFnd:
                if d[0] == smpl: unpaired = d; break
            # ready for job writing
            print unpaired

#    print cfg.get('DATA_QC','memory')

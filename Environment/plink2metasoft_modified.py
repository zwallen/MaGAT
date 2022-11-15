#!/usr/bin/env python

########################################################
# plink2metasoft.py                                    
#   Convert Plink .assoc files to Metasoft input file  
#   Free license -- you are free to use it in any ways 
#   Buhm Han (2012)
########################################################

############################################################
#   Modified by Zachary Wallen (ZW) to be used with        #
#   plink2 result files, and to perform any mathematical   #
#   functions previously performed separately by           #
#   plink2metasoft_subroutine.R script. Modified lines     #
#   will have comments at the end of them to denote what   #
#   was changed from original.                             #
#                                                          #
#   Last updated: 20 May 2021                              #
############################################################

import sys, subprocess, os, math  #### ZW: added 'math' to import list

comple = {'A':'T', 'T':'A', 'G':'C', 'C':'G'}

def flipstrand(string):  #### ZW: added function to use when searching for flipped strands that will handle multi-character alleles
    return ''.join(comple.get(l, l) for l in string)

# PROCESS ARGUMENTS
if len(sys.argv) < 3:
    print_and_exit()
out=sys.argv[1]
files=sys.argv[2:]

# READ FILES
studies=[]
for f in files:
    study={}
    fin=open(f)
    colnames=fin.next().split()
    for line in fin:
        snp={}
        for (x,y) in zip(colnames,line.split()):
            snp[x]=y
        rsid=snp['ID'] #### ZW: changed SNP to ID
        study[rsid]=snp
    fin.close()
    studies.append(study)
    
# UNION OF SNPS
allsnps={}
for study in studies:
    for rsid, snp in study.iteritems():
        allsnps[rsid]=snp

# SORT SNPS
rsids=sorted(allsnps.keys(), 
             key=lambda x:(allsnps[x]['CHROM'],int(allsnps[x]['POS']))) #### ZW: changed CHR and BP to CHROM and POS respectively

# MERGE STUDIES
fout=open(out+'.meta','w') #### ZW: removed '.tmp' from '.meta.tmp', no need since subroutine script not being used
fmap=open(out+'.mmap','w')
flog=open(out+'.log','w')  
for rsid in rsids:
    output=rsid+'\t'
    pivot=-1
    pivotstudyindex=-1
    numstudy=0
    for study in studies:
        studyindex=studies.index(study)+1
        if rsid in study:
            snp=study[rsid]
            if 'OR' in snp:
                beta=str(math.log(float(snp['OR']))) #'log(%s)'%snp['OR'] #### ZW: modified so beta is calculated here and not reliant on subroutine script
            elif 'BETA' in snp:
                beta=snp['BETA']
            else:
                assert 0, 'OR or BETA must be in columns'
            #if 'P' in snp: #### ZW: commented out section because no need to calculate stderr, PLINK already does this
            #    p=snp['P']
            #else:
            #    assert 0, 'P must be in columns' 
            #stderr='abs(%s/qnorm(%s/2))'%(beta,p)
            if 'SE' in snp: #### ZW: created if statement to extract SE column in PLINK2 results if BETA reported
                stderr=snp['SE']
            elif 'LOG(OR)_SE' in snp: #### ZW: created if statement to extract LOG(OR)_SE column in PLINK2 results if OR reported
                stderr=snp['LOG(OR)_SE']
            else:
                assert 0, 'SE or LOG(OR)_SE must be in columns' #### ZW: give error if SE columns are not found
            if pivot == -1:
                pivot=snp # 1ST STUDY's SNP INFO IS PIVOT
                pivotstudyindex=studyindex
            else:
                # CHECK ALLELE TO PIVOT
                if 'A2' in pivot and 'A2' in snp:
                    if pivot['A1'] == snp['A1'] and \
                       pivot['A2'] == snp['A2']:
                        # GOOD
                        pass
                    elif pivot['A1'] == snp['A2'] and \
                         pivot['A2'] == snp['A1']:
                        # SIMPLE FLIP
                        beta=str(-float(beta)) #'-(%s)'%beta  #### ZW: modified so sign of beta is switched here and not reliant on subroutine script
                    elif pivot['A1'] == flipstrand(snp['A1']) and \
                         pivot['A2'] == flipstrand(snp['A2']): #### ZW: changed comple[] to flipstrand()
                        # STRAND INCONSIS., BUT GOOD
                        flog.write('FLIP_STRAND %s in study %d\n'%(rsid,studyindex))
                    elif pivot['A1'] == flipstrand(snp['A2']) and \
                         pivot['A2'] == flipstrand(snp['A1']): #### ZW: changed comple[] to flipstrand()
                        # STRAND INCONSIS., SIMPLE FLIP
                        flog.write('FLIP_STRAND %s in study %d\n'%(rsid,studyindex))
                        beta=str(-float(beta)) #'-(%s)'%beta  #### ZW: modified so sign of beta is switched here and not reliant on subroutine script
                    else:
                        flog.write('EXCLUDE %s due to allele inconsistency: A1:%s A2:%s in study %d but A1:%s A2:%s in study %d\n'
                                   %(rsid, pivot['A1'], pivot['A2'], pivotstudyindex,
                                     snp['A1'], snp['A2'], studyindex))
                        beta='EXCLUDE' #### ZW: added making beta and stderr 'EXCLUDE' if SNP should be excluded
                        stderr='EXCLUDE'
                else: 
                    if pivot['A1'] == snp['A1']:
                        # GOOD
                        pass
                    else:
                        flog.write('EXCLUDE %s due to allele inconsistency: A1:%s in study %d but A1:%s in study %d\n'
                                   %(rsid, pivot['A1'], pivotstudyindex,
                                     snp['A1'], studyindex))
                        beta='EXCLUDE' #### ZW: added making beta and stderr 'EXCLUDE' if SNP should be excluded
                        stderr='EXCLUDE'
                # CHECK CHR & BP TO PIVOT
                if pivot['CHROM'] != snp['CHROM']: #### ZW: changed CHR to CHROM
                    flog.write('WARNING %s has chr inconsistency\n'%rsid)
                if pivot['POS'] != snp['POS']: #### ZW: changed BP to POS
                    flog.write('WARNING %s has basepair inconsistency\n'%rsid)
            output+=beta+' '+stderr+' '
            numstudy+=1
        else:
            output+='NA NA '
    if numstudy == 1:
        flog.write('EXCLUDE %s due to being in single study\n'%rsid)
    elif output.count('EXCLUDE') > 0: #### ZW: added if statement to not write SNPs that need to be excluded
        pass
    else:
        fout.write(output+'\n')
        fmap.write('%s\t%s\t%s\t%s\t%d\n'%
                   (rsid, pivot['CHROM'], pivot['POS'], pivot['A1'], numstudy)) #### ZW: changed CHR and BP to CHROM and POS respectively
fout.close()
fmap.close()
flog.close()

# CALL R TO EVALUATE MATH EXPRESSION        #### ZW: commented this section out, subroutine script not needed with modifications
#subprocess.call(['R --vanilla --slave "--args '+out+'.meta.tmp '+out+'.meta" < plink2metasoft_subroutine.R'],shell=True)
#subprocess.call(['rm '+out+'.meta.tmp'], shell=True)



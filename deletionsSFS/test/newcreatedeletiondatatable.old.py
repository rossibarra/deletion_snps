#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
@author: Edwin
"""

import gzip, argparse, os, glob

parser = argparse.ArgumentParser(description='Process B73 deletion and nondeletion SNP calls in vcf file. Outputs tsv file. Ignores multiallelic sites.')
parser.add_argument('-i', '--vcf', help='VCF input file', type=str, required=True)
parser.add_argument('-b', '--bed', help='BED deletions input directory that contains sorted BED files', type=str, required=True)
parser.add_argument('-o', '--outfile', help='Output filename', type=str, required=False)
parser.add_argument('-O', '--outdir', help='Output directory name', type=str, required=False)
args = parser.parse_args()

#myvcf = 'B73v5.NAM-illumina_filtered-pass-only-two-round-gatk-snps.vcf.gz'
myvcf = args.vcf
compressedfile = False
if myvcf.split('.')[-1] == 'gz':
    compressedfile = True

#
mybeddir = args.bed
if not os.path.isdir(mybeddir):
    print(mybeddir + ' is not accessible, please check your directory location and try again. Also if using relative position for directories, please make sure to prepend with ./\n')
    quit(1)

bedcommonextension = '.deletions.merged.bed'

outFN = args.outfile
if outFN == None:
    outFN = myvcf.replace('.vcf','.datatable.tsv').replace('.gz','')

outFNext = outFN.split('.')[-1]
amboutFN = outFN.replace(outFNext, 'ambiguouscalls.' + outFNext)
junkoutFN = outFN.replace(outFNext, 'junk.' + outFNext)

awdeloutcountsFN = outFN.replace(outFNext, 'awdelcounts.' + outFNext)

outdir = args.outdir
if outdir == None:
    outdir = ''
else:
    #make directory
    if not os.path.exists(outdir):
        os.makedirs(outdir)
        print('Outputting file to: ' + outdir)
    else:
        print('Warning! ' + outdir + ' already exists! Output files will be overwritten if they also exist')

if compressedfile:
    fin = gzip.open(myvcf)
else:
    fin = open(myvcf)

fout = open(outFN, 'w')
ambfout = open(amboutFN, 'w')
junkfout = open(junkoutFN, 'w')

def bedrangebinarysearch(mysearchterm, mysearchnamline, mysearchxsome):
    mysearchranges = mydelrangesdict[mysearchnamline][mysearchxsome]
    mybsearchstart = 0
    myhalflen = 0
    mybsearchend = len(mysearchranges) - 1
    while mybsearchstart <= mybsearchend:
        myhalflen = (mybsearchstart + mybsearchend)//2
        #print(myhalflen)
        mysmallrange = mysearchranges[myhalflen]
        # If found return True. Else traverse the array
        if (mysearchterm == mysmallrange[0]) or (mysearchterm > mysmallrange[0] and mysearchterm < mysmallrange[1]):
            return [True, mybsearchstart]
        # Go Left
        elif mysearchterm < mysmallrange[0]:
            mybsearchend = myhalflen - 1
        # Go Right
        elif mysearchterm >= mysmallrange[1]:
            mybsearchstart = myhalflen + 1
        # If some weird error return not found
        else:
            return [False, -1]
    # if searched every range then return false for not found.
    return [False, -1]

### START NOTES

# for alt allele
# can be A or T or C or G or N or *
# for regular nt it just means that nt i.e. A means A is the alt allele
# for regular nt with commas it means multiple at alleles, i.e. A,T means A is alt1 allele and T is alt2 allele
# for * instead of nt means deletion
# for N instead of nt means ???
# if empty or null then signifies also a deletion: this also applies to ref. ref can be empty/null or A,T,C,G or N.
# if no commas and string length > 1 then short indel


# for genotyping info sorted by info in header for wich variety
# 0/0 means ref allele/ref allele
# 1/1 means alt1 allele/alt1 allele
# 2/2 means alt2 allele/alt2 allele
# k/k means altk allele/altk allele
# ./. means missing allele/missing allele

### END NOTES

# Initialize global vars
print('Initializing variables\n')
#metadatalines = list()
allelecalls = set()
delstr = 'DEL'
missingdatastr = 'MD'
awdelmissingdatastr = 'AMD'
awambiguouscall = 'AMB'
myinfoset = set()
myrefset = set()
myaltset = set()
mysnpcallsset = set()
weirdlist = list()
mypos = 0
mydelrangesdict = dict()
missingdataofnamlineset = set()
mystdalleleset = {'0','1', 'H'}

# take in list of deletion bed files
mydelfiles = glob.glob(mybeddir + '*' + bedcommonextension)

#Create datastructures for the deletion bed files here
print('Reading deletion bed files ending in ' + bedcommonextension + ' in directory ' + mybeddir)
for bed in range(0, len(mydelfiles)):
    mybedFN = mydelfiles[bed]
    mydelID = mydelfiles[bed].split('/')[-1].split('.')[0]
    mydelrangesdict[mydelID] = dict()
    myoldbedstop = 1000000000000
    myoldbedline = 'start'
    myoldchr = 'start'
    print(mydelID)
    for bedline in open(mybedFN):
        # print(bedline)
        bedline = bedline.strip().split()[0:3]
        mybedchr = bedline[0]
        mybedstart = int(bedline[1])
        mybedstop = int(bedline[2])
        if mybedstart >= mybedstop:
            print('Error in bed file where Column 2 is not < Column 3: ' + mybedchr + '\t' + str(mybedstart) + '\t' + str(mybedstop) + '\n')
            #quit(1)
        if mybedchr == myoldchr:
            if mybedstart <= myoldbedstop:
                print('Error in bed file not sorted and/or merged. Please run sort -k1,1v -k2,2n YOURBED, then bedtools merge -i YOURBED. Line is:\n' + myoldbedline + '\n' + mybedchr + '\t' + str(mybedstart) + '\t' + str(mybedstop) + '\n')
                #quit(1)
        if mybedchr not in mydelrangesdict[mydelID].keys():
            mydelrangesdict[mydelID][mybedchr] = list()
        mydelrangesdict[mydelID][mybedchr].append([int(mybedstart), int(mybedstop), 0])
        myoldbedstop = int(bedline[2])
        myoldchr = bedline[0]
        myoldbedline = mybedchr + '\t' + str(mybedstart) + '\t' + str(mybedstop)

# Read first line to test for encoding
print('Begin reading of file ' + myvcf)
myline = fin.readline()
writeheaders = True
readheaders = False
# Check for byte encoding and set file pointer to beginning of file
isbytes = type(myline) == bytes

fin.seek(0)
# Begin iteratring through lines of vcf file
for myline in fin:
    missingdata = False
    awdelflag = False
    referrorflag = False
    ambiguouscall = False
    heterozygouscall = False
    if isbytes:
        # Decode line as utf-8 (unicode)
        myline = myline.decode("utf-8").strip()
    else:
        myline = myline.strip()
    if myline[0] == '#':
        # Read hash tagged header data from vcf file
        myheader = myline.split()
        #metadatalines.append(myline)
        readheaders = True
    else:
        # break
        # Read non has tagged header data from vcf file
        # Also write the headers if writeheaders is true then set to false when done
        if readheaders:
            readheaders = False
            if writeheaders:
                myheader[0] = 'CHROM'
                myhstring = myheader[0] + '\t' + myheader[1] + '\tAWDelID\tMissingSNP_and_SVdel' + '\t' + myheader[3] + '\t' + myheader[4]
                # 
                for h in range(9, len(myheader)):
                    myhstring = myhstring + '\t' + myheader[h] + 'A\t' + myheader[h] + 'B'
                myhstring = myhstring + '\n'
                fout.write(myhstring)
                ambfout.write(myhstring)
                junkfout.write(myhstring)
                writeheaders = False
        # Begin reading of variant calls
        # Fields:
        #   0,    1,   2,     3,         4,         5,          6,             7,                8,                    9:
        # CHROM, POS, ID, REFALLELE, ALTALLELES, PHREDQ, FILTERPASSRESULT, INFOFIELDS, STUFF(GT:AD:DP:GQ:PL), NAMLINE1(0/0:62,0,0:62:99:0,187,2398,187,2524,2790)
        # Read in position and see if it overlaps with a delete in all the deletion bed files.
        myalleles = list()
        myalts = list()
        delstatus = False
        missingdata = False
        #metadata = ''
        myline = myline.split()
        chrom = myline[0]
        pos = myline[1]
        #myinfoset.add(myline[2])
        ref = myline[3]
        myalts = alts = myline[4].split(',')
        myrefset.add(ref)
        for alt in range(0, len(alts)):
            myaltset.add(alts[alt])
            if alts[alt] == ',':
                weirdlist.append(mypos)
            elif alts[alt] == '*':
                delstatus = True
                delpos = alt + 1
            else:
                delpos = -1
        namlines = list()
        # Check if alts has 1 or 2 values. Only accept one alt alleles since having a third is not probable. Also accept *
        # Possible alt values = ['A', 'C', '*', 'T', 'G']
        ### if (not len(myalts) == 0) and (len(myalts) < 3):
        # When more than 1 alt, then check if deletion exists. if not skip.
        if ((len(myalts) == 2) and (myalts [0] == '*' or myalts[1] == '*')) or (len(myalts) ==1):
            # Iterate across nam line SNP calls
            # myheader[9] is the ID for the Ref line. Need to check for homozygous
            i = 9
            mynamline = myheader[i].strip()
            myalleles = myline[i].split(':')[0].split('/')
            for j in range(0, len(myalleles)):
                if not myalleles[j] == '0':
                    referrorflag = True
            if not referrorflag:
                namlines.append(myalleles)
                #metadata = metadata + myline[i] + '\t'
                # myheader[10] is the ID for the first non Ref NAM line. Need to perform normal del check
                for i in range(10, len(myline)):
                    namlinedelstatus = False
                    mynamline = myheader[i].strip()
                    # begin binary search to see if snp call is in AnchorWave deletion
                    myalleles = myline[i].split(':')[0].split('/')
                    if mynamline in mydelrangesdict.keys():
                        if chrom in mydelrangesdict[mynamline].keys():
                            awdelflag, awdelindex = bedrangebinarysearch(int(pos), mynamline, chrom)
                            awdelrange = mydelrangesdict[mynamline][chrom][awdelindex]
                        else:
                            awdelflag = False
                            awdelrange = []
                            awdelindex = -1
                        # Run two seperate loops to test for heterozygosity
                        alleledistance = len(myalleles)
                        if alleledistance > 2:
                            print('More than one allele found on: ' + chrom + '\t' + pos)
                            quit(1)
                        else:
                            j = 0
                            # check homozygous calls
                            if myalleles[j] == myalleles[j+1]:
                                # Check if there is there a deletion at this site
                                if delstatus:
                                    # Check if there is NAMLine overlap with AWDeletion
                                    # Set proper del string and set deletion status for NAMLine 
                                    if myalleles[j] == str(delpos) and awdelflag:
                                        myalleles[j] = myalleles[j+1] = delstr
                                        namlinedelstatus = True
                                    # Since no overlap with AWDel check if there is still a deletion call
                                    # Set * del string and deletion status for NAMLine
                                    elif myalleles[j] == str(delpos):
                                        myalleles[0] = myalleles[1] = '*'
                                        namlinedelstatus = True
                                # Perfrom Independent check for Missing data and if overlap with AWDel
                                # Set proper allele string for overlap or no overlap
                                if myalleles[j] == '.':
                                    missingdata = True
                                    if awdelflag:
                                        myalleles[j] = myalleles[j+1] = awdelmissingdatastr
                                    else:
                                        myalleles[j] = myalleles[j+1] = missingdatastr
                                # Check if there exists an ambiguous SNP call with no NAMDel call
                                # Increment SNPCalls overlapping deletions, set ambiguous call
                                elif awdelflag and myalleles[j] == '1' and not namlinedelstatus:
                                    ambiguouscall = True
                                    myalleles[j] = myalleles[j+1] = awambiguouscall
                                    mydelrangesdict[mynamline][chrom][awdelindex][2] += 1
                                # Report weird sitatuations. Unexpected cases
                                elif not namlinedelstatus and not namlinedelstatus:
                                    print('Weird homozygous call at: ' + chrom + '\t' + pos + '\t' + myalleles[0] + '/' + myalleles[1])
                                else:
                                    print('Woah! homozygous call at: ' + chrom + '\t' + pos + '\t' + myalleles[0] + '/' + myalleles[1] + ' delstatus: ' + str(delstatus))
                            # Now check heterozygous SNP Calls and set Heterozygous calls
                            else:
                                # Check for missing data and set MD flag
                                if myalleles[j] == '.' or myalleles[j+1] == '.':
                                    missingdata = True
                                # Independent check for NAMLine del call and set flag
                                if (myalleles[j] == str(delpos) or myalleles[j+1] == str(delpos)) and delstatus:
                                    namlinedelstatus = True
                                # Independent check for overlap of NAMLine del or MD with AWDel overlap
                                if (namlinedelstatus or missingdata) and awdelflag:
                                    myalleles[0] = myalleles[1] = 'H' + awdelmissingdatastr
                                # if no AWDel overlap but NAMLine del or MD
                                elif namlinedelstatus or missingdata:
                                    myalleles[0] = myalleles[1] = 'H' + missingdatastr
                                # Independent check if expected heterozygous call exists
                                # Check only first char in string and set het flag
                                if myalleles[0][0] in mystdalleleset and myalleles[1][0] in mystdalleleset:
                                    heterozygouscall = True
                                    #print('Heterozygous call at: ' + chrom + '\t' + pos + '\t' + myalleles[0] + '/' + myalleles[1])
                                # Print message if unexpected case for het call
                                else:
                                    print('Weird heterozygous case at: ' + chrom + '\t' + pos + '\t' + myalleles[0] + '/' + myalleles[1])
                        namlines.append(myalleles)
                        #metadata = metadata + myline[i] + '\t'
                    else:
                        # Print error message only once if missing NAMLine AWdeletion data
                        # Set allele call to E
                        if mynamline not in missingdataofnamlineset:
                            missingdataofnamlineset.add(mynamline)
                            print('Data or bed file missing in: ' + mynamline)
                        for j in range(0, len(myalleles)):
                            myalleles[j] = 'E'
                        namlines.append(myalleles)
                        #metadata = metadata + myline[i] + '\t'
                #metadata = metadata[:-1]
                # make a new field where missing data overlaps with a SV deletion
                # fieldshould be ID field idx 2
                # for now heterozygous sites are ignored. NA means ambiguous. Make heterozygous sites ambiguous and report NA?
                # - for deletions, and keep . for missing data
                # what to do when a site is -/INT
                # Fields:
                #   0,    1,   2,     3,         4,         5,          6,             7,                8,                    9:
                # CHROM, POS, ID, REFALLELE, ALTALLELES, PHREDQ, FILTERPASSRESULT, INFOFIELDS, STUFF(GT:AD:DP:GQ:PL), NAMLINE1(0/0:62,0,0:62:99:0,187,2398,187,2524,2790)
                # Place deletion and missingdata info in ID
                # Place our new fields in 9+, then add metadata to end of the line
                # Check if we have an AWDel range, if so set AWDelID, else set ID to None
                if len(awdelrange) > 2:
                    awdelid = chrom + ':' + str(awdelrange[0]) + '-' + str(awdelrange[1])
                else:
                    awdelid = 'None'
                # The StatusStr might be unnecessary
                mystatusstr = str(int(delstatus)) + ',' + str(int(missingdata))
                mystring = chrom + '\t' + pos + '\t' + awdelid + '\t' + mystatusstr + '\t' + ref + '\t'
                for i in range(0, len(alts)):
                    mystring = mystring + alts[i] + ','
                mystring = mystring[:-1] + '\t'
                #mystring = mystring + myline[5] + '\t' + myline[6] + '\t' + myline[7] + '\t' + myline[8] + '\t'
                for i in range(0, len(namlines)):
                    mystring = mystring + str(namlines[i][0]) + '\t' + str(namlines[i][1]) + '\t'
                #mystring = mystring + metadata + '\n'
                mystring = mystring[:-1]
                mystring = mystring + '\n'
                # if call is clean (not ambiguous and not heterozygous biallelic)
                # Multi allelic calls are discarded
                if not ambiguouscall and not heterozygouscall:
                    fout.write(mystring)
                # Ambiguous but not heterozygous
                elif ambiguouscall and not heterozygouscall:
                    ambfout.write(mystring)
                # Biallelic heterozygous or heterozygous and MD or Del
                elif heterozygouscall:
                    junkfout.write(mystring)
                    print('Heterozygous call at: ' + chrom + '\t' + pos)
                # This might output multiallelic calls?
                else:
                    print('Woah! this is weird and not written at: ' + chrom + '\t' + pos)
fin.close()
fout.close()
ambfout.close()
junkfout.close()
fawdelcallsout = open(awdeloutcountsFN, 'w')
fawdelcallsout.write('NAMLine\tAWDelID\tSNPOverlapCounts\n')
for namline in mydelrangesdict.keys():
    for chrom in mydelrangesdict[namline].keys():
        for i in range(0, len(mydelrangesdict[namline][chrom])):
            myawdelrangencount = mydelrangesdict[namline][chrom][i]
            fawdelcallsout.write(namline + '\t' + namline + ':' + chrom + ':' + str(myawdelrangencount[0]) + '-' + str(myawdelrangencount[1]) + '\t' + str(myawdelrangencount[2]) + '\n')
fawdelcallsout.close()

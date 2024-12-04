#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
@author: Edwin
"""

import gzip, argparse

parser = argparse.ArgumentParser(description='Process B73 alignment sam files and output large csv.')
parser.add_argument('-i', '--sam', help='Alignment sorted sam file', type=str, required=True)
parser.add_argument('-o', '--outfile', help='Output filename', type=str, required=False)
parser.add_argument('-r', '--fastaref', help='Input Reference fasta file instead of bed file', type=str, required=False)
parser.add_argument('-q', '--fastaqry', help='Input Query fasta file instead of bed file', type=str, required=False)
parser.add_argument('-O', '--outdir', help='Output directory name', type=str, required=False)
args = parser.parse_args()
mysam = args.sam
refFasta = args.fastaref
qryFasta = args.fastaqry
outFNd = outFNm = args.outfile
if outFNd == None:
    outFNd = mysam.replace('.sam','.deletions.bed').replace('.gz','')
    outFNm = mysam.replace('.sam','.matches.bed').replace('.gz','')
#    outFN = refFasta.split('-')[0] + '-' + refFasta.split('-')[1] + '-' + qryFasta.split('-')[1] + '.bed'
outdir = args.outdir
if outdir == None:
    outdir = ''
else:
    #make directory
    print(outdir)

outputfasta = False

if outputfasta:
    fin = open(b73refchr1)
    myb73chr1seq = fin.readlines()
    myb73chr1seq = myb73chr1seq[1].strip()
    fin.close()

    fin = open(oh7brefchr1)
    oh7bchr1seq = fin.readlines()
    oh7bchr1seq = oh7bchr1seq[1].strip()
    fin.close()

if mysam[-3:] == '.gz':
    fin = gzip.open(mysam)
else:
    fin = open(mysam)
mylines = fin.readlines()
if type(mylines[0]) == bytes:
    for i in range(0,len(mylines)):
        mylines[i] = mylines[i].decode("utf-8").strip().split()
else:
    for i in range(0,len(mylines)):
        mylines[i] = mylines[i].strip().split()

mynumset = {'0', '1', '2', '3', '4', '5', '6', '7', '8', '9'}
mycharset = {'M', 'I', 'D', 'X', '=', 'H'}

foutd = open(outFNd, 'w')
foutm = open(outFNm, 'w')
mybigdeletionslist = list()
mybigmatcheslist = list()
mybiginsertionslist = list()
mybignewinsertionslist = list()
mybighardclippinglist = list()
for j in range(0,len(mylines)):
    if mylines[j][0][0] != '@':
        myqry = mylines[j][0]
        myflag = mylines[j][1]
        myref = mylines[j][2]
        mycigar = mylines[j][5]
        mysequence = mylines[j][9]
        mystartpos = int(mylines[j][3])-1
        if myflag == '0':
            direction = '+'
        elif myflag == '16':
            direction = '-'
        myintstr = ''
        mycigstr = ''
        mycigarsdict = dict()
        mycurrpos = mystartpos
        mycurrstartpos = 0
        mycurrstoppos = 0
        mydelstartpos = 0
        mydeletionslist = list()
        mymatcheslist = list()
        myinsertionslist = list()
        myhardclippinglist = list()
        for i in mycharset:
            mycigarsdict[i] = list()
        for i in range(0,len(mycigar)):
            mychar = mycigar[i]
            if mychar in mynumset:
                myintstr = myintstr + mychar
            else:
                if mychar in mycharset:
                    myint = int(myintstr)
                    #here we detect deletions
                    if mychar == 'D':
                        mycurrpos = mycurrpos + myint
                        mycurrstartpos = mycurrpos - myint
                        mycurrstoppos = mycurrpos
                        mydeletionslist.append([mycurrstartpos, mycurrstoppos])
                        #isprevcigardel = True
                    #here we can detect matches
                    elif mychar == 'M' or mychar == '=' or mychar == 'X':
                        mycurrpos = mycurrpos + myint
                        mycurrstartpos = mycurrpos - myint
                        mycurrstoppos = mycurrpos
                        mymatcheslist.append([mycurrstartpos, mycurrstoppos])
                        #isprevcigardel = False
                    #here we can detect insertions
                    elif mychar == 'I':
                        mycurrstartpos = mycurrpos - myint
                        mycurrstoppos = mycurrpos
                        myinsertionslist.append([mycurrstartpos, mycurrstoppos])
                    elif mychar == 'H':
                        mycurrstartpos = mycurrpos
                        mycurrstoppos = mycurrpos
                        myhardclippinglist.append([mycurrstartpos, mycurrstoppos, myint])
                    mycigarsdict[mychar].append(myint)
                    myintstr = ''
                else:
                    print('New CIGAR type found! :' + mychar + ' with ' + myintstr + ' NTs on ' + myref)
        for i in range(0, len(mydeletionslist)):
            foutd.write(myref + '\t' + str(mydeletionslist[i][0]) + '\t' + str(mydeletionslist[i][1]) + '\t' + myqry + '\t' + '0' + '\t' + direction + '\n')
        newmymatcheslist = list()
        mycount = 0
        for i in range(1, len(mymatcheslist)):
            myprevstartpos = mymatcheslist[i-1][0]
            myprevstoppos = mymatcheslist[i-1][1]
            mycurrstartpos = mymatcheslist[i][0]
            mycurrstoppos = mymatcheslist[i][1]
            myoverlap = myprevstoppos - mycurrstartpos
            if myoverlap <= 0:
                if mycount == 0:
                    myrealstartpos = myprevstartpos
                else:
                    myrealstoppos = mycurrstoppos
                mycount += 1
            else:
                if mycount == 0:
                    newmymatcheslist.append([myprevstartpos, myprevstoppos])
                else:
                    newmymatcheslist.append([myrealstartpos, myrealstoppos])
                mycount = 0
        myrealstoppos = i = mycount = 0
        mynewmatcheslist = list()
        for i in range(1, len(mymatcheslist)):
            myprevstartpos = mymatcheslist[i-1][0]
            myprevstoppos = mymatcheslist[i-1][1]
            mycurrstartpos = mymatcheslist[i][0]
            mycurrstoppos = mymatcheslist[i][1]
            myoverlap = mycurrstartpos - myprevstoppos
            if myoverlap <= 0:
                if mycount == 0:
                    myrealstartpos = myprevstartpos
                else:
                    myrealstoppos = mycurrstoppos
                mycount+=1
            else:
                if mycount == 0:
                    mynewmatcheslist.append([min(myprevstartpos, myprevstoppos), max(myprevstartpos, myprevstoppos)])
                else:
                    mynewmatcheslist.append([min(myrealstartpos, myrealstoppos), max(myrealstartpos, myrealstoppos)])
                mycount = 0
                myrealstartpos = myprevstartpos
            if i == len(mymatcheslist) - 1:
                mynewmatcheslist.append([min(myrealstartpos, mycurrstoppos), max(myrealstartpos, mycurrstoppos)])
        for match in range(0, len(mynewmatcheslist)):
            foutm.write(myref + '\t' + str(mynewmatcheslist[match][0]) + '\t' + str(mynewmatcheslist[match][1]) + '\t' + myqry + '\t' + '0' + '\t' + direction + '\n')
        mybigdeletionslist.append(mydeletionslist)
        mybiginsertionslist.append(myinsertionslist)
        mybighardclippinglist.append(myhardclippinglist)
        mybigmatcheslist.append(mynewmatcheslist)
foutd.close()
foutm.close()
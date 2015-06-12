#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys

fnfull = sys.argv[1]
fid = open(fnfull, 'rU')
raw = fid.readlines()
fid.close()

def convertTimes(strinput):
    tmp = strinput.split(':')
    if len(tmp) == 2:
        # the first one is minutes, second one is seconds
        retval = float(tmp[0]) * 1 + float(tmp[1]) / 60.
        return str(retval)
    elif len(tmp) == 3:
        retval = float(tmp[0]) * 60. + float(tmp[1]) * 1. + float(tmp[2]) / 60.
        return str(retval)


# number of total lines
Nlines = len(raw)
initID  = []
# find the initial index, this is where the data heaader begins
for rawid, rawline in enumerate(raw):
    if 'Plate#1' in rawline:
        initID = rawid

# get metadata from the 4th line
metadat = raw[initID].split('\t')

# Ndata points

Ndat = int(metadat[8])

# Wavelengths read at
WaveRead = metadat[15]
WaveRead = WaveRead.split(' ')
Nwave    = len(WaveRead)

# Acquisition type
AcqType  = raw[initID].split('\t')[5]

# Get Data that is not empty
beginIdx = initID+2 
example  = raw[initID+2].split('\t')
goodCol  = []

for cc in range(len(example)):
    if example[cc] != '' and example[cc] != '\n':
        goodCol.append(int(cc))

# Get column titles
header   = raw[initID+1].split('\t')
header   = [header[idx] for idx in goodCol]
header[0] = 'Time(min)'
header[1] = 'Temperature'

for eachWavelength in range(Nwave):
    
    fnout = '%s_%s%s.%s' % (fnfull[:-4],AcqType,WaveRead[eachWavelength],'tab.txt')
    fid   = open(fnout,'w')
    print 'Wavelength : %d' % eachWavelength

    #Write header
    for hid,hval in enumerate(header):
        fid.write(hval)
        if hid < (len(header)-1):
            fid.write('\t')

    fid.write('\n')

    if eachWavelength > 0:
        spacer = 1
    else:
        spacer = 0
        
    for eachRow in range(Ndat):
        
        currow = beginIdx + Ndat*eachWavelength + eachRow + spacer

        tmp = raw[currow].split('\t')
        goodtmp = [tmp[eachCol] for eachCol in goodCol]
        goodtmp[0] = convertTimes(goodtmp[0])

        for colid,colval in enumerate(goodtmp):
            fid.write(colval)
            if colid < (len(goodtmp)-1):
                fid.write('\t')
                
        fid.write('\n')

    fid.close()
        
            
    
        

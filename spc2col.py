#!/usr/bin/env python

"""
This program converts SPC file from the Fluoromax-3 into a tab-delimited
text file that one can easily import into Excel, Numbers, KaleidaGraph, etc.

Daniel Elnatan (Agard Lab)
- Adapted from tgscpread from MATLAB's Bioinformatics Toolbox

"""

import sys

def frange(start, stop, step):
    
    stepval = float(step)
    npts    = (float(stop)-float(start))
    tmp     = []
    tmp.append(float(start))
    for ii in range(int(npts)):
        tmp.append(tmp[ii] + stepval)
    return tmp

def fread(fid, bytecount, typename):
    import struct
    TYPES = {'int8'  :'b',
            'uint8'  :'B',
            'int16'  :'h',
            'uint16' :'H',
            'int32'  :'i',
            'uint32' :'I',
            'int64'  :'q',
            'uint64' :'Q',
            'float'  :'f',
            'double' :'d',
            'char'   :'s'}
        
    TYPESIZE  = struct.calcsize(TYPES[typename])
    fmt       = '<' + str(bytecount) + TYPES[typename]
    fmt       = str(fmt)
    # enforce little-endian byte ordering
    return struct.unpack(fmt, fid.read(TYPESIZE*bytecount))[0]

for fn_input in sys.argv[1:]:
	#fn_input  = sys.argv[1] # input file
	#fn_output = sys.argv[2] # output file
	fn_output = fn_input[:-4] + '.txt'
	
	fn = fn_input
	
	fid = open(fn, 'rb')
	
	# Begin reading SPC file
	ftflgs = fread(fid, 1, 'uint8')
	fversn = fread(fid, 1, 'uint8')
	
	fexper = fread(fid, 1, 'uint8')
	fexp   = fread(fid, 1, 'uint8')
	fnpts  = fread(fid, 1, 'uint32') # number of sampled wavelengths
	ffirst = fread(fid, 1, 'double') # First wavelength 
	flast  = fread(fid, 1, 'double') # Last wavelength
	fnsub  = fread(fid, 1, 'uint32') # Number of traces
	
	fxtype = fread(fid, 1, 'uint8')
	fytype = fread(fid, 1, 'uint8')
	fztype = fread(fid, 1, 'uint8')
	fpost  = fread(fid, 1, 'uint8')
	
	fdate  = fread(fid, 1, 'uint32')
	fres   = fread(fid, 9, 'uint8')
	fsource = fread(fid, 9, 'uint8')
	fpeakpt = fread(fid, 1, 'uint16')
	fspare  = fread(fid, 8, 'float') # this is actually a single, but python doesn't have the support
	fcmnt   = fread(fid, 130,'uint8')
	fcatxt  = fread(fid, 30, 'uint8')
	flogoff = fread(fid, 1, 'uint32')
	fmods   = fread(fid, 1, 'uint32')
	fprocs  = fread(fid, 1, 'uint8')
	flevel  = fread(fid, 1, 'uint8')
	fsampin = fread(fid, 1, 'uint16')
	ffactor = fread(fid, 1, 'float') # also a single
	fmethod = fread(fid, 48, 'uint8')
	fzinc   = fread(fid, 1, 'float')
	fwplanes = fread(fid, 1, 'uint32')
	fwinc   = fread(fid, 1, 'float')
	fwtype  = fread(fid, 1, 'uint8')
	freserv = fread(fid, 187, 'uint8')
	
	subheader = {'subflgs':[],
	             'subexp' :[],
	             'subindx':[],
	             'subtime':[],
	             'subnext':[],
	             'subnois':[],
	             'subnpts':[],
	             'subscan':[],
	             'subwlevel':[],
	             'subresv':[]}
	
	for traceCount in range(fnsub):
	    
	    subheader['subflgs'].append(fread(fid, 1,'uint8'))                               
	    subheader['subexp'].append(fread(fid, 1, 'uint8'))
	    subheader['subindx'].append(fread(fid, 1, 'uint16'))
	    subheader['subtime'].append(fread(fid, 1, 'float'))
	    subheader['subnext'].append(fread(fid, 1, 'float'))
	    subheader['subnois'].append(fread(fid, 1, 'float'))                           
	    subheader['subnpts'].append(fread(fid, 1, 'uint32'))
	    subheader['subscan'].append(fread(fid, 1, 'uint32'))
	    subheader['subwlevel'].append(fread(fid, 1, 'float'))
	    subheader['subresv'].append(fread(fid, 4, 'uint8'))
	
	Ydata = []
	
	for ii in range(fnpts):
	    Ydata.append(fread(fid,1,'uint32'))
	
	fid.close()
	
	Xdata = frange(ffirst,flast, (flast-ffirst+1)/fnpts)
	Ydata = [(2**fexp*y / 2**32) for y in Ydata]
	
	fid = open(fn_output,'wt')
	fid.write('X\tY\n')
	
	for rr in range(fnpts):
	    fid.write(str(Xdata[rr]) + '\t' + str(Ydata[rr]))
	    if rr < fnpts:
	        fid.write('\n')
	
	fid.close()



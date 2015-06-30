#===============================================================================
#
# author: Miguel Betegon
# last edited: 150529
# version 0.4
#
# NOTE: ALL positions used in these classes are 0-offset. Positions parsed from
# files with a different offset will be adjusted to be 0-offset.
#
#===============================================================================

import cPickle
from string import maketrans
from itertools import groupby



class Obj:
	"""General object that can be customized on the fly"""
	def __init__(self):
		pass
		
	def __getattribute__(self, name):
		if name in self.__dict__[name]:
			return self.__dict__[name]
		return None
	
	def __setattr__(self, name, value):
		self.__dict__[name] = value


class Graph:
	
	def __init__(self):
		pass
	
	def bar(self, xValues, yValues, height=25, char='*'):
		maxVal = float(max(yValues))
		print maxVal
		for i in range(height, 0, -1):
			row = []
			for j in yValues:
				if j/maxVal*height > i-1:
					row.append(char)
				else:
					row.append(' ')
			print ''.join(row)
		xAxis = []
		for x in range(0, len(xValues), 10):
			xAxis.append(str(xValues[x])+' '*(10-len(str(xValues[x]))))
		print ''.join(xAxis)
			
	


class File:
	"""The File class opens and closes files and provides iterators for normal text, fa, fq, bed (4 or 6 columns) and bowtie files"""
	def __init__(self, fileName, mode='r'):
		"""Initialise the File object and open the file"""
		self.fileName = fileName
		self.fh = None
		self.entries = -1
		self.mode = mode
		self.fh = open(self.fileName, mode)
	
	def numEntries(self):
		"""Return a quick rough estimate of the number of entries (fq records, fa records or lines, depending on the file)"""
		if self.entries == -1:
			self.entries = self._entryCount()
		return self.entries
	
	def __iter__(self):
		"""Line iterator. Yields lines from the file"""
		for line in self.fh:
			yield line
	
	def close(self):
		"""Close the File object"""
		if self.fh != None:
			self.fh.close()
			self.fh = None
	
	def write(self, strg):
		"""Write string strg to the file"""
		self.fh.write(strg)
	
	def _entryCount(self):
		"""Calculate a quick rough estimate of the number of entries (fq records, fa records or lines, depending on the file)"""
		count = 0
		chars = 0.0
		for line in self.fh:
			if count == 0:
				if line[0] == '@':
					linesPerEntry = 4
				elif line[0] == '>':
					linesPerEntry = 2
				else:
					linesPerEntry = 1
			chars += len(line)
			count += 1
			if count == 1000000:
				break
		if count < 1000000:
			return int(count)/linesPerEntry
		charsPerLine = chars/count
		self.fh.seek(0,2)
		totalBytes = float(self.fh.tell())
		self.fh.seek(0,0)
		totalEntries = totalBytes / (charsPerLine * linesPerEntry)
		return (int(totalEntries)/1000000)*1000000
	
	def faIterOneLine(self):
		"""Fast fa iterator, ONLY for fa files with one line of sequence per read.\nYields FaRecord objects"""
		for line in self.fh:
			yield FaRecord(line.rstrip()[1:], self.fh.next().rstrip())
	
	def fqIter(self):
		"""Fast fq iterator, ONLY for fa files with one line of sequence per read\nYields FqRecord objects"""
		for line in self.fh:
			seq = self.fh.next().rstrip()
			self.fh.next()
			score = self.fh.next().rstrip()
			yield FqRecord(line.rstrip()[1:], seq, score)
	
	def faIter(self):
		"""Standard fa iterator for records with more than one line of sequence.\nYields FaRecord objects"""
		groups = (x[1] for x in groupby(self.fh, lambda line: line[0] == ">"))
		for group in groups:
			title = group.next()[1:].strip()
			seq = "".join(s.strip() for s in groups.next())
			yield FaRecord(title, seq)
	
	def bwtIter(self):
		"""Bowtie 1.x output iterator.\nYields BowtieAlignment objects"""
		for line in self.fh:
			yield BowtieAlignment(line.rstrip())
			
	def bedIter(self):
		"""BED file iterator (accepts 3 and 6 column files).\nYields Feature objects"""
		for line in self.fh:
			fields = line.strip().split('\t')
			if len(fields) == 3:
				yield Feature(fields[0], '', fields[1], fields[2])
			else:
				yield Feature(fields[0], fields[3], fields[1], fields[2], fields[5])
	
	def blastIter(self):
		"""BLAST -m8 output iterator.\nYields BlastHit objects"""
		for line in self.fh:
			yield BlastHit(line)
				

class FaRecord:
	"""FaRecord objects capture fasta records (title and sequence) and provide functionality"""
	
	def __init__(self, title, seq):
		"""Initialize the object"""
		self.title = title
		self.seq = Sequence(seq)
		
	def __str__(self):
		"""Return the text representation of the record in the form:\n>title\nsequence"""
		return '>%s\n%s\n' % (self.title, str(self.seq.wrap()))
	
	def __len__(self):
		"""Return the length of the sequence in the record"""
		return len(self.seq)
	
	def __getslice__(self, start, end):
		"""Return a new FaRecord object where the sequence is trimmed like sequence[start:end]"""
		return FaRecord(self.title, self.seq[start:end])
	
	
	
class FqRecord(FaRecord):
	"""FqRecord objects capture fastq records (title, sequence and score) and provide functionality"""
	
	def __init__(self, title, seq, score):
		"""Initialize the object"""
		FaRecord.__init__(self, title, seq)
		self.score = score
		
	def __getslice__(self, start, end):
		"""Return a new FqRecord object where the sequence and score are trimmed [start:end]"""
		return FqRecord(self.title, self.seq[start:end], self.score[start:end])
	
	def __str__(self):
		"""Return the text representation of the record in the form:\n@title\nsequence\n+\nscore"""
		return '@%s\n%s\n+\n%s\n' % (self.title, self.seq, self.score)

				
class BlastHit:
	
	def __init__(self, line):
		fields = line.strip().split('\t')
		self.queryId = fields[0]
		self.subjectId = fields[1]
		self.percentId = float(fields[2])
		self.alignmentLength = int(fields[3])
		self.numMismatches = int(fields[4])
		self.numGaps = int(fields[5])
		self.queryStart = int(fields[6])-1 # Positions start at 1 in blast output
		self.queryEnd = int(fields[7])-1
		self.subjectStart = int(fields[8])-1 # Positions start at 1 in blast output
		self.subjectEnd = int(fields[9])-1
		self.eValue = float(fields[10])
		self.bitScore = float(fields[11])
	
	

class gPickle:
	""" Wrapper for python's cPickle. Now we can (un)serialize in just one line """
	
	@staticmethod
	def serialize(obj, fileName):
		"""Save object obj in file fileName"""
		pickleFh = open(fileName, 'w')
		cPickle.dump(obj, pickleFh)
		pickleFh.close()
	
	@staticmethod
	def unserialize(fileName):
		"""Get object from file fileName"""
		try:
			pickleFh = open(fileName, 'r')
		except Exception:
			print "Error: pickle file %s not found." % fileName
			exit(1)
		obj = cPickle.load(pickleFh)
		pickleFh.close()
		return obj



class Sequence:
	
	_transTable = maketrans('ACGTNXacgtnx-', 'TGCANXtgcanx-')
	_toDNATable = maketrans('Uu', 'Tt')
	_toRNATable = maketrans('Tt', 'Uu')
	aminoacids = {
		'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
		'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
		'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
		'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
		'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
		'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
		'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
		'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
		'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
		'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
		'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
		'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
		'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
		'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
		'TAC':'Y', 'TAT':'Y', 'TAA':'*', 'TAG':'*',
		'TGC':'C', 'TGT':'C', 'TGA':'*', 'TGG':'W'}
	
	""" Private method, don't call from outside class """
	@staticmethod
	def _translateCodon(codon):
		if codon.upper() not in Sequence.aminoacids:
			return '?'
		return Sequence.aminoacids[codon.upper()]
	
	""" Return string with translated sequence """
	def protein(self):
		return Sequence(''.join([Sequence._translateCodon(self.sequence[3*i:(3*i)+3]) for i in range(0,len(self.sequence)/3)]))
	
	def __init__(self, sequence):
		self.sequence = sequence
	
	def __str__(self):
		return self.sequence
	
	def __len__(self):
		return len(self.sequence)
	
	def __getitem__(self, i):
		return self.sequence[i]

	def __getslice__(self, start, stop):
		return Sequence(self.sequence[start:stop])
	
	""" Return sequence string with \n every 'wrap' characters """
	def wrap(self, wrapAt=80):
		return Sequence('\n'.join([self.sequence[x: x+wrapAt] for x in xrange(0, self.length(), wrapAt)]))
	
	""" Return the number of ocurrences of c in the sequence """
	def count(self, c):
		return self.sequence.count(c)
	
	def toDna(self):
		return Sequence(self.sequence.translate(Sequence._toDNATable))
	
	def toRna(self):
		return Sequence(self.sequence.translate(Sequence._toRNATable))
	
	def gcContent(self):
		return (self.sequence.count('G')+self.sequence.count('C')+self.sequence.count('g')+self.sequence.count('c'))/float(len(self.sequence))
	
	""" Return the lzw compression score of the sequence """
	def lzw(self):
		output = []
		table = dict(dict((chr(i), i) for i in range(256)))
		s = ''
		for ch in self.sequence:
			it = s + ch
			if it in table:
				s = it
			else:
				output.append(table[s])
				table[it] = len(table)
				s = ch
		output.append(table[s])
		return len(output)
	
	""" Return the reverse complement of the sequence """
	def reverseCompl(self):
		return Sequence(self.sequence.translate(Sequence._transTable)[::-1])
		
	""" Return a dictionary with the nucleotide frequencies of the sequence. {'A':X, 'C':X, 'G':X, 'T':X, 'N':X} """
	def nucleotideFreq(self):
		return dict([(nucl, self.sequence.count(nucl)/float(len(self.sequence))) for nucl in 'ACGTN'])
	
	""" 2 sequence objects are equal if their sequence strings are the same """
	def __eq__(self, sequence):
		if self.get().upper() == sequence.get().upper():
			return True
		return False
	
	def getOrfs(self):
		for frame in (0,1,2):
			inOrf = False
			orfStarts = []
			for i in range(0, len(self.sequence)-2, 3):
				i += frame
				codon = self.sequence[i:i+3]
				if codon == 'ATG':
					orfStarts.append(i)
					inOrf = True
				if codon  in ('TAA', 'TAG', 'TGA') and inOrf:
					orfEnd = i+3
					inOrf = False
					for orfStart in orfStarts:
						yield (orfStart, orfEnd)
					orfStarts = []
			


class Read(Feature):
	
	def __init__(self, chrm, name, start, end, strand, sequence, score=None):
		Feature.__init__(self, chrm, name, start, end, strand)
		self.sequence = Sequence(sequence)
		self.score = score
	
	def getSequence(self):
		return self.sequence
	
	def __str__(self):
		return str(self.sequence)
	
	def __len__(self):
		return len(self.sequence)



# Positions in a bowtie file are 0-offset
class BowtieAlignment(Read):
	
	def __init__(self, line):
		self.line = line.strip()
		fields = self.line.split('\t')
		title = fields[0]
		strand = fields[1]
		chrm = fields[2]
		start = int(fields[3])
		seq = fields[4]
		end = start + len(seq)
		Read.__init__(self, chrm, title, start, end, strand, seq)
		self.col7 = fields[6]
		
	def getScore(self):
		return self.score
	
	def __str__(self):
		return self.line
	
	def __len__(self):
		return len(self.seq)

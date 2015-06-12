import sys

if len(sys.argv) < 2:
    print "Usage: platereader.py input.txt"
    exit(0)

inputFile = sys.argv[1]
inputFh = open(inputFile)
lines = inputFh.readlines()[0].split('\r')
inputFh.close()

outFile = inputFile[:-4] + "_csv.txt"
outFh = open(outFile, 'w')

for line in lines:
    if line.startswith('#'):
        continue
    elif line.startswith('Plate:'):
        continue
    elif line.startswith('Time'):
        continue
    else:
        fields = line.split('\t')
        if len(fields) < 3:
            continue
        time = fields[0]
        temperature = fields[1]
        values = [x  for x in fields[2:] if x != '']
    outFh.write("%s\t%s\n" % (time, '\t'.join(values)))

outFh.close()

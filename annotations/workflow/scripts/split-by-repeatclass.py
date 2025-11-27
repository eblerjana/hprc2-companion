import sys
import subprocess

def name_to_class(name):
	if "DNA" in name:
		return "DNA"
	elif "LINE" in name:
		return "LINE"
	elif "LTR" in name:
		return "LTR"
	elif "Satellite" in name:
		return "Satellite"
	elif "SINE" in name:
		return "SINE"
	elif "SegDup" in name:
		return "SegDup"
	else:
		return "OTHER"


outname = sys.argv[1]


files = {
	"DNA": open(outname + '_DNA.bed', 'w'),
	"LINE":  open(outname + '_LINE.bed', 'w'), 
	"SINE":  open(outname + '_SINE.bed', 'w'),
	"Satellite":  open(outname + '_Satellite.bed', 'w'),
	"LTR":  open(outname + '_LTR.bed', 'w'),
	"SegDup": open(outname + '_SegDup.bed', 'w'),
	"OTHER": open(outname + '_OTHER.bed', 'w')
}

for line in sys.stdin:
	if line.startswith('#'):
		continue
	fields = line.strip().split()
	repeatclass = name_to_class(fields[4])
	files[repeatclass].write('\t'.join([fields[0], fields[1], fields[2]]) + '\n')

for repeatclass in files.keys():
	files[repeatclass].close()

	cmd = "bedtools merge -i " + outname + '_' + repeatclass + '.bed' +  " > " + outname + '_' + repeatclass + '_merged.bed'
	subprocess.call(cmd, shell=True)



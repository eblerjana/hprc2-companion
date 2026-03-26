import sys
import subprocess
import gzip
from collections import defaultdict

def name_to_class(name):
	if "AMPL" in name:
		return "AMPL"
	elif "CEN" in name:
		return "CEN"
	elif "DYZ" in name:
		return "DYZ"
	elif "HET" in name:
		return "HET"
	elif "SAT" in name:
		return "SAT"
	elif "TELO" in name:
		return "TELO"
#	elif "ERRBASE" in name:
#		return "ERRBASE"
	elif "PAR" in name:
		return "PAR"
	elif "XDR" in name:
		return "XDR"
	elif "XTR" in name:
		return "XTR"
	else:
		return None

def parse_fai(filename):
	chr_to_len = {}
	for line in open(filename, 'r'):
		fields = line.strip().split()
		chr_to_len[fields[0]] = int(fields[1])
	return chr_to_len


bed_file = sys.argv[1]
vcf_file = sys.argv[2]
fai_file = sys.argv[3]
outname = sys.argv[4]

repeatclasses = ["AMPL", "CEN", "DYZ", "HET", "SAT", "TELO", "PAR", "XDR", "XTR"]


## step 1: create separate BED file for each repeat class

files = {}

for r in repeatclasses:
	files[r] = open(outname + "_" + r + ".bed", "w")

for line in open(bed_file, 'r'):
	if line.startswith('#'):
		continue
	fields = line.strip().split()
	repeatclass = name_to_class(fields[3])
	if repeatclass is not None:
		files[repeatclass].write('\t'.join([fields[0], fields[1], fields[2], repeatclass]) + '\n')

for repeatclass in files.keys():
	files[repeatclass].close()
	cmd = "bedtools merge -i " + outname + '_' + repeatclass + '.bed' +  " > " + outname + '_' + repeatclass + '_merged.bed'
	subprocess.call(cmd, shell=True)


## step 2: create point BED for VCF variants

with open(outname + '_variants.bed', 'w') as outbed:
	for line in gzip.open(vcf_file, 'rt'):
		if line.startswith('#'):
			continue
		fields = line.strip().split()
		outbed.write('\t'.join([fields[0], fields[1], str(int(fields[1]) + 1)]) + '\n')


## step 3: annotate variants

annotation_beds = [outname + "_" + r + "_merged.bed" for r in repeatclasses]
annotation_names = repeatclasses

cmd = "bedtools annotate -i " + outname + "_variants.bed  -files " + " ".join(annotation_beds)  + " -names " +  " ".join(annotation_names) + " > " + outname + "_annotated.bed"
subprocess.call(cmd, shell=True)

## step 4: create variant to repeatclass mapping

pos_to_region = {}

for line in open(outname + "_annotated.bed", 'r'):
	if line.startswith('#'):
		header = line.strip().split()[1:]
		continue
	fields = line.strip().split()
	pos = (fields[0], fields[1])
	assigned = False
	for i,overlap in enumerate(fields[3:]):
		if float(overlap) > 0.0:
			assert not assigned
			pos_to_region[pos] = header[i] + '_' + fields[0]	
			assigned = True
	if not assigned:
		pos_to_region[pos] = "OTHER_" + fields[0]


## step 5: write file with all interval sizes

region_to_length = defaultdict(lambda: 0)
chr_to_len = parse_fai(fai_file)

for filename, region in zip(annotation_beds, repeatclasses):
	for line in open(filename, 'r'):
		fields = line.strip().split()
		length = int(fields[2]) - int(fields[1])
		region_to_length[region + '_' + fields[0]] += length
		chr_to_len[fields[0]] -= length

with open(outname + "_mapping.tsv", 'w') as mapping:
	for r,l in region_to_length.items():
		mapping.write("# " + str(r) + "\t" + str(l) + "\n")

	for r,l in chr_to_len.items():
		mapping.write("# OTHER_" + str(r) + "\t" + str(l) + "\n")

	for k,v in pos_to_region.items():
		mapping.write("\t".join([k[0], k[1], v]) + "\n")

import sys

genome_index = sys.argv[1]

chrom_to_length = {}
for line in open(genome_index, 'r'):
	fields = line.split()
	chrom_to_length[fields[0]] = fields[1]

counter = 0

for line in sys.stdin:
	if line.startswith('##'):
		if line.startswith('##contig=<'):
			continue
		print(line.strip())
		continue
	if line.startswith('#'):
		for chrom,length in chrom_to_length.items():
			print("##contig=<ID=" + chrom + ",length=" + length + ">")
		print(line.strip())
		continue
	fields = line.split()
	if fields[0] not in chrom_to_length:
		continue
	print(line.strip())

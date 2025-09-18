import sys

min_bubble_size = int(sys.argv[1])

for line in sys.stdin:
	if line.startswith('#'):
		continue
	fields = line.strip().split()
	chromosome = fields[0]
	allele_lengths = [len(fields[3])] + [len(a) for a in fields[4].split(',')]
	max_allele_length = max(allele_lengths)
	if max_allele_length < min_bubble_size:
		continue
	start_pos = int(fields[1]) - 1
	end_pos = start_pos + len(fields[3])
	print("\t".join([chromosome, str(start_pos), str(end_pos)]))

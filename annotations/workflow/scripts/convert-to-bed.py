import sys

print('#chrom\tstart\tend\tQV')

for line in sys.stdin:
	if line.startswith('ref'):
		continue
	fields = line.strip().split()
	chrom = fields[0].split('_')[-1].split(':')[0]
	start = fields[0].split('_')[-1].split(':')[-1].split('-')[0]
	end = fields[0].split('_')[-1].split(':')[-1].split('-')[-1]

	qv = fields[16]
	print('\t'.join([chrom, start, end, qv]))	

import sys

for line in sys.stdin:
	fields = line.strip().split()
	print('\t'.join([fields[0], fields[1], fields[2], "SegDup"]))
	print('\t'.join([fields[3], fields[4], fields[5], "SegDup"]))

import argparse

def add_columns(first, second):
	pos_to_data = {}
	header = []
	nr_new_fields = 0
	nr_second = 0
	nr_printed = 0
	for line in open(second, 'r'):
		fields = line.strip().split()
		if line.startswith('#'):
			header = fields[2:]
			nr_new_fields = len(header)
			continue
		k = (fields[0], fields[1])
		assert not k in pos_to_data
		pos_to_data[k] = fields[2:]
		nr_second += 1

	for line in open(first, 'r'):
		fields = line.strip().split()
		if line.startswith('#'):
			header = fields + header
			print('\t'.join(header))
			continue
		k = (fields[0], fields[1])
		if k in pos_to_data:
			print('\t'.join(fields + pos_to_data[k]))
			nr_printed += 1
		else:
			print('\t'.join(fields + ['None'] * nr_new_fields))
	
	assert nr_printed == nr_second


parser = argparse.ArgumentParser(prog='merge_vcfs.py', description=__doc__)
parser.add_argument('file1', metavar='FILE1', help='First TSV')
parser.add_argument('file2', metavar='FILE2', help='Second TSV')
args = parser.parse_args()

add_columns(args.file1, args.file2)

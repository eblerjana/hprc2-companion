import sys
import argparse

if __name__ == '__main__':
	parser = argparse.ArgumentParser(prog='determine-contig-names.py', description=__doc__)
	parser.add_argument('--exclude', nargs='+', default=[])
	args = parser.parse_args()

	for line in sys.stdin:
		if line.startswith('>'):
			name = line.strip().split()[0][1:]
			if not name in args.exclude:
				print(name)

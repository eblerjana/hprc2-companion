import sys
import argparse

if __name__ == '__main__':
	parser = argparse.ArgumentParser(prog='determine-contig-names.py', description=__doc__)
	parser.add_argument('--exclude', nargs='+', default=[])
	args = parser.parse_args()

	for line in sys.stdin:
		if line.startswith('>'):
			name = line.strip().split()[0][1:]
			include_contig = True
			for f in args.exclude:
				if name == f:
					include_contig = False
				if name.startswith(f + "_"):
					include_contig = False
			if include_contig:
				print(name)

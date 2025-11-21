import sys

for line in sys.stdin:
	if line.startswith("##"):
		print(line.strip())
		continue
	fields = line.strip().split()[:10]
	if line.startswith("#"):
		print("\t".join(fields[:10]))
		continue
	fields[9] = "1/1"
	print("\t".join(fields))

import sys

sample = sys.argv[1]

for line in sys.stdin:
	if line.startswith("##"):
		print(line.strip())
		continue
	if line.startswith("#"):
		fields = line.strip().split()
		print("\t".join(fields[:9]) + "\t" + sample)
		continue
	fields = line.strip().split() [:8]
	fields.append("GT")
	fields.append("./.")
	print("\t".join(fields))

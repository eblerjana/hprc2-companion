import sys

for line in sys.stdin:
	if line.startswith("##"):
		print(line.strip())
		continue
	elif line.startswith("#"):
		print("##INFO=<ID=ID,Number=A,Type=String,Description=\"Variant IDs.\">")
		print(line.strip())
		continue
	else:
		print(line.strip())

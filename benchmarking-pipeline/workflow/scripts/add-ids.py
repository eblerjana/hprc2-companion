import sys
import gzip

input_vcf = sys.argv[1]
input_sample = sys.argv[2]

var_to_id = {}

for line in gzip.open(input_vcf, 'rt'):
	if line.startswith('#'):
		continue
	fields = line.split()
	info_fields = {s.split('=')[0] : s.split('=')[1] for s in fields[7].split(';') if '=' in s}
	assert 'ID' in info_fields
	alleles = fields[4].split(',')
	# assume biallelic VCF
	assert len(alleles) == 1
	variant = (fields[0], fields[1], fields[3], alleles[0])
	assert variant not in var_to_id
	var_to_id[variant] = info_fields['ID']


for line in sys.stdin:
	if line.startswith('##'):
		print(line.strip())
		continue
	if line.startswith('#'):
		fields = line.strip().split()
		fields[-1] = input_sample
		print("\t".join(fields))
		continue
	fields = line.strip().split()
	info_fields = {s.split('=')[0] : s.split('=')[1] for s in fields[7].split(';') if '=' in s}
	assert not 'ID' in info_fields
	alleles = fields[4].split(',')
	assert len(alleles) == 1
	variant = (fields[0], fields[1], fields[3], alleles[0])
	assert variant in var_to_id
	info_fields['ID'] = var_to_id[variant]
	if fields[7] == '.':
		updated = 'ID=' + var_to_id[variant]
	else:
		updated = fields[7] + ';ID=' + var_to_id[variant]
	fields[7] = updated
	print('\t'.join(fields))

import sys
import argparse
import pandas as pd

def get_alleles(fields, present_in_ht = False):
	"""
	Get list of REF/ALT alleles.
	"""
	all_alleles = [fields[3]] + [a for a in fields[4].split(',')]
	if present_in_ht:
		filtered_alleles = []
		# only keep alleles covered by haplotypes
		for gt in fields[9:]:
			a = gt.split('|')
			for allele in a:
				if allele != '.':
					filtered_alleles.append(all_alleles[int(allele)])
		alleles = filtered_alleles
	else:
		alleles = all_alleles
	return alleles


def get_af(fields):
	"""
	Extract allele specific allele frequencies from INFO field.
	"""
	allele_to_af = {}
	info_fields = {k.split('=')[0] : k.split('=')[1] for k in fields[7].split(';') if '=' in k}
	assert 'AF' in info_fields
	total = 0.0
	for allele, freq in zip(fields[4].split(','), info_fields['AF'].split(',')):
		allele_to_af[allele] = float(freq)
		total += float(freq)
	assert total <= 1.00001
	if total > 1.0:
		# handle numerical inaccuracies
		total = 1.0
	assert total >= 0.0
	allele_to_af[fields[3]] = 1.0 - total
	return allele_to_af


def get_bubble_stats(fields):
	"""
	Get length of bubble and the total number of bubble alleles.
	"""
	all_alleles = [fields[3]] + [a for a in fields[4].split(',')]
	bubble_length = max([len(a) for a in all_alleles])
	return len(all_alleles), bubble_length


def gt_category(fields):
	"""
	Determine genotype status (assuming one sample).
	"""
	alleles = [a.strip() for a in fields[-1].split('|')]
	assert len(alleles) == 2
	if fields[-1].strip() == '0|0':
		return 'absent'
	elif alleles[0] != alleles[1]:
		return 'heterozygous'
	else:
		assert alleles[0] == alleles[1]
		return 'homozygous'


def get_info_field(fields, tag):
	"""
	Get given info field value.
	"""
	info_fields = {k.split('=')[0] : k.split('=')[1] for k in fields[7].split(';') if '=' in k}
	if tag in info_fields:
		return info_fields[tag]
	else:
		return None


def write_bubble_stats(truthfile, panelfile, outname):
	"""
	Write per bubble statistics.
	"""
	total_variants = 0
	skipped_variants = 0
	with open(outname, 'w') as outfile: 
		true_alleles = {}
		gt_cat = {}
		for line in open(truthfile, 'r'):
			if line.startswith('#'):
				continue
			fields = line.strip().split()
			pos = (fields[0], int(fields[1]))
			alleles = get_alleles(fields, True)
			total_alleles, bubble_len = get_bubble_stats(fields)
			assert len(alleles) < 3
			true_alleles[pos] = alleles
			allele_id_field = get_info_field(fields, 'ID')
			allele_ids = set([])
			for a in allele_id_field.split(','):
				for i in a.split(':'):
					allele_ids.add(i)
			if len(alleles) != 0:
				gt_cat[pos] = (gt_category(fields), total_alleles, bubble_len, ','.join(list(allele_ids)))

		outfile.write('\t'.join(['#chromosome', 'position', 'bubble_len', 'nr_bubble_alleles', 'nr_unique_kmers', 'true_genotype', 'nr_alleles_in_truth', 'nr_alleles_in_panel', 'missed_alleles_AF', 'allele_ids']) + '\n')
		nr_missed = 0
		nr_considered = 0
		for line in open(panelfile, 'r'):
			if line.startswith('#'):
				continue
			fields = line.strip().split()
			pos = (fields[0], int(fields[1]))
			alleles = get_alleles(fields, True)
			assert pos in true_alleles
			# check if both alleles are in panel
			if len(true_alleles[pos]) < 2:
				skipped_variants += 1
				continue

			total_variants += 1
			first_allele_present = true_alleles[pos][0] in alleles
			second_allele_present = true_alleles[pos][1] in alleles

			allele_frequencies = get_af(fields)
			consider_first = allele_frequencies[true_alleles[pos][0]] > 0.0
			consider_second = allele_frequencies[true_alleles[pos][1]] > 0.0

			missed_af = []
			nr_alleles_in_truth = 0
			nr_alleles_in_panel = 0

			if gt_cat[pos][0] == 'heterozygous':
				if consider_first:
					nr_alleles_in_truth += 1
					if first_allele_present:
						nr_alleles_in_panel += 1
					else:
						missed_af.append(allele_frequencies[true_alleles[pos][0]])
				else:
					nr_missed += 1

				if consider_second:
					nr_alleles_in_truth += 1
					if second_allele_present:
						nr_alleles_in_panel += 1
					else:
						missed_af.append(allele_frequencies[true_alleles[pos][1]])
				else:
					nr_missed += 1

			else:
				assert consider_first == consider_second
				assert first_allele_present == second_allele_present
				assert true_alleles[pos][0] == true_alleles[pos][1]
				if consider_first:
					nr_alleles_in_truth += 1
					if first_allele_present:
						nr_alleles_in_panel += 1
					else:
						missed_af.append(allele_frequencies[true_alleles[pos][0]])
				else:
					nr_missed += 1


			nr_considered += nr_alleles_in_truth
			missed_afs = ','.join([str(a) for a in missed_af]) if missed_af else "None"
			nr_unique_kmers = str(get_info_field(fields, 'UK'))
			outfile.write('\t'.join([pos[0], str(pos[1]), str(gt_cat[pos][2]), str(gt_cat[pos][1]), nr_unique_kmers, gt_cat[pos][0], str(nr_alleles_in_truth), str(nr_alleles_in_panel), missed_afs, gt_cat[pos][3]]) + '\n')
	print('total variants: ' + str(total_variants))
	print('skipped variants: ' + str(skipped_variants))
	print('Alleles unique to sample: ' + str(nr_missed))
	print('Alleles considered: ' + str(nr_considered))


def count(df, key, counts):
	total_p = 0
	total_t = 0
	for p,t in zip(df['nr_alleles_in_panel'], df['nr_alleles_in_truth']):
		total_p += p
		total_t += t
	counts[key][0] = total_p
	counts[key][1] = total_t

	
def write_summary_stats(filename, outname, sample, size):
	df = pd.read_csv(filename, sep='\t')
	chromosomes = ['whole_genome'] + list(df['#chromosome'].unique())

	print("Processing chromosomes: " + ','.join(chromosomes))
	with open(outname, 'w') as outfile:
		outfile.write('\t'.join(['sample', 'sampling_size', 'chromosome', 'true_genotype', 'region', 'total_alleles', 'covered_alleles', 'covered_alleles[%]']) + '\n')
		for chrom in chromosomes:
			for i, genotype in enumerate(['all', 'absent', 'heterozygous', 'homozygous']):
				category_to_counts = {}
				category_to_counts['all_bubbles'] = [0,0]
				category_to_counts['bubbles>=50bp'] = [0,0]
				category_to_counts['bubbles<50bp'] = [0,0]
				category_to_counts['multiallelic_bubbles'] = [0,0]
				category_to_counts['biallelic_bubbles'] = [0,0]
				category_to_counts['UK>0'] = [0,0]
				category_to_counts['UK=0'] = [0,0]

				df_region = df if chrom == 'whole_genome' else df[df['#chromosome'] == chrom]
				df_plot = df_region if genotype == 'all' else df_region[df_region['true_genotype'] == genotype]

				# stats for all bubbles
				df_subset = df_plot
				count(df_subset, 'all_bubbles', category_to_counts)

				# stats for bubbles >= 50bp
				df_subset = df_plot[df_plot['bubble_len'] >= 50]
				count(df_subset, 'bubbles>=50bp', category_to_counts)

				# stats for bubbles < 50bp
				df_subset = df_plot[df_plot['bubble_len'] < 50]
				count(df_subset, 'bubbles<50bp', category_to_counts)


				# stats for multiallelic bubbles
				df_subset = df_plot[df_plot['nr_bubble_alleles'] > 2]
				count(df_subset, 'multiallelic_bubbles', category_to_counts)

				# stats for biallelic bubbles
				df_subset = df_plot[df_plot['nr_bubble_alleles'] == 2]
				count(df_subset, 'biallelic_bubbles', category_to_counts)


				# stats for bubbles with unique kmers
				df_subset = df_plot[ (df_plot['nr_unique_kmers'] != 0) & (df_plot['nr_unique_kmers'] != "None") ]
				count(df_subset, 'UK>0', category_to_counts)

				# stats for bubbles with no unique kmers
				df_subset = df_plot[df_plot['nr_unique_kmers'] == 0]
				count(df_subset, 'UK=0', category_to_counts)

				categories = [k for k in category_to_counts.keys()]
				for c in categories:
					# print ratios
					ratio = category_to_counts[c][0] / category_to_counts[c][1] if category_to_counts[c][1] > 0 else 0.0
					outfile.write('\t'.join([sample, size, chrom, genotype, c, str(category_to_counts[c][1]), str(category_to_counts[c][0]), str(ratio * 100.0)]) + '\n')
	


if __name__ == '__main__':
	parser = argparse.ArgumentParser(prog='evaluate_panel.py', description=__doc__)
	parser.add_argument('-panel', metavar='PANEL', help='Sampled panel VCF.', required=True)
	parser.add_argument('-truth', metavar='TRUTH', help='Ground truth haplotypes.', required=True)
	parser.add_argument('-outname', metavar='PREFIX', help='Prefix of output files.', required=True)
	parser.add_argument('-sample', metavar='SAMPLE', help='Sample name.', required=True)
	parser.add_argument('-size', metavar='SIZE', help='Sampling size.', required=True)
	args = parser.parse_args()

	bubble_outfile = args.outname + '_bubble_' + args.sample + '_' + args.size + '.tsv'
	summary_outfile = args.outname + '_summary_' + args.sample + '_' + args.size + '.tsv'

	# write per bubble statistics
	write_bubble_stats(args.truth, args.panel, bubble_outfile)

	# write summary statistics
	write_summary_stats(bubble_outfile, summary_outfile, args.sample, args.size)


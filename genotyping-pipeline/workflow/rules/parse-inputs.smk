import os.path

# set panel. If PanGenie-ready panels are given, use these. Otherwise,
# create them from provided MC graph

PANEL_BI = config["panel_bi"]
PANEL_MULTI = config["panel_multi"]
MC_GFA = config["mc_gfa"]
MC_VCF = config["mc_vcf"]


if not PANEL_BI or not PANEL_MULTI:
	# make sure MC data is provided
	assert MC_GFA and MC_VCF, "No panels and no MC graph provided."
	PANEL_BI = None # TODO results of vcf preparation step
	PANEL_MULTI = None # TODO results of vcf preparation step

REFERENCE = config["reference"]
assert os.path.isfile(REFERENCE), "File " + REFERENCE + " does not exist."

ILLUMINA = {}
ONT = {}
assert os.path.isfile(config["sample-sheet"]), "File " + config["sample-sheet"] + " does not exist."

for line in open(config["sample-sheet"], 'r'):
	if line.startswith("#"):
		continue
	fields = line.strip().split()
	ILLUMINA[fields[1]] = fields[7]
	assert os.path.isfile(fields[7]), "File " + fields[7] + " does not exist."
	ONT[fields[1]] = fields[8]
	assert os.path.isfile(fields[8]), "File " + fields[8] + " does not exist."

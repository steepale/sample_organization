import sys
import re

# infile
infile = sys.argv[1]

# outfile
outfile = open(sys.argv[2], 'w')

# Reference files
dnaseq_list_file = '/Users/Alec/Documents/Bioinformatics/MDV_Project/databases/samples/all_samples_dnaseq_2014_017NNN-N.txt'
rnaseq_list_file = '/Users/Alec/Documents/Bioinformatics/MDV_Project/databases/samples/all_samples_rnaseq_2014_017NNN-N.txt'
microarray_list_file = '/Users/Alec/Documents/Bioinformatics/MDV_Project/databases/samples/all_samples_microarrays_2014_017NNN-N.txt'
cyto_list_file = '/Users/Alec/Documents/Bioinformatics/MDV_Project/databases/samples/all_samples_cytogenetics_2014_017NNN-N.txt'
immunohisto_list_file = '/Users/Alec/Documents/Bioinformatics/MDV_Project/databases/samples/all_samples_immunohistochem_2017.txt'

# Create a dictionary for birds with mutliple high quality tumors
birds2tumor_scores = {}
# Iterate through file to contruct dictionary
for line in open(infile):
	if line[0] != '#':
		line = line.rstrip()
		col = line.split('\t')
		bird_id = col[4]
		tumor_score = col[11]
		if bird_id not in birds2tumor_scores.keys():
			birds2tumor_scores[bird_id] = [tumor_score]
		else:
			birds2tumor_scores[bird_id].append(tumor_score)

# Read samples into set that have undergone each specific molecular test
# make all sets into dictionary
test_dic = {}
for test_file in [dnaseq_list_file, rnaseq_list_file, microarray_list_file, cyto_list_file, immunohisto_list_file]:
	for line in open(test_file):
		if line[0] != '#':
			line = line.rstrip()
			sample = line
			for test in ['dnaseq', 'rnaseq', 'microarrays', 'cytogenetics', 'immunohistochem']:
				if re.search(test, test_file):
					if test not in test_dic.keys():
						test_dic[test] = set()
						test_dic[test].add(sample)
					else:
						test_dic[test].add(sample)

# Iterate through infile and collect stats for filters
for line in open(infile):
	# Write a header
	if line [0] == '#' and line[1] != '#':
		outfile.write('##SAMPLE_ID'+'\t'+'Collected tissue sample identification'+'\n')
		outfile.write('##TUBE_LABEL'+'\t'+'How sample was labeled on storage tube'+'\n')
		outfile.write('##SEX'+'\t'+'The sex of the bird'+'\n')
		outfile.write('##GERM_TUM'+'\t'+'Whether sample is characterized as germline or somatic'+'\n')
		outfile.write('##TISSUE'+'\t'+'Type of collected tissue'+'\n')
		outfile.write('##TUMOR_NUM'+'\t'+'Number of COLLECTED tumors from bird'+'\n')
		outfile.write('##TUMOR_SCORE'+'\t'+'Grade of tumor; Scale of 1-3, 3 represents most homologous largest tumor in appearence, interpretation by eye'+'\n')
		outfile.write('##COHORT'+'\t'+'Experimental cohort'+'\n')
		outfile.write('##STORAGE'+'\t'+'The storage strategy of tissue based on potential for biological sampling; DNA, RNA, cytogenetics, liquid_nitrogen, immunohistochem'+'\n')
		outfile.write('##SAMPLING'+'\t'+'The actual biological tests that have been performed from sample; dna_seq, affy_array, rna_seq, cytogenetics, immunohistochem'+'\n')
		outfile.write('##FILTER'+'\t'+'Justification for priority score'+'\n')
		outfile.write('###FILTERS:'+'\n')
		outfile.write('###filter_1:'+'\t'+'All samples that underwent DNA-sequencing from Summer 2014 cohort (plus matching germlines)'+'\n')
		outfile.write('###filter_2:'+'\t'+'All samples that underwent RNA-sequencing from Summer 2014 cohort'+'\n')
		outfile.write('###filter_3:'+'\t'+'Gonadal tumors with tumor score greater than or equal to 2'+'\n')
		outfile.write('###filter_4:'+'\t'+'The tumors from birds that possess multiple high quality tumors'+'\n')
		outfile.write('###filter_5:'+'\t'+'Tumors that have a tumor score of 3 regardless of tissue'+'\n')
		outfile.write('###filter_6:'+'\t'+'Tumors that have a tumor score of 2 regardless of tissue'+'\n')
		outfile.write('###filter_7:'+'\t'+'Remaining gonadal tumor samples'+'\n')
		outfile.write('###filter_8:'+'\t'+'Remaining tumor samples'+'\n')
		outfile.write('###filter_9:'+'\t'+'Remaining germline samples'+'\n')
		outfile.write('#SAMPLE_ID'+'\t'+'TUBE_LABEL'+'\t'+'SEX'+'\t'+'GERM_TUM'+'\t'+'TISSUE'+'\t'+'TUMOR_NUM'+'\t'+'TUMOR_SCORE'+'\t'+'COHORT'+'\t'+'STORAGE'+'\t'+'SAMPLING'+'\t'+'FILTER'+'\n')
	elif line[0] != '#':
		line = line.rstrip()
		col = line.split('\t')
		sample_id = col[0]
		tube_label = col[1]
		germ_tum = col[2]
		tissue = col[3]
		bird_id = col[4]
		pen_id = col[5]
		dod = col[6]
		mort_type = col[7]
		sex = col[8]
		tumor_num = col[9]
		tumor_score = col[11]
		cohort = col[12]
		# Adjust storage
		storage = set(col[10].split(';'))
		if 'immunohistochem' in storage:
			storage.discard('immunohistochem')
		if len(storage) > 0:
			storage = ';'.join(map(str,storage))
			if storage == '':
				storage = 'na'
		# Create a set for sampling
		sampling = set()
		for test in ['dnaseq', 'rnaseq', 'microarrays', 'cytogenetics', 'immunohistochem']:
			if sample_id in test_dic[test]:
				sampling.add(test)
		if len(sampling) > 0:
			sampling = ';'.join(map(str,sampling))
		else:
			sampling = 'na'
		# Create an if then for tumor score exception
		if tumor_score == 'na':
			tumor_score_com = '0'
		else:
			tumor_score_com = tumor_score
		# Create a field for birds with multiple high quality tumors
		score_sheet = []
		if len(birds2tumor_scores[bird_id]) >= 3:
			for score in birds2tumor_scores[bird_id]:
				if score != 'na':
					if int(score) >= 2:
						score_sheet.append(score)
		if len(score_sheet) >= 3:
			heavy_hit_bird = 'yes'
		else:
			heavy_hit_bird = 'no'
		# Create conditions for filters
		# Samples already output
		output_set = set()
		# filter_1: All samples that underwent DNA-sequencing from Summer 2014 cohort (plus matching germlines)
		if 'dnaseq' in sampling.split(';') and sample_id not in output_set:
			output_set.add(sample_id)
			filter_status = '1'
		# filter_2: All samples that underwent RNA-sequencing from Summer 2014 cohort 
		elif 'rnaseq' in sampling.split(';') and germ_tum != 'germline' and sample_id not in output_set:
			output_set.add(sample_id)
			filter_status = '2'
		# filter_3: Gonadal tumors with tumor score greater than or equal to 2 
		elif tissue == 'gonad' and germ_tum != 'germline' and int(tumor_score_com) >= 2 and sample_id not in output_set:
			output_set.add(sample_id)
			filter_status = '3'
		# filter_4: The tumors from birds that possess multiple high quality tumors 
		elif heavy_hit_bird == 'yes' and germ_tum != 'germline' and sample_id not in output_set:
			output_set.add(sample_id)
			filter_status = '4'
		# filter_5: Tumors that have a tumor score of 3 regardless of tissue 
		elif tumor_score == '3' and germ_tum != 'germline' and sample_id not in output_set:
			output_set.add(sample_id)
			filter_status = '5'
		# filter_6: Tumors that have a tumor score of 2 regardless of tissue 
		elif tumor_score == '2' and germ_tum != 'germline' and sample_id not in output_set:
			output_set.add(sample_id)
			filter_status = '6'
		# filter_7: Remaining gonadal tumors 
		elif tissue == 'gonad' and germ_tum != 'germline' and sample_id not in output_set:
			output_set.add(sample_id)
			filter_status = '7'	
		# filter_8: Remaining tumors samples 
		elif germ_tum != 'germline' and sample_id not in output_set:
			output_set.add(sample_id)
			filter_status = '8'
		# filter_9: Remaining samples 
		elif sample_id not in output_set:
			output_set.add(sample_id)
			filter_status = '9'
		# Write to outfile
		outfile.write(sample_id+'\t'+tube_label+'\t'+sex+'\t'+germ_tum+'\t'+tissue+'\t'+tumor_num+'\t'+tumor_score+'\t'+cohort+'\t'+storage+'\t'+sampling+'\t'+filter_status+'\n')
#Close outfile
outfile.close()

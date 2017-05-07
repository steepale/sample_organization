import sys
import re

# infile
infile = sys.argv[1]

# outfile
outfile = open(sys.argv[2], 'w')

# Iterate through the infile, collect stats, and reorganize data
for line in open(infile):
	# Write output header
	if line[0] == '#':
		outfile.write('##SAMPLE_ID'+'\t'+'Collected tissue sample identification'+'\n')
		outfile.write('##TUBE_LABEL'+'\t'+'How sample was labeled on storage tube'+'\n')
		outfile.write('##GERM_TUM'+'\t'+'Whether sample is characterized as germline or somatic'+'\n')
		outfile.write('##TISSUE'+'\t'+'Type of collected tissue'+'\n')
		outfile.write('##BIRD_ID'+'\t'+'Identification of bird, from which tissue was collected'+'\n')
		outfile.write('##PEN_ID'+'\t'+'Identification of housing unit from which bird resided during experiment'+'\n')
		outfile.write('##DOD'+'\t'+'Date Of Death of bird'+'\n')
		outfile.write('##MORT_TYPE'+'\t'+'Mortality type, whether bird died from Mareks Disease or was euthanized'+'\n')
		outfile.write('##SEX'+'\t'+'The sex of the bird'+'\n')
		outfile.write('##TUMOR_NUM'+'\t'+'Number of COLLECTED tumors from bird'+'\n')
		outfile.write('##SAMPLING'+'\t'+'The storage strategy of tissue based on potential for biological sampling; DNA, RNA, cytogenetics'+'\n')
		outfile.write('##TUMOR_SCORE'+'\t'+'Grade of tumor; Scale of 1-3, 3 represents most homologous largest tumor in appearence, interpretation by eye'+'\n')
		outfile.write('#SAMPLE_ID'+'\t'+'TUBE_LABEL'+'\t'+'GERM_TUM'+'\t'+'TISSUE'+'\t'+'BIRD_ID'+'\t'+'PEN_ID'+'\t'+'DOD'+'\t'+'MORT_TYPE'+'\t'+'SEX'+'\t'+'TUMOR_NUM'+'\t'+'SAMPLING'+'\t'+'TUMOR_SCORE'+'\n')
	# collect stats if not header
	elif line[0] != '#':
		line = line.rstrip()
		col = line.split('\t')
		in_sample = col[0]
		dod = col[1]
		in_sex = col[2]
		if in_sex == 'F':
			sex = 'female'
		elif in_sex == 'M':
			sex = 'male'
		in_tumor = col[3]
		in_tumor_score = col[4]
		in_normal = col[5]
		in_oct = col[6]
		in_notes = col[7]
		bird_id = 'S'+ in_sample.zfill(3)
		tube_label = 'na'
		pen_id = 'na'
		mort_type = 'euthanized'
		# Create lists of all delimited fields
		tumor_list = in_tumor.split(';')
		tumor_score_list = in_tumor_score.split(';')
		normal_list = in_normal.split(';')
		# continue variable creation
		tumor_num = str(len(tumor_list) + 1)
		# Iterate through the germline samples
		for g_i, g in enumerate(normal_list):
			sample_id = bird_id + '-N_' + str(g_i)
			germ_tum = 'germline'
			tissue = g.lower()
			sampling_list = ['liquid_nitrogen']
			if in_oct != 'NA' and re.search(g, in_oct) and len(sampling_list) != 2:
				sampling_list.append('immunohistochem')
			tumor_score = 'na'
			sampling = ';'.join(map(str,sampling_list))
			# Write the data to outfile
			outfile.write(sample_id+'\t'+tube_label+'\t'+germ_tum+'\t'+tissue+'\t'+bird_id+'\t'+pen_id+'\t'+dod+'\t'+mort_type+'\t'+sex+'\t'+tumor_num+'\t'+sampling+'\t'+tumor_score+'\n')
		# Iterate through the tumor samples
		for t_i, t in enumerate(tumor_list):
			sample_id = bird_id + '-T_' + str(t_i + 1)
			germ_tum = 'tumor'
			tissue = t.lower()
			sampling_list = ['liquid_nitrogen']
			if in_oct != 'NA' and re.search(t, in_oct) and len(sampling_list) != 2:
				sampling_list.append('immunohistochem')
			if tumor_score_list[t_i] == 'NA':
				tumor_score_list[t_i] == 'na'
			else:
				tumor_score = tumor_score_list[t_i]
			sampling = ';'.join(map(str,sampling_list))
			# Write the data to outfile
			outfile.write(sample_id+'\t'+tube_label+'\t'+germ_tum+'\t'+tissue+'\t'+bird_id+'\t'+pen_id+'\t'+dod+'\t'+mort_type+'\t'+sex+'\t'+tumor_num+'\t'+sampling+'\t'+tumor_score+'\n')
# Close outfile
outfile.close()

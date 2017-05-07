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
		outfile.write('##SAMPLING'+'\t'+'The biological data that was collected from tissue; DNA, RNA, cytogenetics'+'\n')
		outfile.write('##TUMOR_SCORE'+'\t'+'Grade of tumor; Scale of 1-3, 3 represents most homologous largest tumor in appearence, interpretation by eye'+'\n')
		outfile.write('#SAMPLE_ID'+'\t'+'TUBE_LABEL'+'\t'+'GERM_TUM'+'\t'+'TISSUE'+'\t'+'BIRD_ID'+'\t'+'PEN_ID'+'\t'+'DOD'+'\t'+'MORT_TYPE'+'\t'+'SEX'+'\t'+'TUMOR_NUM'+'\t'+'SAMPLING'+'\t'+'TUMOR_SCORE'+'\n')
	# collect stats if not header
	elif line[0] != '#':
		line = line.rstrip()
		col = line.split('\t')
		bird_id = col[0]
		pen_id = col[1]
		dod = col[2]
		mort_type = col[3]
		sex = col[4].lower()
		tumor_pres = col[5]
		c1 = col[6]
		c1_dna = col[7]
		c1_rna = col[8]
		t1 = col[9]
		t1_dna = col[10]
		t1_rna = col[11]
		t1_cyto = col[12]
		t1_score = col[13]
		t2 = col[14]
		t2_dna = col[15]
		t2_rna = col[16]
		t2_cyto = col[17]
		t2_score = col[18]
		t3 = col[19]
		t3_dna = col[20]
		t3_rna = col[21]
		t3_cyto = col[22]
		t3_score = col[23]
		t4 = col[24]
		t4_dna = col[25]
		t4_rna = col[26]
		t4_cyto = col[27]
		t4_score = col[28]
		t5 = col[29]
		t5_dna = col[30]
		t5_rna = col[31]
		t5_cyto = col[32]
		t5_score = col[33]
		t6 = col[34]
		t6_dna = col[35]
		t6_rna = col[36]
		t6_cyto = col[37]
		t6_score = col[38]
		t7 = col[39]
		t7_dna = col[40]
		t7_rna = col[41]
		t7_cyto = col[42]
		t7_score = col[43]
		# Create a list of tumor calls
		t_list = [t1, t2, t3, t4, t5, t6, t7]
		#Create a list of tumors
		tumors = []
		for t in t_list:
			if t != 'NA':
				tumors.append(t)
		# Create a list of tumor dna sampling
		tum_dna = [t1_dna, t2_dna, t3_dna, t4_dna, t5_dna, t6_dna, t7_dna]
		# Create a list of tumor rna sampling
		tum_rna = [t1_rna, t2_rna, t3_rna, t4_rna, t5_rna, t6_rna, t7_rna]
		# Create a list of tumor cyto sampling
		tum_cyto = [t1_cyto, t2_cyto, t3_cyto, t4_cyto, t5_cyto, t6_cyto, t7_cyto]
		# Create a list of tumor scores
		tum_score = [t1_score, t2_score, t3_score, t4_score, t5_score, t6_score, t7_score]
		# Outfile fields
		tumor_num = len(tumors)
		# Iterate through controls and write them to output file
		for c_i, control in enumerate([c1]):
			sample_id = bird_id + '-' + str(c_i)
			# Filter out samples with no controls
			if control != 'NA':	
				# Tube label
				if re.search(';', control):
					label_end = control.split(';')[1]
					tissue = control.split(';')[0].lower()
					tube_label = bird_id + '-' + label_end
				else:
					label_end = ''
					tissue = control.lower()
					tube_label = bird_id
				germ_tum = 'germline'
				# Biological sampling strategy
				bio_data = []
				if c1_dna == 'X':
					bio_data.append('dna')
				if c1_rna == 'X':
					bio_data.append('rna')
				bio_data = ';'.join(map(str,bio_data))
				tumor_score = 'na'
				# Mortality type
				if mort_type == 'K':
					mort_type = 'euthanized'
				elif mort_type == 'D':
					mort_type = 'disease'
				# Write germline outputs
				outfile.write(sample_id+'\t'+tube_label+'\t'+germ_tum+'\t'+tissue+'\t'+bird_id+'\t'+pen_id+'\t'+dod+'\t'+mort_type+'\t'+sex+'\t'+str(tumor_num)+'\t'+bio_data+'\t'+tumor_score+'\n')
		# Iterate through tumors and write them to output file
		for t_i, tumor in enumerate(t_list):
			sample_id = bird_id + '-' + str(t_i + 1)
			# Filter out samples with no controls
			if tumor != 'NA':	
				# Tube label
				if re.search(';', tumor):
					label_end = tumor.split(';')[1]
					tissue = tumor.split(';')[0].lower()
					tube_label = bird_id + '-' + label_end
				else:
					label_end = str(t_i + 1)
					tissue = tumor.lower()
					tube_label = bird_id + '-' + label_end
				germ_tum = 'tumor'
				# Biological sampling strategy
				bio_data = {}
				for dna_i, t_dna in enumerate(tum_dna): 
					if t_dna == 'X':
						if dna_i in bio_data.keys():
							bio_data[dna_i].append('dna')
						elif dna_i not in bio_data.keys():
							bio_data[dna_i] = ['dna']
				for rna_i, t_rna in enumerate(tum_rna): 
					if t_rna == 'X':
						if rna_i in bio_data.keys():
							bio_data[rna_i].append('rna')
						elif rna_i not in bio_data.keys():
							bio_data[rna_i] = ['rna']
				for cyto_i, t_cyto in enumerate(tum_cyto): 
					if t_cyto == 'X':
						if cyto_i in bio_data.keys():
							bio_data[cyto_i].append('cyto')
						elif cyto_i not in bio_data.keys():
							bio_data[cyto_i] = ['cyto']
				# Extract the sampling info from the bio_data dictionary based on tumor sample
				if t_i in bio_data.keys():
					bio_data = ';'.join(map(str,bio_data[t_i]))
				else:
					bio_data = 'na'
				# Create dictionary for tumor scores
				score_dic = {}
				for score_i, t_score in enumerate(tum_score):
					if t_score != 'NA':
						if score_i in score_dic.keys():
							score_dic[score_i].append(t_score)
						elif score_i not in score_dic.keys():
							score_dic[score_i] = t_score
				if t_i in score_dic.keys():
					tumor_score = score_dic[t_i]
				else:
					tumor_score = 'na'
				# Mortality type
				if mort_type == 'K':
					mort_type = 'euthanized'
				elif mort_type == 'D':
					mort_type = 'disease'
				# Write germline outputs
				outfile.write(sample_id+'\t'+tube_label+'\t'+germ_tum+'\t'+tissue+'\t'+bird_id+'\t'+pen_id+'\t'+dod+'\t'+mort_type+'\t'+sex+'\t'+str(tumor_num)+'\t'+bio_data+'\t'+tumor_score+'\n')


# Close outfile
outfile.close()

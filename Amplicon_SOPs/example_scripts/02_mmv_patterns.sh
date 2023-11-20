#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#        Step 3:   Rename files
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

#Generate the following list based on the last 2 columns in your samples.txt file.
#Replace all '\t' with '.[12].fastq.gz' '
#Replace all '\n' with '-R#1.fastq.gz' \nmmv ''
#---------------------------------------------------
#Golay_L1-Golay_R1.[12].fastq.gz	GUT073-R#1.fastq.gz
#Golay_L1-Golay_R2.[12].fastq.gz	GUT056-R#1.fastq.gz
#Golay_L1-Golay_R3.[12].fastq.gz	GUT053-R#1.fastq.gz
#Golay_L1-Golay_R4.[12].fastq.gz	GUT012-R#1.fastq.gz
#Golay_L1-Golay_R5.[12].fastq.gz	GUT043-R#1.fastq.gz
#...etc
#---------------------------------------------------

mmv 'Golay_L1-Golay_R13.[12].fastq.gz' 'sampleID1-R#1.fastq.gz'
mmv 'Golay_L2-Golay_R2.[12].fastq.gz' 'sampleID2-R#1.fastq.gz'
mmv 'Golay_L3-Golay_R14.[12].fastq.gz' 'sampleID3-R#1.fastq.gz'
mmv 'Golay_L1-Golay_R4.[12].fastq.gz' 'sampleID4-R#1.fastq.gz'
mmv 'Golay_L3-Golay_R4.[12].fastq.gz' 'sampleID5-R#1.fastq.gz'

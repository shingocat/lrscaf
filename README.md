LRScaf: improving draft genomes using long noisy reads

Hybrid assembly strategy is a reasonable and promising approach to utilize strengths and settle weaknesses in Next-Generation Sequencing (NGS) and Third-Generation Sequencing (TGS) technologies. According to this principle, we here present a new toolkit named LRScaf (Long Reads Scaffolder) by applied TGS data to improve draft genome assembly. The main features are: short running time, accuracy, and being contiguity. To scaffold rice genome, it could be done in 20 mins with minimap mapper.

############

Installig requirements:

############

java version: 1.8+.
the related libarary under lib folder.
ant.
############

To build hass project:

############

downlaod the source code and the required libaray in lib. Then use maven to build this project.

############

LRScaf Parameter:

############

Parameter	Abbreviation	XML Code	Details

xml	x	NA	The parameter XML file! All command-line parameters would be omitted if set.

contig	c	contig	The file of pre-assembled contigs in fasta format.

m5	m5	m5	The aligned file in -m 5 format of BLASR.

m4	m4	m4	The aligned file in -m 4 format of BLASR.

sam	sam	sam	The aligned file in sam format.

mm	mm	mm	The aligned file in Minimap format

output	o	output	The scaffolding output folder.

help	h	NA	The help information.

miniCntLen	micl	min_contig_length	The minimum contigs length for scaffolding! Default: <3000>.

miniPBLen	mipl	min_pacbio_length	The minimum PacBio long read's length for scaffolding! Default: <5000>.

identity	i	identity	The identity threshold for BLASR alignment! Default: <0.8>.

miniOLLen	mioll	min_overlap_length	The minimum overlap length threshold for BLASR alignment! Default: <2400>.

miniOLRatio	miolr	min_overlap_ratio	The minimum overlap ratio threshold for BLASR alignment! Default: <0.8>.

maOHLen	maohl	max_overhang_length	The maximum overhang length threshold for BLASR alignment! Default: <300>.

maOHRatio	maohr	max_overhang_ratio	The maximum overhang ratio threshold for BLASR alignment! Default: <0.1>.

maELen	mael	max_end_length	The maximum ending length of PacBio® Long Read! Default: <300>.

maERatio	maer	max_end_ratio	The maximum ending ratio of PacBio® Long Read! Default: <0.1>.

miSLN	misl	min_supported_links	The minimum support links number! Default: <1>.

doll	doll	use_overlap_link	The indicator for using overlap relationship of contig! Default: , discard overlap relationship!

ratio	r	ratio	The ratio for deleting error prone edges! Default: <0.2>.

mr	mr	repeat_mask	The indicator for masking repeats! Default: .

tiplength	tl	tip_length	The maximum tip length.

iqrtime	iqrt	iqr_time	The IQR time for defined repeat outlier.

mmcm	mmcm	mmcm	The parameter to remove Minimap aligned records. Default: <8>.

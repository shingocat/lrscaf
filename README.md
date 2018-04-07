<h1>LRScaf: improving draft genomes using long noisy reads</h1>

<p>Hybrid assembly strategy is a reasonable and promising approach to utilize strengths and settle weaknesses in Next-Generation Sequencing (NGS) and Third-Generation Sequencing (TGS) technologies. According to this principle, we here present a new toolkit named LRScaf (Long Reads Scaffolder) by applied TGS data to improve draft genome assembly. The main features are: short running time, accuracy, and being contiguity. To scaffold rice genome, it could be done in 20 mins with minimap mapper.</p>

#######################################################<br>
Requirements:
#######################################################<br>
Java version: 1.8+.

#######################################################<br>
To build LRScaf:
#######################################################<br>
1. There is a jar package named LRScaf-<version>.jar under target folder. User could download it and run it with command "java -jar LRScafp-<version>.jar -x <configure.xml>"
2. If you need to compile downlaod the source code and the required library in lib. Then use maven to build this project.

#######################################################<br>
Parameters of LRScaf
#######################################################<br>
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

#######################################################<br>
License
#######################################################<br>
LRScaf is free software: you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by the
Free Software Foundation, either version 3 of the License, or (at your
option) any later version.

This program is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
General Public License for more details.

You should have received a copy of the GNU General Public License along
with this program. If not, see <http://www.gnu.org/licenses/>.

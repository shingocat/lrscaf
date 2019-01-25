<b>LRScaf: improving draft genomes using long noisy reads</b>
<br>
<p>
Hybrid assembly strategy is a reasonable and promising approach to utilize strengths and settle weaknesses in Next-Generation Sequencing (NGS) and Third-Generation Sequencing (TGS) technologies. According to this principle, we here present a new toolkit named LRScaf (Long Reads Scaffolder) by applied TGS data to improve draft genome assembly. The main features are: short running time, accuracy, and being contiguity. To scaffold rice genome, it could be done in 20 mins with minimap mapper. In human, LRScaf could improve the draft assembly NG50 from 127.5 Kb to 10.4 Mb on 20x PacBio CHM1 dataset and NG50 from 115.7 Kb to 17.4 Mb on ~35x Nanopore NA12878 dataset.
</p>
<br>
################################################################################<br>
<b>Requirements</b><br>
################################################################################<br>
Java version: 1.8+.<br>
<br>
################################################################################<br>
<b>Building LRScaf project</b><br>
################################################################################<br>
There are two ways to build and run this project:
<li>There is a jar package named LRScaf-&ltversion&gt.jar under target folder in releases. User could run it with command: <i>"java -jar LRScaf-&ltversion&gt.jar -x &ltconfigure.xml&gt"</i>. The details of configuration XML file are described below.</li>
<li>If you want to compile the source code by yourself, you could download the source code and then compile and build this project by maven &lthttps://maven.apache.org/&gt in following steps:</li>
# 1. download the latest releases version and unzip the package<br>
<i>unzip lrscaf-&ltversion&gt.zip</i><br>
# 2. change the working folder<br>
<i>cd lrscaf-&ltversion&gt</i><br>
# 3. complie source code and package the project, and a jar package named LRScaf-&ltversion&gt.jar would be under the target folder.<br>
<i>mvn package</i><br>
<br>
################################################################################<br>
<b>Quick starting</b><br>
################################################################################<br>
# XML configuration style<br>
<i>java -jar LRScaf-&ltversion&gt.jar -x &ltconfigure.xml&gt</i><br>
# or command-line in short style<br>
<i>java -jar LRScaf-&ltversion&gt.jar -c &ltdraft_assembly.fasta&gt -a &ltalignment.m4&gt -t &ltm4&gt -o &ltoutput_foloder&gt [options]</i><br>
# or command-line in long style<br>
<i>java -jar LRScaf-&ltversion&gt.jar --contig &ltdraft_assembly.fasta&gt --alignedFile &ltalignment.m4&gt -t &ltm4&gt --output &ltoutput_foloder&gt [options]</i><br>
<br>
################################################################################<br>
<b>Parameters of LRScaf</b><br>
################################################################################<br>
<p>LRScaf supports parameters set by XML confiuration file or command-line. It recommends to use XML configuration file. There is a template configuration file of XML format, named "scafconf.xml", in the project. In command-line, LRScaf supports long (dash-dash) and short (dash) style of GNU like options. And the following table would show each parameter meaning and default value if available.</p>
<p>The first and second columns are the command-line paremeters in long and its coressponding short style.</p>
<p>The third column is the code in XML configuration file. NA is not available in XML configuration file.</p>
<p>The fourth column is the details and default value of this option if available.</p>
<table>
<tr><th>Parameter</th><th>Abbreviation</th><th>XML Code</th><th>Details</th></tr>
<tr><td>xml</td><td>x</td><td>NA</td><td>The XML configuration file. All command-line parameters would be omitted if this is set.</td></tr>
<tr><td>contig</td><td>c</td><td>contig</td><td>The contigs file of draft assembly in fasta format.</td></tr>
<tr><td>m5</td><td>m5</td><td>m5</td><td>The alignment file in -m 5 format of BLASR.</td></tr>
<tr><td>m4</td><td>m4</td><td>m4</td><td>The alignment file in -m 4 format of BLASR.</td></tr>
<tr><td>sam</td><td>sam</td><td>sam</td><td>The alignment file in sam format of BLASR.</td></tr>
<tr><td>mm</td><td>mm</td><td>mm</td><td>The alignment file in PAF format of Minimap.</td></tr>
<tr><td>output</td><td>o</td><td>output</td><td>The output folder.</td></tr>
<tr><td>miniCntLen</td><td>micl</td><td>min_contig_length</td><td>The minimum contigs length to be included for scaffolding. Default: &lt500&gt bp.</td></tr>
<tr><td>identity</td><td>i</td><td>identity</td><td>The identity threshold for filtering invalid alignment. Default: &lt0.8&gt.<br>This value <b>must be</b> modify according to the mapper. <br>For the BLASR alignment file, the higher value means the higher identity. <br>For the Minimap alignment file, the value should not be larger than 0.3 and the value could be set to 0.1.</td></tr>
<tr><td>miniOLLen</td><td>mioll</td><td>min_overlap_length</td><td>The minimum overlap length of contig. Default: &lt400&gt bp.</td></tr>
<tr><td>miniOLRatio</td><td>miolr</td><td>min_overlap_ratio</td><td>The minimum overlap length ratio of contig. Default: &lt0.8&gt.<br>If the overlap length is large than the miniOLLen, <br>it will compute the ratio of overlap length which is overlap_length/contig_length.</td></tr>
<tr><td>maOHLen</td><td>maohl</td><td>max_overhang_length</td><td>The maximum overhang length of contig. Default: &lt500&gt bp.</td></tr>
<tr><td>maOHRatio</td><td>maohr</td><td>max_overhang_ratio</td><td>The maximum overhang ratio of contig. Default: &lt0.1&gt. <br>If the overhang length is less than the maohl, <br>it will compute the ratio of overhang length which is overhang_lenght/contig_length.</td></tr>
<tr><td>maELen</td><td>mael</td><td>max_end_length</td><td>The maximum ending length of long read. Default: &lt500&gt bp.</td></tr>
<tr><td>maERatio</td><td>maer</td><td>max_end_ratio</td><td>The maximum ending ratio of long read. Default: &lt0.1&gt.<br>It will compute the ending length (ending_len) by long_read_length * maer, <br>then def_ending_len = (mael >= ending_len ? ending_len : mael).</td></tr>
<tr><td>miSLN</td><td>misl</td><td>min_supported_links</td><td>The minimum support links. Default: &lt2&gt. <br>If the depth of long reads less than 10x, the misl could be set to 1.</td></tr>
<tr><td>ratio</td><td>r</td><td>ratio</td><td>The ratio for deleting error prone edges in divergence nodes. Default: &lt0.2&gt.</td></tr>
<tr><td>mr</td><td>mr</td><td>repeat_mask</td><td>The indicator for masking repeats. Default: &lttrue&gt. It recommends to be true.</td></tr>
<tr><td>tiplength</td><td>tl</td><td>tip_length</td><td>The maximum tip length. Default: &lt1500&gt bp.</td></tr>
<tr><td>iqrtime</td><td>iqrt</td><td>iqr_time</td><td>The IQR times for setting contigs as repeats by their coverages. Default: &lt3&gt.</td></tr>
<tr><td>mmcm</td><td>mmcm</td><td>mmcm</td><td>The parameter to filter invalid Minimap alignments. Default: &lt8&gt. <b>Only for Minimap alignment</b>.</td></tr>
<tr><td>process</td><td>p</td><td>process</td><td>The multi-threads settings. Default:&lt4&gt.</td></tr>
<tr><td>help</td><td>h</td><td>NA</td><td>Print this help information.</td></tr>
</table>
<br>
################################################################################<br>
<b>XML Configuration File Content</b><br>
################################################################################<br>
&lt?xml version="1.0" encoding="UTF-8"?&gt<br>
&ltscaffold&gt<br>
&nbsp;  &lt!--The input file for scaffolding, including contigs and aligned files (i.e. m5, m4 or mm file) --&gt<br>
&nbsp; 	&ltinput&gt<br>
&nbsp; &nbsp; &ltcontig&gtDraft assembly in fasta format.&lt/contig&gt<br>
&nbsp; &nbsp; &ltm4&gtThe aligned file in BLASR -m 4 format.&lt/m4&gt<br>
&nbsp; 	&lt/input&gt<br>
&nbsp; 	&lt!-- The output folder for scaffolding --&gt<br>
&nbsp; 	&ltoutput&gtThe output folder.&lt/output&gt<br>
&nbsp; 	&lt!-- The parameters for scaffolding--&gt<br>
&nbsp; 	&ltparas&gt<br>
&nbsp; &nbsp; &lt!--More details are showed in README.md--&gt<br>
&nbsp; &nbsp; &ltmin_contig_length&gt500&lt/min_contig_length&gt<br>
&nbsp; &nbsp; &ltidentity&gt0.8&lt/identity&gt<br>
&nbsp; &nbsp; &ltmin_overlap_length&gt400&lt/min_overlap_length&gt<br>
&nbsp; &nbsp; &ltmin_overlap_ratio&gt0.8&lt/min_overlap_ratio&gt<br>
&nbsp; &nbsp; &ltmax_overhang_length&gt500&lt/max_overhang_length&gt<br>
&nbsp; &nbsp; &ltmax_overhang_ratio&gt0.1&lt/max_overhang_ratio&gt<br>
&nbsp; &nbsp; &ltmax_end_length&gt500&lt/max_end_length&gt<br>
&nbsp; &nbsp; &ltmax_end_ratio&gt0.1&lt/max_end_ratio&gt<br>
&nbsp; &nbsp; &ltmin_supported_links&gt2&lt/min_supported_links&gt<br>
&nbsp; &nbsp; &lttips_length&gt1500&lt/tips_length&gt<br>
&nbsp; &nbsp; &ltratio&gt0.2&lt/ratio&gt<br>
&nbsp; &nbsp; &ltrepeat_mask&gttrue&lt/repeat_mask&gt<br>
&nbsp; &nbsp; &ltiqr_time&gt3&lt/iqr_time&gt<br>
&nbsp; &nbsp; &ltmmcm&gt8&lt/mmcm&gt &lt!--only for Minimap Alignment.--&gt<br>
&nbsp; &nbsp; &ltprocess&gt4&lt/process&gt<br>
&nbsp; 	&lt/paras&gt<br>
&lt/scaffold&gt<br>
<br>
################################################################################<br>
<b>License</b><br>
################################################################################<br>
<p>LRScaf is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
<br>
This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
<br>
You should have received a copy of the GNU General Public License along with this program. If not, see &lthttp://www.gnu.org/licenses/&gt. </p>
<br>
If you have any questions, please feel free to contact me &ltqinmao@caas.cn&gt.
<br>
<br>
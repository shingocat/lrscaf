<b>LRScaf: improving draft genomes using long noisy reads</b>
<br>
<p>
Hybrid assembly strategy is a reasonable and promising approach to utilize strengths and settle weaknesses in Next-Generation Sequencing (NGS) and Third-Generation Sequencing (TGS) technologies. According to this principle, we here present a new toolkit named LRScaf (Long Reads Scaffolder) by applied TGS data to improve draft genome assembly. The main features are: short running time, accuracy, and being contiguity. To scaffold rice genome, it could be done in 20 mins with minimap mapper.
</p>
<br>
################################################################################<br>
<b>Requirements:</b><br>
################################################################################<br>
Java version: 1.8+.<br>
<br>
################################################################################<br>
<b>Building LRScaf project:</b><br>
################################################################################<br>
There are two ways to build and run this project:
<li>There is a jar package named LRScaf-&ltversion&gt.jar under target folder in releases. User could run it with command: <i>"java -jar LRScaf-&ltversion&gt.jar -x &ltconfigure.xml&gt"</i>. The details of configure XML file are described below.</li>
<li>If you want to compile the source code by yourself, you could download the source code and then use command: <i>"maven build"</i>  to build this project.</li>
<br>
################################################################################<br>
<b>Quick starting:</b><br>
################################################################################<br>
<i>java -jar LRScaf-&ltversion&gt.jar -x &ltconfigure.xml&gt</i><br>
<br>
################################################################################<br>
<b>Parameters of LRScaf</b><br>
################################################################################<br>
<p>LRScaf supports parameters set by XML confiuration file or command-line. It recommends to use XML configuration file. There is a template configuration file of XML format, named "scafconf.xml", in the project. And the following table would show each parameter meaning and default value if available.</p>
<table>
<tr><th>Parameter</th><th>Abbreviation</th><th>XML Code</th><th>Details</th></tr>
<tr><td>xml</td><td>x</td><td>NA</td><td>XML configure file! All command-line parameters would be omitted if this is set.</td></tr>
<tr><td>contig</td><td>c</td><td>contig</td><td>The contigs file for draft assembly in fasta format.</td></tr>
<tr><td>m5</td><td>m5</td><td>m5</td><td>The aligned file in -m 5 format of BLASR.</td></tr>
<tr><td>m4</td><td>m4</td><td>m4</td><td>The aligned file in -m 4 format of BLASR.</td></tr>
<tr><td>sam</td><td>sam</td><td>sam</td><td>The aligned file in sam format of BLASR.</td></tr>
<tr><td>mm</td><td>mm</td><td>mm</td><td>The aligned file in Minimap defualt output format.</td></tr>
<tr><td>output</td><td>o</td><td>output</td><td>The output folder.</td></tr>
<tr><td>miniCntLen</td><td>micl</td><td>min_contig_length</td><td>The minimum contigs length for scaffolding. Default: &lt500&gt.</td></tr>
<tr><td>identity</td><td>i</td><td>identity</td><td>The identity threshold for BLASR alignment! Default: &lt0.8&gt.</td></tr>
<tr><td>miniOLLen</td><td>mioll</td><td>min_overlap_length</td><td>The minimum overlap length threshold for BLASR alignment! Default: &lt400&gt.</td></tr>
<tr><td>miniOLRatio</td><td>miolr</td><td>min_overlap_ratio</td><td>The minimum overlap ratio threshold for BLASR alignment! Default: &lt0.8&gt.</td></tr>
<tr><td>maOHLen</td><td>maohl</td><td>max_overhang_length</td><td>The maximum overhang length threshold for BLASR alignment! Default: &lt300&gt.</td></tr>
<tr><td>maOHRatio</td><td>maohr</td><td>max_overhang_ratio</td><td>The maximum overhang ratio threshold for BLASR alignment! Default: &lt0.1&tg.</td></tr>
<tr><td>maELen</td><td>mael</td><td>max_end_length</td><td>The maximum ending length of PacBio® Long Read! Default: &lt300&gt.</td></tr>
<tr><td>maERatio</td><td>maer</td><td>max_end_ratio</td><td>The maximum ending ratio of PacBio® Long Read! Default: &lt0.1&gt.</td></tr>
<tr><td>miSLN</td><td>misl</td><td>min_supported_links</td><td>The minimum support links number! Default: &lt2&gt.</td></tr>
<tr><td>doll</td><td>doll</td><td>use_overlap_link</td><td>The indicator for using overlap relationship of contig! Default: &lttrue&gt, discard overlap relationship!</td></tr>
<tr><td>ratio</td><td>r</td><td>ratio</td><td>The ratio for deleting error prone edges! Default: &lt0.2&gt.</td></tr>
<tr><td>mr</td><td>mr</td><td>repeat_mask</td><td>The indicator for masking repeats! Default: &lttrue&gt.</td></tr>
<tr><td>tiplength</td><td>tl</td><td>tip_length</td><td>The maximum tip length. Default: &lt1500&gt</td></tr>
<tr><td>iqrtime</td><td>iqrt</td><td>iqr_time</td><td>The IQR time for defined repeat outlier. Default: &lt3&gt</td></tr>
<tr><td>mmcm</td><td>mmcm</td><td>mmcm</td><td>The parameter to remove Minimap aligned records. Default: &lt8&gt.</td></tr>
<tr><td>help</td><td>h</td><td>NA</td><td>The help information.</td></tr>
</table>
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
If you have any questions, please feel free to contact me &ltmqin@outlook.com&gt.
<br>
<br>
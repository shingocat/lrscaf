/*
*File: agis.ps.util2.EdgeBunlder2.java
*User: mqin
*Email: mqin@ymail.com
*Date: 2016年7月17日
*/
package agis.ps.util;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.List;
import java.util.Map;
import java.util.Vector;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import agis.ps.file.PBLinkReader;
import agis.ps.file.PBLinkWriter;
import agis.ps.file.TriadLinkWriter;
import agis.ps.link.Edge;
import agis.ps.link.M5Record;
import agis.ps.link.MRecord;
import agis.ps.link.PBLink;
import agis.ps.link.PBLinkM;
import agis.ps.seqs.Contig;

public class M5EdgeBundler {

	private static Logger logger = LoggerFactory.getLogger(M5EdgeBundler.class);
	private Parameter paras = null;
	private int minPBLen = 0;
	private int minCNTLen = 0;
	private double identity = 0.0d;
	private Map<String, Contig> contigs;
	private LinkBuilder linkBuilder = null;
	private PBLinkWriter linkWriter = null;
	private TriadLinkWriter tlWriter = null;

	public M5EdgeBundler(Parameter paras) {
		this.paras = paras;
		this.minCNTLen = paras.getMinContLen();
		this.minPBLen = paras.getMinPBLen();
		this.identity = paras.getIdentity();
		this.linkBuilder = new LinkBuilder(paras);
		this.linkWriter = new PBLinkWriter(paras);
		this.linkWriter.init();
		this.tlWriter = new TriadLinkWriter(paras);
		this.tlWriter.init();
	}

	public M5EdgeBundler(Parameter paras, Map<String, Contig> contigs) {
		this(paras);
		this.contigs = contigs;
	}

	public List<Edge> building() {
		List<Edge> edges = null;
		this.buildingLinks();
		PBLinkReader linkReader = new PBLinkReader(paras);
		List<PBLink> links = linkReader.read();
		EdgeBundler eb = new EdgeBundler(paras);
		edges = eb.pbLink2Edges(links);
		return edges;
	}

	private void buildingLinks() {
		File file = null;
		FileReader fr = null;
		BufferedReader br = null;
		String alnFile = paras.getAlgFile();
		int index = 0;
		try {
			file = new File(alnFile);
			if (!file.exists()) {
				logger.error(this.getClass().getName() + "\t" + "The m5 file do not exist!");
				return;
			}
			fr = new FileReader(file);
			br = new BufferedReader(fr);
			String line = null;
			String[] arrs = null;
			List<MRecord> records = new Vector<MRecord>(10);
			String id = null; // pacbio long read id
			MRecord m = null;
			while (true) {
				index++;
				line = br.readLine();
				if (line != null) {
					line.trim();
					line = line.replaceAll(System.getProperty("line.separator"), "");
					arrs = line.split("\\s+");
					if (arrs[0].equalsIgnoreCase("qName") && arrs[1].equalsIgnoreCase("qLength"))
						continue;
					if (id == null) {
						id = arrs[0];
						m = this.validate(arrs);
						if (m != null)
							records.add(m);
					} else { // id != null
						if (id.equals(arrs[0])) {
							m = this.validate(arrs);
							if (m != null)
								records.add(m);
						} else {
							if (records.size() >= 2) {
								List<PBLinkM> links = linkBuilder.mRecord2Link(records);
								if (links != null && links.size() > 0)
									linkWriter.write(links);
							}
							records.clear();
							id = arrs[0];
							m = this.validate(arrs);
							if (m != null)
								records.add(m);
						}
					}
				} else {
					// if line == null; ending of file
					if (records.size() >= 2) {
						List<PBLinkM> links = linkBuilder.mRecord2Link(records);
						if (links != null && links.size() > 0)
							linkWriter.write(links);
					}
					break;
				}
			}
		} catch (IOException e) {
			logger.error(this.getClass().getName() + "\t" + e.getMessage());
			logger.debug("index " + index);
		} catch (Exception e) {
			logger.error(this.getClass().getName() + "\t" + e.getMessage());
			logger.debug("index " + index);
		} finally {
			try {
				if (br != null)
					br.close();
			} catch (Exception e) {
				logger.error(this.getClass().getName() + "\t" + e.getMessage());
			}
		}
		// write triadlink to file
		tlWriter.write2(linkBuilder.getTriadLinks());
		tlWriter.close();
		linkWriter.close();
		return;
	}

	private MRecord validate(String[] arrs) {
		// if less than minimum pacbio length
		if (Integer.valueOf(arrs[1]) < minPBLen)
			return null;
		// if less tan minimum contig length
		if (Integer.valueOf(arrs[6]) < minCNTLen)
			return null;
		// if the identity less than specified value;
		double sum = Double.valueOf(arrs[11]) + Double.valueOf(arrs[12]) + Double.valueOf(arrs[13])
				+ Double.valueOf(arrs[14]);
		double value = Double.valueOf(arrs[11]) / sum;
		if (value < identity)
			return null;
		M5Record m5 = new M5Record();
		m5.setqName(arrs[0]);
		m5.setqLength(Integer.valueOf(arrs[1]));
		m5.setqStart(Integer.valueOf(arrs[2]));
		m5.setqEnd(Integer.valueOf(arrs[3]));
		m5.setqStrand(arrs[4].equals("+") ? Strand.FORWARD : Strand.REVERSE);
		m5.settName(arrs[5]);
		m5.settLength(Integer.valueOf(arrs[6]));
		m5.settStart(Integer.valueOf(arrs[7]));
		m5.settEnd(Integer.valueOf(arrs[8]));
		m5.settStrand(arrs[9].equals("+") ? Strand.FORWARD : Strand.REVERSE);
		m5.setScore(Integer.valueOf(arrs[10]));
		m5.setNumMatch(Integer.valueOf(arrs[11]));
		m5.setNumMismatch(Integer.valueOf(arrs[12]));
		m5.setNumIns(Integer.valueOf(arrs[13]));
		m5.setNumDel(Integer.valueOf(arrs[14]));
		m5.setMapQV(Integer.valueOf(arrs[15]));
		// m5.setqAlignedSeq(arrs[16]);
		// m5.setMatchPattern(arrs[17]);
		// m5.settAlignedSeq(arrs[18]);
		// setting identity
		m5.setIdentity(value);
		return m5;
	}
}

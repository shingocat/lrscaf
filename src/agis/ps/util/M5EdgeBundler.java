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
import java.util.HashMap;
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

public class M5EdgeBundler {

	private static Logger logger = LoggerFactory.getLogger(M5EdgeBundler.class);
	private Parameter paras = null;
	private int minPBLen = 0;
	private int minCNTLen = 0;
	private double identity = 0.0d;
	private LinkBuilder linkBuilder = null;
	private PBLinkWriter linkWriter = null;
	private TriadLinkWriter tlWriter = null;
	private List<Edge> edges = null;
	// for finding repeat contigs
	private Map<String, List<MRecord>> cntMs = new HashMap<String, List<MRecord>>();
	private List<String> repeats;

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

	public List<Edge> building() {
		this.findingRepeats();
		this.buildingLinks();
		PBLinkReader linkReader = new PBLinkReader(paras);
		List<PBLink> links = linkReader.read();
		logger.info(this.getClass().getName() + "\tFlitered links: " + links.size());
		EdgeBundler eb = new EdgeBundler(paras);
		edges = eb.pbLink2Edges(links);
		return edges;
	}
	
	private void findingRepeats()
	{
		RepeatFinder rf = new RepeatFinder(paras);
		repeats = rf.findRepeat2();
		rf = null;
	}

	private void buildingLinks() {
		long start = System.currentTimeMillis();
		File file = null;
		FileReader fr = null;
		BufferedReader br = null;
		String alnFile = paras.getAlgFile();
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
				line = br.readLine();
				if (line != null) {
					arrs = line.split("\\s+");
					if (arrs[0].equalsIgnoreCase("qName") && arrs[1].equalsIgnoreCase("qLength"))
						continue;
//					if(arrs[0].equalsIgnoreCase("m130605_032054_42207_c100515142550000001823076608221373_s1_p0/147441/0_4277"))
//						logger.debug("breakpoint");
					if (id == null) {
						id = arrs[0];
						m = MRecordValidator.validate(arrs, paras);
						if (m != null)
						{
							records.add(m);
						}
					} else { // id != null
						if (id.equals(arrs[0])) {
							m = MRecordValidator.validate(arrs, paras);
							if (m != null)
							{
								records.add(m);
							}
						} else {
							if (records.size() >= 2) {
								List<PBLinkM> links = linkBuilder.mRecord2Link(records, repeats);
								if (links != null && links.size() > 0)
									linkWriter.write(links);
							}
							records.clear();
							id = arrs[0];
							m = MRecordValidator.validate(arrs, paras);
							if (m != null)
							{
								records.add(m);
							}
						}
					}
				} else {
					// if line == null; ending of file
					if (records.size() >= 2) {
						List<PBLinkM> links = linkBuilder.mRecord2Link(records, repeats);
						if (links != null && links.size() > 0)
							linkWriter.write(links);
					}
					break;
				}
			}
		} catch (IOException e) {
			logger.error(this.getClass().getName() + "\t" + e.getMessage());
		} catch (Exception e) {
			logger.error(this.getClass().getName() + "\t" + e.getMessage());
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
		// Repeats finder;
		long end = System.currentTimeMillis();
		logger.info("Links building, erase time : " + (end - start) + " ms");
		return;
	}
}

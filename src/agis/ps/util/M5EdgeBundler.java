/*
*File: agis.ps.util2.EdgeBunlder2.java
*User: mqin
*Email: mqin@ymail.com
*Date: 2016年7月17日
*/
package agis.ps.util;

import java.util.List;
import java.util.Vector;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import agis.ps.file.M5Reader;
import agis.ps.file.PBLinkReader;
import agis.ps.file.PBLinkWriter;
import agis.ps.file.SimilarityCntWriter;
import agis.ps.file.TriadLinkWriter;
import agis.ps.link.CntFileEncapsulate;
import agis.ps.link.Edge;
import agis.ps.link.M5FileEncapsulate;
import agis.ps.link.M5Record;
import agis.ps.link.MRecord;
import agis.ps.link.PBLink;
import agis.ps.link.PBLinkM;

public class M5EdgeBundler {

	private static Logger logger = LoggerFactory.getLogger(M5EdgeBundler.class);
	private Parameter paras = null;
	private LinkBuilder linkBuilder = null;
	private PBLinkWriter linkWriter = null;
	private TriadLinkWriter tlWriter = null;
	private List<Edge> edges = null;
	// for finding repeat contigs
	private List<String> repeats;
	private M5FileEncapsulate m5file;
	private CntFileEncapsulate cntfile;

	public M5EdgeBundler(Parameter paras) {
		this.paras = paras;
		this.linkBuilder = new LinkBuilder(paras);
		this.linkWriter = new PBLinkWriter(paras);
//		this.linkWriter.init();
		this.tlWriter = new TriadLinkWriter(paras);
//		this.tlWriter.init();
	}

	public List<Edge> building() {
		this.readAlignFile();
		this.findingRepeats();
//		this.buildingLinks();
		// clear raw data;
		m5file = null;
		// link to edges;
		PBLinkReader linkReader = new PBLinkReader(paras);
		List<PBLink> links = linkReader.read();
		logger.info(this.getClass().getName() + "\tFlitered links: " + links.size());
		EdgeBundler eb = new EdgeBundler(paras);
		edges = eb.pbLink2Edges(links, cntfile);
		return edges;
	}
	
	private void readAlignFile()
	{
		M5Reader reader = new M5Reader(paras);
//		m5file = reader.readByArrayCopy();
//		m5file = reader.readByLineNumberReader();
		cntfile = m5file.getCntFileEncapsulate();
	}
	
	private void findingRepeats()
	{
//		RepeatFinder rf = new RepeatFinder(paras);
////		repeats = rf.findRepeat2();
//		repeats = rf.findByFileEncapsulate(m5file);
//		rf = null;
	}
	
//	
//	private void buildingLinks()
//	{
//		long start = System.currentTimeMillis();
//		try {
//
//			List<MRecord> records = new Vector<MRecord>(10);
//			String id = null; // pacbio long read id
//			MRecord m = null;
//			M5Record m5 = null;
//			int index = 0;
//			while (true) {
//				m5 = m5file.getM5Record(index);
////				if(m5.getqName().equals("10639"))
////					logger.debug("breakpoint");
//				if(m5 != null)
//				{
//					if (id == null) {
//						id = m5.getqName();
//						m = MRecordValidator.validate(m5, paras);
//						if (m != null)
//						{
//							records.add(m);
//						}
//					} else { // id != null
//						if (id.equals(m5.getqName())) {
//							m = MRecordValidator.validate(m5, paras);
//							if (m != null)
//							{
//								records.add(m);
//							}
//						} else {
//							if (records.size() >= 2) {
//								List<PBLinkM> links = linkBuilder.mRecord2Link(records, repeats);
//								if (links != null && links.size() > 0)
//									linkWriter.write(links);
//							}
//							records.clear();
//							id = m5.getqName();
//							m = MRecordValidator.validate(m5, paras);
//							if (m != null)
//							{
//								records.add(m);
//							}
//						}
//					}
//				} else {
//					if (records.size() >= 2) {
//						List<PBLinkM> links = linkBuilder.mRecord2Link(records, repeats);
//						if (links != null && links.size() > 0)
//							linkWriter.write(links);
//					}
//					break;
//				}					
//				index++;			
//			}
//		} catch (Exception e) {
//			logger.error(this.getClass().getName() + "\t" + e.getMessage());
//		} 
//		// write triadlink to file
//		tlWriter.write2(linkBuilder.getTriadLinks());
//		tlWriter.close();
//		linkWriter.close();
//		// write similarity contigs into file;
////		SimilarityCntWriter scw = new SimilarityCntWriter(paras);
////		scw.write(linkBuilder.getSimCnts());
//		// Repeats finder;
//		long end = System.currentTimeMillis();
//		logger.info("Links building, erase time : " + (end - start) + " ms");
//		return;
//	}
	// return the ContigFileEncapsulate
	public CntFileEncapsulate getCntFile()
	{
		return this.cntfile;
	}
	
}

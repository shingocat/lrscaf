/*
*File: agis.ps.util.EdgeBundler.java
*User: mqin
*Email: mqin@ymail.com
*Date: 2016年1月13日
*/
package agis.ps.util;

import java.io.File;
import java.io.IOException;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Vector;

import org.apache.lucene.analysis.Analyzer;
import org.apache.lucene.analysis.standard.StandardAnalyzer;
import org.apache.lucene.document.Document;
import org.apache.lucene.index.DirectoryReader;
import org.apache.lucene.index.IndexReader;
import org.apache.lucene.queryparser.classic.QueryParser;
import org.apache.lucene.search.IndexSearcher;
import org.apache.lucene.search.Query;
import org.apache.lucene.search.ScoreDoc;
import org.apache.lucene.search.TopDocs;
import org.apache.lucene.store.Directory;
import org.apache.lucene.store.SimpleFSDirectory;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import agis.ps.file.GapRecordWriter;
import agis.ps.link.Edge;
import agis.ps.link.PBLink;
import agis.ps.link.PBLinkM;
import agis.ps.seqs.Contig;
import agis.ps.seqs.PBGapSeq;

public class EdgeBundler {
	// private static int MIN_LINK_NUM = 3;
	static final Logger logger = LoggerFactory.getLogger(EdgeBundler.class);
	private List<PBLinkM> links;
	private Parameter paras;
	private List<Edge> edges = null;
	private Map<String, Contig> contigs = null;
	private List<GapRecord> gaps = new Vector<GapRecord>(50);
	private Directory directory = null;
	private IndexReader reader = null;
	private IndexSearcher searcher = null;
	private Analyzer analyzer = null;
	private QueryParser parser = null;
	private static int index = 0;
	
	public EdgeBundler(Parameter paras)
	{
		this.paras = paras;
	}

	public EdgeBundler(List<PBLinkM> links, Parameter paras, Map<String, Contig> contigs) {
		this.links = links;
		this.paras = paras;
		this.contigs = contigs;
	}
	
	public List<Edge> pbLink2Edges(List<PBLink> links)
	{
		long start = System.currentTimeMillis();
		if (links == null || links.size() == 0)
			throw new IllegalArgumentException(
					this.getClass().getName() + "The Links is empty when passed to EdgeBundler");
		if (edges != null)
			edges = null;
		this.initCntIndexer();
		edges = new Vector<Edge>(links.size()/5);
		boolean isUseOLLink = paras.isUseOLLink();
		int minSupLink = paras.getMinSupLinks();
		logger.info(this.getClass().getName() + "\tMinimum supported links: " + minSupLink);
		logger.info(this.getClass().getName() + "\tUsed overlap link: " + isUseOLLink);
		// the all the same origin and terminus to a hash map;
		// both direction
		Map<String, List<PBLink>> temp = new HashMap<String, List<PBLink>>();
		for (PBLink pb : links) {
			String id1 = pb.getOrigin() + ":->:" + pb.getTerminus();
			String id2 = pb.getTerminus() + ":->:" + pb.getOrigin();
			if (temp.containsKey(id1)) {
				temp.get(id1).add(pb);
				id1 = null;
				id2 = null;
				pb = null;
			} else if (temp.containsKey(id2)) {
				temp.get(id2).add(pb);
				id1 = null;
				id2 = null;
				pb = null;
			} else {
				List<PBLink> lins = new Vector<PBLink>(10);
				lins.add(pb);
				temp.put(id1, lins);
				lins = null;
				id1 = null;
				id2 = null;
				pb = null;
			}
		}

		// loop in temp hash map
		for (String s : temp.keySet()) {
			try {
				index++;
				// divide into two set, i.e. A->B; B->A
				List<PBLink> all = temp.get(s);
				List<PBLink> s1 = new Vector<PBLink>(10);
				List<PBLink> s2 = new Vector<PBLink>(10);
				String[] ids = s.split(":->:");
				for (PBLink p : all) {
					if (ids[0].equalsIgnoreCase(p.getOrigin()))
						s1.add(p);
					else
						s2.add(p);
				}

				// compute the most type in s1 and s2;
				// A for + +; B for + -; C for - -; D for - +;
				// FOR S1 SET
				List<PBLink> s1As = new Vector<PBLink>(10);
				List<PBLink> s1Bs = new Vector<PBLink>(10);
				List<PBLink> s1Cs = new Vector<PBLink>(10);
				List<PBLink> s1Ds = new Vector<PBLink>(10);
				for (PBLink p : s1) {
					Strand oStrand = p.getOStrand();
					Strand tStrand = p.getTStrand();
					if (oStrand.equals(Strand.FORWARD)) {
						if (tStrand.equals(Strand.FORWARD)) {
							s1As.add(p);
						} else {
							s1Bs.add(p);
						}
					} else {
						if (tStrand.equals(Strand.REVERSE)) {
							s1Cs.add(p);
						} else {
							s1Ds.add(p);
						}
					}
					oStrand = null;
					tStrand = null;
				}
				// FOR S2 SET
				List<PBLink> s2As = new Vector<PBLink>(10);
				List<PBLink> s2Bs = new Vector<PBLink>(10);
				List<PBLink> s2Cs = new Vector<PBLink>(10);
				List<PBLink> s2Ds = new Vector<PBLink>(10);
				for (PBLink p : s2) {
					Strand oStrand = p.getOStrand();
					Strand tStrand = p.getTStrand();
					if (oStrand.equals(Strand.FORWARD)) {
						if (tStrand.equals(Strand.FORWARD)) {
							s2As.add(p);
						} else {
							s2Bs.add(p);
						}
					} else {
						if (tStrand.equals(Strand.REVERSE)) {
							s2Cs.add(p);
						} else {
							s2Ds.add(p);
						}
					}
					oStrand = null;
					tStrand = null;
				}

				int sA = s1As.size() + s2Cs.size(); // A(+)->B(+):B(-)->A(-)
				int sB = s1Bs.size() + s2Bs.size(); // A(+)->B(-):B(+)->A(-)
				int sC = s1Cs.size() + s2As.size(); // A(-)->B(-):B(+)->A(+)
				int sD = s1Ds.size() + s2Ds.size(); // A(-)->B(+):B(-)->A(+)

				int max = MathTool.max(sA, sB, sC, sD);
				if (max == sA) // for the A statement;
				{
					List<Integer> dists = new Vector<Integer>(10);
					for (PBLink p : s1As)
						dists.add(p.getDistance());
					for (PBLink p : s2Cs)
						dists.add(p.getDistance());
					int mean = MathTool.mean(dists);
					int sd = MathTool.sd(dists);
					int upper = mean + 3 * sd;
					int low = mean - 3 * sd;
					// if the distance larger than mean + 3 * sd, then remove;
					Iterator<PBLink> it1 = s1As.iterator();
					while (it1.hasNext()) {
						PBLink pb = it1.next();
						if (pb.getDistance() > upper || pb.getDistance() < low)
							it1.remove();
					}
					Iterator<PBLink> it2 = s2Cs.iterator();
					while (it2.hasNext()) {
						PBLink pb = it2.next();
						if (pb.getDistance() > upper || pb.getDistance() < low)
							it2.remove();
					}
					// if the supported link count less then user specified
					int supLink = s1As.size() + s2Cs.size();
					// logger.debug(this.getClass().getName() + "\tSupported
					// link: "
					// + supLink);
					if (supLink < minSupLink)
						continue;
					// recompute mean and sd;
					dists = null;
					dists = new Vector<Integer>();
					for (PBLink p : s1As)
						dists.add(p.getDistance());
					for (PBLink p : s2Cs)
						dists.add(p.getDistance());
					mean = MathTool.mean(dists);
					sd = MathTool.sd(dists);

					// if the mean is less tan zero and the user specified do
					// not
					// use overlap link;
					if (!isUseOLLink && mean < 0)
						continue;
					// there are very similarity of two contigs with high
					// identity
					// and PB reads also have high error ratio, so it will lead
					// to
					// pb read links two high identity contigs by overlap two
					// large
					// defining a standard when they over larger or less than
					// this
					// filtering parameter, it will omit this edges;
					if (mean < 0) {
						low = Math.abs(mean) - 4 * sd;
						upper = Math.abs(mean) + 4 * sd;
//						int c1Len = contigs.get(ids[0]).getLength();
//						int c2Len = contigs.get(ids[1]).getLength();
						int c1Len = this.indexCntLength(ids[0]);
						int c2Len = this.indexCntLength(ids[1]);
						if (c1Len >= low && c1Len <= upper)
							continue;
						if (c2Len >= low && c2Len <= upper)
							continue;
					}
					// for store gap record information
					if (mean > 0) {
						GapRecord gr = new GapRecord();
						gr.setStart(ids[0]);
						gr.setEnd(ids[1]);
						for (PBLink p : s1As) {
							if (p.getDistance() > 0) {
								PBGapSeq seq = new PBGapSeq();
//								seq.setId(p.getOrigin().getqName());
//								int p1 = p.getOrigin().getqEnd();
//								int p2 = p.getTerminus().getqStart();
								seq.setId(p.getPbId());
								int p1 = p.getOPEnd();
								int p2 = p.getTPStart();
								if (p1 <= p2) {
									seq.setStart(p1);
									seq.setEnd(p2);
								} else {
									seq.setStart(p2);
									seq.setEnd(p1);
								}
								seq.setStrand(Strand.FORWARD);
								gr.addSeq(seq);
							}
						}
						for (PBLink p : s2Cs) {
							if (p.getDistance() > 0) {
								PBGapSeq seq = new PBGapSeq();
//								seq.setId(p.getOrigin().getqName());
//								int p1 = p.getOrigin().getqEnd();
//								int p2 = p.getTerminus().getqStart();
								seq.setId(p.getPbId());
								int p1 = p.getOPEnd();
								int p2 = p.getTPStart();
								if (p1 <= p2) {
									seq.setStart(p1);
									seq.setEnd(p2);
								} else {
									seq.setStart(p2);
									seq.setEnd(p1);
								}
								seq.setStrand(Strand.REVERSE);
								gr.addSeq(seq);
							}
						}
						gaps.add(gr);
					}
					// edge + +;
					Edge edge = new Edge();
//					Contig origin = contigs.get(ids[0]);
//					Contig terminus = contigs.get(ids[1]);
					Contig origin = new Contig();
					Contig terminus = new Contig();
					origin.setID(ids[0]);
//					origin.setLength(this.indexCntLength(ids[0]));
					terminus.setID(ids[1]);
//					terminus.setLength(this.indexCntLength(ids[1]));
					edge.setOrigin(origin);
					edge.setTerminus(terminus);
					edge.setDistMean(mean);
					edge.setDistSd(sd);
					edge.setFake(false);
					// if the mean is minus, it means that edge is overlap, set
					// it
					// to be overlap and not valid;
					if (mean < 0)
						edge.setOL(true);
					if (isUseOLLink)
						edge.setValid(true);
					else
						edge.setValid(false);
					edge.setoStrand(Strand.FORWARD);
					edge.settStrand(Strand.FORWARD);
					edge.setLinkNum(dists.size());
					edges.add(edge);
					// edge - -;
					edge = null;
					// origin = null;
					// terminus = null;
					edge = new Edge();
					Contig tC = origin;
					origin = terminus;
					terminus = tC;
					tC = null;
					// origin = new Contig();
					// terminus = new Contig();
					// origin.setID(ids[1]);
					// origin.setLength(contigs.get(ids[1]).getLength());
					// terminus.setID(ids[0]);
					// terminus.setLength(contigs.get(ids[0]).getLength());
					edge.setOrigin(origin);
					edge.setTerminus(terminus);
					edge.setDistMean(mean);
					edge.setDistSd(sd);
					edge.setFake(false);
					// if the mean is minus, it means that edge is overlap, set
					// it
					// to be overlap and not valid;
					if (mean < 0)
						edge.setOL(true);
					if (isUseOLLink)
						edge.setValid(true);
					else
						edge.setValid(false);
					edge.setoStrand(Strand.REVERSE);
					edge.settStrand(Strand.REVERSE);
					edge.setLinkNum(dists.size());
					edges.add(edge);
					edge = null;
					origin = null;
					terminus = null;
				} else if (max == sB) // for the B statement
				{
					List<Integer> dists = new Vector<Integer>(10);
					for (PBLink p : s1Bs)
						dists.add(p.getDistance());
					for (PBLink p : s2Bs)
						dists.add(p.getDistance());
					int mean = MathTool.mean(dists);
					int sd = MathTool.sd(dists);
					int upper = mean + 3 * sd;
					int low = mean - 3 * sd;
					// if the distance larger than mean + 3 * sd, then remove;
					Iterator<PBLink> it1 = s1Bs.iterator();
					while (it1.hasNext()) {
						PBLink pb = it1.next();
						if (pb.getDistance() > upper || pb.getDistance() < low)
							it1.remove();
					}
					Iterator<PBLink> it2 = s2Bs.iterator();
					while (it2.hasNext()) {
						PBLink pb = it2.next();
						if (pb.getDistance() > upper || pb.getDistance() < low)
							it2.remove();
					}
					// if the supported link count less then user specified
					int supLink = s1Bs.size() + s2Bs.size();
					// logger.debug(this.getClass().getName() + "\tSupported
					// link: "
					// + supLink);
					if (supLink < minSupLink)
						continue;
					// recompute mean and sd;
					dists = null;
					dists = new Vector<Integer>();
					for (PBLink p : s1Bs)
						dists.add(p.getDistance());
					for (PBLink p : s2Bs)
						dists.add(p.getDistance());
					mean = MathTool.mean(dists);
					sd = MathTool.sd(dists);
					// if the mean is less tan zero and the user specified do
					// not
					// use overlap link;
					if (!isUseOLLink && mean < 0)
						continue;
					// there are very similarity of two contigs with high
					// identity
					// and PB reads also have high error ratio, so it will lead
					// to
					// pb read links two high identity contigs by overlap two
					// large
					// defining a standard when they over larger or less than
					// this
					// filtering parameter, it will omit this edges;
					if (mean < 0) {
						low = Math.abs(mean) - 4 * sd;
						upper = Math.abs(mean) + 4 * sd;
//						int c1Len = contigs.get(ids[0]).getLength();
//						int c2Len = contigs.get(ids[1]).getLength();
						int c1Len = this.indexCntLength(ids[0]);
						int c2Len = this.indexCntLength(ids[1]);
						if (c1Len >= low && c1Len <= upper)
							continue;
						if (c2Len >= low && c2Len <= upper)
							continue;
					}
					// for store gap record information
					if (mean > 0) {
						GapRecord gr = new GapRecord();
						gr.setStart(ids[0]);
						gr.setEnd(ids[1]);
						for (PBLink p : s1Bs) {
							if (p.getDistance() > 0) {
								PBGapSeq seq = new PBGapSeq();
//								seq.setId(p.getOrigin().getqName());
//								int p1 = p.getOrigin().getqEnd();
//								int p2 = p.getTerminus().getqStart();
								seq.setId(p.getPbId());
								int p1 = p.getOPEnd();
								int p2 = p.getTPStart();
								if (p1 <= p2) {
									seq.setStart(p1);
									seq.setEnd(p2);
								} else {
									seq.setStart(p2);
									seq.setEnd(p1);
								}
								seq.setStrand(Strand.FORWARD);
								gr.addSeq(seq);
							}
						}
						for (PBLink p : s2Bs) {
							if (p.getDistance() > 0) {
								PBGapSeq seq = new PBGapSeq();
//								seq.setId(p.getOrigin().getqName());
//								int p1 = p.getOrigin().getqEnd();
//								int p2 = p.getTerminus().getqStart();
								seq.setId(p.getPbId());
								int p1 = p.getOPEnd();
								int p2 = p.getTPStart();
								if (p1 <= p2) {
									seq.setStart(p1);
									seq.setEnd(p2);
								} else {
									seq.setStart(p2);
									seq.setEnd(p1);
								}
								seq.setStrand(Strand.REVERSE);
								gr.addSeq(seq);
							}
						}
						gaps.add(gr);
					}
					// edge + -;
					Edge edge = new Edge();
//					Contig origin = contigs.get(ids[0]);
//					Contig terminus = contigs.get(ids[1]);
					 Contig origin = new Contig();
					 Contig terminus = new Contig();
					 origin.setID(ids[0]);
					// origin.setLength(contigs.get(ids[0]).getLength());
					 terminus.setID(ids[1]);
					// terminus.setLength(contigs.get(ids[1]).getLength());
					edge.setOrigin(origin);
					edge.setTerminus(terminus);
					edge.setDistMean(mean);
					edge.setDistSd(sd);
					edge.setFake(false);
					// if the mean is minus, it means that edge is overlap, set
					// it
					// to be overlap and not valid;
					if (mean < 0)
						edge.setOL(true);
					if (isUseOLLink)
						edge.setValid(true);
					else
						edge.setValid(false);
					edge.setoStrand(Strand.FORWARD);
					edge.settStrand(Strand.REVERSE);
					edge.setLinkNum(dists.size());
					edges.add(edge);
					// edge - +;
					edge = null;
					// origin = null;
					// terminus = null;
					edge = new Edge();
					Contig tC = origin;
					origin = terminus;
					terminus = tC;
					tC = null;
					// origin = new Contig();
					// terminus = new Contig();
					// origin.setID(ids[1]);
					// origin.setLength(contigs.get(ids[1]).getLength());
					// terminus.setID(ids[0]);
					// terminus.setLength(contigs.get(ids[0]).getLength());
					edge.setOrigin(origin);
					edge.setTerminus(terminus);
					edge.setDistMean(mean);
					edge.setDistSd(sd);
					edge.setFake(false);
					// if the mean is minus, it means that edge is overlap, set
					// it
					// to be overlap and not valid;
					if (mean < 0)
						edge.setOL(true);
					if (isUseOLLink)
						edge.setValid(true);
					else
						edge.setValid(false);
					edge.setoStrand(Strand.FORWARD);
					edge.settStrand(Strand.REVERSE);
					edge.setLinkNum(dists.size());
					edges.add(edge);
					edge = null;
					origin = null;
					terminus = null;
				} else if (max == sC) // for the C statement
				{
					List<Integer> dists = new Vector<Integer>(10);
					for (PBLink p : s1Cs)
						dists.add(p.getDistance());
					for (PBLink p : s2As)
						dists.add(p.getDistance());
					int mean = MathTool.mean(dists);
					int sd = MathTool.sd(dists);
					int upper = mean + 3 * sd;
					int low = mean - 3 * sd;
					// if the distance larger than mean + 3 * sd, then remove;
					Iterator<PBLink> it1 = s1Cs.iterator();
					while (it1.hasNext()) {
						PBLink pb = it1.next();
						if (pb.getDistance() > upper || pb.getDistance() < low)
							it1.remove();
					}
					Iterator<PBLink> it2 = s2As.iterator();
					while (it2.hasNext()) {
						PBLink pb = it2.next();
						if (pb.getDistance() > upper || pb.getDistance() < low)
							it2.remove();
					}
					// if the supported link count less then user specified
					int supLink = s1Cs.size() + s2As.size();
					// logger.debug(this.getClass().getName() + "\tSupported
					// link: "
					// + supLink);
					if (supLink < minSupLink)
						continue;
					// recompute mean and sd;
					dists = null;
					dists = new Vector<Integer>();
					for (PBLink p : s1Cs)
						dists.add(p.getDistance());
					for (PBLink p : s2As)
						dists.add(p.getDistance());
					mean = MathTool.mean(dists);
					sd = MathTool.sd(dists);
					// if the mean is less tan zero and the user specified do
					// not
					// use overlap link;
					if (!isUseOLLink && mean < 0)
						continue;
					// there are very similarity of two contigs with high
					// identity
					// and PB reads also have high error ratio, so it will lead
					// to
					// pb read links two high identity contigs by overlap two
					// large
					// defining a standard when they over larger or less than
					// this
					// filtering parameter, it will omit this edges;
					if (mean < 0) {
						low = Math.abs(mean) - 4 * sd;
						upper = Math.abs(mean) + 4 * sd;
//						int c1Len = contigs.get(ids[0]).getLength();
//						int c2Len = contigs.get(ids[1]).getLength();
						int c1Len = this.indexCntLength(ids[0]);
						int c2Len = this.indexCntLength(ids[1]);
						if (c1Len >= low && c1Len <= upper)
							continue;
						if (c2Len >= low && c2Len <= upper)
							continue;
					}
					// for store gap record information
					if (mean > 0) {
						GapRecord gr = new GapRecord();
						gr.setStart(ids[0]);
						gr.setEnd(ids[1]);
						for (PBLink p : s1Cs) {
							if (p.getDistance() > 0) {
								PBGapSeq seq = new PBGapSeq();
//								seq.setId(p.getOrigin().getqName());
//								int p1 = p.getOrigin().getqEnd();
//								int p2 = p.getTerminus().getqStart();
								seq.setId(p.getPbId());
								int p1 = p.getOPEnd();
								int p2 = p.getTPStart();
								if (p1 <= p2) {
									seq.setStart(p1);
									seq.setEnd(p2);
								} else {
									seq.setStart(p2);
									seq.setEnd(p1);
								}
								seq.setStrand(Strand.FORWARD);
								gr.addSeq(seq);
							}
						}
						for (PBLink p : s2As) {
							if (p.getDistance() > 0) {
								PBGapSeq seq = new PBGapSeq();
//								seq.setId(p.getOrigin().getqName());
//								int p1 = p.getOrigin().getqEnd();
//								int p2 = p.getTerminus().getqStart();
								seq.setId(p.getPbId());
								int p1 = p.getOPEnd();
								int p2 = p.getTPStart();
								if (p1 <= p2) {
									seq.setStart(p1);
									seq.setEnd(p2);
								} else {
									seq.setStart(p2);
									seq.setEnd(p1);
								}
								seq.setStrand(Strand.REVERSE);
								gr.addSeq(seq);
							}
						}
						gaps.add(gr);
					}

					// edge - -;
					Edge edge = new Edge();
//					Contig origin = contigs.get(ids[0]);
//					Contig terminus = contigs.get(ids[1]);
					 Contig origin = new Contig();
					 Contig terminus = new Contig();
					 origin.setID(ids[0]);
					// origin.setLength(contigs.get(ids[0]).getLength());
					 terminus.setID(ids[1]);
					// terminus.setLength(contigs.get(ids[1]).getLength());
					edge.setOrigin(origin);
					edge.setTerminus(terminus);
					edge.setDistMean(mean);
					edge.setDistSd(sd);
					edge.setFake(false);
					// if the mean is minus, it means that edge is overlap, set
					// it
					// to be overlap and not valid;
					if (mean < 0)
						edge.setOL(true);
					if (isUseOLLink)
						edge.setValid(true);
					else
						edge.setValid(false);
					edge.setoStrand(Strand.REVERSE);
					edge.settStrand(Strand.REVERSE);
					edge.setLinkNum(dists.size());
					edges.add(edge);
					// edge + +;
					edge = null;
					// origin = null;
					// terminus = null;
					edge = new Edge();
					Contig tC = origin;
					origin = terminus;
					terminus = tC;
					tC = null;
					// origin = new Contig();
					// terminus = new Contig();
					// origin.setID(ids[1]);
					// origin.setLength(contigs.get(ids[1]).getLength());
					// terminus.setID(ids[0]);
					// terminus.setLength(contigs.get(ids[0]).getLength());
					edge.setOrigin(origin);
					edge.setTerminus(terminus);
					edge.setDistMean(mean);
					edge.setDistSd(sd);
					edge.setFake(false);
					// if the mean is minus, it means that edge is overlap, set
					// it
					// to be overlap and not valid;
					if (mean < 0)
						edge.setOL(true);
					if (isUseOLLink)
						edge.setValid(true);
					else
						edge.setValid(false);
					edge.setoStrand(Strand.FORWARD);
					edge.settStrand(Strand.FORWARD);
					edge.setLinkNum(dists.size());
					edges.add(edge);
					edge = null;
					origin = null;
					terminus = null;
				} else if (max == sD) // for the D statement
				{
					List<Integer> dists = new Vector<Integer>(10);
					for (PBLink p : s1Ds)
						dists.add(p.getDistance());
					for (PBLink p : s2Ds)
						dists.add(p.getDistance());
					int mean = MathTool.mean(dists);
					int sd = MathTool.sd(dists);
					int upper = mean + 3 * sd;
					int low = mean - 3 * sd;
					// if the distance larger than mean + 3 * sd, then remove;
					Iterator<PBLink> it1 = s1Ds.iterator();
					while (it1.hasNext()) {
						PBLink pb = it1.next();
						if (pb.getDistance() > upper || pb.getDistance() < low)
							it1.remove();
					}
					Iterator<PBLink> it2 = s2Ds.iterator();
					while (it2.hasNext()) {
						PBLink pb = it2.next();
						if (pb.getDistance() > upper || pb.getDistance() < low)
							it2.remove();
					}
					// if the supported link count less then user specified
					int supLink = s1Ds.size() + s2Ds.size();
					// logger.debug(this.getClass().getName() + "\tSupported
					// link: "
					// + supLink);
					if (supLink < minSupLink)
						continue;
					// recompute mean and sd;
					dists = null;
					dists = new Vector<Integer>();
					for (PBLink p : s1Ds)
						dists.add(p.getDistance());
					for (PBLink p : s2Ds)
						dists.add(p.getDistance());
					mean = MathTool.mean(dists);
					sd = MathTool.sd(dists);
					// if the mean is less tan zero and the user specified do
					// not
					// use overlap link;
					if (!isUseOLLink && mean < 0)
						continue;
					// there are very similarity of two contigs with high
					// identity
					// and PB reads also have high error ratio, so it will lead
					// to
					// pb read links two high identity contigs by overlap two
					// large
					// defining a standard when they over larger or less than
					// this
					// filtering parameter, it will omit this edges;
					if (mean < 0) {
						low = Math.abs(mean) - 4 * sd;
						upper = Math.abs(mean) + 4 * sd;
//						int c1Len = contigs.get(ids[0]).getLength();
//						int c2Len = contigs.get(ids[1]).getLength();
						int c1Len = this.indexCntLength(ids[0]);
						int c2Len = this.indexCntLength(ids[1]);
						if (c1Len >= low && c1Len <= upper)
							continue;
						if (c2Len >= low && c2Len <= upper)
							continue;
					}
					// for store gap record information
					if (mean > 0) {
						GapRecord gr = new GapRecord();
						gr.setStart(ids[0]);
						gr.setEnd(ids[1]);
						for (PBLink p : s1Ds) {
							if (p.getDistance() > 0) {
								PBGapSeq seq = new PBGapSeq();
//								seq.setId(p.getOrigin().getqName());
//								int p1 = p.getOrigin().getqEnd();
//								int p2 = p.getTerminus().getqStart();
								seq.setId(p.getPbId());
								int p1 = p.getOPEnd();
								int p2 = p.getTPStart();
								if (p1 <= p2) {
									seq.setStart(p1);
									seq.setEnd(p2);
								} else {
									seq.setStart(p2);
									seq.setEnd(p1);
								}
								seq.setStrand(Strand.FORWARD);
								gr.addSeq(seq);
							}
						}
						for (PBLink p : s2Ds) {
							if (p.getDistance() > 0) {
								PBGapSeq seq = new PBGapSeq();
//								seq.setId(p.getOrigin().getqName());
//								int p1 = p.getOrigin().getqEnd();
//								int p2 = p.getTerminus().getqStart();
								seq.setId(p.getPbId());
								int p1 = p.getOPEnd();
								int p2 = p.getTPStart();
								if (p1 <= p2) {
									seq.setStart(p1);
									seq.setEnd(p2);
								} else {
									seq.setStart(p2);
									seq.setEnd(p1);
								}
								seq.setStrand(Strand.REVERSE);
								gr.addSeq(seq);
							}
						}
						gaps.add(gr);
					}
					// edge - +;
					Edge edge = new Edge();
//					Contig origin = contigs.get(ids[0]);
//					Contig terminus = contigs.get(ids[1]);
					 Contig origin = new Contig();
					 Contig terminus = new Contig();
					 origin.setID(ids[0]);
					// origin.setLength(contigs.get(ids[0]).getLength());
					 terminus.setID(ids[1]);
					// terminus.setLength(contigs.get(ids[1]).getLength());
					edge.setOrigin(origin);
					edge.setTerminus(terminus);
					edge.setDistMean(mean);
					edge.setDistSd(sd);
					edge.setFake(false);
					// if the mean is minus, it means that edge is overlap, set
					// it
					// to be overlap and not valid;
					if (mean < 0)
						edge.setOL(true);
					if (isUseOLLink)
						edge.setValid(true);
					else
						edge.setValid(false);
					edge.setoStrand(Strand.REVERSE);
					edge.settStrand(Strand.FORWARD);
					edge.setLinkNum(dists.size());
					edges.add(edge);
					// edge + -;
					edge = null;
					// origin = null;
					// terminus = null;
					edge = new Edge();
					Contig tC = origin;
					origin = terminus;
					terminus = tC;
					tC = null;
					// origin = new Contig();
					// terminus = new Contig();
					// origin.setID(ids[1]);
					// origin.setLength(contigs.get(ids[1]).getLength());
					// terminus.setID(ids[0]);
					// terminus.setLength(contigs.get(ids[0]).getLength());
					edge.setOrigin(origin);
					edge.setTerminus(terminus);
					edge.setDistMean(mean);
					edge.setDistSd(sd);
					edge.setFake(false);
					// if the mean is minus, it means that edge is overlap, set
					// it
					// to be overlap and not valid;
					if (mean < 0)
						edge.setOL(true);
					if (isUseOLLink)
						edge.setValid(true);
					else
						edge.setValid(false);
					edge.setoStrand(Strand.REVERSE);
					edge.settStrand(Strand.FORWARD);
					edge.setLinkNum(dists.size());
					edges.add(edge);
					edge = null;
					origin = null;
					terminus = null;
				}

				// gc the variables
				all = null;
				s1 = null;
				s2 = null;
				s1As = null;
				s1Bs = null;
				s1Cs = null;
				s1Ds = null;
				s2As = null;
				s2Bs = null;
				s2Cs = null;
				s2Ds = null;
			} catch (Exception e) {
				logger.error(this.getClass().getName() + "\t" + e.getMessage());
			}
		}

		GapRecordWriter grw = new GapRecordWriter(paras, gaps);
		grw.write();
		gaps = null;
		long end = System.currentTimeMillis();
		logger.info("Edge Bundling, erase time: " + (end - start) + " ms");
		return edges;
	}

	public List<Edge> pbLinkM5Bundling() {
		if (links == null || links.size() == 0)
			throw new IllegalArgumentException(this.getClass().getName() + "The PBLinkM5 could not be empty!");
		return this.pbLinkM5Bundling(links, paras);
	}

	// the first implemented
	private List<Edge> pbLinkM5Bundling(List<PBLinkM> links, Parameter paras) {
		long start = System.currentTimeMillis();
		if (links == null || links.size() == 0)
			throw new IllegalArgumentException(
					this.getClass().getName() + "The Links is empty when passed to EdgeBundler");
		boolean isUseOLLink = paras.isUseOLLink();
		if (edges != null)
			edges = null;
		edges = new Vector<Edge>(200);
		int minSupLink = paras.getMinSupLinks();
		logger.info(this.getClass().getName() + "\tMinimum supported links:" + minSupLink);
		logger.info(this.getClass().getName() + "\tUsed overlap link:" + isUseOLLink);
		// the all the same origin and terminus to a hash map;
		// both direction
		Map<String, List<PBLinkM>> temp = new HashMap<String, List<PBLinkM>>();
		for (PBLinkM pb : links) {
			String id1 = pb.getOrigin().gettName() + ":->:" + pb.getTerminus().gettName();
			String id2 = pb.getTerminus().gettName() + ":->:" + pb.getOrigin().gettName();
			if (temp.containsKey(id1)) {
				temp.get(id1).add(pb);
				id1 = null;
				id2 = null;
				pb = null;
			} else if (temp.containsKey(id2)) {
				temp.get(id2).add(pb);
				id1 = null;
				id2 = null;
				pb = null;
			} else {
				List<PBLinkM> lins = new Vector<PBLinkM>(10);
				lins.add(pb);
				temp.put(id1, lins);
				lins = null;
				id1 = null;
				id2 = null;
				pb = null;
			}
		}

		// loop in temp hash map

		for (String s : temp.keySet()) {
			try {
				index++;
				// divide into two set, i.e. A->B; B->A
				List<PBLinkM> all = temp.get(s);
				List<PBLinkM> s1 = new Vector<PBLinkM>(10);
				List<PBLinkM> s2 = new Vector<PBLinkM>(10);
				String[] ids = s.split(":->:");
				for (PBLinkM p : all) {
					if (ids[0].equalsIgnoreCase(p.getOrigin().gettName()))
						s1.add(p);
					else
						s2.add(p);
				}

				// compute the most type in s1 and s2;
				// A for + +; B for + -; C for - -; D for - +;
				// FOR S1 SET
				List<PBLinkM> s1As = new Vector<PBLinkM>(10);
				List<PBLinkM> s1Bs = new Vector<PBLinkM>(10);
				List<PBLinkM> s1Cs = new Vector<PBLinkM>(10);
				List<PBLinkM> s1Ds = new Vector<PBLinkM>(10);
				for (PBLinkM p : s1) {
					Strand oStrand = p.getOrigin().gettStrand();
					Strand tStrand = p.getTerminus().gettStrand();
					if (oStrand.equals(Strand.FORWARD)) {
						if (tStrand.equals(Strand.FORWARD)) {
							s1As.add(p);
						} else {
							s1Bs.add(p);
						}
					} else {
						if (tStrand.equals(Strand.REVERSE)) {
							s1Cs.add(p);
						} else {
							s1Ds.add(p);
						}
					}
					oStrand = null;
					tStrand = null;
				}
				// FOR S2 SET
				List<PBLinkM> s2As = new Vector<PBLinkM>(10);
				List<PBLinkM> s2Bs = new Vector<PBLinkM>(10);
				List<PBLinkM> s2Cs = new Vector<PBLinkM>(10);
				List<PBLinkM> s2Ds = new Vector<PBLinkM>(10);
				for (PBLinkM p : s2) {
					Strand oStrand = p.getOrigin().gettStrand();
					Strand tStrand = p.getTerminus().gettStrand();
					if (oStrand.equals(Strand.FORWARD)) {
						if (tStrand.equals(Strand.FORWARD)) {
							s2As.add(p);
						} else {
							s2Bs.add(p);
						}
					} else {
						if (tStrand.equals(Strand.REVERSE)) {
							s2Cs.add(p);
						} else {
							s2Ds.add(p);
						}
					}
					oStrand = null;
					tStrand = null;
				}

				int sA = s1As.size() + s2Cs.size(); // A(+)->B(+):B(-)->A(-)
				int sB = s1Bs.size() + s2Bs.size(); // A(+)->B(-):B(+)->A(-)
				int sC = s1Cs.size() + s2As.size(); // A(-)->B(-):B(+)->A(+)
				int sD = s1Ds.size() + s2Ds.size(); // A(-)->B(+):B(-)->A(+)

				int max = MathTool.max(sA, sB, sC, sD);
				if (max == sA) // for the A statement;
				{
					List<Integer> dists = new Vector<Integer>(10);
					for (PBLinkM p : s1As)
						dists.add(p.getDistance());
					for (PBLinkM p : s2Cs)
						dists.add(p.getDistance());
					int mean = MathTool.mean(dists);
					int sd = MathTool.sd(dists);
					int upper = mean + 3 * sd;
					int low = mean - 3 * sd;
					// if the distance larger than mean + 3 * sd, then remove;
					Iterator<PBLinkM> it1 = s1As.iterator();
					while (it1.hasNext()) {
						PBLinkM pb = it1.next();
						if (pb.getDistance() > upper || pb.getDistance() < low)
							it1.remove();
					}
					Iterator<PBLinkM> it2 = s2Cs.iterator();
					while (it2.hasNext()) {
						PBLinkM pb = it2.next();
						if (pb.getDistance() > upper || pb.getDistance() < low)
							it2.remove();
					}
					// if the supported link count less then user specified
					int supLink = s1As.size() + s2Cs.size();
					// logger.debug(this.getClass().getName() + "\tSupported
					// link: "
					// + supLink);
					if (supLink < minSupLink)
						continue;
					// recompute mean and sd;
					dists = null;
					dists = new Vector<Integer>();
					for (PBLinkM p : s1As)
						dists.add(p.getDistance());
					for (PBLinkM p : s2Cs)
						dists.add(p.getDistance());
					mean = MathTool.mean(dists);
					sd = MathTool.sd(dists);

					// if the mean is less tan zero and the user specified do
					// not
					// use overlap link;
					if (!isUseOLLink && mean < 0)
						continue;
					// there are very similarity of two contigs with high
					// identity
					// and PB reads also have high error ratio, so it will lead
					// to
					// pb read links two high identity contigs by overlap two
					// large
					// defining a standard when they over larger or less than
					// this
					// filtering parameter, it will omit this edges;
					if (mean < 0) {
						low = Math.abs(mean) - 4 * sd;
						upper = Math.abs(mean) + 4 * sd;
						int c1Len = contigs.get(ids[0]).getLength();
						int c2Len = contigs.get(ids[1]).getLength();
						if (c1Len >= low && c1Len <= upper)
							continue;
						if (c2Len >= low && c2Len <= upper)
							continue;
					}
					// for store gap record information
					if (mean > 0) {
						GapRecord gr = new GapRecord();
						gr.setStart(ids[0]);
						gr.setEnd(ids[1]);
						for (PBLinkM p : s1As) {
							if (p.getDistance() > 0) {
								PBGapSeq seq = new PBGapSeq();
								seq.setId(p.getOrigin().getqName());
								int p1 = p.getOrigin().getqEnd();
								int p2 = p.getTerminus().getqStart();
								if (p1 <= p2) {
									seq.setStart(p1);
									seq.setEnd(p2);
								} else {
									seq.setStart(p2);
									seq.setEnd(p1);
								}
								seq.setStrand(Strand.FORWARD);
								gr.addSeq(seq);
							}
						}
						for (PBLinkM p : s2Cs) {
							if (p.getDistance() > 0) {
								PBGapSeq seq = new PBGapSeq();
								seq.setId(p.getOrigin().getqName());
								int p1 = p.getOrigin().getqEnd();
								int p2 = p.getTerminus().getqStart();
								if (p1 <= p2) {
									seq.setStart(p1);
									seq.setEnd(p2);
								} else {
									seq.setStart(p2);
									seq.setEnd(p1);
								}
								seq.setStrand(Strand.REVERSE);
								gr.addSeq(seq);
							}
						}
						gaps.add(gr);
					}
					// edge + +;
					Edge edge = new Edge();
					Contig origin = contigs.get(ids[0]);
					Contig terminus = contigs.get(ids[1]);
					// Contig origin = new Contig();
					// Contig terminus = new Contig();
					// origin.setID(ids[0]);
					// origin.setLength(contigs.get(ids[0]).getLength());
					// terminus.setID(ids[1]);
					// terminus.setLength(contigs.get(ids[1]).getLength());
					edge.setOrigin(origin);
					edge.setTerminus(terminus);
					edge.setDistMean(mean);
					edge.setDistSd(sd);
					edge.setFake(false);
					// if the mean is minus, it means that edge is overlap, set
					// it
					// to be overlap and not valid;
					if (mean < 0)
						edge.setOL(true);
					if (isUseOLLink)
						edge.setValid(true);
					else
						edge.setValid(false);
					edge.setoStrand(Strand.FORWARD);
					edge.settStrand(Strand.FORWARD);
					edge.setLinkNum(dists.size());
					edges.add(edge);
					// edge - -;
					edge = null;
					// origin = null;
					// terminus = null;
					edge = new Edge();
					Contig tC = origin;
					origin = terminus;
					terminus = tC;
					tC = null;
					// origin = new Contig();
					// terminus = new Contig();
					// origin.setID(ids[1]);
					// origin.setLength(contigs.get(ids[1]).getLength());
					// terminus.setID(ids[0]);
					// terminus.setLength(contigs.get(ids[0]).getLength());
					edge.setOrigin(origin);
					edge.setTerminus(terminus);
					edge.setDistMean(mean);
					edge.setDistSd(sd);
					edge.setFake(false);
					// if the mean is minus, it means that edge is overlap, set
					// it
					// to be overlap and not valid;
					if (mean < 0)
						edge.setOL(true);
					if (isUseOLLink)
						edge.setValid(true);
					else
						edge.setValid(false);
					edge.setoStrand(Strand.REVERSE);
					edge.settStrand(Strand.REVERSE);
					edge.setLinkNum(dists.size());
					edges.add(edge);
					edge = null;
					origin = null;
					terminus = null;
				} else if (max == sB) // for the B statement
				{
					List<Integer> dists = new Vector<Integer>(10);
					for (PBLinkM p : s1Bs)
						dists.add(p.getDistance());
					for (PBLinkM p : s2Bs)
						dists.add(p.getDistance());
					int mean = MathTool.mean(dists);
					int sd = MathTool.sd(dists);
					int upper = mean + 3 * sd;
					int low = mean - 3 * sd;
					// if the distance larger than mean + 3 * sd, then remove;
					Iterator<PBLinkM> it1 = s1Bs.iterator();
					while (it1.hasNext()) {
						PBLinkM pb = it1.next();
						if (pb.getDistance() > upper || pb.getDistance() < low)
							it1.remove();
					}
					Iterator<PBLinkM> it2 = s2Bs.iterator();
					while (it2.hasNext()) {
						PBLinkM pb = it2.next();
						if (pb.getDistance() > upper || pb.getDistance() < low)
							it2.remove();
					}
					// if the supported link count less then user specified
					int supLink = s1Bs.size() + s2Bs.size();
					// logger.debug(this.getClass().getName() + "\tSupported
					// link: "
					// + supLink);
					if (supLink < minSupLink)
						continue;
					// recompute mean and sd;
					dists = null;
					dists = new Vector<Integer>();
					for (PBLinkM p : s1Bs)
						dists.add(p.getDistance());
					for (PBLinkM p : s2Bs)
						dists.add(p.getDistance());
					mean = MathTool.mean(dists);
					sd = MathTool.sd(dists);
					// if the mean is less tan zero and the user specified do
					// not
					// use overlap link;
					if (!isUseOLLink && mean < 0)
						continue;
					// there are very similarity of two contigs with high
					// identity
					// and PB reads also have high error ratio, so it will lead
					// to
					// pb read links two high identity contigs by overlap two
					// large
					// defining a standard when they over larger or less than
					// this
					// filtering parameter, it will omit this edges;
					if (mean < 0) {
						low = Math.abs(mean) - 4 * sd;
						upper = Math.abs(mean) + 4 * sd;
						int c1Len = contigs.get(ids[0]).getLength();
						int c2Len = contigs.get(ids[1]).getLength();
						if (c1Len >= low && c1Len <= upper)
							continue;
						if (c2Len >= low && c2Len <= upper)
							continue;
					}
					// for store gap record information
					if (mean > 0) {
						GapRecord gr = new GapRecord();
						gr.setStart(ids[0]);
						gr.setEnd(ids[1]);
						for (PBLinkM p : s1Bs) {
							if (p.getDistance() > 0) {
								PBGapSeq seq = new PBGapSeq();
								seq.setId(p.getOrigin().getqName());
								int p1 = p.getOrigin().getqEnd();
								int p2 = p.getTerminus().getqStart();
								if (p1 <= p2) {
									seq.setStart(p1);
									seq.setEnd(p2);
								} else {
									seq.setStart(p2);
									seq.setEnd(p1);
								}
								seq.setStrand(Strand.FORWARD);
								gr.addSeq(seq);
							}
						}
						for (PBLinkM p : s2Bs) {
							if (p.getDistance() > 0) {
								PBGapSeq seq = new PBGapSeq();
								seq.setId(p.getOrigin().getqName());
								int p1 = p.getOrigin().getqEnd();
								int p2 = p.getTerminus().getqStart();
								if (p1 <= p2) {
									seq.setStart(p1);
									seq.setEnd(p2);
								} else {
									seq.setStart(p2);
									seq.setEnd(p1);
								}
								seq.setStrand(Strand.REVERSE);
								gr.addSeq(seq);
							}
						}
						gaps.add(gr);
					}
					// edge + -;
					Edge edge = new Edge();
					Contig origin = contigs.get(ids[0]);
					Contig terminus = contigs.get(ids[1]);
					// Contig origin = new Contig();
					// Contig terminus = new Contig();
					// origin.setID(ids[0]);
					// origin.setLength(contigs.get(ids[0]).getLength());
					// terminus.setID(ids[1]);
					// terminus.setLength(contigs.get(ids[1]).getLength());
					edge.setOrigin(origin);
					edge.setTerminus(terminus);
					edge.setDistMean(mean);
					edge.setDistSd(sd);
					edge.setFake(false);
					// if the mean is minus, it means that edge is overlap, set
					// it
					// to be overlap and not valid;
					if (mean < 0)
						edge.setOL(true);
					if (isUseOLLink)
						edge.setValid(true);
					else
						edge.setValid(false);
					edge.setoStrand(Strand.FORWARD);
					edge.settStrand(Strand.REVERSE);
					edge.setLinkNum(dists.size());
					edges.add(edge);
					// edge - +;
					edge = null;
					// origin = null;
					// terminus = null;
					edge = new Edge();
					Contig tC = origin;
					origin = terminus;
					terminus = tC;
					tC = null;
					// origin = new Contig();
					// terminus = new Contig();
					// origin.setID(ids[1]);
					// origin.setLength(contigs.get(ids[1]).getLength());
					// terminus.setID(ids[0]);
					// terminus.setLength(contigs.get(ids[0]).getLength());
					edge.setOrigin(origin);
					edge.setTerminus(terminus);
					edge.setDistMean(mean);
					edge.setDistSd(sd);
					edge.setFake(false);
					// if the mean is minus, it means that edge is overlap, set
					// it
					// to be overlap and not valid;
					if (mean < 0)
						edge.setOL(true);
					if (isUseOLLink)
						edge.setValid(true);
					else
						edge.setValid(false);
					edge.setoStrand(Strand.FORWARD);
					edge.settStrand(Strand.REVERSE);
					edge.setLinkNum(dists.size());
					edges.add(edge);
					edge = null;
					origin = null;
					terminus = null;
				} else if (max == sC) // for the C statement
				{
					List<Integer> dists = new Vector<Integer>(10);
					for (PBLinkM p : s1Cs)
						dists.add(p.getDistance());
					for (PBLinkM p : s2As)
						dists.add(p.getDistance());
					int mean = MathTool.mean(dists);
					int sd = MathTool.sd(dists);
					int upper = mean + 3 * sd;
					int low = mean - 3 * sd;
					// if the distance larger than mean + 3 * sd, then remove;
					Iterator<PBLinkM> it1 = s1Cs.iterator();
					while (it1.hasNext()) {
						PBLinkM pb = it1.next();
						if (pb.getDistance() > upper || pb.getDistance() < low)
							it1.remove();
					}
					Iterator<PBLinkM> it2 = s2As.iterator();
					while (it2.hasNext()) {
						PBLinkM pb = it2.next();
						if (pb.getDistance() > upper || pb.getDistance() < low)
							it2.remove();
					}
					// if the supported link count less then user specified
					int supLink = s1Cs.size() + s2As.size();
					// logger.debug(this.getClass().getName() + "\tSupported
					// link: "
					// + supLink);
					if (supLink < minSupLink)
						continue;
					// recompute mean and sd;
					dists = null;
					dists = new Vector<Integer>();
					for (PBLinkM p : s1Cs)
						dists.add(p.getDistance());
					for (PBLinkM p : s2As)
						dists.add(p.getDistance());
					mean = MathTool.mean(dists);
					sd = MathTool.sd(dists);
					// if the mean is less tan zero and the user specified do
					// not
					// use overlap link;
					if (!isUseOLLink && mean < 0)
						continue;
					// there are very similarity of two contigs with high
					// identity
					// and PB reads also have high error ratio, so it will lead
					// to
					// pb read links two high identity contigs by overlap two
					// large
					// defining a standard when they over larger or less than
					// this
					// filtering parameter, it will omit this edges;
					if (mean < 0) {
						low = Math.abs(mean) - 4 * sd;
						upper = Math.abs(mean) + 4 * sd;
						int c1Len = contigs.get(ids[0]).getLength();
						int c2Len = contigs.get(ids[1]).getLength();
						if (c1Len >= low && c1Len <= upper)
							continue;
						if (c2Len >= low && c2Len <= upper)
							continue;
					}
					// for store gap record information
					if (mean > 0) {
						GapRecord gr = new GapRecord();
						gr.setStart(ids[0]);
						gr.setEnd(ids[1]);
						for (PBLinkM p : s1Cs) {
							if (p.getDistance() > 0) {
								PBGapSeq seq = new PBGapSeq();
								seq.setId(p.getOrigin().getqName());
								int p1 = p.getOrigin().getqEnd();
								int p2 = p.getTerminus().getqStart();
								if (p1 <= p2) {
									seq.setStart(p1);
									seq.setEnd(p2);
								} else {
									seq.setStart(p2);
									seq.setEnd(p1);
								}
								seq.setStrand(Strand.FORWARD);
								gr.addSeq(seq);
							}
						}
						for (PBLinkM p : s2As) {
							if (p.getDistance() > 0) {
								PBGapSeq seq = new PBGapSeq();
								seq.setId(p.getOrigin().getqName());
								int p1 = p.getOrigin().getqEnd();
								int p2 = p.getTerminus().getqStart();
								if (p1 <= p2) {
									seq.setStart(p1);
									seq.setEnd(p2);
								} else {
									seq.setStart(p2);
									seq.setEnd(p1);
								}
								seq.setStrand(Strand.REVERSE);
								gr.addSeq(seq);
							}
						}
						gaps.add(gr);
					}

					// edge - -;
					Edge edge = new Edge();
					Contig origin = contigs.get(ids[0]);
					Contig terminus = contigs.get(ids[1]);
					// Contig origin = new Contig();
					// Contig terminus = new Contig();
					// origin.setID(ids[0]);
					// origin.setLength(contigs.get(ids[0]).getLength());
					// terminus.setID(ids[1]);
					// terminus.setLength(contigs.get(ids[1]).getLength());
					edge.setOrigin(origin);
					edge.setTerminus(terminus);
					edge.setDistMean(mean);
					edge.setDistSd(sd);
					edge.setFake(false);
					// if the mean is minus, it means that edge is overlap, set
					// it
					// to be overlap and not valid;
					if (mean < 0)
						edge.setOL(true);
					if (isUseOLLink)
						edge.setValid(true);
					else
						edge.setValid(false);
					edge.setoStrand(Strand.REVERSE);
					edge.settStrand(Strand.REVERSE);
					edge.setLinkNum(dists.size());
					edges.add(edge);
					// edge + +;
					edge = null;
					// origin = null;
					// terminus = null;
					edge = new Edge();
					Contig tC = origin;
					origin = terminus;
					terminus = tC;
					tC = null;
					// origin = new Contig();
					// terminus = new Contig();
					// origin.setID(ids[1]);
					// origin.setLength(contigs.get(ids[1]).getLength());
					// terminus.setID(ids[0]);
					// terminus.setLength(contigs.get(ids[0]).getLength());
					edge.setOrigin(origin);
					edge.setTerminus(terminus);
					edge.setDistMean(mean);
					edge.setDistSd(sd);
					edge.setFake(false);
					// if the mean is minus, it means that edge is overlap, set
					// it
					// to be overlap and not valid;
					if (mean < 0)
						edge.setOL(true);
					if (isUseOLLink)
						edge.setValid(true);
					else
						edge.setValid(false);
					edge.setoStrand(Strand.FORWARD);
					edge.settStrand(Strand.FORWARD);
					edge.setLinkNum(dists.size());
					edges.add(edge);
					edge = null;
					origin = null;
					terminus = null;
				} else if (max == sD) // for the D statement
				{
					List<Integer> dists = new Vector<Integer>(10);
					for (PBLinkM p : s1Ds)
						dists.add(p.getDistance());
					for (PBLinkM p : s2Ds)
						dists.add(p.getDistance());
					int mean = MathTool.mean(dists);
					int sd = MathTool.sd(dists);
					int upper = mean + 3 * sd;
					int low = mean - 3 * sd;
					// if the distance larger than mean + 3 * sd, then remove;
					Iterator<PBLinkM> it1 = s1Ds.iterator();
					while (it1.hasNext()) {
						PBLinkM pb = it1.next();
						if (pb.getDistance() > upper || pb.getDistance() < low)
							it1.remove();
					}
					Iterator<PBLinkM> it2 = s2Ds.iterator();
					while (it2.hasNext()) {
						PBLinkM pb = it2.next();
						if (pb.getDistance() > upper || pb.getDistance() < low)
							it2.remove();
					}
					// if the supported link count less then user specified
					int supLink = s1Ds.size() + s2Ds.size();
					// logger.debug(this.getClass().getName() + "\tSupported
					// link: "
					// + supLink);
					if (supLink < minSupLink)
						continue;
					// recompute mean and sd;
					dists = null;
					dists = new Vector<Integer>();
					for (PBLinkM p : s1Ds)
						dists.add(p.getDistance());
					for (PBLinkM p : s2Ds)
						dists.add(p.getDistance());
					mean = MathTool.mean(dists);
					sd = MathTool.sd(dists);
					// if the mean is less tan zero and the user specified do
					// not
					// use overlap link;
					if (!isUseOLLink && mean < 0)
						continue;
					// there are very similarity of two contigs with high
					// identity
					// and PB reads also have high error ratio, so it will lead
					// to
					// pb read links two high identity contigs by overlap two
					// large
					// defining a standard when they over larger or less than
					// this
					// filtering parameter, it will omit this edges;
					if (mean < 0) {
						low = Math.abs(mean) - 4 * sd;
						upper = Math.abs(mean) + 4 * sd;
						int c1Len = contigs.get(ids[0]).getLength();
						int c2Len = contigs.get(ids[1]).getLength();
						if (c1Len >= low && c1Len <= upper)
							continue;
						if (c2Len >= low && c2Len <= upper)
							continue;
					}
					// for store gap record information
					if (mean > 0) {
						GapRecord gr = new GapRecord();
						gr.setStart(ids[0]);
						gr.setEnd(ids[1]);
						for (PBLinkM p : s1Ds) {
							if (p.getDistance() > 0) {
								PBGapSeq seq = new PBGapSeq();
								seq.setId(p.getOrigin().getqName());
								int p1 = p.getOrigin().getqEnd();
								int p2 = p.getTerminus().getqStart();
								if (p1 <= p2) {
									seq.setStart(p1);
									seq.setEnd(p2);
								} else {
									seq.setStart(p2);
									seq.setEnd(p1);
								}
								seq.setStrand(Strand.FORWARD);
								gr.addSeq(seq);
							}
						}
						for (PBLinkM p : s2Ds) {
							if (p.getDistance() > 0) {
								PBGapSeq seq = new PBGapSeq();
								seq.setId(p.getOrigin().getqName());
								int p1 = p.getOrigin().getqEnd();
								int p2 = p.getTerminus().getqStart();
								if (p1 <= p2) {
									seq.setStart(p1);
									seq.setEnd(p2);
								} else {
									seq.setStart(p2);
									seq.setEnd(p1);
								}
								seq.setStrand(Strand.REVERSE);
								gr.addSeq(seq);
							}
						}
						gaps.add(gr);
					}
					// edge - +;
					Edge edge = new Edge();
					Contig origin = contigs.get(ids[0]);
					Contig terminus = contigs.get(ids[1]);
					// Contig origin = new Contig();
					// Contig terminus = new Contig();
					// origin.setID(ids[0]);
					// origin.setLength(contigs.get(ids[0]).getLength());
					// terminus.setID(ids[1]);
					// terminus.setLength(contigs.get(ids[1]).getLength());
					edge.setOrigin(origin);
					edge.setTerminus(terminus);
					edge.setDistMean(mean);
					edge.setDistSd(sd);
					edge.setFake(false);
					// if the mean is minus, it means that edge is overlap, set
					// it
					// to be overlap and not valid;
					if (mean < 0)
						edge.setOL(true);
					if (isUseOLLink)
						edge.setValid(true);
					else
						edge.setValid(false);
					edge.setoStrand(Strand.REVERSE);
					edge.settStrand(Strand.FORWARD);
					edge.setLinkNum(dists.size());
					edges.add(edge);
					// edge + -;
					edge = null;
					// origin = null;
					// terminus = null;
					edge = new Edge();
					Contig tC = origin;
					origin = terminus;
					terminus = tC;
					tC = null;
					// origin = new Contig();
					// terminus = new Contig();
					// origin.setID(ids[1]);
					// origin.setLength(contigs.get(ids[1]).getLength());
					// terminus.setID(ids[0]);
					// terminus.setLength(contigs.get(ids[0]).getLength());
					edge.setOrigin(origin);
					edge.setTerminus(terminus);
					edge.setDistMean(mean);
					edge.setDistSd(sd);
					edge.setFake(false);
					// if the mean is minus, it means that edge is overlap, set
					// it
					// to be overlap and not valid;
					if (mean < 0)
						edge.setOL(true);
					if (isUseOLLink)
						edge.setValid(true);
					else
						edge.setValid(false);
					edge.setoStrand(Strand.REVERSE);
					edge.settStrand(Strand.FORWARD);
					edge.setLinkNum(dists.size());
					edges.add(edge);
					edge = null;
					origin = null;
					terminus = null;
				}

				// gc the variables
				all = null;
				s1 = null;
				s2 = null;
				s1As = null;
				s1Bs = null;
				s1Cs = null;
				s1Ds = null;
				s2As = null;
				s2Bs = null;
				s2Cs = null;
				s2Ds = null;
			} catch (Exception e) {
				logger.debug(this.getClass().getName() + "\tindex = " + index);
				logger.debug(this.getClass().getName() + "\t" + e.getMessage());
				logger.error(this.getClass().getName() + "\t" + e.getMessage());
			}
		}

		GapRecordWriter grw = new GapRecordWriter(paras, gaps);
		grw.write();
		gaps = null;
		long end = System.currentTimeMillis();
		logger.info("Edge Bundling, erase time: " + (end - start) + " ms");
		return edges;
	}

	// the second implemented
	private List<Edge> pbLinkM5Bundling2(List<PBLinkM> links, Parameter paras) {
		if (links == null || links.size() == 0)
			throw new IllegalArgumentException(
					this.getClass().getName() + "The Links is empty when passed to EdgeBundler");
		boolean isUseOLLink = paras.isUseOLLink();
		if (edges != null)
			edges = null;
		edges = new Vector<Edge>(6);
		int minSupLink = paras.getMinSupLinks();
		logger.info(this.getClass().getName() + "\tMinimum supported links:" + minSupLink);
		logger.info(this.getClass().getName() + "\tUsed overlap link:" + isUseOLLink);

		// the all the same origin and terminus to a hash map;
		// both direction
		Map<String, List<PBLinkM>> temp = new HashMap<String, List<PBLinkM>>();
		for (PBLinkM pb : links) {
			String id1 = pb.getOrigin().gettName() + ":->:" + pb.getTerminus().gettName();
			String id2 = pb.getTerminus().gettName() + ":->:" + pb.getOrigin().gettName();
			if (temp.containsKey(id1)) {
				temp.get(id1).add(pb);
				id1 = null;
				id2 = null;
				pb = null;
			} else if (temp.containsKey(id2)) {
				temp.get(id2).add(pb);
				id1 = null;
				id2 = null;
				pb = null;
			} else {
				List<PBLinkM> lins = new Vector<PBLinkM>(10);
				lins.add(pb);
				temp.put(id1, lins);
				lins = null;
				id1 = null;
				id2 = null;
				pb = null;
			}
		}

		// loop in temp hash map

		for (String s : temp.keySet()) {
			try {
				index++;
				// divide into two set, i.e. A->B; B->A
				List<PBLinkM> all = temp.get(s);
				List<PBLinkM> s1 = new Vector<PBLinkM>(10);
				List<PBLinkM> s2 = new Vector<PBLinkM>(10);
				String[] ids = s.split(":->:");
				for (PBLinkM p : all) {
					if (ids[0].equalsIgnoreCase(p.getOrigin().gettName()))
						s1.add(p);
					else
						s2.add(p);
				}

				// compute the most type in s1 and s2;
				// A for + +; B for + -; C for - -; D for - +;
				// FOR S1 SET
				List<PBLinkM> s1As = new Vector<PBLinkM>(10);
				List<PBLinkM> s1Bs = new Vector<PBLinkM>(10);
				List<PBLinkM> s1Cs = new Vector<PBLinkM>(10);
				List<PBLinkM> s1Ds = new Vector<PBLinkM>(10);
				for (PBLinkM p : s1) {
					Strand oStrand = p.getOrigin().gettStrand();
					Strand tStrand = p.getTerminus().gettStrand();
					if (oStrand.equals(Strand.FORWARD)) {
						if (tStrand.equals(Strand.FORWARD)) {
							s1As.add(p);
						} else {
							s1Bs.add(p);
						}
					} else {
						if (tStrand.equals(Strand.REVERSE)) {
							s1Cs.add(p);
						} else {
							s1Ds.add(p);
						}
					}
					oStrand = null;
					tStrand = null;
				}
				// FOR S2 SET
				List<PBLinkM> s2As = new Vector<PBLinkM>(10);
				List<PBLinkM> s2Bs = new Vector<PBLinkM>(10);
				List<PBLinkM> s2Cs = new Vector<PBLinkM>(10);
				List<PBLinkM> s2Ds = new Vector<PBLinkM>(10);
				for (PBLinkM p : s2) {
					Strand oStrand = p.getOrigin().gettStrand();
					Strand tStrand = p.getTerminus().gettStrand();
					if (oStrand.equals(Strand.FORWARD)) {
						if (tStrand.equals(Strand.FORWARD)) {
							s2As.add(p);
						} else {
							s2Bs.add(p);
						}
					} else {
						if (tStrand.equals(Strand.REVERSE)) {
							s2Cs.add(p);
						} else {
							s2Ds.add(p);
						}
					}
					oStrand = null;
					tStrand = null;
				}

				int sA = s1As.size() + s2Cs.size(); // A(+)->B(+):B(-)->A(-)
				int sB = s1Bs.size() + s2Bs.size(); // A(+)->B(-):B(+)->A(-)
				int sC = s1Cs.size() + s2As.size(); // A(-)->B(-):B(+)->A(+)
				int sD = s1Ds.size() + s2Ds.size(); // A(-)->B(+):B(-)->A(+)

				int max = MathTool.max(sA, sB, sC, sD);
				if (max == sA) // for the A statement;
				{
					List<Integer> dists = new Vector<Integer>(10);
					for (PBLinkM p : s1As)
						dists.add(p.getDistance());
					for (PBLinkM p : s2Cs)
						dists.add(p.getDistance());
					int mean = MathTool.mean(dists);
					int sd = MathTool.sd(dists);
					int upper = mean + 3 * sd;
					int low = mean - 3 * sd;
					// if the distance larger than mean + 3 * sd, then remove;
					Iterator<PBLinkM> it1 = s1As.iterator();
					while (it1.hasNext()) {
						PBLinkM pb = it1.next();
						if (pb.getDistance() > upper || pb.getDistance() < low)
							it1.remove();
					}
					Iterator<PBLinkM> it2 = s2Cs.iterator();
					while (it2.hasNext()) {
						PBLinkM pb = it2.next();
						if (pb.getDistance() > upper || pb.getDistance() < low)
							it2.remove();
					}
					// if the supported link count less then user specified
					int supLink = s1As.size() + s2Cs.size();
					// logger.debug(this.getClass().getName() + "\tSupported
					// link: "
					// + supLink);
					if (supLink < minSupLink)
						continue;
					// recompute mean and sd;
					dists = null;
					dists = new Vector<Integer>();
					for (PBLinkM p : s1As)
						dists.add(p.getDistance());
					for (PBLinkM p : s2Cs)
						dists.add(p.getDistance());
					mean = MathTool.mean(dists);
					sd = MathTool.sd(dists);

					// if the mean is less tan zero and the user specified do
					// not
					// use overlap link;
					if (!isUseOLLink && mean < 0)
						continue;
					// there are very similarity of two contigs with high
					// identity
					// and PB reads also have high error ratio, so it will lead
					// to
					// pb read links two high identity contigs by overlap two
					// large
					// defining a standard when they over larger or less than
					// this
					// filtering parameter, it will omit this edges;
					if (mean < 0) {
						low = Math.abs(mean) - 4 * sd;
						upper = Math.abs(mean) + 4 * sd;
						int c1Len = contigs.get(ids[0]).getLength();
						int c2Len = contigs.get(ids[1]).getLength();
						if (c1Len >= low && c1Len <= upper)
							continue;
						if (c2Len >= low && c2Len <= upper)
							continue;
					}
					// for store gap record information
					if (mean > 0) {
						GapRecord gr = new GapRecord();
						gr.setStart(ids[0]);
						gr.setEnd(ids[1]);
						for (PBLinkM p : s1As) {
							if (p.getDistance() > 0) {
								PBGapSeq seq = new PBGapSeq();
								seq.setId(p.getOrigin().getqName());
								int p1 = p.getOrigin().getqEnd();
								int p2 = p.getTerminus().getqStart();
								if (p1 <= p2) {
									seq.setStart(p1);
									seq.setEnd(p2);
								} else {
									seq.setStart(p2);
									seq.setEnd(p1);
								}
								seq.setStrand(Strand.FORWARD);
								gr.addSeq(seq);
							}
						}
						for (PBLinkM p : s2Cs) {
							if (p.getDistance() > 0) {
								PBGapSeq seq = new PBGapSeq();
								seq.setId(p.getOrigin().getqName());
								int p1 = p.getOrigin().getqEnd();
								int p2 = p.getTerminus().getqStart();
								if (p1 <= p2) {
									seq.setStart(p1);
									seq.setEnd(p2);
								} else {
									seq.setStart(p2);
									seq.setEnd(p1);
								}
								seq.setStrand(Strand.REVERSE);
								gr.addSeq(seq);
							}
						}
						gaps.add(gr);
					}
					// edge + +;
					Edge edge = new Edge();
					Contig origin = contigs.get(ids[0]);
					Contig terminus = contigs.get(ids[1]);
					// Contig origin = new Contig();
					// Contig terminus = new Contig();
					// origin.setID(ids[0]);
					// origin.setLength(contigs.get(ids[0]).getLength());
					// terminus.setID(ids[1]);
					// terminus.setLength(contigs.get(ids[1]).getLength());
					edge.setOrigin(origin);
					edge.setTerminus(terminus);
					edge.setDistMean(mean);
					edge.setDistSd(sd);
					edge.setFake(false);
					// if the mean is minus, it means that edge is overlap, set
					// it
					// to be overlap and not valid;
					if (mean < 0)
						edge.setOL(true);
					if (isUseOLLink)
						edge.setValid(true);
					else
						edge.setValid(false);
					edge.setoStrand(Strand.FORWARD);
					edge.settStrand(Strand.FORWARD);
					edge.setLinkNum(dists.size());
					edges.add(edge);
					// edge - -;
					edge = null;
					// origin = null;
					// terminus = null;
					edge = new Edge();
					Contig tC = origin;
					origin = terminus;
					terminus = tC;
					tC = null;
					// origin = new Contig();
					// terminus = new Contig();
					// origin.setID(ids[1]);
					// origin.setLength(contigs.get(ids[1]).getLength());
					// terminus.setID(ids[0]);
					// terminus.setLength(contigs.get(ids[0]).getLength());
					edge.setOrigin(origin);
					edge.setTerminus(terminus);
					edge.setDistMean(mean);
					edge.setDistSd(sd);
					edge.setFake(false);
					// if the mean is minus, it means that edge is overlap, set
					// it
					// to be overlap and not valid;
					if (mean < 0)
						edge.setOL(true);
					if (isUseOLLink)
						edge.setValid(true);
					else
						edge.setValid(false);
					edge.setoStrand(Strand.REVERSE);
					edge.settStrand(Strand.REVERSE);
					edge.setLinkNum(dists.size());
					edges.add(edge);
					edge = null;
					origin = null;
					terminus = null;
				} else if (max == sB) // for the B statement
				{
					List<Integer> dists = new Vector<Integer>(10);
					for (PBLinkM p : s1Bs)
						dists.add(p.getDistance());
					for (PBLinkM p : s2Bs)
						dists.add(p.getDistance());
					int mean = MathTool.mean(dists);
					int sd = MathTool.sd(dists);
					int upper = mean + 3 * sd;
					int low = mean - 3 * sd;
					// if the distance larger than mean + 3 * sd, then remove;
					Iterator<PBLinkM> it1 = s1Bs.iterator();
					while (it1.hasNext()) {
						PBLinkM pb = it1.next();
						if (pb.getDistance() > upper || pb.getDistance() < low)
							it1.remove();
					}
					Iterator<PBLinkM> it2 = s2Bs.iterator();
					while (it2.hasNext()) {
						PBLinkM pb = it2.next();
						if (pb.getDistance() > upper || pb.getDistance() < low)
							it2.remove();
					}
					// if the supported link count less then user specified
					int supLink = s1Bs.size() + s2Bs.size();
					// logger.debug(this.getClass().getName() + "\tSupported
					// link: "
					// + supLink);
					if (supLink < minSupLink)
						continue;
					// recompute mean and sd;
					dists = null;
					dists = new Vector<Integer>();
					for (PBLinkM p : s1Bs)
						dists.add(p.getDistance());
					for (PBLinkM p : s2Bs)
						dists.add(p.getDistance());
					mean = MathTool.mean(dists);
					sd = MathTool.sd(dists);
					// if the mean is less tan zero and the user specified do
					// not
					// use overlap link;
					if (!isUseOLLink && mean < 0)
						continue;
					// there are very similarity of two contigs with high
					// identity
					// and PB reads also have high error ratio, so it will lead
					// to
					// pb read links two high identity contigs by overlap two
					// large
					// defining a standard when they over larger or less than
					// this
					// filtering parameter, it will omit this edges;
					if (mean < 0) {
						low = Math.abs(mean) - 4 * sd;
						upper = Math.abs(mean) + 4 * sd;
						int c1Len = contigs.get(ids[0]).getLength();
						int c2Len = contigs.get(ids[1]).getLength();
						if (c1Len >= low && c1Len <= upper)
							continue;
						if (c2Len >= low && c2Len <= upper)
							continue;
					}
					// for store gap record information
					if (mean > 0) {
						GapRecord gr = new GapRecord();
						gr.setStart(ids[0]);
						gr.setEnd(ids[1]);
						for (PBLinkM p : s1Bs) {
							if (p.getDistance() > 0) {
								PBGapSeq seq = new PBGapSeq();
								seq.setId(p.getOrigin().getqName());
								int p1 = p.getOrigin().getqEnd();
								int p2 = p.getTerminus().getqStart();
								if (p1 <= p2) {
									seq.setStart(p1);
									seq.setEnd(p2);
								} else {
									seq.setStart(p2);
									seq.setEnd(p1);
								}
								seq.setStrand(Strand.FORWARD);
								gr.addSeq(seq);
							}
						}
						for (PBLinkM p : s2Bs) {
							if (p.getDistance() > 0) {
								PBGapSeq seq = new PBGapSeq();
								seq.setId(p.getOrigin().getqName());
								int p1 = p.getOrigin().getqEnd();
								int p2 = p.getTerminus().getqStart();
								if (p1 <= p2) {
									seq.setStart(p1);
									seq.setEnd(p2);
								} else {
									seq.setStart(p2);
									seq.setEnd(p1);
								}
								seq.setStrand(Strand.REVERSE);
								gr.addSeq(seq);
							}
						}
						gaps.add(gr);
					}
					// edge + -;
					Edge edge = new Edge();
					Contig origin = contigs.get(ids[0]);
					Contig terminus = contigs.get(ids[1]);
					// Contig origin = new Contig();
					// Contig terminus = new Contig();
					// origin.setID(ids[0]);
					// origin.setLength(contigs.get(ids[0]).getLength());
					// terminus.setID(ids[1]);
					// terminus.setLength(contigs.get(ids[1]).getLength());
					edge.setOrigin(origin);
					edge.setTerminus(terminus);
					edge.setDistMean(mean);
					edge.setDistSd(sd);
					edge.setFake(false);
					// if the mean is minus, it means that edge is overlap, set
					// it
					// to be overlap and not valid;
					if (mean < 0)
						edge.setOL(true);
					if (isUseOLLink)
						edge.setValid(true);
					else
						edge.setValid(false);
					edge.setoStrand(Strand.FORWARD);
					edge.settStrand(Strand.REVERSE);
					edge.setLinkNum(dists.size());
					edges.add(edge);
					// edge - +;
					edge = null;
					// origin = null;
					// terminus = null;
					edge = new Edge();
					Contig tC = origin;
					origin = terminus;
					terminus = tC;
					tC = null;
					// origin = new Contig();
					// terminus = new Contig();
					// origin.setID(ids[1]);
					// origin.setLength(contigs.get(ids[1]).getLength());
					// terminus.setID(ids[0]);
					// terminus.setLength(contigs.get(ids[0]).getLength());
					edge.setOrigin(origin);
					edge.setTerminus(terminus);
					edge.setDistMean(mean);
					edge.setDistSd(sd);
					edge.setFake(false);
					// if the mean is minus, it means that edge is overlap, set
					// it
					// to be overlap and not valid;
					if (mean < 0)
						edge.setOL(true);
					if (isUseOLLink)
						edge.setValid(true);
					else
						edge.setValid(false);
					edge.setoStrand(Strand.FORWARD);
					edge.settStrand(Strand.REVERSE);
					edge.setLinkNum(dists.size());
					edges.add(edge);
					edge = null;
					origin = null;
					terminus = null;
				} else if (max == sC) // for the C statement
				{
					List<Integer> dists = new Vector<Integer>(10);
					for (PBLinkM p : s1Cs)
						dists.add(p.getDistance());
					for (PBLinkM p : s2As)
						dists.add(p.getDistance());
					int mean = MathTool.mean(dists);
					int sd = MathTool.sd(dists);
					int upper = mean + 3 * sd;
					int low = mean - 3 * sd;
					// if the distance larger than mean + 3 * sd, then remove;
					Iterator<PBLinkM> it1 = s1Cs.iterator();
					while (it1.hasNext()) {
						PBLinkM pb = it1.next();
						if (pb.getDistance() > upper || pb.getDistance() < low)
							it1.remove();
					}
					Iterator<PBLinkM> it2 = s2As.iterator();
					while (it2.hasNext()) {
						PBLinkM pb = it2.next();
						if (pb.getDistance() > upper || pb.getDistance() < low)
							it2.remove();
					}
					// if the supported link count less then user specified
					int supLink = s1Cs.size() + s2As.size();
					// logger.debug(this.getClass().getName() + "\tSupported
					// link: "
					// + supLink);
					if (supLink < minSupLink)
						continue;
					// recompute mean and sd;
					dists = null;
					dists = new Vector<Integer>();
					for (PBLinkM p : s1Cs)
						dists.add(p.getDistance());
					for (PBLinkM p : s2As)
						dists.add(p.getDistance());
					mean = MathTool.mean(dists);
					sd = MathTool.sd(dists);
					// if the mean is less tan zero and the user specified do
					// not
					// use overlap link;
					if (!isUseOLLink && mean < 0)
						continue;
					// there are very similarity of two contigs with high
					// identity
					// and PB reads also have high error ratio, so it will lead
					// to
					// pb read links two high identity contigs by overlap two
					// large
					// defining a standard when they over larger or less than
					// this
					// filtering parameter, it will omit this edges;
					if (mean < 0) {
						low = Math.abs(mean) - 4 * sd;
						upper = Math.abs(mean) + 4 * sd;
						int c1Len = contigs.get(ids[0]).getLength();
						int c2Len = contigs.get(ids[1]).getLength();
						if (c1Len >= low && c1Len <= upper)
							continue;
						if (c2Len >= low && c2Len <= upper)
							continue;
					}
					// for store gap record information
					if (mean > 0) {
						GapRecord gr = new GapRecord();
						gr.setStart(ids[0]);
						gr.setEnd(ids[1]);
						for (PBLinkM p : s1Cs) {
							if (p.getDistance() > 0) {
								PBGapSeq seq = new PBGapSeq();
								seq.setId(p.getOrigin().getqName());
								int p1 = p.getOrigin().getqEnd();
								int p2 = p.getTerminus().getqStart();
								if (p1 <= p2) {
									seq.setStart(p1);
									seq.setEnd(p2);
								} else {
									seq.setStart(p2);
									seq.setEnd(p1);
								}
								seq.setStrand(Strand.FORWARD);
								gr.addSeq(seq);
							}
						}
						for (PBLinkM p : s2As) {
							if (p.getDistance() > 0) {
								PBGapSeq seq = new PBGapSeq();
								seq.setId(p.getOrigin().getqName());
								int p1 = p.getOrigin().getqEnd();
								int p2 = p.getTerminus().getqStart();
								if (p1 <= p2) {
									seq.setStart(p1);
									seq.setEnd(p2);
								} else {
									seq.setStart(p2);
									seq.setEnd(p1);
								}
								seq.setStrand(Strand.REVERSE);
								gr.addSeq(seq);
							}
						}
						gaps.add(gr);
					}

					// edge - -;
					Edge edge = new Edge();
					Contig origin = contigs.get(ids[0]);
					Contig terminus = contigs.get(ids[1]);
					// Contig origin = new Contig();
					// Contig terminus = new Contig();
					// origin.setID(ids[0]);
					// origin.setLength(contigs.get(ids[0]).getLength());
					// terminus.setID(ids[1]);
					// terminus.setLength(contigs.get(ids[1]).getLength());
					edge.setOrigin(origin);
					edge.setTerminus(terminus);
					edge.setDistMean(mean);
					edge.setDistSd(sd);
					edge.setFake(false);
					// if the mean is minus, it means that edge is overlap, set
					// it
					// to be overlap and not valid;
					if (mean < 0)
						edge.setOL(true);
					if (isUseOLLink)
						edge.setValid(true);
					else
						edge.setValid(false);
					edge.setoStrand(Strand.REVERSE);
					edge.settStrand(Strand.REVERSE);
					edge.setLinkNum(dists.size());
					edges.add(edge);
					// edge + +;
					edge = null;
					// origin = null;
					// terminus = null;
					edge = new Edge();
					Contig tC = origin;
					origin = terminus;
					terminus = tC;
					tC = null;
					// origin = new Contig();
					// terminus = new Contig();
					// origin.setID(ids[1]);
					// origin.setLength(contigs.get(ids[1]).getLength());
					// terminus.setID(ids[0]);
					// terminus.setLength(contigs.get(ids[0]).getLength());
					edge.setOrigin(origin);
					edge.setTerminus(terminus);
					edge.setDistMean(mean);
					edge.setDistSd(sd);
					edge.setFake(false);
					// if the mean is minus, it means that edge is overlap, set
					// it
					// to be overlap and not valid;
					if (mean < 0)
						edge.setOL(true);
					if (isUseOLLink)
						edge.setValid(true);
					else
						edge.setValid(false);
					edge.setoStrand(Strand.FORWARD);
					edge.settStrand(Strand.FORWARD);
					edge.setLinkNum(dists.size());
					edges.add(edge);
					edge = null;
					origin = null;
					terminus = null;
				} else if (max == sD) // for the D statement
				{
					List<Integer> dists = new Vector<Integer>(10);
					for (PBLinkM p : s1Ds)
						dists.add(p.getDistance());
					for (PBLinkM p : s2Ds)
						dists.add(p.getDistance());
					int mean = MathTool.mean(dists);
					int sd = MathTool.sd(dists);
					int upper = mean + 3 * sd;
					int low = mean - 3 * sd;
					// if the distance larger than mean + 3 * sd, then remove;
					Iterator<PBLinkM> it1 = s1Ds.iterator();
					while (it1.hasNext()) {
						PBLinkM pb = it1.next();
						if (pb.getDistance() > upper || pb.getDistance() < low)
							it1.remove();
					}
					Iterator<PBLinkM> it2 = s2Ds.iterator();
					while (it2.hasNext()) {
						PBLinkM pb = it2.next();
						if (pb.getDistance() > upper || pb.getDistance() < low)
							it2.remove();
					}
					// if the supported link count less then user specified
					int supLink = s1Ds.size() + s2Ds.size();
					// logger.debug(this.getClass().getName() + "\tSupported
					// link: "
					// + supLink);
					if (supLink < minSupLink)
						continue;
					// recompute mean and sd;
					dists = null;
					dists = new Vector<Integer>();
					for (PBLinkM p : s1Ds)
						dists.add(p.getDistance());
					for (PBLinkM p : s2Ds)
						dists.add(p.getDistance());
					mean = MathTool.mean(dists);
					sd = MathTool.sd(dists);
					// if the mean is less tan zero and the user specified do
					// not
					// use overlap link;
					if (!isUseOLLink && mean < 0)
						continue;
					// there are very similarity of two contigs with high
					// identity
					// and PB reads also have high error ratio, so it will lead
					// to
					// pb read links two high identity contigs by overlap two
					// large
					// defining a standard when they over larger or less than
					// this
					// filtering parameter, it will omit this edges;
					if (mean < 0) {
						low = Math.abs(mean) - 4 * sd;
						upper = Math.abs(mean) + 4 * sd;
						int c1Len = contigs.get(ids[0]).getLength();
						int c2Len = contigs.get(ids[1]).getLength();
						if (c1Len >= low && c1Len <= upper)
							continue;
						if (c2Len >= low && c2Len <= upper)
							continue;
					}
					// for store gap record information
					if (mean > 0) {
						GapRecord gr = new GapRecord();
						gr.setStart(ids[0]);
						gr.setEnd(ids[1]);
						for (PBLinkM p : s1Ds) {
							if (p.getDistance() > 0) {
								PBGapSeq seq = new PBGapSeq();
								seq.setId(p.getOrigin().getqName());
								int p1 = p.getOrigin().getqEnd();
								int p2 = p.getTerminus().getqStart();
								if (p1 <= p2) {
									seq.setStart(p1);
									seq.setEnd(p2);
								} else {
									seq.setStart(p2);
									seq.setEnd(p1);
								}
								seq.setStrand(Strand.FORWARD);
								gr.addSeq(seq);
							}
						}
						for (PBLinkM p : s2Ds) {
							if (p.getDistance() > 0) {
								PBGapSeq seq = new PBGapSeq();
								seq.setId(p.getOrigin().getqName());
								int p1 = p.getOrigin().getqEnd();
								int p2 = p.getTerminus().getqStart();
								if (p1 <= p2) {
									seq.setStart(p1);
									seq.setEnd(p2);
								} else {
									seq.setStart(p2);
									seq.setEnd(p1);
								}
								seq.setStrand(Strand.REVERSE);
								gr.addSeq(seq);
							}
						}
						gaps.add(gr);
					}
					// edge - +;
					Edge edge = new Edge();
					Contig origin = contigs.get(ids[0]);
					Contig terminus = contigs.get(ids[1]);
					// Contig origin = new Contig();
					// Contig terminus = new Contig();
					// origin.setID(ids[0]);
					// origin.setLength(contigs.get(ids[0]).getLength());
					// terminus.setID(ids[1]);
					// terminus.setLength(contigs.get(ids[1]).getLength());
					edge.setOrigin(origin);
					edge.setTerminus(terminus);
					edge.setDistMean(mean);
					edge.setDistSd(sd);
					edge.setFake(false);
					// if the mean is minus, it means that edge is overlap, set
					// it
					// to be overlap and not valid;
					if (mean < 0)
						edge.setOL(true);
					if (isUseOLLink)
						edge.setValid(true);
					else
						edge.setValid(false);
					edge.setoStrand(Strand.REVERSE);
					edge.settStrand(Strand.FORWARD);
					edge.setLinkNum(dists.size());
					edges.add(edge);
					// edge + -;
					edge = null;
					// origin = null;
					// terminus = null;
					edge = new Edge();
					Contig tC = origin;
					origin = terminus;
					terminus = tC;
					tC = null;
					// origin = new Contig();
					// terminus = new Contig();
					// origin.setID(ids[1]);
					// origin.setLength(contigs.get(ids[1]).getLength());
					// terminus.setID(ids[0]);
					// terminus.setLength(contigs.get(ids[0]).getLength());
					edge.setOrigin(origin);
					edge.setTerminus(terminus);
					edge.setDistMean(mean);
					edge.setDistSd(sd);
					edge.setFake(false);
					// if the mean is minus, it means that edge is overlap, set
					// it
					// to be overlap and not valid;
					if (mean < 0)
						edge.setOL(true);
					if (isUseOLLink)
						edge.setValid(true);
					else
						edge.setValid(false);
					edge.setoStrand(Strand.REVERSE);
					edge.settStrand(Strand.FORWARD);
					edge.setLinkNum(dists.size());
					edges.add(edge);
					edge = null;
					origin = null;
					terminus = null;
				}

				// gc the variables
				all = null;
				s1 = null;
				s2 = null;
				s1As = null;
				s1Bs = null;
				s1Cs = null;
				s1Ds = null;
				s2As = null;
				s2Bs = null;
				s2Cs = null;
				s2Ds = null;
			} catch (Exception e) {
				logger.debug(this.getClass().getName() + "\tindex = " + index);
				logger.debug(this.getClass().getName() + "\t" + e.getMessage());
				logger.error(this.getClass().getName() + "\t" + e.getMessage());
			}
		}

		GapRecordWriter grw = new GapRecordWriter(paras, gaps);
		grw.write();
		gaps = null;
		return edges;
	}

	public List<Edge> pesudoEdging() {
		if (edges == null || edges.size() == 0)
			throw new IllegalArgumentException(
					this.getClass().getName() + "The Edges is empty, could not be pesudo edging!");
		// define the original size of new Edge list is to be 1.5 folds;
		List<Edge> values = new Vector<Edge>(edges.size() / 2 + edges.size());
		for (Edge e : edges) {
			values.add(e);
			Edge tEdge = new Edge();
			tEdge.setOrigin(e.getTerminus());
			tEdge.setTerminus(e.getOrigin());
			if (!edges.contains(tEdge)) {
				tEdge.setFake(true);
				tEdge.setDistMean(e.getDistMean());
				tEdge.setDistSd(e.getDistSd());
				tEdge.setLinkNum(e.getLinkNum());
				tEdge.setOL(e.isOL());
				tEdge.setValid(e.isValid());
				if (e.getoStrand().equals(Strand.FORWARD))
					tEdge.settStrand(Strand.REVERSE);
				else
					tEdge.settStrand(Strand.FORWARD);
				if (e.gettStrand().equals(Strand.FORWARD))
					tEdge.setoStrand(Strand.REVERSE);
				else
					tEdge.setoStrand(Strand.FORWARD);
				values.add(tEdge);
			}
		}
		return values;
	}
	
	private void initCntIndexer()
	{
		try {
			String path = paras.getOutFolder() + System.getProperty("file.separator") + "cnt.index";
			directory = new SimpleFSDirectory(new File(path).toPath());
			reader = DirectoryReader.open(directory);
			searcher = new IndexSearcher(reader);
			analyzer = new StandardAnalyzer();
			parser = new QueryParser("id", analyzer);
		} catch (IOException e) {
			logger.error(this.getClass().getName() + "\t" + e.getMessage());
		}
	}
	
	private int indexCntLength(String id)
	{
		int len = 0;
		try {
			Query query = parser.parse(id);
			TopDocs tds = searcher.search(query, 10);
			for (ScoreDoc sd : tds.scoreDocs) {
				Document doc = searcher.doc(sd.doc);
				len = Integer.valueOf(doc.get("len"));
			}
		} catch (Exception e) {
			logger.error(this.getClass().getName() + "\t" + e.getMessage());
		}
		return len;
	}
	
	private String indexCntSeq(String id)
	{
		String seq = "";
		try {
			Query query = parser.parse(id);
			TopDocs tds = searcher.search(query, 10);
			for (ScoreDoc sd : tds.scoreDocs) {
				Document doc = searcher.doc(sd.doc);
				seq = doc.get("seq");
			}
		} catch (Exception e) {
			logger.error(this.getClass().getName() + "\t" + e.getMessage());
		}
		return seq;
	}
}

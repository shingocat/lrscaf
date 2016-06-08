/*
*File: agis.ps.util.EdgeBundler.java
*User: mqin
*Email: mqin@ymail.com
*Date: 2016年1月13日
*/
package agis.ps.util;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Vector;

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

	public EdgeBundler() {
		// do nothing;
	}

	public EdgeBundler(List<PBLinkM> links, Parameter paras) {
		this.links = links;
		this.paras = paras;
	}
	
	public EdgeBundler(List<PBLinkM> links, Parameter paras, Map<String, Contig> contigs)
	{
		this(links, paras);
		this.contigs = contigs;
	}

	public List<Edge> pbLinkM5Bundling() {
		if (links == null || links.size() == 0)
			throw new IllegalArgumentException(this.getClass().getName() + "The PBLinkM5 could not be empty!");
		return this.pbLinkM5Bundling2(links, paras);
	}
	
	// the second implemented
	public List<Edge> pbLinkM5Bundling2(List<PBLinkM> links, Parameter paras)
	{
		if(links == null || links.size() == 0)
			throw new IllegalArgumentException(
					this.getClass().getName() + "The Links is empty when passed to EdgeBundler");
		boolean isUseOLLink = paras.isUseOLLink();
		if(edges != null)
			edges = null;
		edges = new Vector<Edge>(200);
		int minSupLink = paras.getMinSupLinks();
		logger.debug(this.getClass().getName() + "\tMinimum supported links:" + minSupLink);
		logger.debug(this.getClass().getName() + "\tUsed overlap link:" + isUseOLLink);
		// the all the same origin and terminus to a hash map;
		// both direction
		Map<String, List<PBLinkM>> temp = new HashMap<String, List<PBLinkM>>();
		for(PBLinkM pb : links)
		{
			String id1 = pb.getOrigin().gettName() + ":->:" + pb.getTerminus().gettName();
			String id2 = pb.getTerminus().gettName() + ":->:" + pb.getOrigin().gettName();
			if(temp.containsKey(id1))
			{
				temp.get(id1).add(pb);
				id1 = null;
				id2 = null;
				pb = null;
			} else if(temp.containsKey(id2))
			{
				temp.get(id2).add(pb);
				id1 = null;
				id2 = null;
				pb = null;
			} else
			{
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
		for(String s : temp.keySet())
		{
			// divide into two set, i.e. A->B; B->A
			List<PBLinkM> all = temp.get(s);
			List<PBLinkM> s1 = new Vector<PBLinkM>(10);
			List<PBLinkM> s2 = new Vector<PBLinkM>(10);
			String [] ids = s.split(":->:");
			for(PBLinkM p : all)
			{
				if(ids[0].equalsIgnoreCase(p.getOrigin().gettName()))
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
			for(PBLinkM p : s1)
			{
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
//			FOR S2 SET
			List<PBLinkM> s2As = new Vector<PBLinkM>(10);
			List<PBLinkM> s2Bs = new Vector<PBLinkM>(10);
			List<PBLinkM> s2Cs = new Vector<PBLinkM>(10);
			List<PBLinkM> s2Ds = new Vector<PBLinkM>(10);
			for(PBLinkM p : s2)
			{
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
			if(max == sA) // for the A statement;
			{
				List<Integer> dists = new Vector<Integer>(10);
				for(PBLinkM p : s1As)
					dists.add(p.getDistance());
				for(PBLinkM p : s2Cs)
					dists.add(p.getDistance());
				int mean = MathTool.mean(dists);
				int sd = MathTool.sd(dists);
				int upper = mean + 3 * sd;
				int low = mean - 3 * sd;
				// if the distance larger than mean + 3 * sd, then remove;
				try {
					for (int i = 0; i < s1As.size(); i++) {
						PBLinkM pb = s1As.get(i);
						if (pb.getDistance() > upper || pb.getDistance() < low)
							s1As.remove(pb);
					}
					for (int i = 0; i < s2Cs.size(); i++) {
						PBLinkM pb = s2Cs.get(i);
						if (pb.getDistance() > upper || pb.getDistance() < low)
							s2Cs.remove(pb);
					}
				} catch (Exception e) {
					logger.debug(this.getClass().getName() + e.getMessage() + "\t" + e.getClass().getName());
					logger.error(this.getClass().getName() + e.getMessage() + "\t" + e.getClass().getName());
					logger.info(this.getClass().getName() + e.getMessage() + "\t" + e.getClass().getName());
					continue;
				}
				// if the supported link count less then user specified
				int supLink = s1As.size() + s2Cs.size();
//				logger.debug(this.getClass().getName() + "\tSupported link: " + supLink);
				if( supLink < minSupLink)
					continue;
				// recompute mean and sd;
				dists = null;
				dists = new Vector<Integer>();
				for(PBLinkM p : s1As)
					dists.add(p.getDistance());
				for(PBLinkM p : s2Cs)
					dists.add(p.getDistance());
				mean = MathTool.mean(dists);
				sd = MathTool.sd(dists);
				// if the mean is less tan zero and the user specified do not use overlap link;
				if(!isUseOLLink && mean < 0)
					continue;
				// for store gap record information
				if(mean > 0)
				{
					GapRecord gr = new GapRecord();
					gr.setStart(ids[0]);
					gr.setEnd(ids[1]);
					for(PBLinkM p : s1As)
					{
						if(p.getDistance() > 0)
						{
							PBGapSeq seq = new PBGapSeq();
							seq.setId(p.getOrigin().getqName());
							int p1 = p.getOrigin().getqEnd();
							int p2 = p.getTerminus().getqStart();
							if(p1 <= p2)
							{
								seq.setStart(p1);
								seq.setEnd(p2);
							} else
							{
								seq.setStart(p2);
								seq.setEnd(p1);
							}
							seq.setStrand(Strand.FORWARD);
							gr.addSeq(seq);
						}
					}
					for(PBLinkM p : s2Cs)
					{
						if(p.getDistance() > 0)
						{
							PBGapSeq seq = new PBGapSeq();
							seq.setId(p.getOrigin().getqName());
							int p1 = p.getOrigin().getqEnd();
							int p2 = p.getTerminus().getqStart();
							if(p1 <= p2)
							{
								seq.setStart(p1);
								seq.setEnd(p2);
							} else
							{
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
				Contig origin = new Contig();
				Contig terminus = new Contig();
				origin.setID(ids[0]);
				origin.setLength(contigs.get(ids[0]).getLength());
				terminus.setID(ids[1]);
				terminus.setLength(contigs.get(ids[1]).getLength());
				edge.setOrigin(origin);
				edge.setTerminus(terminus);
				edge.setDistMean(mean);
				edge.setDistSd(sd);
				edge.setFake(false);
				// if the mean is minus, it means that edge is overlap, set it
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
//				 edge - -;
				edge = null;
				origin = null;
				terminus = null;
				edge = new Edge();
				origin = new Contig();
				terminus = new Contig();
				origin.setID(ids[1]);
				origin.setLength(contigs.get(ids[1]).getLength());
				terminus.setID(ids[0]);
				terminus.setLength(contigs.get(ids[0]).getLength());
				edge.setOrigin(origin);
				edge.setTerminus(terminus);
				edge.setDistMean(mean);
				edge.setDistSd(sd);
				edge.setFake(false);
				// if the mean is minus, it means that edge is overlap, set it
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
			} else if(max == sB) // for the B statement
			{
				List<Integer> dists = new Vector<Integer>(10);
				for(PBLinkM p : s1Bs)
					dists.add(p.getDistance());
				for(PBLinkM p : s2Bs)
					dists.add(p.getDistance());
				int mean = MathTool.mean(dists);
				int sd = MathTool.sd(dists);
				int upper = mean + 3 * sd;
				int low = mean - 3 * sd;
				// if the distance larger than mean + 3 * sd, then remove;
				try {
					for (int i = 0; i < s1Bs.size(); i++) {
						PBLinkM pb = s1Bs.get(i);
						if (pb.getDistance() > upper || pb.getDistance() < low)
							s1Bs.remove(pb);
					}
					for (int i = 0; i < s2Bs.size(); i++) {
						PBLinkM pb = s2Bs.get(i);
						if (pb.getDistance() > upper || pb.getDistance() < low)
							s2Bs.remove(pb);
					}
				} catch (Exception e) {
					logger.debug(this.getClass().getName() + e.getMessage() + "\t" + e.getClass().getName());
					logger.error(this.getClass().getName() + e.getMessage() + "\t" + e.getClass().getName());
					logger.info(this.getClass().getName() + e.getMessage() + "\t" + e.getClass().getName());
					continue;
				}
				// if the supported link count less then user specified
				int supLink = s1Bs.size() + s2Bs.size();
//				logger.debug(this.getClass().getName() + "\tSupported link: " + supLink);
				if( supLink < minSupLink)
					continue;
				// recompute mean and sd;
				dists = null;
				dists = new Vector<Integer>();
				for(PBLinkM p : s1Bs)
					dists.add(p.getDistance());
				for(PBLinkM p : s2Bs)
					dists.add(p.getDistance());
				mean = MathTool.mean(dists);
				sd = MathTool.sd(dists);
				// if the mean is less tan zero and the user specified do not use overlap link;
				if(!isUseOLLink && mean < 0)
					continue;
				// for store gap record information
				if(mean > 0)
				{
					GapRecord gr = new GapRecord();
					gr.setStart(ids[0]);
					gr.setEnd(ids[1]);
					for(PBLinkM p : s1Bs)
					{
						if(p.getDistance() > 0)
						{
							PBGapSeq seq = new PBGapSeq();
							seq.setId(p.getOrigin().getqName());
							int p1 = p.getOrigin().getqEnd();
							int p2 = p.getTerminus().getqStart();
							if(p1 <= p2)
							{
								seq.setStart(p1);
								seq.setEnd(p2);
							} else
							{
								seq.setStart(p2);
								seq.setEnd(p1);
							}
							seq.setStrand(Strand.FORWARD);
							gr.addSeq(seq);
						}
					}
					for(PBLinkM p : s2Bs)
					{
						if(p.getDistance() > 0)
						{
							PBGapSeq seq = new PBGapSeq();
							seq.setId(p.getOrigin().getqName());
							int p1 = p.getOrigin().getqEnd();
							int p2 = p.getTerminus().getqStart();
							if(p1 <= p2)
							{
								seq.setStart(p1);
								seq.setEnd(p2);
							} else
							{
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
				Contig origin = new Contig();
				Contig terminus = new Contig();
				origin.setID(ids[0]);
				origin.setLength(contigs.get(ids[0]).getLength());
				terminus.setID(ids[1]);
				terminus.setLength(contigs.get(ids[1]).getLength());
				edge.setOrigin(origin);
				edge.setTerminus(terminus);
				edge.setDistMean(mean);
				edge.setDistSd(sd);
				edge.setFake(false);
				// if the mean is minus, it means that edge is overlap, set it
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
//				 edge - +;
				edge = null;
				origin = null;
				terminus = null;
				edge = new Edge();
				origin = new Contig();
				terminus = new Contig();
				origin.setID(ids[1]);
				origin.setLength(contigs.get(ids[1]).getLength());
				terminus.setID(ids[0]);
				terminus.setLength(contigs.get(ids[0]).getLength());
				edge.setOrigin(origin);
				edge.setTerminus(terminus);
				edge.setDistMean(mean);
				edge.setDistSd(sd);
				edge.setFake(false);
				// if the mean is minus, it means that edge is overlap, set it
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
			} else if(max == sC) // for the C statement
			{
				List<Integer> dists = new Vector<Integer>(10);
				for(PBLinkM p : s1Cs)
					dists.add(p.getDistance());
				for(PBLinkM p : s2As)
					dists.add(p.getDistance());
				int mean = MathTool.mean(dists);
				int sd = MathTool.sd(dists);
				int upper = mean + 3 * sd;
				int low = mean - 3 * sd;
				// if the distance larger than mean + 3 * sd, then remove;
				try {
					for (int i = 0; i < s1Cs.size(); i++) {
						PBLinkM pb = s1Cs.get(i);
						if (pb.getDistance() > upper || pb.getDistance() < low)
							s1Cs.remove(pb);
					}
					for (int i = 0; i < s2As.size(); i++) {
						PBLinkM pb = s2As.get(i);
						if (pb.getDistance() > upper || pb.getDistance() < low)
							s2As.remove(pb);
					}
				} catch (Exception e) {
					logger.debug(this.getClass().getName() + e.getMessage() + "\t" + e.getClass().getName());
					logger.error(this.getClass().getName() + e.getMessage() + "\t" + e.getClass().getName());
					logger.info(this.getClass().getName() + e.getMessage() + "\t" + e.getClass().getName());
					continue;
				}
				// if the supported link count less then user specified
				int supLink = s1Cs.size() + s2As.size();
//				logger.debug(this.getClass().getName() + "\tSupported link: " + supLink);
				if( supLink < minSupLink)
					continue;
				// recompute mean and sd;
				dists = null;
				dists = new Vector<Integer>();
				for(PBLinkM p : s1Cs)
					dists.add(p.getDistance());
				for(PBLinkM p : s2As)
					dists.add(p.getDistance());
				mean = MathTool.mean(dists);
				sd = MathTool.sd(dists);
				// if the mean is less tan zero and the user specified do not use overlap link;
				if(!isUseOLLink && mean < 0)
					continue;
				// for store gap record information
				if(mean > 0)
				{
					GapRecord gr = new GapRecord();
					gr.setStart(ids[0]);
					gr.setEnd(ids[1]);
					for(PBLinkM p : s1Cs)
					{
						if(p.getDistance() > 0)
						{
							PBGapSeq seq = new PBGapSeq();
							seq.setId(p.getOrigin().getqName());
							int p1 = p.getOrigin().getqEnd();
							int p2 = p.getTerminus().getqStart();
							if(p1 <= p2)
							{
								seq.setStart(p1);
								seq.setEnd(p2);
							} else
							{
								seq.setStart(p2);
								seq.setEnd(p1);
							}
							seq.setStrand(Strand.FORWARD);
							gr.addSeq(seq);
						}
					}
					for(PBLinkM p : s2As)
					{
						if(p.getDistance() > 0)
						{
							PBGapSeq seq = new PBGapSeq();
							seq.setId(p.getOrigin().getqName());
							int p1 = p.getOrigin().getqEnd();
							int p2 = p.getTerminus().getqStart();
							if(p1 <= p2)
							{
								seq.setStart(p1);
								seq.setEnd(p2);
							} else
							{
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
				Contig origin = new Contig();
				Contig terminus = new Contig();
				origin.setID(ids[0]);
				origin.setLength(contigs.get(ids[0]).getLength());
				terminus.setID(ids[1]);
				terminus.setLength(contigs.get(ids[1]).getLength());
				edge.setOrigin(origin);
				edge.setTerminus(terminus);
				edge.setDistMean(mean);
				edge.setDistSd(sd);
				edge.setFake(false);
				// if the mean is minus, it means that edge is overlap, set it
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
//				 edge + +;
				edge = null;
				origin = null;
				terminus = null;
				edge = new Edge();
				origin = new Contig();
				terminus = new Contig();
				origin.setID(ids[1]);
				origin.setLength(contigs.get(ids[1]).getLength());
				terminus.setID(ids[0]);
				terminus.setLength(contigs.get(ids[0]).getLength());
				edge.setOrigin(origin);
				edge.setTerminus(terminus);
				edge.setDistMean(mean);
				edge.setDistSd(sd);
				edge.setFake(false);
				// if the mean is minus, it means that edge is overlap, set it
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
			} else if(max == sD) // for the D statement
			{
				List<Integer> dists = new Vector<Integer>(10);
				for(PBLinkM p : s1Ds)
					dists.add(p.getDistance());
				for(PBLinkM p : s2Ds)
					dists.add(p.getDistance());
				int mean = MathTool.mean(dists);
				int sd = MathTool.sd(dists);
				int upper = mean + 3 * sd;
				int low = mean - 3 * sd;
				// if the distance larger than mean + 3 * sd, then remove;
				try {
					for (int i = 0; i < s1Ds.size(); i++) {
						PBLinkM pb = s1Ds.get(i);
						if (pb.getDistance() > upper || pb.getDistance() < low)
							s1Ds.remove(pb);
					}
					for (int i = 0; i < s2Ds.size(); i++) {
						PBLinkM pb = s2Ds.get(i);
						if (pb.getDistance() > upper || pb.getDistance() < low)
							s2Ds.remove(pb);
					}
				} catch (Exception e) {
					logger.debug(this.getClass().getName() + e.getMessage() + "\t" + e.getClass().getName());
					logger.error(this.getClass().getName() + e.getMessage() + "\t" + e.getClass().getName());
					logger.info(this.getClass().getName() + e.getMessage() + "\t" + e.getClass().getName());
					continue;
				}
				// if the supported link count less then user specified
				int supLink = s1Ds.size() + s2Ds.size();
//				logger.debug(this.getClass().getName() + "\tSupported link: " + supLink);
				if( supLink < minSupLink)
					continue;
				// recompute mean and sd;
				dists = null;
				dists = new Vector<Integer>();
				for(PBLinkM p : s1Ds)
					dists.add(p.getDistance());
				for(PBLinkM p : s2Ds)
					dists.add(p.getDistance());
				mean = MathTool.mean(dists);
				sd = MathTool.sd(dists);
				// if the mean is less tan zero and the user specified do not use overlap link;
				if(!isUseOLLink && mean < 0)
					continue;
				// for store gap record information
				if(mean > 0)
				{
					GapRecord gr = new GapRecord();
					gr.setStart(ids[0]);
					gr.setEnd(ids[1]);
					for(PBLinkM p : s1Ds)
					{
						if(p.getDistance() > 0)
						{
							PBGapSeq seq = new PBGapSeq();
							seq.setId(p.getOrigin().getqName());
							int p1 = p.getOrigin().getqEnd();
							int p2 = p.getTerminus().getqStart();
							if(p1 <= p2)
							{
								seq.setStart(p1);
								seq.setEnd(p2);
							} else
							{
								seq.setStart(p2);
								seq.setEnd(p1);
							}
							seq.setStrand(Strand.FORWARD);
							gr.addSeq(seq);
						}
					}
					for(PBLinkM p : s2Ds)
					{
						if(p.getDistance() > 0)
						{
							PBGapSeq seq = new PBGapSeq();
							seq.setId(p.getOrigin().getqName());
							int p1 = p.getOrigin().getqEnd();
							int p2 = p.getTerminus().getqStart();
							if(p1 <= p2)
							{
								seq.setStart(p1);
								seq.setEnd(p2);
							} else
							{
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
				Contig origin = new Contig();
				Contig terminus = new Contig();
				origin.setID(ids[0]);
				origin.setLength(contigs.get(ids[0]).getLength());
				terminus.setID(ids[1]);
				terminus.setLength(contigs.get(ids[1]).getLength());
				edge.setOrigin(origin);
				edge.setTerminus(terminus);
				edge.setDistMean(mean);
				edge.setDistSd(sd);
				edge.setFake(false);
				// if the mean is minus, it means that edge is overlap, set it
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
//				 edge + -;
				edge = null;
				origin = null;
				terminus = null;
				edge = new Edge();
				origin = new Contig();
				terminus = new Contig();
				origin.setID(ids[1]);
				origin.setLength(contigs.get(ids[1]).getLength());
				terminus.setID(ids[0]);
				terminus.setLength(contigs.get(ids[0]).getLength());
				edge.setOrigin(origin);
				edge.setTerminus(terminus);
				edge.setDistMean(mean);
				edge.setDistSd(sd);
				edge.setFake(false);
				// if the mean is minus, it means that edge is overlap, set it
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
		}
		GapRecordWriter grw = new GapRecordWriter(paras, gaps);
		grw.write();
		gaps = null;
		return edges;
	}
	

	public List<Edge> pbLinkM5Bundling(List<PBLinkM> links, Parameter paras) {
		if (links == null || links.size() == 0)
			throw new IllegalArgumentException(
					this.getClass().getName() + "The Links is empty when passed to EdgeBundler");
		boolean isUseOLLink = paras.isUseOLLink();
		if (edges != null)
			edges = null;
		// setting default size of vector is 200 elements;
		edges = new Vector<Edge>(200);
		// storing all the same origin and terminus to a hash map;
		Map<String, List<PBLinkM>> temp = new HashMap<String, List<PBLinkM>>();
		for (PBLinkM pb : links) {
			// not include overlap links
			// if (pb.isOverLap())
			// continue;
			String id = pb.getOrigin().gettName() + ":->:" + pb.getTerminus().gettName();
			if (temp.containsKey(id)) {
				temp.get(id).add(pb);
				id = null;
				pb = null;
			} else {
				List<PBLinkM> tLinks = new Vector<PBLinkM>();
				tLinks.add(pb);
				temp.put(id, tLinks);
				tLinks = null;
				id = null;
				pb = null;
			}
		}

		// build the edge for fitted for criterion;
		int count = 0;
		for (String s : temp.keySet()) {
			// if the support links less than minimum and larger than maximum
			// support links, it will omitted;
			int numLinks = temp.get(s).size();
			if (numLinks >= paras.getMinSupLinks() && numLinks <= paras.getMaxSupLinks()) {
				count += 1;
				// statistical analysis contig pairs distance
				List<PBLinkM> tlinks = temp.get(s);
//				List<Integer> dists = new Vector<Integer>();
//				for (PBLinkM5 pb : tlinks) {
//					dists.add(pb.getDistance());
//				}
//
//				
//				int	mean = MathTool.mean(dists);
//				int	sd = MathTool.sd(dists);
//				
//				int upper = mean + 3 * sd;
//				int low = mean - 3 * sd;
//				// if the distance larger than mean + 2 * sd, then remove;
//				try {
//					// This will show iteration error when remove element
//					// for (PBLink pb : tlinks) {
//					// if (pb.getDistance() > upper || pb.getDistance() < low)
//					// tlinks.remove(pb);
//					// }
//					for (int i = 0; i < tlinks.size(); i++) {
//						PBLinkM5 pb = tlinks.get(i);
//						if (pb.getDistance() > upper || pb.getDistance() < low)
//							tlinks.remove(pb);
//					}
//				} catch (Exception e) {
//					logger.debug(this.getClass().getName() + e.getMessage() + "\t" + e.getClass().getName());
//					logger.error(this.getClass().getName() + e.getMessage() + "\t" + e.getClass().getName());
//					logger.info(this.getClass().getName() + e.getMessage() + "\t" + e.getClass().getName());
//					continue;
//				}
				// checking the most frequencies contig pair type, and retain
				// the most pair type; and recompute the mean and sd after
				// remove outlier
				int typeA = 0; // + +;
				List<Integer> typeADists = new Vector<Integer>(10);
				List<PBLinkM> typeALinks = new Vector<PBLinkM>(10);
				int typeB = 0; // + -;
				List<Integer> typeBDists = new Vector<Integer>(10); 
				List<PBLinkM> typeBLinks = new Vector<PBLinkM>(10);
				int typeC = 0; // - -;
				List<Integer> typeCDists = new Vector<Integer>(10);
				List<PBLinkM> typeCLinks = new Vector<PBLinkM>(10);
				int typeD = 0; // - +;
				List<Integer> typeDDists = new Vector<Integer>(10);
				List<PBLinkM> typeDLinks = new Vector<PBLinkM>(10);

//				dists = null;
//				dists = new Vector<Integer>();
				Integer oLen = null; // original contig length;
				Integer tLen = null; // terminus contig length;
				for (PBLinkM pb : tlinks) {
//					dists.add(pb.getDistance());
					if (oLen == null)
						oLen = pb.getOrigin().gettLength();
					if (tLen == null)
						tLen = pb.getTerminus().gettLength();
					Strand oStrand = pb.getOrigin().gettStrand();
					Strand tStrand = pb.getTerminus().gettStrand();
					if (oStrand.equals(Strand.FORWARD)) {
						if (tStrand.equals(Strand.FORWARD)) {
							typeA += 1;
							typeALinks.add(pb);
							typeADists.add(pb.getDistance());
						} else {
							typeB += 1;
							typeBLinks.add(pb);
							typeBDists.add(pb.getDistance());
						}
					} else {
						if (tStrand.equals(Strand.REVERSE)) {
							typeC += 1;
							typeCLinks.add(pb);
							typeCDists.add(pb.getDistance());
						} else {
							typeD += 1;
							typeDLinks.add(pb);
							typeDDists.add(pb.getDistance());
						}
					}
				}
				int max = MathTool.max(typeA, typeB, typeC, typeD);
				// recompute the mean and sd after remove outlier
				// dists.clear();
				// for (PBLinkM5 pb : tlinks) {
				// dists.add(pb.getDistance());
				// }
				int mean = 0;
				int sd = 0;
				if(max == typeA)
				{
					mean = MathTool.mean(typeADists);
					sd = MathTool.sd(typeADists);
					
					int upper = mean + 3 * sd;
					int low = mean - 3 * sd;
					// if the distance larger than mean + 2 * sd, then remove;
					try {
						// This will show iteration error when remove element
						// for (PBLink pb : tlinks) {
						// if (pb.getDistance() > upper || pb.getDistance() < low)
						// tlinks.remove(pb);
						// }
						for (int i = 0; i < typeALinks.size(); i++) {
							PBLinkM pb = typeALinks.get(i);
							if (pb.getDistance() > upper || pb.getDistance() < low)
								typeALinks.remove(pb);
						}
					} catch (Exception e) {
						logger.debug(this.getClass().getName() + e.getMessage() + "\t" + e.getClass().getName());
						logger.error(this.getClass().getName() + e.getMessage() + "\t" + e.getClass().getName());
						logger.info(this.getClass().getName() + e.getMessage() + "\t" + e.getClass().getName());
						continue;
					}
					
					// recompute mean and sd
					typeADists = null;
					typeADists = new Vector<Integer>(10);
					for(PBLinkM p : typeALinks)
					{
						typeADists.add(p.getDistance());
					}
					mean = MathTool.mean(typeADists);
					sd = MathTool.sd(typeADists);
					
					// initiated the edge;
					Edge edge = new Edge();
					Contig origin = new Contig();
					Contig terminus = new Contig();
					String[] arr = s.split(":->:");
					origin.setID(arr[0]);
					origin.setLength(oLen);
					terminus.setID(arr[1]);
					terminus.setLength(tLen);
					edge.setOrigin(origin);
					edge.setTerminus(terminus);
					edge.setDistMean(mean);
					edge.setDistSd(sd);
					edge.setFake(false);
					// if the mean is minus, it means that edge is overlap, set it
					// to be overlap and not valid;
					if (mean < 0)
						edge.setOL(true);
					if (isUseOLLink)
						edge.setValid(true);
					else
						edge.setValid(false);
					edge.setoStrand(Strand.FORWARD);
					edge.settStrand(Strand.FORWARD);
					edge.setLinkNum(typeA);
					edges.add(edge);
				} else if(max == typeB)
				{
					mean = MathTool.mean(typeBDists);
					sd = MathTool.sd(typeBDists);
					
					int upper = mean + 3 * sd;
					int low = mean - 3 * sd;
					// if the distance larger than mean + 2 * sd, then remove;
					try {
						// This will show iteration error when remove element
						// for (PBLink pb : tlinks) {
						// if (pb.getDistance() > upper || pb.getDistance() < low)
						// tlinks.remove(pb);
						// }
						for (int i = 0; i < typeBLinks.size(); i++) {
							PBLinkM pb = typeBLinks.get(i);
							if (pb.getDistance() > upper || pb.getDistance() < low)
								typeBLinks.remove(pb);
						}
					} catch (Exception e) {
						logger.debug(this.getClass().getName() + e.getMessage() + "\t" + e.getClass().getName());
						logger.error(this.getClass().getName() + e.getMessage() + "\t" + e.getClass().getName());
						logger.info(this.getClass().getName() + e.getMessage() + "\t" + e.getClass().getName());
						continue;
					}
					// recompute mean and sd
					typeBDists = null;
					typeBDists = new Vector<Integer>(10);
					for(PBLinkM p : typeBLinks)
					{
						typeBDists.add(p.getDistance());
					}
					mean = MathTool.mean(typeBDists);
					sd = MathTool.sd(typeBDists);
					
					// initiated the edge;
					Edge edge = new Edge();
					Contig origin = new Contig();
					Contig terminus = new Contig();
					String[] arr = s.split(":->:");
					origin.setID(arr[0]);
					origin.setLength(oLen);
					terminus.setID(arr[1]);
					terminus.setLength(tLen);
					edge.setOrigin(origin);
					edge.setTerminus(terminus);
					edge.setDistMean(mean);
					edge.setDistSd(sd);
					edge.setFake(false);
					// if the mean is minus, it means that edge is overlap, set it
					// to be overlap and not valid;
					if (mean < 0)
						edge.setOL(true);
					if (isUseOLLink)
						edge.setValid(true);
					else
						edge.setValid(false);
					edge.setoStrand(Strand.FORWARD);
					edge.settStrand(Strand.REVERSE);
					edge.setLinkNum(typeB);
					edges.add(edge);
				} else if(max == typeC)
				{
					mean = MathTool.mean(typeCDists);
					sd = MathTool.sd(typeCDists);
					
					int upper = mean + 3 * sd;
					int low = mean - 3 * sd;
					// if the distance larger than mean + 2 * sd, then remove;
					try {
						// This will show iteration error when remove element
						// for (PBLink pb : tlinks) {
						// if (pb.getDistance() > upper || pb.getDistance() < low)
						// tlinks.remove(pb);
						// }
						for (int i = 0; i < typeCLinks.size(); i++) {
							PBLinkM pb = typeCLinks.get(i);
							if (pb.getDistance() > upper || pb.getDistance() < low)
								typeCLinks.remove(pb);
						}
					} catch (Exception e) {
						logger.debug(this.getClass().getName() + e.getMessage() + "\t" + e.getClass().getName());
						logger.error(this.getClass().getName() + e.getMessage() + "\t" + e.getClass().getName());
						logger.info(this.getClass().getName() + e.getMessage() + "\t" + e.getClass().getName());
						continue;
					}
					// recompute mean and sd
					typeCDists = null;
					typeCDists = new Vector<Integer>(10);
					for(PBLinkM p : typeCLinks)
					{
						typeCDists.add(p.getDistance());
					}
					mean = MathTool.mean(typeCDists);
					sd = MathTool.sd(typeCDists);
					
					// initiated the edge;
					Edge edge = new Edge();
					Contig origin = new Contig();
					Contig terminus = new Contig();
					String[] arr = s.split(":->:");
					origin.setID(arr[0]);
					origin.setLength(oLen);
					terminus.setID(arr[1]);
					terminus.setLength(tLen);
					edge.setOrigin(origin);
					edge.setTerminus(terminus);
					edge.setDistMean(mean);
					edge.setDistSd(sd);
					edge.setFake(false);
					// if the mean is minus, it means that edge is overlap, set it
					// to be overlap and not valid;
					if (mean < 0)
						edge.setOL(true);
					if (isUseOLLink)
						edge.setValid(true);
					else
						edge.setValid(false);
					edge.setoStrand(Strand.REVERSE);
					edge.settStrand(Strand.REVERSE);
					edge.setLinkNum(typeC);
					edges.add(edge);
				} else if(max == typeD)
				{
					mean = MathTool.mean(typeDDists);
					sd = MathTool.sd(typeDDists);
					
					int upper = mean + 3 * sd;
					int low = mean - 3 * sd;
					// if the distance larger than mean + 2 * sd, then remove;
					try {
						// This will show iteration error when remove element
						// for (PBLink pb : tlinks) {
						// if (pb.getDistance() > upper || pb.getDistance() < low)
						// tlinks.remove(pb);
						// }
						for (int i = 0; i < typeDLinks.size(); i++) {
							PBLinkM pb = typeDLinks.get(i);
							if (pb.getDistance() > upper || pb.getDistance() < low)
								typeDLinks.remove(pb);
						}
					} catch (Exception e) {
						logger.debug(this.getClass().getName() + e.getMessage() + "\t" + e.getClass().getName());
						logger.error(this.getClass().getName() + e.getMessage() + "\t" + e.getClass().getName());
						logger.info(this.getClass().getName() + e.getMessage() + "\t" + e.getClass().getName());
						continue;
					}
					// recompute mean and sd
					typeDDists = null;
					typeDDists = new Vector<Integer>(10);
					for(PBLinkM p : typeDLinks)
					{
						typeDDists.add(p.getDistance());
					}
					mean = MathTool.mean(typeDDists);
					sd = MathTool.sd(typeDDists);
					
					// initiated the edge;
					Edge edge = new Edge();
					Contig origin = new Contig();
					Contig terminus = new Contig();
					String[] arr = s.split(":->:");
					origin.setID(arr[0]);
					origin.setLength(oLen);
					terminus.setID(arr[1]);
					terminus.setLength(tLen);
					edge.setOrigin(origin);
					edge.setTerminus(terminus);
					edge.setDistMean(mean);
					edge.setDistSd(sd);
					edge.setFake(false);
					// if the mean is minus, it means that edge is overlap, set it
					// to be overlap and not valid;
					if (mean < 0)
						edge.setOL(true);
					if (isUseOLLink)
						edge.setValid(true);
					else
						edge.setValid(false);
					edge.setoStrand(Strand.REVERSE);
					edge.settStrand(Strand.FORWARD);
					edge.setLinkNum(typeD);
					edges.add(edge);
				}
				
//				int mean = MathTool.mean(typeADists);
//				int sd = MathTool.sd(typeADists);
//				
//				// initiated the edge;
//				Edge edge = new Edge();
//				Contig origin = new Contig();
//				Contig terminus = new Contig();
//				String[] arr = s.split(":->:");
//				origin.setID(arr[0]);
//				origin.setLength(oLen);
//				terminus.setID(arr[1]);
//				terminus.setLength(tLen);
//				edge.setOrigin(origin);
//				edge.setTerminus(terminus);
//				edge.setDistMean(mean);
//				edge.setDistSd(sd);
//				edge.setFake(false);
//				// if the mean is minus, it means that edge is overlap, set it
//				// to be overlap and not valid;
//				if (mean < 0)
//					edge.setOL(true);
//				if (isUseOLLink)
//					edge.setValid(true);
//				else
//					edge.setValid(false);
//				// edge.setOrigin(temp.get(s).get(0).getOrigin());
//				// edge.setTerminus(temp.get(s).get(0).getTerminus());
//				if (max == typeA) {
//					edge.setoStrand(Strand.FORWARD);
//					edge.settStrand(Strand.FORWARD);
//					edge.setLinkNum(typeA);
//				} else if (max == typeB) {
//					edge.setoStrand(Strand.FORWARD);
//					edge.settStrand(Strand.REVERSE);
//					edge.setLinkNum(typeB);
//				} else if (max == typeC) {
//					edge.setoStrand(Strand.REVERSE);
//					edge.settStrand(Strand.REVERSE);
//					edge.setLinkNum(typeC);
//				} else if (max == typeD) {
//					edge.setoStrand(Strand.REVERSE);
//					edge.settStrand(Strand.FORWARD);
//					edge.setLinkNum(typeD);
//				}
//				edges.add(edge);
			}
		}

		// merge two opposite direction edge to only one edge base
		// not implemented now
		return edges;
	}

	public List<Edge> pbLinkBundling(List<PBLink> links, Map<String, String> paras) throws Exception {
		if (edges != null)
			edges = null;
		edges = new Vector<Edge>();
		// storing all the same origin and terminus to a hash map;
		Map<String, List<PBLink>> temp = new HashMap<String, List<PBLink>>();
		for (PBLink pb : links) {
			String id = pb.getOrigin().getID() + ":" + pb.getTerminus().getID();
			if (temp.containsKey(id)) {
				temp.get(id).add(pb);
			} else {
				List<PBLink> tLinks = new Vector<PBLink>();
				tLinks.add(pb);
				temp.put(id, tLinks);
			}
		}
		// build the edge for fitted for criterion;
		int count = 0;
		for (String s : temp.keySet()) {
			// if (temp.get(s).size() >= EdgeBundler.MIN_LINK_NUM) {
			if (temp.get(s).size() >= Integer.valueOf(paras.get("MIN_LINK_NUM"))) {
				count += 1;
				// statistical analysis contig pairs distance
				List<PBLink> tlinks = temp.get(s);
				List<Integer> dists = new Vector<Integer>();
				for (PBLink pb : tlinks) {
					dists.add(pb.getDistance());
				}
				int mean = 0;
				int sd = 0;
				try {
					mean = MathTool.mean(dists);
					sd = MathTool.sd(dists);
				} catch (Exception e) {
					logger.debug(this.getClass().getName() + e.getMessage() + "\t" + e.getClass().getName());
					logger.error(this.getClass().getName() + e.getMessage() + "\t" + e.getClass().getName());
					logger.info(this.getClass().getName() + e.getMessage() + "\t" + e.getClass().getName());
					continue;
				}
				int upper = mean + 2 * sd;
				int low = mean - 2 * sd;
				// if the distance larger than mean + 2 * sd, then remove;
				try {
					// This will show iteraction error when remove element
					// for (PBLink pb : tlinks) {
					// if (pb.getDistance() > upper || pb.getDistance() < low)
					// tlinks.remove(pb);
					// }
					for (int i = 0; i < tlinks.size(); i++) {
						PBLink pb = tlinks.get(i);
						if (pb.getDistance() > upper || pb.getDistance() < low)
							tlinks.remove(pb);
					}
				} catch (Exception e) {
					logger.debug(String.valueOf(count));
					logger.debug(this.getClass().getName() + e.getMessage() + "\t" + e.getClass().getName());
					logger.error(this.getClass().getName() + e.getMessage() + "\t" + e.getClass().getName());
					logger.info(this.getClass().getName() + e.getMessage() + "\t" + e.getClass().getName());
					continue;
				}
				// checking the most frequencies contig pair type
				int typeA = 0; // + +;
				int typeB = 0; // + -;
				int typeC = 0; // - -;
				int typeD = 0; // - +;
				for (PBLink pb : tlinks) {
					Strand oStrand = pb.getoStrand();
					Strand tStrand = pb.gettStrand();
					if (oStrand.equals(Strand.FORWARD)) {
						if (tStrand.equals(Strand.FORWARD)) {
							typeA += 1;
						} else {
							typeB += 1;
						}
					} else {
						if (tStrand.equals(Strand.REVERSE)) {
							typeC += 1;
						} else {
							typeD += 1;
						}
					}
				}
				int max = MathTool.max(typeA, typeB, typeC, typeD);

				Edge edge = new Edge();
				edge.setOrigin(temp.get(s).get(0).getOrigin());
				edge.setTerminus(temp.get(s).get(0).getTerminus());
				if (max == typeA) {
					edge.setoStrand(Strand.FORWARD);
					edge.settStrand(Strand.FORWARD);
					edge.setLinkNum(typeA);
				} else if (max == typeB) {
					edge.setoStrand(Strand.FORWARD);
					edge.settStrand(Strand.REVERSE);
					edge.setLinkNum(typeB);
				} else if (max == typeC) {
					edge.setoStrand(Strand.REVERSE);
					edge.settStrand(Strand.REVERSE);
					edge.setLinkNum(typeC);
				} else if (max == typeD) {
					edge.setoStrand(Strand.REVERSE);
					edge.settStrand(Strand.FORWARD);
					edge.setLinkNum(typeD);
				}
				edges.add(edge);
			}
		}

		// merge two opposite direction edge to only one edge base
		// not implemented now
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
}

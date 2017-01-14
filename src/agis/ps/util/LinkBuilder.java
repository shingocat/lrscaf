/*
*File: agis.ps.util.LinkBuilder.java
*User: mqin
*Email: mqin@ymail.com
*Date: 2016年1月12日
*/
package agis.ps.util;

import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Vector;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import agis.ps.link.MRecord;
import agis.ps.link.PBLink;
import agis.ps.link.TriadLink;
import agis.ps.seqs.Contig;

public class LinkBuilder {
	private static Logger logger = LoggerFactory.getLogger(LinkBuilder.class);
	private LinkedList<TriadLink> triadlinks = new LinkedList<TriadLink>();
	private List<PBLink> links = null;
	// private LinkedList<String> simcnts = new LinkedList<String>();
	private Parameter paras;
	private int minOLLen;
	private double minOLRatio;
	private int maxOHLen;
	private double maxOHRatio;
	private int maxEndLen;
	private double maxEndRatio;
	private double olRatio = 0.5;
	private double olweight = 0.6;
	private double identweight = 0.4;

	public LinkBuilder(Parameter paras) {
		this.paras = paras;
		minOLLen = paras.getMinOLLen();
		minOLRatio = paras.getMinOLRatio();
		maxOHLen = paras.getMaxOHLen();
		maxOHRatio = paras.getMaxOHRatio();
		maxEndLen = paras.getMaxEndLen();
		maxEndRatio = paras.getMaxEndRatio();
	}

	public List<PBLink> mRecords2Links(Map<String, List<MRecord>> records, List<String> repeats) {
		long start = System.currentTimeMillis();
		if (links == null)
			links = new Vector<PBLink>();
		links.clear();
		int bug = 0;
		try {
			for (String s : records.keySet()) {
//				if(s.equals("PB3577"))
//					logger.debug("breakpoint");
//				if(bug == 346)
//					logger.debug("breakpoint");
				List<MRecord> rs = records.get(s);
				List<PBLink> temp = this.mRecord2PBLink2(rs, repeats);
				if (temp != null)
					links.addAll(temp);
				bug++;
			}
		} catch (Exception e) {
			logger.error(bug + this.getClass().getName() + "\t" + e.getMessage());
		}
		long end = System.currentTimeMillis();
		logger.info("Valid Links Acount: " + links.size());
		logger.info("Building Link, erase time : " + (end - start) + " ms");
		return links;
	}
	
	private List<PBLink> mRecord2PBLink2(List<MRecord> records, List<String> repeats)
	{
		LinkedList<PBLink> links = new LinkedList<PBLink>();
		LinkedList<MRecord> valids = new LinkedList<MRecord>();
		int size = valids.size();
		MRecord former = null;
		MRecord current = null;
		Iterator<MRecord> it = records.iterator();
		Collections.sort(records, new ByLocOrderComparator());
		int index = 0;
		boolean isRemove = false;
		while(it.hasNext())
		{
			current = it.next();
			String currentCNTId = current.gettName();
			// check whether is repeats;
			if(repeats.contains(currentCNTId))
				continue;
			String currentPBId = current.getqName();
			int currentPBLen = current.getqLength();
			int currentPBStart = current.getqStart();
			int currentPBEnd = current.getqEnd() - 1;
			int currentCNTLen = current.gettLength();
			int currentCNTStart = current.gettStart();
			int currentCNTEnd = current.gettEnd() - 1;
			Strand currentStrand = current.gettStrand();
			boolean isInner = false;
			int endLen = maxEndLen;
			int endDefLen = (int) Math.round(currentPBLen * maxEndRatio);
			// if max end length larger than the endDefLen, used the endDefLen
			// as threshold
			if (endDefLen < maxEndLen)
				endLen = endDefLen;
			if (currentPBStart >= endLen && currentPBStart <= (currentPBLen - endLen)) {
				if (currentPBEnd <= (currentPBLen - endLen))
					isInner = true;
			}
			// the overlap length of contig
			int ol_len = currentCNTEnd - currentCNTStart + 1;
			double ratio = (double) ol_len / currentCNTLen;
			// the overhang length
			int contLeftLen = currentCNTStart;
			int contRigthLen = currentCNTLen - currentCNTEnd - 1;
			// for the reversed strand, BLASR define the coordinate according to
			// aligned seq
			if (currentStrand.equals(Strand.REVERSE)) {
				int temp = contLeftLen;
				contLeftLen = contRigthLen;
				contRigthLen = temp;
			}
			int ohLen = maxOHLen;
			int ohDefLen = (int) Math.round(currentCNTLen * maxOHRatio);
			if (ohDefLen < maxOHLen)
				ohLen = ohDefLen;
			if (isInner) {
				// the overlap length less than specified value, next;
				if (ol_len < minOLLen)
					continue;
				// if the overlap length enough for specified value, the ratio
				// is less than specified value, also next;
				if (ratio < minOLRatio)
					continue;
				// for two end length of contig not allow larger than maxOHLen
				if (contLeftLen > ohLen || contRigthLen > ohLen)
					continue;
			} else {
				if (currentPBStart <= endLen && currentPBEnd <= endLen) {
					// checking the right side, for the end point in p1-p2
					if (contRigthLen > ohLen)
						continue;
				} else if (currentPBStart <= endLen && currentPBEnd <= (currentPBLen - endLen)) {
					if (ol_len < minOLLen)
						continue;
					if (contRigthLen > ohLen)
						continue;
				} else if (currentPBStart <= endLen && currentPBEnd >= (currentPBLen - endLen)) {
					// start point in p1-p2 and end point in p3-p4
					// do no afford info, discard
					continue;
				} else if (currentPBStart >= endLen && currentPBStart <= (currentPBLen - endLen) && currentPBEnd >= (currentPBLen - endLen)) {
					// start point in p2-p3 and end point in p3-p4
					if (ol_len < minOLLen)
						continue;
					if (contLeftLen > ohLen)
						continue;
				} else if (currentPBStart >= (currentPBLen - endLen) && currentPBEnd >= (currentPBLen - endLen)) {
					if (contLeftLen > ohLen)
						continue;
				}
			}
			
			// building PBLink
			if(former == null)
			{
				valids.add(current);
				former = current;
				continue;
			} else
			{
				int dist = this.getDistance(former, current);
				int formerPBStart = former.getqStart();
				int formerPBEnd = former.getqEnd();
				if(dist >= 0)
				{
					PBLink link = new PBLink();
					link.setPbId(currentPBId);
					link.setDist(dist);
					link.setOrigin(former.gettName());
					link.setOStrand(former.gettStrand());
					link.setOPStart(formerPBStart);
					link.setOPEnd(formerPBEnd);
					link.setTerminus(currentCNTId);
					link.setTStrand(currentStrand);
					link.setTPStart(currentPBStart);
					link.setTPEnd(currentPBEnd);
					links.add(link);
					valids.add(current);
					former = current;
				} else
				{
					int formerPBOLLen = formerPBEnd - formerPBStart;
					int currentPBOLLen = currentPBEnd - currentPBStart;
					double formerRatio  = (double)(Math.abs(dist))/formerPBOLLen;
					double currentRatio = (double)(Math.abs(dist))/currentPBOLLen;
					if(formerRatio >= olRatio || currentRatio >= olRatio)
					{
						double formerIdentity = former.getIdentity();
						double currentIdentity = current.getIdentity();
						int formerScore = (int)Math.round(formerPBOLLen * olweight + 
								formerIdentity * 1000 * identweight);
						int currentScore = (int)Math.round(currentPBOLLen * olweight +
								currentIdentity * 1000 * identweight);
						if(formerScore >= currentScore)
						{
							continue;
						} else
						{
							isRemove = true;
							if(!links.isEmpty())
								links.removeLast(); // delete the last link
							if(!valids.isEmpty())
								valids.removeLast(); // delete the former from valids;
							if(valids.isEmpty())
							{
								valids.add(current);
								former = current;
								continue;
							} else
							{
								former = valids.getLast();
								dist = this.getDistance(former, current);
								PBLink link = new PBLink();
								link.setPbId(currentPBId);
								link.setDist(dist);
								link.setOrigin(former.gettName());
								link.setOStrand(former.gettStrand());
								link.setOPStart(formerPBStart);
								link.setOPEnd(formerPBEnd);
								link.setTerminus(currentCNTId);
								link.setTStrand(currentStrand);
								link.setTPStart(currentPBStart);
								link.setTPEnd(currentPBEnd);
								links.add(link);
								valids.add(current);
								former = current;
							}
						}
					} else
					{
						PBLink link = new PBLink();
						link.setPbId(currentPBId);
						link.setDist(dist);
						link.setOrigin(former.gettName());
						link.setOStrand(former.gettStrand());
						link.setOPStart(formerPBStart);
						link.setOPEnd(formerPBEnd);
						link.setTerminus(currentCNTId);
						link.setTStrand(currentStrand);
						link.setTPStart(currentPBStart);
						link.setTPEnd(currentPBEnd);
						links.add(link);
						valids.add(current);
						former = current;
					}
				}
			}
			
			// building triadlink
			size = valids.size();
			if(size >= 3)
			{
				// delete the corresponding triadlink;
				if(isRemove)
				{ 
					for(int i = 0; i <= index; i++)
						if(!triadlinks.isEmpty())
							triadlinks.removeLast();
					index--;
					isRemove = false;
				}
				for(int i = 0; i <= index; i++)
				{
					TriadLink tl = new TriadLink();
					Contig pre = new Contig();
					pre.setID(valids.get(i).gettName());
					tl.setPrevious(pre);
					Contig lst = new Contig();
					lst.setID(valids.get(size - 1).gettName());
					tl.setLast(lst);
					if((i - index) == 0)
					{
						Contig mid = new Contig();
						mid.setID(valids.get(size - 2).gettName());
						tl.setMiddle(mid);
					}
					tl.setSupLinks(1);
					if(triadlinks.contains(tl))
					{
						int loc = triadlinks.indexOf(tl);
						int supLink = triadlinks.get(loc).getSupLinks();
						supLink += 1;
						triadlinks.get(loc).setSupLinks(supLink);
					} else
					{
						triadlinks.add(tl);
					}
				}
				index++;
			}
		}
		if(links.isEmpty())
			return null;
		return links;
	}
	

	private List<PBLink> mRecord2PBLink(List<MRecord> records, List<String> repeats) {
		int size = records.size();
		if (records == null || size <= 1)
			return null;
		if (size >= 3)
			this.pesudoTriadLink(records);
		Iterator<MRecord> it = records.iterator();
		List<PBLink> pbLinks = new Vector<PBLink>();
		List<MRecord> valids = new Vector<MRecord>(records.size());
		while (it.hasNext()) {
			MRecord m = it.next();
			int pbLen = m.getqLength();
			int pbStart = m.getqStart();
			int pbEnd = m.getqEnd() - 1;
			int contLen = m.gettLength();
			int contStart = m.gettStart();
			int contEnd = m.gettEnd() - 1;
			Strand tStrand = m.gettStrand();
			boolean isInner = false;
			int endLen = maxEndLen;
			int endDefLen = (int) Math.round(pbLen * maxEndRatio);
			// if max end length larger than the endDefLen, used the endDefLen
			// as threshold
			if (endDefLen < maxEndLen)
				endLen = endDefLen;
			if (pbStart >= endLen && pbStart <= (pbLen - endLen)) {
				if (pbEnd <= (pbLen - endLen))
					isInner = true;
			}
			// the overlap length of contig
			int ol_len = contEnd - contStart + 1;
			double ratio = (double) ol_len / contLen;
			// the overhang length
			int contLeftLen = contStart;
			int contRigthLen = contLen - contEnd - 1;
			// for the reversed strand, BLASR define the coordinate according to
			// aligned seq
			if (tStrand.equals(Strand.REVERSE)) {
				int temp = contLeftLen;
				contLeftLen = contRigthLen;
				contRigthLen = temp;
			}
			int ohLen = maxOHLen;
			int ohDefLen = (int) Math.round(contLen * maxOHRatio);
			if (ohDefLen < maxOHLen)
				ohLen = ohDefLen;
			if (isInner) {
				// the overlap length less than specified value, next;
				if (ol_len < minOLLen)
					continue;
				// if the overlap length enough for specified value, the ratio
				// is less than specified value, also next;
				if (ratio < minOLRatio)
					continue;
				// for two end length of contig not allow larger than maxOHLen
				if (contLeftLen > ohLen || contRigthLen > ohLen)
					continue;
			} else {
				if (pbStart <= endLen && pbEnd <= endLen) {
					// checking the right side, for the end point in p1-p2
					if (contRigthLen > ohLen)
						continue;
				} else if (pbStart <= endLen && pbEnd <= (pbLen - endLen)) {
					if (ol_len < minOLLen)
						continue;
					if (contRigthLen > ohLen)
						continue;
				} else if (pbStart <= endLen && pbEnd >= (pbLen - endLen)) {
					// start point in p1-p2 and end point in p3-p4
					// do no afford info, discard
					continue;
				} else if (pbStart >= endLen && pbStart <= (pbLen - endLen) && pbEnd >= (pbLen - endLen)) {
					// start point in p2-p3 and end point in p3-p4
					if (ol_len < minOLLen)
						continue;
					if (contLeftLen > ohLen)
						continue;
				} else if (pbStart >= (pbLen - endLen) && pbEnd >= (pbLen - endLen)) {
					if (contLeftLen > ohLen)
						continue;
				}
			}
			// if the mrecord is valid;
			valids.add(m);
		}
		
		size = valids.size();
		if (size < 2)
			return null;

		// deleting repeats;
		if (paras.isRepMask() && !repeats.isEmpty()) {
			Iterator<MRecord> iterator = valids.iterator();
			while (iterator.hasNext()) {
				MRecord m = iterator.next();
				String cntId = m.gettName();
				if (repeats.contains(cntId))
					iterator.remove();
			}
		}

		// sorting the contig_pairs;
		Collections.sort(valids, new ByLocOrderComparator());
		// checking the similarity contigs
		// valids = validateContigPair(valids);
//		valids = this.validateSimilarityContigs(valids);
//		valids = validateOverlapContigs(valids);
		valids = this.deleteOverlaps(valids);
		size = valids.size();

		if (size < 2)
			return null;

		// build only the successive link; A->B->C, it will build A->B
		// and B->C, omitted A->C
		for (int i = 0; i <= size - 2; i++) {
			MRecord m1 = valids.get(i);
			MRecord m2 = valids.get(i + 1);

			if (m1.gettName().equalsIgnoreCase(m2.gettName())) {
				m1 = null;
				m2 = null;
				continue;
			}
			// // do not considering contain case
			String pId = m1.getqName();
			String oId = m1.gettName();
			int oPBS = m1.getqStart(); // origin pacbio start point;
//			int oCntS = m1.gettStart(); // origin contig start point;
			int oPBE = m1.getqEnd(); // origin pacbio end point;
//			int oCntE = m1.gettEnd(); // origin contig end point;
//			int oCntLen = m1.gettLength();// origin contig length;
			Strand oStrand = m1.gettStrand();
			String tId = m2.gettName();
			int tPBS = m2.getqStart(); // terminus pacbio start point;
//			int tCntS = m2.gettStart(); // terminus contig start point;
			int tPBE = m2.getqEnd(); // terminus pacbio end point;
//			int tCntE = m2.gettEnd(); // terminus contig end point;
//			int tCntLen = m2.gettLength(); // terminus contig length
			Strand tStrand = m2.gettStrand();
			int dist = this.getDistance(m1, m2);

			PBLink p = new PBLink();
			p.setPbId(pId);
			p.setOrigin(oId);
			p.setOPStart(oPBS);
			p.setOPEnd(oPBE);
			p.setOStrand(oStrand);
			p.setTerminus(tId);
			p.setTPStart(tPBS);
			p.setTPEnd(tPBE);
			p.setTStrand(tStrand);
			p.setDist(dist);
			pbLinks.add(p);
			p = null;
			m1 = null;
			m2 = null;
		}
		// build TriadLink
		if (size >= 3) {
			// List<TriadLink> triad = new Vector<TriadLink>();
			for (int i = 0; i <= size - 3; i++) {
				MRecord m1 = valids.get(i);
				MRecord m2 = valids.get(i + 1);
				MRecord m3 = valids.get(i + 2);
				Contig pre = new Contig();
				pre.setID(m1.gettName());
				// pre.setLength(m1.gettLength());
				Contig mid = new Contig();
				mid.setID(m2.gettName());
				// mid.setLength(m2.gettLength());
				Contig lst = new Contig();
				lst.setID(m3.gettName());
				// lst.setLength(m3.gettLength());
				TriadLink tl = new TriadLink(pre, mid, lst);
				tl.setSupLinks(1);
				if (triadlinks.contains(tl)) {
					int index = triadlinks.indexOf(tl);
					int supLink = triadlinks.get(index).getSupLinks();
					supLink += 1;
					triadlinks.get(index).setSupLinks(supLink);
					m1 = null;
					m2 = null;
					m3 = null;
					pre = null;
					mid = null;
					lst = null;
					tl = null;
				} else {
					triadlinks.add(tl);
					m1 = null;
					m2 = null;
					m3 = null;
					pre = null;
					mid = null;
					lst = null;
					tl = null;
				}
			}
		}
		// build crossing links
		if (size >= 4) {
			MRecord m1 = valids.get(0);
			MRecord m2 = valids.get(size - 1);
			Contig first = new Contig();
			first.setID(m1.gettName());
			Contig last = new Contig();
			last.setID(m2.gettName());
			TriadLink tl = new TriadLink();
			tl.setPrevious(first);
			tl.setLast(last);
			tl.setSupLinks(1);
			if (triadlinks.contains(tl)) {
				int index = triadlinks.indexOf(tl);
				int supLink = triadlinks.get(index).getSupLinks();
				supLink += 1;
				triadlinks.get(index).setSupLinks(supLink);
				m1 = null;
				m2 = null;
				first = null;
				last = null;
				tl = null;
			} else {
				triadlinks.add(tl);
				m1 = null;
				m2 = null;
				first = null;
				last = null;
				tl = null;
			}
		}
		return pbLinks;
	}
	
	private List<MRecord> deleteOverlaps(List<MRecord> data)
	{
//		List<MRecord> valids = new LinkedList<MRecord>();
		int size = data.size();
		List<MRecord> removes = new Vector<MRecord>(size);
		Iterator<MRecord> it = data.iterator();
		MRecord current = null;
		MRecord next = null;
		boolean isFirst = true;
		double olRatio = 0.5; // larger than this value, it will remove;
//		int constant = 100; // 100 bp to defined whether contigs is similarity
		double olweight = 0.6;
		double identweight = 0.4;
		while(it.hasNext())
		{
			if(isFirst)
			{
				current = it.next();
				isFirst = false;
			} else
			{
				next = it.next();
				int dist = this.getDistance(current, next);
				if(dist >= 0)
				{
					current = next;
					continue;
				}
				// it the distance is large than one of the cnt length, remove it;
//				int currentCNTLen = current.gettLength();
//				int nextCNTLen = next.gettLength();
				int currentPS = current.getqStart();
				int currentPE = current.getqEnd();
				int currentPBOLLen = currentPE - currentPS;
				int nextPS = next.getqStart();
				int nextPE = next.getqEnd();
				int nextPBOLLen = nextPE - nextPS;
//				int ol = currentPE - currentPS;
				double currentRatio = (double)(Math.abs(dist))/currentPBOLLen;
				double nextRatio = (double)(Math.abs(dist))/nextPBOLLen;
				if(currentRatio >= olRatio || nextRatio >= olRatio)
				{
					double currentIdentity = current.getIdentity();
					double nextIdentity = next.getIdentity();
					int currentScore = (int) Math.round(currentPBOLLen * olweight + currentIdentity * 1000 * identweight);
					int nextScore = (int) Math.round(nextPBOLLen * olweight + nextIdentity * 1000 * identweight);
					if(currentScore >= nextScore)
					{
						removes.add(next);
						continue;
					} else
					{
						removes.add(current);
						current = next;
						continue;
					}
				}
				current = next;
			} //  end if(isFirst)
		} // end while(it.hasNext)
		if(!removes.isEmpty())
			data.removeAll(removes);
		return data;
	}

	// public List<PBLinkM> mRecord2Link(List<MRecord> records, List<String>
	// repeats)
	// {
	// if(records.size() >= 3)
	// this.pesudoTriadLink(records);
	// Iterator<MRecord> it = records.iterator();
	// List<PBLinkM> pbLinks = new Vector<PBLinkM>();
	// List<MRecord> valids = new Vector<MRecord>(records.size());
	// while (it.hasNext()) {
	// MRecord m = it.next();
	// int minOLLen = paras.getMinOLLen();
	// double minOLRatio = paras.getMinOLRatio();
	// int maxOHLen = paras.getMaxOHLen();
	// double maxOHRatio = paras.getMaxOHRatio();
	// int maxEndLen = paras.getMaxEndLen();
	// double maxEndRatio = paras.getMaxEndRatio();
	//
	// int pbLen = m.getqLength();
	// int pbStart = m.getqStart();
	// int pbEnd = m.getqEnd() - 1;
	// int contLen = m.gettLength();
	// int contStart = m.gettStart();
	// int contEnd = m.gettEnd() - 1;
	// Strand tStrand = m.gettStrand();
	//
	// boolean isInner = false;
	// int endDefLen = (int) Math.round(pbLen * maxEndRatio);
	// // if max end length larger than the endDefLen, used the endDefLen
	// // as threshold
	// if (endDefLen < maxEndLen)
	// maxEndLen = endDefLen;
	// if (pbStart >= maxEndLen && pbStart <= (pbLen - maxEndLen)) {
	// if (pbEnd <= (pbLen - maxEndLen))
	// isInner = true;
	// }
	// // the overlap length of contig
	// int ol_len = contEnd - contStart + 1;
	// double ratio = (double) ol_len / contLen;
	// // the overhang length
	// int contLeftLen = contStart;
	// int contRigthLen = contLen - contEnd - 1;
	// // for the reversed strand, BLASR define the coordinate according to
	// // aligned seq
	// if (tStrand.equals(Strand.REVERSE)) {
	// int temp = contLeftLen;
	// contLeftLen = contRigthLen;
	// contRigthLen = temp;
	// }
	// int ohDefLen = (int) Math.round(contLen * maxOHRatio);
	// if (ohDefLen < maxOHLen)
	// maxOHLen = ohDefLen;
	// if (isInner) {
	// // the overlap length less than specified value, next;
	// if (ol_len < minOLLen)
	// continue;
	// // if the overlap length enough for specified value, the ratio
	// // is less than specified value, also next;
	// if (ratio < minOLRatio)
	// continue;
	// // for two end length of contig not allow larger than maxOHLen
	// if (contLeftLen > maxOHLen || contRigthLen > maxOHLen)
	// continue;
	// } else {
	// if(pbStart <= maxEndLen && pbEnd <= maxEndLen)
	// { // checking the right side, for the end point in p1-p2
	// if(contRigthLen > maxOHLen)
	// continue;
	// } else if(pbStart <= maxEndLen && pbEnd <= (pbLen - maxEndLen))
	// {
	// if(ol_len < minOLLen)
	// continue;
	// if(contRigthLen > maxOHLen)
	// continue;
	// } else if(pbStart <= maxEndLen && pbEnd >= (pbLen - maxEndLen))
	// {
	// // start point in p1-p2 and end point in p3-p4
	// // do no afford info, discard
	// continue;
	// } else if(pbStart >= maxEndLen && pbStart <= (pbLen - maxEndLen)
	// && pbEnd >= (pbLen - maxEndLen))
	// { // start point in p2-p3 and end point in p3-p4
	// if(ol_len < minOLLen)
	// continue;
	// if(contLeftLen > maxOHLen)
	// continue;
	// } else if(pbStart >= (pbLen - maxEndLen) && pbEnd >= (pbLen - maxEndLen))
	// {
	// if(contLeftLen > maxOHLen)
	// continue;
	// }
	// }
	// // if the mrecord is valid;
	// valids.add(m);
	// }
	//
	// // deleting repeats;
	// if(paras.isRepMask())
	// {
	// Iterator<MRecord> iterator = valids.iterator();
	// while(iterator.hasNext())
	// {
	// MRecord m = iterator.next();
	// String cntId = m.gettName();
	// if(repeats.contains(cntId))
	// iterator.remove();
	// }
	// }
	//
	// // iterator valid mrecords
	// if (valids.size() < 2) {
	// return null;
	// } else {
	// // sorting the contig_pairs;
	// Collections.sort(valids, new ByLocOrderComparator());
	// // checking the similarity contigs
	//// valids = validateContigPair(valids);
	// valids = this.validateSimilarityContigs(valids);
	// valids = validateOverlapContigs(valids);
	// int cpSize = valids.size();
	//
	// if(cpSize < 2)
	// return null;
	//
	// // build only the successive link; A->B->C, it will build A->B
	// // and B->C, omitted A->C
	// for (int i = 0; i <= cpSize - 2; i++) {
	// MRecord m1 = valids.get(i);
	// MRecord m2 = valids.get(i + 1);
	//
	// if (m1.gettName().equalsIgnoreCase(m2.gettName())) {
	// m1 = null;
	// m2 = null;
	// continue;
	// }
	//// // do not considering contain case
	//// int m1PS = m1.getqStart();
	//// int m1PE = m1.getqEnd();
	//// int m2PS = m2.getqStart();
	//// int m2PE = m2.getqEnd();
	//// if(m1PS <= m2PS && m1PE >= m2PE)
	//// {
	//// m1 = null;
	//// m2 = null;
	//// continue;
	//// }
	// PBLinkM p = new PBLinkM();
	// p.setOrigin(m1);
	// p.setTerminus(m2);
	// p.setId(m1.getqName());
	//// int distance = p.getDistance();
	//// if(distance < 0)
	//// {
	//// int olength = m1.gettLength();
	//// int tlength = m2.gettLength();
	//// distance = Math.abs(distance);
	//// if(distance >= olength || distance >= tlength)
	//// continue;
	//// }
	// pbLinks.add(p);
	// p = null;
	// m1 = null;
	// m2 = null;
	// }
	// // build TriadLink
	// if (cpSize >= 3) {
	//// List<TriadLink> triad = new Vector<TriadLink>();
	// for (int i = 0; i <= cpSize - 3; i++) {
	// MRecord m1 = valids.get(i);
	// MRecord m2 = valids.get(i + 1);
	// MRecord m3 = valids.get(i + 2);
	// Contig pre = new Contig();
	// pre.setID(m1.gettName());
	// // pre.setLength(m1.gettLength());
	// Contig mid = new Contig();
	// mid.setID(m2.gettName());
	// // mid.setLength(m2.gettLength());
	// Contig lst = new Contig();
	// lst.setID(m3.gettName());
	// // lst.setLength(m3.gettLength());
	// TriadLink tl = new TriadLink(pre, mid, lst);
	// tl.setSupLinks(1);
	// if (tls.contains(tl)) {
	// int index = tls.indexOf(tl);
	// int supLink = tls.get(index).getSupLinks();
	// supLink += 1;
	// tls.get(index).setSupLinks(supLink);
	// m1 = null;
	// m2 = null;
	// m3 = null;
	// pre = null;
	// mid = null;
	// lst = null;
	// tl = null;
	// } else {
	// tls.add(tl);
	// m1 = null;
	// m2 = null;
	// m3 = null;
	// pre = null;
	// mid = null;
	// lst = null;
	// tl = null;
	// }
	// }
	//// if(triads.size() > 0)
	//// this.tlWriter.write2(triads);
	// }
	// // build crossing links
	// if(cpSize >= 4)
	// {
	// MRecord m1 = valids.get(0);
	// MRecord m2 = valids.get(cpSize - 1);
	// Contig first = new Contig();
	// first.setID(m1.gettName());
	// Contig last = new Contig();
	// last.setID(m2.gettName());
	// TriadLink tl = new TriadLink();
	// tl.setPrevious(first);
	// tl.setLast(last);
	// tl.setSupLinks(1);
	// if (tls.contains(tl)) {
	// int index = tls.indexOf(tl);
	// int supLink = tls.get(index).getSupLinks();
	// supLink += 1;
	// tls.get(index).setSupLinks(supLink);
	// m1 = null;
	// m2 = null;
	// first = null;
	// last = null;
	// tl = null;
	// } else {
	// tls.add(tl);
	// m1 = null;
	// m2 = null;
	// first = null;
	// last = null;
	// tl = null;
	// }
	// }
	//
	// return pbLinks;
	// }
	// }
	//
//	private List<MRecord> validateSimilarityContigs(List<MRecord> data) {
//		Map<String, Vector<MRecord>> sims = this.findSimilarityCnts(data);
//		int size = sims.size();
//		if (size != 0) {
//			for (int i = 1; i <= size; i++) {
//				List<MRecord> ms = sims.get(String.valueOf(i));
//				// data.remove(ms); // try remove all similarity contigs
//				List<Integer> scores = new Vector<Integer>(ms.size());
//				Iterator<MRecord> it = ms.iterator();
//				double olweight = 0.6;
//				double identweight = 0.4;
//				while (it.hasNext()) {
//					MRecord m = it.next();
//					int pS = m.getqStart();
//					int pE = m.getqEnd();
//					int ol = pE - pS;
//					int ident = (int) Math.round(m.getIdentity() * 1000);
//					int score = (int) Math.round(ol * olweight + ident * identweight);
//					scores.add(score);
//				}
//				int max = MathTool.max(scores);
//				it = ms.iterator();
//				while (it.hasNext()) {
//					MRecord m = it.next();
//					int pS = m.getqStart();
//					int pE = m.getqEnd();
//					int ol = pE - pS;
//					int ident = (int) Math.round(m.getIdentity() * 1000);
//					int score = (int) Math.round(ol * olweight + ident * identweight);
//					if (score != max) {
//						if (data.contains(m))
//							data.remove(m);
//					}
//				}
//			}
//		}
//		return data;
//	}

//	private List<MRecord> validateOverlapContigs(List<MRecord> data) {
//		List<MRecord> removes = new Vector<MRecord>(data.size());
//		Iterator<MRecord> it = data.iterator();
//		boolean isFirst = true;
//		MRecord former = null;
//		MRecord current = null;
//		double olRatio = 0.9; // larger than this value, it will remove;
//		while (it.hasNext()) {
//			if (isFirst) {
//				former = it.next();
//				isFirst = false;
//				continue;
//			} else {
//				// if(former.getqName().equalsIgnoreCase("m130605_032054_42207_c100515142550000001823076608221373_s1_p0/147280/0_8554"))
//				// logger.debug("breakpoint");
//				current = it.next();
//				int dist = this.getDistance(former, current);
//				if (dist >= 0) {
//					former = current;
//					continue;
//				}
//
//				// if the distance between contigs is large than one of the
//				// length, then
//				// remove the overlap one;
//				int fCL = former.gettLength();
//				int cCL = current.gettLength();
//				if (Math.abs(dist) >= fCL) {
//					if (fCL < cCL) {
//						removes.add(former);
//						former = current;
//						continue;
//					} else {
//						removes.add(current);
//						continue;
//					}
//				}
//				if (Math.abs(dist) >= cCL) {
//					if (cCL < fCL) {
//						removes.add(current);
//						continue;
//					} else {
//						removes.add(former);
//						former = current;
//						continue;
//					}
//				}
//
//				// former code
//				int fPS = former.getqStart();
//				int fPE = former.getqEnd();
//				int fl = fPE - fPS;
//				int cPS = current.getqStart();
//				int cPE = current.getqEnd();
//				int cl = cPE - cPS;
//				int ol = fPE - cPS;
//				double fOLRatio = (double) ol / fl;
//				double cOLRatio = (double) ol / cl;
//				if (fOLRatio >= olRatio || cOLRatio >= olRatio) {
//					if (fOLRatio >= cOLRatio) {
//						removes.add(former);
//						former = current;
//						continue;
//					} else {
//						removes.add(current);
//						continue;
//					}
//				} else {
//					former = current;
//				}
//			}
//		}
//		data.removeAll(removes);
//		return data;
//	}

	private int getDistance(MRecord origin, MRecord terminus) {
		int dist = 0;
		// int oPBS = origin.getqStart(); // origin pacbio start point;
		int oCntS = origin.gettStart(); // origin contig start point;
		int oPBE = origin.getqEnd(); // origin pacbio end point;
		int oCntE = origin.gettEnd(); // origin contig end point;
		int tPBS = terminus.getqStart(); // terminus pacbio start point;
		int tCntS = terminus.gettStart(); // terminus contig start point;
		// int tPBE = terminus.getqEnd(); // terminus pacbio end point;
		int tCntE = terminus.gettEnd(); // terminus contig end point;
		int oCntLen = origin.gettLength();// origin contig length;
		int tCntLen = terminus.gettLength(); // terminus contig length
		if (origin.gettStrand().equals(Strand.FORWARD) && terminus.gettStrand().equals(Strand.FORWARD)) { // +
																											// +
			int oRightLen = oCntLen - oCntE;
			int tLeftLen = tCntS;
			dist = tPBS - oPBE - oRightLen - tLeftLen;
		} else if (origin.gettStrand().equals(Strand.REVERSE) && terminus.gettStrand().equals(Strand.REVERSE)) { // -
																													// -
			int oLeftLen = oCntS;
			int tRightLen = tCntLen - tCntE;
			dist = tPBS - oPBE - oLeftLen - tRightLen;
		} else if (origin.gettStrand().equals(Strand.FORWARD) && terminus.gettStrand().equals(Strand.REVERSE)) { // +
																													// -
			int oRightLen = oCntLen - oCntE;
			int tRightLen = tCntLen - tCntE;
			dist = tPBS - oPBE - oRightLen - tRightLen;
		} else if (origin.gettStrand().equals(Strand.REVERSE) && terminus.gettStrand().equals(Strand.FORWARD)) { // -
																													// +
			int oLeftLen = oCntS;
			int tLeftLen = tCntS;
			dist = tPBS - oPBE - oLeftLen - tLeftLen;
		}
		return dist;
	}

	// test in yeast for some pacbio read do not valid by the filter parameter
	// but it provide some link info!
	private void pesudoTriadLink(List<MRecord> ms) {
		if (triadlinks == null)
			triadlinks = new LinkedList<TriadLink>();
		Collections.sort(ms, new ByLocOrderComparator());
		int cpSize = ms.size();
		MRecord m1 = ms.get(0);
		MRecord m2 = ms.get(cpSize - 1);
		Contig first = new Contig();
		first.setID(m1.gettName());
		Contig last = new Contig();
		last.setID(m2.gettName());
		TriadLink tl = new TriadLink();
		tl.setPrevious(first);
		tl.setLast(last);
		tl.setSupLinks(1);
		if (triadlinks.contains(tl)) {
			int index = triadlinks.indexOf(tl);
			int supLink = triadlinks.get(index).getSupLinks();
			supLink += 1;
			triadlinks.get(index).setSupLinks(supLink);
			m1 = null;
			m2 = null;
			first = null;
			last = null;
			tl = null;
		} else {
			triadlinks.add(tl);
			m1 = null;
			m2 = null;
			first = null;
			last = null;
			tl = null;
		}
	}

//	private Map<String, Vector<MRecord>> findSimilarityCnts(List<MRecord> data) {
//		Map<String, Vector<MRecord>> sims = new HashMap<String, Vector<MRecord>>();
//		Vector<MRecord> tempSims = new Vector<MRecord>(5);
//		Iterator<MRecord> it = data.iterator();
//		int count = 1;
//		int preStart = 0;
//		int preEnd = 0;
//		int constant = 100; // 100 bp to defined whether contigs is similarity
//		boolean isFirst = true;
//		while (it.hasNext()) {
//			MRecord m = it.next();
//			if (isFirst) { // the first time
//				preStart = m.getqStart();
//				preEnd = m.getqEnd();
//				tempSims.add(m);
//				isFirst = false;
//				continue;
//			} else {
//				int curStart = m.getqStart();
//				int curEnd = m.getqEnd();
//
//				int left = curStart - preStart;
//				int right = curEnd - preEnd;
//				left = Math.abs(left);
//				right = Math.abs(right);
//
//				if (left <= constant && right <= constant) { // similarity
//																// contigs
//					if (!tempSims.contains(m))
//						tempSims.add(m);
//					preStart = curStart;
//					preEnd = curEnd;
//				} else {
//					preStart = curStart;
//					preEnd = curEnd;
//					if (tempSims.size() <= 1) {
//						tempSims.clear();
//						if (!tempSims.contains(m))
//							tempSims.add(m);
//					} else {
//						sims.put(String.valueOf(count), tempSims);
//						tempSims = new Vector<MRecord>(5);
//						count++;
//					}
//				}
//			}
//		}
//		if (!tempSims.isEmpty() && tempSims.size() >= 2) {
//			sims.put(String.valueOf(count), tempSims);
//			tempSims = null;
//		}
//		return sims;
//	}

	public List<TriadLink> getTriadLinks() {
		return triadlinks;
	}
}

class ByLocOrderComparator implements Comparator<Object> {
	@Override
	public int compare(Object cFirst, Object cSecond) {
		MRecord m1 = (MRecord) cFirst;
		MRecord m2 = (MRecord) cSecond;
		int diff = Integer.valueOf(m1.getqStart()) - Integer.valueOf(m2.getqStart());
		if (diff > 0)
			return 1;
		if (diff < 0)
			return -1;
		else
			return 0;
	}
}

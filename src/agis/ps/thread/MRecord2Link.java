/** 
** Usage: TODO
** Author: mqin
** Email: mqin@outlook.com
** Date: 2019年1月23日
*/
package agis.ps.thread;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Hashtable;
//import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Vector;
import java.util.concurrent.Callable;
//import java.util.concurrent.CountDownLatch;

import agis.ps.link.ByLocOrderComparator;
import agis.ps.link.ILink;
import agis.ps.link.Link;
import agis.ps.link.MRecord;
import agis.ps.link.TriadLink;
//import agis.ps.seqs.Contig;
//import agis.ps.seqs.Scaffold;
import agis.ps.seqs.Sequence;
import agis.ps.util.Parameter;
import agis.ps.util.Strand;

public class MRecord2Link implements Callable<Map<String, List<ILink>>> {
	private List<ILink> triadlinks = null;
	private List<ILink> links = null;
	private Parameter paras;
	private List<List<MRecord>> records;
	private List<String> repeats;
	private int minOLLen;
	private double minOLRatio;
	private int maxOHLen;
	private double maxOHRatio;
	private int maxEndLen;
	private double maxEndRatio;
	// default overlap 100 bp will consider which is the best.
	private int olLength = -1500;
	private double olRatio = 0.6;
	private double olweight = 0.6;
	private double identweight = 0.4;
	// private Map<String, Contig> cnts;
	private Map<String, Sequence> seqs;

	// public MRecord2Link(Parameter paras, List<List<MRecord>> records,
	// List<String> repeats, Map<String, Contig> cnts
	// ) {
	// this.paras = paras;
	// this.records = records;
	// this.repeats = repeats;
	// this.minOLLen = paras.getMinOLLen();
	// this.minOLRatio = paras.getMinOLRatio();
	// this.maxOHLen = paras.getMaxOHLen();
	// this.maxOHRatio = paras.getMaxOHRatio();
	// this.maxEndLen = paras.getMaxEndLen();
	// this.maxEndRatio = paras.getMaxEndRatio();
	// this.cnts = cnts;
	// this.links = new Vector<ILink>();
	// this.triadlinks = new Vector<ILink>();
	// }

	public MRecord2Link(Parameter paras, List<List<MRecord>> records, List<String> repeats,
			Map<String, Sequence> seqs) {
		this.paras = paras;
		this.records = records;
		this.repeats = repeats;
		this.seqs = seqs;
		this.minOLLen = paras.getMinOLLen();
		this.minOLRatio = paras.getMinOLRatio();
		this.maxOHLen = paras.getMaxOHLen();
		this.maxOHRatio = paras.getMaxOHRatio();
		this.maxEndLen = paras.getMaxEndLen();
		this.maxEndRatio = paras.getMaxEndRatio();
		this.links = new Vector<ILink>();
		this.triadlinks = new Vector<ILink>();
	}

	@Override
	public Map<String, List<ILink>> call() {
		int size = records.size();
		Map<String, List<ILink>> map = new Hashtable<String, List<ILink>>();
		for (int i = 0; i < size; i++) {
			this.mrecord2links(records.get(i));
		}
		map.put("LINKS", links);
		map.put("TRIADLINKS", triadlinks);
		return map;
	}

	private void mrecord2links(List<MRecord> records) {
		List<MRecord> valids = this.validMecordsIsRepeatOrOverhangAndOverlap(records);
		if (valids.size() <= 1)
			return;
		List<MRecord> ms = this.validBestMecords(valids);
		valids = null;
		if (ms.size() <= 1)
			return;
		// building links
		this.buildLinks(ms);
		// building triadlinks
		this.buildTriadLinks(ms);
	}

	private void buildLinks(List<MRecord> ms) {
		MRecord former = null;
		MRecord current = null;
		int size = ms.size();
		for (int i = 0; i < size - 1; i++) {
			former = ms.get(i);
			current = ms.get(i + 1);
			int dist = this.getDistance(former, current);
			Link link = new Link();
			link.setLrId(former.getqName());
			link.setDist(dist);
			link.setOriginal(this.seqs.get(former.gettName()));
			link.setoStrand(former.gettStrand());
			link.setoStart(former.getqStart());
			link.setoEnd(former.getqEnd());
			link.setTerminus(this.seqs.get(current.gettName()));
			link.settStrand(current.gettStrand());
			link.settStart(current.getqStart());
			link.settEnd(current.getqEnd());
			links.add(link);
		}
	}

	private void buildTriadLinks(List<MRecord> ms) {
		int size = ms.size();
		if (size >= 3) {
			int index = 0;
			for (int i = 0; i < size - 2; i++) {
				TriadLink tl = new TriadLink();
				MRecord pre = ms.get(i);
				MRecord mid = ms.get(i + 1);
				index = i + 2;
				while ((size - index) > 0) {
					MRecord lst = ms.get(index);
					Sequence preCnt = new Sequence();
					Sequence midCnt = new Sequence();
					Sequence lstCnt = new Sequence();
					preCnt.setId(pre.gettName());
					midCnt.setId(mid.gettName());
					lstCnt.setId(lst.gettName());
					tl.setPrevious(preCnt);
					tl.setMiddle(midCnt);
					tl.setLast(lstCnt);
					tl.setSupLinks(1);
					if (triadlinks.contains(tl)) {
						int loc = triadlinks.indexOf(tl);
						int supLink = ((TriadLink) triadlinks.get(loc)).getSupLinks();
						supLink += 1;
						((TriadLink) triadlinks.get(loc)).setSupLinks(supLink);
					} else {
						triadlinks.add(tl);
					}
					index++;
				}
			}
		}
	}

	private List<MRecord> validBestMecords(List<MRecord> records) {
		// remove overlap contigs large than threshold;
		// They are always similarity contigs;
		int size = records.size();
		Collections.sort(records, new ByLocOrderComparator());
		ArrayList<MRecord> valids = new ArrayList<MRecord>(records.size());
		MRecord former = records.get(0);
		MRecord current = null;
		String formerSeqId = null;
		String currentSeqId = null;
		// boolean isAdded = false;
		for (int i = 1; i < size; i++) {
			current = records.get(i);
			formerSeqId = former.gettName();
			currentSeqId = current.gettName();
			if (formerSeqId.equals(currentSeqId))
				continue;
			int dist = this.getDistance(former, current);
			if (dist >= 0) {
				valids.add(former);
				if (i == (size - 1))
					valids.add(current);
				else
					former = current;
			} else {
				int formerPBStart = former.getqStart();
				int formerPBEnd = former.getqEnd();
				int currentPBStart = current.getqStart();
				int currentPBEnd = current.getqEnd();
				// checking score;
				int formerPBOLLen = formerPBEnd - formerPBStart;
				int currentPBOLLen = currentPBEnd - currentPBStart;
				double formerRatio = (double) (Math.abs(dist)) / formerPBOLLen;
				double currentRatio = (double) (Math.abs(dist)) / currentPBOLLen;
				// should be more generalization
				if ((formerRatio >= olRatio || currentRatio >= olRatio) || (dist <= olLength)) {
					double formerIdentity = former.getIdentity();
					double currentIdentity = current.getIdentity();
					int formerScore = (int) Math.round(formerPBOLLen * olweight + formerIdentity * 1000 * identweight);
					int currentScore = (int) Math
							.round(currentPBOLLen * olweight + currentIdentity * 1000 * identweight);
					if (formerScore >= currentScore) {
						if (i == (size - 1))
							valids.add(former);
						continue;
					} else {
						if (i == (size - 1)) {
							valids.add(current);
							continue;
						} else {
							former = current;
						}
					}
				} else {
					valids.add(former);
					if (i == (size - 1))
						valids.add(current);
					else
						former = current;
				}
			}
		}
		return valids;
	}

	private List<MRecord> validMecordsIsRepeatOrOverhangAndOverlap(List<MRecord> records) {
		List<MRecord> valids = new ArrayList<MRecord>(records.size());
		// remove repeat;
		int size = records.size();
		for (int i = 0; i < size; i++) {
			MRecord record = records.get(i);
			String tId = record.gettName();
			if (this.paras.isRepMask()) {
				if (this.repeats.contains(tId)) {
					continue;
				}
			}
			// if(!repeats.contains(tId))
			int currentPBLen = record.getqLength();
			int currentPBStart = record.getqStart();
			int currentPBEnd = record.getqEnd() - 1;
			int currentCNTLen = record.gettLength();
			int currentCNTStart = record.gettStart();
			int currentCNTEnd = record.gettEnd() - 1;
			Strand currentStrand = record.gettStrand();
			boolean isInner = false;
			int endLen = this.maxEndLen;
			int endDefLen = (int) Math.round(currentPBLen * this.maxEndRatio);
			// if max end length larger than the endDefLen, used the endDefLen
			// as threshold
			if (endDefLen < this.maxEndLen)
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
			int ohLen = this.maxOHLen;
			int ohDefLen = (int) Math.round(currentCNTLen * this.maxOHRatio);
			if (ohDefLen < this.maxOHLen)
				ohLen = ohDefLen;
			if (isInner) {
				// the overlap length less than specified value, next;
				if (ol_len < this.minOLLen)
					continue;
				// if the overlap length enough for specified value, the ratio
				// is less than specified value, also next;
				if (ratio < this.minOLRatio)
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
					// checking ol length, for the end point in p2, p3;
					if (ol_len < this.minOLLen)
						continue;
					if (contRigthLen > ohLen)
						continue;
				} else if (currentPBStart <= endLen && currentPBEnd >= (currentPBLen - endLen)) {
					// start point in p1-p2 and end point in p3-p4
					// do no afford info, discard
					if (contLeftLen > ohLen && contRigthLen > ohLen)
						continue;
				} else if (currentPBStart >= endLen && currentPBStart <= (currentPBLen - endLen)
						&& currentPBEnd >= (currentPBLen - endLen)) {
					// start point in p2-p3 and end point in p3-p4
					if (ol_len < this.minOLLen)
						continue;
					if (contLeftLen > ohLen)
						continue;
				} else if (currentPBStart >= (currentPBLen - endLen) && currentPBEnd >= (currentPBLen - endLen)) {
					if (contLeftLen > ohLen)
						continue;
				}
			}
			valids.add(record);
		}
		return valids;
	}

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

		if (origin.gettStrand().equals(Strand.FORWARD) && terminus.gettStrand().equals(Strand.FORWARD)) {
			// + +
			int oRightLen = oCntLen - oCntE;
			int tLeftLen = tCntS;
			dist = tPBS - oPBE - oRightLen - tLeftLen;
		} else if (origin.gettStrand().equals(Strand.REVERSE) && terminus.gettStrand().equals(Strand.REVERSE)) {
			// - -
			int oLeftLen = oCntS;
			int tRightLen = tCntLen - tCntE;
			dist = tPBS - oPBE - oLeftLen - tRightLen;
		} else if (origin.gettStrand().equals(Strand.FORWARD) && terminus.gettStrand().equals(Strand.REVERSE)) {
			// + -
			int oRightLen = oCntLen - oCntE;
			int tRightLen = tCntLen - tCntE;
			dist = tPBS - oPBE - oRightLen - tRightLen;
		} else if (origin.gettStrand().equals(Strand.REVERSE) && terminus.gettStrand().equals(Strand.FORWARD)) {
			// - +
			int oLeftLen = oCntS;
			int tLeftLen = tCntS;
			dist = tPBS - oPBE - oLeftLen - tLeftLen;
		}
		return dist;
	}

}

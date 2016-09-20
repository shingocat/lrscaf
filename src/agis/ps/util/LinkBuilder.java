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
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Vector;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import agis.ps.file.TriadLinkWriter;
import agis.ps.link.MRecord;
import agis.ps.link.PBLinkM;
import agis.ps.link.TriadLink;
import agis.ps.seqs.Contig;

public class LinkBuilder {
	private static Logger logger = LoggerFactory.getLogger(LinkBuilder.class);
	private List<TriadLink> tls = new Vector<TriadLink>(100); 
	private List<MRecord> ms;
	private Parameter paras;
	
	public LinkBuilder(Parameter paras)
	{
		this.paras = paras;
	}

	public LinkBuilder(List<MRecord> ms, Parameter paras) {
		this.ms = ms;
		this.paras = paras;
	}

	public List<PBLinkM> mRecord2Link() {
		if (ms == null || ms.size() == 0)
			throw new IllegalArgumentException("LinkBuilder: The MRecords could not be empty!");
		return this.mRecord2Link(ms, paras);
	}
	
	public List<PBLinkM> mRecord2Link(List<MRecord> records, List<String> repeats)
	{
		if(records.size() >= 3)
			this.pesudoTriadLink(records);
		Iterator<MRecord> it = records.iterator();
		List<PBLinkM> pbLinks = new Vector<PBLinkM>();
		List<MRecord> valids = new Vector<MRecord>(records.size());
		while (it.hasNext()) {
			MRecord m = it.next();
			int minOLLen = paras.getMinOLLen();
			double minOLRatio = paras.getMinOLRatio();
			int maxOHLen = paras.getMaxOHLen();
			double maxOHRatio = paras.getMaxOHRatio();
			int maxEndLen = paras.getMaxEndLen();
			double maxEndRatio = paras.getMaxEndRatio();

			int pbLen = m.getqLength();
			int pbStart = m.getqStart();
			int pbEnd = m.getqEnd() - 1;
			int contLen = m.gettLength();
			int contStart = m.gettStart();
			int contEnd = m.gettEnd() - 1;
			Strand tStrand = m.gettStrand();

			boolean isInner = false;
			int endDefLen = (int) Math.round(pbLen * maxEndRatio);
			// if max end length larger than the endDefLen, used the endDefLen
			// as threshold
			if (endDefLen < maxEndLen)
				maxEndLen = endDefLen;
			if (pbStart >= maxEndLen && pbStart <= (pbLen - maxEndLen)) {
				if (pbEnd <= (pbLen - maxEndLen))
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
			int ohDefLen = (int) Math.round(contLen * maxOHRatio);
			if (ohDefLen < maxOHLen)
				maxOHLen = ohDefLen;
			if (isInner) {
				// the overlap length less than specified value, next;
				if (ol_len < minOLLen)
					continue;
				// if the overlap length enough for specified value, the ratio
				// is less than specified value, also next;
				if (ratio < minOLRatio)
					continue;
				// for two end length of contig not allow larger than maxOHLen
				if (contLeftLen > maxOHLen || contRigthLen > maxOHLen)
					continue;
			} else {
				if(pbStart <= maxEndLen && pbEnd <= maxEndLen)
				{ // checking the right side, for the end point in p1-p2
					if(contRigthLen > maxOHLen)
						continue;
				} else if(pbStart <= maxEndLen && pbEnd <= (pbLen - maxEndLen))
				{
					if(ol_len < minOLLen)
						continue;
					if(contRigthLen > maxOHLen)
						continue;
				} else if(pbStart <= maxEndLen && pbEnd >= (pbLen - maxEndLen))
				{
					// start point in p1-p2 and end point in p3-p4
					// do no afford info, discard
					continue;
				} else if(pbStart >= maxEndLen && pbStart <= (pbLen - maxEndLen) 
						&& pbEnd >= (pbLen - maxEndLen))
				{ // start point in p2-p3 and end point in p3-p4
					if(ol_len < minOLLen)
						continue;
					if(contLeftLen > maxOHLen)
						continue;
				} else if(pbStart >= (pbLen - maxEndLen) && pbEnd >= (pbLen - maxEndLen))
				{
					if(contLeftLen > maxOHLen)
						continue;
				}
			}
			// if the mrecord is valid;
			valids.add(m);
		}
		
		// deleting repeats;
		if(paras.isRepMask())
		{
			Iterator<MRecord> iterator = valids.iterator();
			while(iterator.hasNext())
			{
				MRecord m = iterator.next();
				String cntId = m.gettName();
				if(repeats.contains(cntId))
					iterator.remove();
			}
		}

		// iterator valid mrecords
		if (valids.size() < 2) {
			return null;
		} else {
			// sorting the contig_pairs;
			Collections.sort(valids, new ByLocOrderComparator());
			// checking the similarity contigs
//			valids = validateContigPair(valids);
			valids = this.validateSimilarityContigs(valids);
			valids = validateOverlapContigs(valids);
			int cpSize = valids.size();

			// build only the successive link; A->B->C, it will build A->B
			// and B->C, omitted A->C
			for (int i = 0; i <= cpSize - 2; i++) {
				MRecord m1 = valids.get(i);
				MRecord m2 = valids.get(i + 1);

				if (m1.gettName().equalsIgnoreCase(m2.gettName())) {
					m1 = null;
					m2 = null;
					continue;
				}
//				// do not considering contain case
//				int m1PS = m1.getqStart();
//				int m1PE = m1.getqEnd();
//				int m2PS = m2.getqStart();
//				int m2PE = m2.getqEnd();
//				if(m1PS <= m2PS && m1PE >= m2PE)
//				{
//					m1 = null;
//					m2 = null;
//					continue;
//				}
				PBLinkM p = new PBLinkM();
				p.setOrigin(m1);
				p.setTerminus(m2);
				p.setId(m1.getqName());
//				int distance = p.getDistance();
//				if(distance < 0)
//				{
//					int olength = m1.gettLength();
//					int tlength = m2.gettLength();
//					distance = Math.abs(distance);
//					if(distance >= olength || distance >= tlength)
//						continue;
//				}
				pbLinks.add(p);
				p = null;
				m1 = null;
				m2 = null;
			}
			// build TriadLink
			if (cpSize >= 3) {
//				List<TriadLink> triad  = new Vector<TriadLink>();
				for (int i = 0; i <= cpSize - 3; i++) {
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
					if (tls.contains(tl)) {
						int index = tls.indexOf(tl);
						int supLink = tls.get(index).getSupLinks();
						supLink += 1;
						tls.get(index).setSupLinks(supLink);
						m1 = null;
						m2 = null;
						m3 = null;
						pre = null;
						mid = null;
						lst = null;
						tl = null;
					} else {
						tls.add(tl);
						m1 = null;
						m2 = null;
						m3 = null;
						pre = null;
						mid = null;
						lst = null;
						tl = null;
					}
				}
//				if(triads.size() > 0)
//					this.tlWriter.write2(triads);
			}
			// build crossing links
			if(cpSize >= 4)
			{
				MRecord m1 = valids.get(0);
				MRecord m2 = valids.get(cpSize - 1);
				Contig first = new Contig();
				first.setID(m1.gettName());
				Contig last = new Contig();
				last.setID(m2.gettName());
				TriadLink tl = new TriadLink();
				tl.setPrevious(first);
				tl.setLast(last);
				tl.setSupLinks(1);
				if (tls.contains(tl)) {
					int index = tls.indexOf(tl);
					int supLink = tls.get(index).getSupLinks();
					supLink += 1;
					tls.get(index).setSupLinks(supLink);
					m1 = null;
					m2 = null;
					first = null;
					last = null;
					tl = null;
				} else {
					tls.add(tl);
					m1 = null;
					m2 = null;
					first = null;
					last = null;
					tl = null;
				}
			}

			return pbLinks;
		}
	}
	
	public List<PBLinkM> mRecord2Link(List<MRecord> records)
	{
		Iterator<MRecord> it = records.iterator();
		List<PBLinkM> pbLinks = new Vector<PBLinkM>();
		List<MRecord> valids = new Vector<MRecord>(records.size());
		while (it.hasNext()) {
			MRecord m = it.next();
			int minOLLen = paras.getMinOLLen();
			double minOLRatio = paras.getMinOLRatio();
			int maxOHLen = paras.getMaxOHLen();
			double maxOHRatio = paras.getMaxOHRatio();
			int maxEndLen = paras.getMaxEndLen();
			double maxEndRatio = paras.getMaxEndRatio();

			int pbLen = m.getqLength();
			int pbStart = m.getqStart();
			int pbEnd = m.getqEnd() - 1;
			int contLen = m.gettLength();
			int contStart = m.gettStart();
			int contEnd = m.gettEnd() - 1;
			Strand tStrand = m.gettStrand();

			boolean isInner = false;
			int endDefLen = (int) Math.round(pbLen * maxEndRatio);
			// if max end length larger than the endDefLen, used the endDefLen
			// as threshold
			if (endDefLen < maxEndLen)
				maxEndLen = endDefLen;
			if (pbStart >= maxEndLen && pbStart <= (pbLen - maxEndLen)) {
				if (pbEnd <= (pbLen - maxEndLen))
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
			int ohDefLen = (int) Math.round(contLen * maxOHRatio);
			if (ohDefLen < maxOHLen)
				maxOHLen = ohDefLen;
			if (isInner) {
				// the overlap length less than specified value, next;
				if (ol_len < minOLLen)
					continue;
				// if the overlap length enough for specified value, the ratio
				// is less than specified value, also next;
				if (ratio < minOLRatio)
					continue;
				// for two end length of contig not allow larger than maxOHLen
				if (contLeftLen > maxOHLen || contRigthLen > maxOHLen)
					continue;
			} else {
				if(pbStart <= maxEndLen && pbEnd <= maxEndLen)
				{ // checking the right side, for the end point in p1-p2
					if(contRigthLen > maxOHLen)
						continue;
				} else if(pbStart <= maxEndLen && pbEnd <= (pbLen - maxEndLen))
				{
					if(ol_len < minOLLen)
						continue;
					if(contRigthLen > maxOHLen)
						continue;
				} else if(pbStart <= maxEndLen && pbEnd >= (pbLen - maxEndLen))
				{
					// start point in p1-p2 and end point in p3-p4
					// do no afford info, discard
					continue;
				} else if(pbStart >= maxEndLen && pbStart <= (pbLen - maxEndLen) 
						&& pbEnd >= (pbLen - maxEndLen))
				{ // start point in p2-p3 and end point in p3-p4
					if(ol_len < minOLLen)
						continue;
					if(contLeftLen > maxOHLen)
						continue;
				} else if(pbStart >= (pbLen - maxEndLen) && pbEnd >= (pbLen - maxEndLen))
				{
					if(contLeftLen > maxOHLen)
						continue;
				}
//				if (pbStart <= maxEndLen && pbEnd <= (pbLen - maxEndLen)) {
//					// check right side, for the end point in p1-p2 and p2-p3
//					if (ol_len < minOLLen)
//						continue;
//					if (contRigthLen > maxOHLen)
//						continue;
//				} else if (pbStart <= maxEndLen && pbEnd >= (pbLen - maxEndLen)) {
//					// for start point in p1-p2 and end point in p3-p4, outer
//					// case!
//					// do no afford info, next directly;
//					// if (ol_len < minOLLen)
//					// continue;
//					continue;
//				} else if (pbStart >= maxEndLen && pbEnd >= (pbLen - maxEndLen)) {
//					// check left side, for start point in p2-p3 and end point
//					// in p3-p4,
//					// and start point in p3-p4 and end point in p3-p4;
//					if (ol_len < minOLLen)
//						continue;
//					if (contLeftLen > maxOHLen)
//						continue;
//				}
			}
			// if the mrecord is valid;
			valids.add(m);
		}

		// iterator valid mrecords
		if (valids.size() < 2) {
			return null;
		} else {
			// sorting the contig_pairs;
			Collections.sort(valids, new ByLocOrderComparator());
			// checking the similarity contigs
//			valids = validateContigPair(valids);
			valids = this.validateSimilarityContigs(valids);
			int cpSize = valids.size();

			// build only the successive link; A->B->C, it will build A->B
			// and B->C, omitted A->C
			for (int i = 0; i <= cpSize - 2; i++) {
				MRecord m1 = valids.get(i);
				MRecord m2 = valids.get(i + 1);

				if (m1.gettName().equalsIgnoreCase(m2.gettName())) {
					m1 = null;
					m2 = null;
					continue;
				}
				PBLinkM p = new PBLinkM();
				p.setOrigin(m1);
				p.setTerminus(m2);
				p.setId(m1.getqName());
				pbLinks.add(p);
				p = null;
				m1 = null;
				m2 = null;
			}
			// build TriadLink
			if (cpSize >= 3) {
//				List<TriadLink> triads = new Vector<TriadLink>();
				for (int i = 0; i <= cpSize - 3; i++) {
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
					if (tls.contains(tl)) {
						int index = tls.indexOf(tl);
						int supLink = tls.get(index).getSupLinks();
						supLink += 1;
						tls.get(index).setSupLinks(supLink);
						m1 = null;
						m2 = null;
						m3 = null;
						pre = null;
						mid = null;
						lst = null;
						tl = null;
					} else {
						tls.add(tl);
						m1 = null;
						m2 = null;
						m3 = null;
						pre = null;
						mid = null;
						lst = null;
						tl = null;
					}
				}
//				if(triads.size() > 0)
//					this.tlWriter.write2(triads);
			}
			return pbLinks;
		}
	}

	// method to change valid m5record to link two contigs;
	private List<PBLinkM> mRecord2Link(List<MRecord> mRecords, Parameter paras) {
		long start = System.currentTimeMillis();
		List<PBLinkM> pbLinks = new Vector<PBLinkM>();
		HashMap<String, List<MRecord>> pSet = new HashMap<String, List<MRecord>>();
		// get all valid M5Record, and storing in a hash map by the pacbio read
		// id;
		for (MRecord m : mRecords) {
			int minOLLen = paras.getMinOLLen();
			double minOLRatio = paras.getMinOLRatio();
			int maxOHLen = paras.getMaxOHLen();
			double maxOHRatio = paras.getMaxOHRatio();
			int maxEndLen = paras.getMaxEndLen();
			double maxEndRatio = paras.getMaxEndRatio();

			int pbLen = m.getqLength();
			int pbStart = m.getqStart();
			int pbEnd = m.getqEnd() - 1;
			int contLen = m.gettLength();
			int contStart = m.gettStart();
			int contEnd = m.gettEnd() - 1;
			Strand tStrand = m.gettStrand();

			// all the basic test condition is implemented in M5Reader class
			// with PB and CNT length and identity
			// pacbio length less than specified value, next;
			// if (pbLen < minPBLen.intValue())
			// continue;
			// // contig length less than specified value, next;
			// if (contLen < minContLen.intValue())
			// continue;
			// // if the contig and pacbio identity less than specified value,
			// next;
			// if (m5.getIdentity() < identity)
			// continue;
			// check the contig inner or outer of PacBio read, four points
			// values defined it;
			// pb read model: p1----p2-------p3----p4,
			// where p1 denoted as read zero point;
			// p2 denoted as (pb_read_len * max_end_ratio) if not larger than
			// 500bp, else 500bp;
			// p3 denoted as (pb_read_len * max_end_ratio) if not larger than
			// 500bp, else 500bp;
			// p4 denoted as read terminus point;
			// if overlap start point in [p1,p2], and end point in [p1,p2],
			// defined contig is outer;
			// if overlap start point in [p1,p2], and end point in [p2,p3],
			// defined contig is outer;
			// if overlap start point in [p1,p2], and end point in [p3,p4],
			// defined contig is outer;
			// if overlap start point in [p2,p3], and end point in [p2,p3],
			// defined contig is inner;
			// if overlap start point in [p2,p3], and end point in [p3,p4],
			// defined contig is outer;
			// if overlap start point in [p3,p4], and end point in [p3,p4],
			// defined contig is outer;
			boolean isInner = false;
			int endDefLen = (int) Math.round(pbLen * maxEndRatio);
			// if max end length larger than the endDefLen, used the endDefLen
			// as threshold
			if (endDefLen < maxEndLen)
				maxEndLen = endDefLen;
			if (pbStart >= maxEndLen && pbStart <= (pbLen - maxEndLen)) {
				if (pbEnd <= (pbLen - maxEndLen))
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
			int ohDefLen = (int) Math.round(contLen * maxOHRatio);
			if (ohDefLen < maxOHLen)
				maxOHLen = ohDefLen;
			if (isInner) {
				// the overlap length less than specified value, next;
				if (ol_len < minOLLen)
					continue;
				// if the overlap length enough for specified value, the ratio
				// is less than specified value, also next;
				if (ratio < minOLRatio)
					continue;
				// for two end length of contig not allow larger than maxOHLen
				if (contLeftLen > maxOHLen || contRigthLen > maxOHLen)
					continue;
			} else {
				if (pbStart <= maxEndLen && pbEnd <= (pbLen - maxEndLen)) {
					// check right side, for the end point in p1-p2 and p2-p3
					if (ol_len < minOLLen)
						continue;
					if (contRigthLen > maxOHLen)
						continue;
				} else if (pbStart <= maxEndLen && pbEnd >= (pbLen - maxEndLen)) {
					// for start point in p1-p2 and end point in p3-p4, outer
					// case!
					// do no afford info, next directly;
					// if (ol_len < minOLLen)
					// continue;
					continue;
				} else if (pbStart >= maxEndLen && pbEnd >= (pbLen - maxEndLen)) {
					// check left side, for start point in p2-p3 and end point
					// in p3-p4,
					// and start point in p3-p4 and end point in p3-p4;
					if (ol_len < minOLLen)
						continue;
					if (contLeftLen > maxOHLen)
						continue;
				}
			}
			// put all valid records to a hashmap by key as the read id;
			if (pSet.containsKey(m.getqName())) {
				pSet.get(m.getqName()).add(m);
			} else {
				List<MRecord> records = new Vector<MRecord>();
				records.add(m);
				pSet.put(m.getqName(), records);
				records = null;
			}
		}

		logger.debug(this.getClass().getName() + "\t" + "Valid link: " + pSet.size());
		List<String> repeats = null;
		if (paras.isRepMask()) {
			RepeatFinder rf = new RepeatFinder(paras);
			repeats = rf.findRepeat(pSet);
		}
		// transform valid M5Record to Pacbio links
		int count = 0;
		List<TriadLink> triads = new Vector<TriadLink>();
		for (String s : pSet.keySet()) {
			count++;
			List<MRecord> contig_pairs = pSet.get(s);
			// remove the repeats
			if (paras.isRepMask()) {
				List<MRecord> temp = new Vector<MRecord>(contig_pairs.size());
				try {
					for (MRecord m : contig_pairs) {
						String cName = m.gettName();
						if (!repeats.contains(cName))
							temp.add(m);
					}
				} catch (Exception e) {
					logger.debug(
							this.getClass().getName() + "\t" + count + e.getMessage() + "\t" + e.getClass().getName());
				}
				contig_pairs = null;
				contig_pairs = temp;
				temp = null;
			}
			// sorting the contig_pairs;
			Collections.sort(contig_pairs, new ByLocOrderComparator());
			// checking the similarity contigs
			if (contig_pairs.size() >= 2)
				contig_pairs = validateContigPair(contig_pairs);
			int cpSize = contig_pairs.size();

			// at least having two contigs under the same pacbio read;
			if (cpSize > 1) {
				// build only the successive link; A->B->C, it will build A->B
				// and B->C, omitted A->C
				for (int i = 0; i <= cpSize - 2; i++) {
					MRecord m1 = contig_pairs.get(i);
					MRecord m2 = contig_pairs.get(i + 1);

					if (m1.gettName().equalsIgnoreCase(m2.gettName())) {
						m1 = null;
						m2 = null;
						continue;
					}
					PBLinkM p = new PBLinkM();
					p.setOrigin(m1);
					p.setTerminus(m2);
					p.setId(s);
					pbLinks.add(p);
					p = null;
					m1 = null;
					m2 = null;
				}
				// build TriadLink
				if (cpSize >= 3) {
					for (int i = 0; i <= cpSize - 3; i++) {
						MRecord m1 = contig_pairs.get(i);
						MRecord m2 = contig_pairs.get(i + 1);
						MRecord m3 = contig_pairs.get(i + 2);
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
						if (triads.contains(tl)) {
							int index = triads.indexOf(tl);
							int supLink = triads.get(index).getSupLinks();
							supLink += 1;
							triads.get(index).setSupLinks(supLink);
							m1 = null;
							m2 = null;
							m3 = null;
							pre = null;
							mid = null;
							lst = null;
							tl = null;
						} else {
							triads.add(tl);
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
				// build all the combination link;
				/*
				 * for(int i = 0; i <= cpSize - 2; i++) { for(int j = i + 1 ; j
				 * <= cpSize - 1; j++) { M5Record m1 = contig_pairs.get(i);
				 * M5Record m2 = contig_pairs.get(j);
				 * if(m1.gettName().equalsIgnoreCase(m2.gettName())) { m1 =
				 * null; m2 = null; continue; } PBLinkM5 p = new PBLinkM5();
				 * p.setOrigin(m1); p.setTerminus(m2); p.setId(s);
				 * pbLinks.add(p); p = null; m1 = null; m2 = null; } }
				 */
			}
		}
		TriadLinkWriter tlw = new TriadLinkWriter(paras, triads);
		tlw.write();
		pSet = null;
		long end = System.currentTimeMillis();
		logger.info("Link Building, erase time: " + (end - start) + " ms");
		return pbLinks;
	}

	

	// this method used to checking almost identical contigs linking by same
	// pacbio reads
	private List<MRecord> validateContigPair(List<MRecord> data) {
		List<MRecord> simCNTs = new Vector<MRecord>(5);
		List<MRecord> delCNTs = new Vector<MRecord>(5);
		// int minOLLen = paras.getMinOLLen();// defined by parameter;
		int minOLLen = 300; // defined by myself;
		int start = 0;
		int end = 0;
		for (int i = 0; i < data.size(); i++) {
			MRecord m = data.get(i);
			if (i == 0) {
				start = m.getqStart();
				end = m.getqEnd();
				continue;
			}
			int tStart = m.getqStart();
			int tEnd = m.getqEnd();
			int ol = 0;
			if (end > tEnd) {
				ol = tEnd - tStart;
			} else {
				ol = end - tStart;
			}
			if (ol > minOLLen) {
				if (!simCNTs.contains(data.get(i - 1)))
					simCNTs.add(data.get(i - 1));
				if (!simCNTs.contains(m))
					simCNTs.add(m);
			} else {
				start = tStart;
				end = tEnd;
			}
		}
		if (simCNTs.size() != 0) {
			double olweight = 0.6;
			double identweight = 0.4;
			long score = 0;
			// logger.debug("Similairty contigs:");
			for (int i = 0; i < simCNTs.size(); i++) {
				MRecord m = simCNTs.get(i);
				// logger.debug(m.gettName());
				int pS = m.getqStart();
				int pE = m.getqEnd();
				int ol = pE - pS;
				int ident = (int) Math.round(m.getIdentity() * 1000);
				if (i == 0) {
					score = Math.round(ol * olweight + ident * identweight);
					continue;
				}
				long temp = Math.round(ol * olweight + ident * identweight);
				if (score <= temp) {
					delCNTs.add(simCNTs.get(i - 1));
					score = temp;
				} else {
					delCNTs.add(m);
				}
			}
			data.removeAll(delCNTs);
		}
		return data;
	}
	
	private List<MRecord> validateSimilarityContigs(List<MRecord> data)
	{
		Map<String, List<MRecord>> sims = this.findSimilarityCnts(data);
		int size = sims.size();
		if(size != 0)
		{
			for(int i = 1; i <= size; i++)
			{
				List<MRecord> ms = sims.get(String.valueOf(i));
//				data.remove(ms); // try remove all similarity contigs
				List<Integer> scores = new Vector<Integer>(ms.size());
				Iterator<MRecord> it = ms.iterator();
				double olweight = 0.6;
				double identweight = 0.4;
				while(it.hasNext())
				{
					MRecord m = it.next();
					int pS = m.getqStart();
					int pE = m.getqEnd();
					int ol = pE - pS;
					int ident = (int) Math.round(m.getIdentity() * 1000);
					int score = (int) Math.round(ol * olweight + ident * identweight);
					scores.add(score);
				}
				int max = MathTool.max(scores);
				it = ms.iterator();
				while(it.hasNext())
				{
					MRecord m = it.next();
					int pS = m.getqStart();
					int pE = m.getqEnd();
					int ol = pE - pS;
					int ident = (int) Math.round(m.getIdentity() * 1000);
					int score = (int) Math.round(ol * olweight + ident * identweight);
					if(score != max)
					{
						if(data.contains(m))
							data.remove(m);
					}
				}
			}
		}
		return data;
	}
	
	private List<MRecord> validateOverlapContigs(List<MRecord> data)
	{
		List<MRecord> removes = new Vector<MRecord>(data.size());
		Iterator<MRecord> it = data.iterator();
		boolean isFirst = true;
		MRecord former = null;
		MRecord current = null;
		double olRatio = 0.9; // larger than this value, it will remove;
		while(it.hasNext())
		{
			if(isFirst)
			{
				former = it.next();
				isFirst = false;
				continue;
			} else
			{
//				if(former.getqName().equalsIgnoreCase("m130605_032054_42207_c100515142550000001823076608221373_s1_p0/147280/0_8554"))
//					logger.debug("breakpoint");
				current = it.next();
				int dist = this.getDistance(former, current);
				if(dist >= 0)
				{
					former = current;
					continue;
				}
				
				// former code
				int fPS = former.getqStart();
				int fPE = former.getqEnd();
				int fl = fPE - fPS;
				int cPS = current.getqStart();
				int cPE = current.getqEnd();
				int cl = cPE - cPS;
				int ol = fPE - cPS;
				double fOLRatio = (double)ol / fl;
				double cOLRatio = (double)ol / cl;
				if(fOLRatio >= olRatio || cOLRatio >= olRatio)
				{
					if(fOLRatio >= cOLRatio)
					{
						removes.add(former);
						former = current;
						continue;
					} else
					{
						removes.add(current);
						continue;
					}
//					if((fOLRatio >= olRatio))
//					{
//						logger.debug("remove former" + former.getqName());
//						logger.debug("remove former " + former.gettName());
//						removes.add(former);
//						former = current;
//						continue;
//					}
//					if((cOLRatio >= olRatio))
//					{
//						logger.debug("remover current " + current.getqName());
//						logger.debug("remover current " + current.gettName());
//						removes.add(current);
//						continue;
//					}
				} else
				{
					former = current;
				}
			}
		}
		data.removeAll(removes);
		return data;
	}
	
	private int getDistance(MRecord origin, MRecord terminus)
	{
		int dist = 0;
		int oPBS = origin.getqStart(); // origin pacbio start point;
		int oCntS = origin.gettStart(); // origin contig start point;
		int oPBE = origin.getqEnd(); // origin pacbio end point;
		int oCntE = origin.gettEnd(); // origin contig end point;
		int tPBS = terminus.getqStart(); // terminus pacbio start point;
		int tCntS = terminus.gettStart(); // terminus contig start point;
		int tPBE = terminus.getqEnd(); // terminus pacbio end point;
		int tCntE = terminus.gettEnd(); // terminus contig end point;
		int oCntLen = origin.gettLength();// origin contig length;
		int tCntLen = terminus.gettLength(); // terminus contig length
		if(origin.gettStrand().equals(Strand.FORWARD) && terminus.gettStrand().equals(Strand.FORWARD))
		{ // + +
			int oRightLen = oCntLen - oCntE;
			int tLeftLen = tCntS;
			dist = tPBS - oPBE - oRightLen - tLeftLen;
		} else if(origin.gettStrand().equals(Strand.REVERSE) && terminus.gettStrand().equals(Strand.REVERSE))
		{ // - -
			int oLeftLen = oCntS;
			int tRightLen = tCntLen -tCntE;
			dist = tPBS - oPBE - oLeftLen - tRightLen;
		} else if(origin.gettStrand().equals(Strand.FORWARD) && terminus.gettStrand().equals(Strand.REVERSE))
		{ // + -
			int oRightLen = oCntLen - oCntE;
			int tRightLen = tCntLen -tCntE;
			dist = tPBS - oPBE - oRightLen - tRightLen;
		} else if(origin.gettStrand().equals(Strand.REVERSE) && terminus.gettStrand().equals(Strand.FORWARD))
		{ // - +
			int oLeftLen = oCntS;
			int tLeftLen = tCntS;
			dist = tPBS - oPBE - oLeftLen - tLeftLen;
		}
		return dist;
	}
	
	// test in yeast for some pacbio read do not valid by the filter parameter
	// but it provide some link info!
	private void pesudoTriadLink(List<MRecord> ms)
	{
		if(tls == null)
			tls = new Vector<TriadLink>(100);
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
		if (tls.contains(tl)) {
			int index = tls.indexOf(tl);
			int supLink = tls.get(index).getSupLinks();
			supLink += 1;
			tls.get(index).setSupLinks(supLink);
			m1 = null;
			m2 = null;
			first = null;
			last = null;
			tl = null;
		} else {
			tls.add(tl);
			m1 = null;
			m2 = null;
			first = null;
			last = null;
			tl = null;
		}
	}
	
	
	private Map<String, List<MRecord>> findSimilarityCnts(List<MRecord> data)
	{
		Map<String, List<MRecord>> sims = new HashMap<String, List<MRecord>>();
		List<MRecord> tempSims = new Vector<MRecord>(5);
		Iterator<MRecord> it = data.iterator();
		int count = 1;
		int preStart = 0;
		int preEnd = 0;
		int constant = 100; // 100 bp to defined whether contigs is similarity
		boolean isFirst = true;
		while(it.hasNext())
		{
			MRecord m = it.next();
			if(isFirst)
			{ // the first time
				preStart = m.getqStart();
				preEnd = m.getqEnd();
				tempSims.add(m);
				isFirst = false;
				continue;
			} else
			{
				int curStart = m.getqStart();
				int curEnd = m.getqEnd();
				
				int left = curStart - preStart;
				int right = curEnd - preEnd;
				left = Math.abs(left);
				right = Math.abs(right);
				
				if(left <= constant && right <= constant)
				{ // similarity contigs
					if(!tempSims.contains(m))
						tempSims.add(m);
					preStart = curStart;
					preEnd = curEnd;
				} else
				{
					preStart = curStart;
					preEnd = curEnd;
					if(tempSims.size() <= 1)
					{
						tempSims.clear();
						if(!tempSims.contains(m))
							tempSims.add(m);
					} else
					{
						sims.put(String.valueOf(count), tempSims);
						tempSims = new Vector<MRecord>(5);
						count++;
					}
				}
			}			
		}
		if(!tempSims.isEmpty() && tempSims.size() >= 2)
		{
			sims.put(String.valueOf(count), tempSims);
			tempSims = null;
		}
		return sims;
	}

	public List<TriadLink> getTriadLinks()
	{
		return tls;
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
		// String[] ar1 = ((String) cFirst).split("==");
		// String[] ar2 = ((String) cSecond).split("==");
		// int diff = Integer.valueOf(ar1[1]) - Integer.valueOf(ar2[1]);
		// if (diff > 0)
		// return 1;
		// if (diff < 0)
		// return -1;
		// else
		// return 0;
	}

}

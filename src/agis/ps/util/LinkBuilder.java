/*
*File: agis.ps.util.LinkBuilder.java
*User: mqin
*Email: mqin@ymail.com
*Date: 2016年1月12日
*/
package agis.ps.util;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Vector;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import agis.ps.SimplePath;
import agis.ps.file.TriadLinkWriter;
import agis.ps.link.M5Record;
import agis.ps.link.MRecord;
import agis.ps.link.PBLink;
import agis.ps.link.PBLinkM;
import agis.ps.link.TriadLink;
import agis.ps.seqs.Contig;

public class LinkBuilder {
	private static Logger logger = LoggerFactory.getLogger(LinkBuilder.class);
	private List<MRecord> ms;
	private Parameter paras;

	public LinkBuilder() {
		// do nothing;
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
	
	// directly building pb links, do not consideration of overlap length and overhang;
	private List<PBLinkM> mRecord2Link2(List<MRecord> records, Parameter paras)
	{
		long start = System.currentTimeMillis();
		// initiate the pblinks based on half of records size
		List<PBLinkM> pbLinks = new Vector<PBLinkM>(records.size()/2);
		Map<String, List<MRecord>> pSet = new HashMap<String, List<MRecord>>(records.size()/2);
		// with some filtering setting;
		// considering overhang of contig;
		for(MRecord m : records)
		{

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
			maxEndLen = 100;
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
			// for the reversed strand, BLASR define the coordinate according to aligned seq 
			if(tStrand.equals(Strand.REVERSE))
			{
				int temp = contLeftLen;
				contLeftLen = contRigthLen;
				contRigthLen = temp;
			}
			int ohDefLen = (int) Math.round(contLen * maxOHRatio);
//			if (ohDefLen < maxOHLen)
//				maxOHLen = ohDefLen;
			maxOHLen = 50;
			if (isInner) {
				// the overlap length less than specified value, next;
//				if (ol_len < minOLLen)
//					continue;
//				// if the overlap length enough for specified value, the ratio
//				// is less than specified value, also next;
//				if (ratio < minOLRatio)
//					continue;
				// for two end length of contig not allow larger than maxOHLen
				if (contLeftLen > maxOHLen || contRigthLen > maxOHLen)
					continue;
			} else {
				if (pbStart <= maxEndLen && pbEnd <= (pbLen - maxEndLen)) {
					// check right side, for the end point in p1-p2 and p2-p3
//					if (ol_len < minOLLen)
//						continue;
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
//					if (ol_len < minOLLen)
//						continue;
					if (contLeftLen > maxOHLen)
						continue;
				}
			}
			
			if (pSet.containsKey(m.getqName())) {
				pSet.get(m.getqName()).add(m);
			} else {
				List<MRecord> rs = new Vector<MRecord>();
				rs.add(m);
				pSet.put(m.getqName(), rs);
				rs = null;
			}
		}
		logger.debug(this.getClass().getName() + "\t" + "Valid link: " + pSet.size());

		List<String> repeats = null;
		if(paras.isRepMask())
		{
			RepeatFinder rf = new RepeatFinder(paras);
			repeats =  rf.findRepeat(pSet);
		}
		// transform valid M5Record to Pacbio links
		int count = 0;
		List<TriadLink> triads = new Vector<TriadLink>();
		for (String s : pSet.keySet()) {
			count++;
			List<MRecord> contig_pairs = pSet.get(s);
			// remove the repeats
			if(paras.isRepMask())
			{
				Iterator<MRecord> it = contig_pairs.iterator();
				while(it.hasNext())
				{
					MRecord m = it.next();
					if(repeats.contains(m.gettName()))
					{
						it.remove();
					}
				}
			}
			int cpSize = contig_pairs.size();

			// at least having two contigs under the same pacbio read;
			if (cpSize > 1) {
				// sorting the contig_pairs;
				Collections.sort(contig_pairs, new ByLocOrderComparator());
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
				if(cpSize >= 3)
				{
					for(int i = 0; i <= cpSize - 3; i++)
					{
						MRecord m1 = contig_pairs.get(i);
						MRecord m2 = contig_pairs.get(i + 1);
						MRecord m3 = contig_pairs.get(i + 2);
						Contig pre = new Contig();
						pre.setID(m1.gettName());
//						pre.setLength(m1.gettLength());
						Contig mid = new Contig();
						mid.setID(m2.gettName());
//						mid.setLength(m2.gettLength());
						Contig lst = new Contig();
						lst.setID(m3.gettName());
//						lst.setLength(m3.gettLength());
						TriadLink tl = new TriadLink(pre, mid, lst);
						tl.setSupLinks(1);
						if(triads.contains(tl))
						{
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
						} else
						{
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
			}
		}
		TriadLinkWriter tlw = new TriadLinkWriter(paras, triads);
		tlw.write();
		pSet = null;
		
		long end = System.currentTimeMillis();
		logger.info("Link building earasing time : " + (end - start) + " ms");
		return pbLinks;
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
			// for the reversed strand, BLASR define the coordinate according to aligned seq 
			if(tStrand.equals(Strand.REVERSE))
			{
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
		if(paras.isRepMask())
		{
			RepeatFinder rf = new RepeatFinder(paras);
			repeats =  rf.findRepeat(pSet);
		}
		// transform valid M5Record to Pacbio links
		// String [] repeats = new String[]{"1035","1045","1049","1059"};
//		List<String> repeats = new Vector<String>(4);
//		repeats.add("1035");
//		repeats.add("1045");
//		repeats.add("1049");
//		repeats.add("1059");
//		repeats.add("1025");
//		repeats.add("1027");
//		repeats.add("1037");
//		repeats.add("1013");
//		repeats.add("1015");
//		repeats.add("1019");
		int count = 0;
//		Map<TriadLink, Integer> triads = new HashMap<TriadLink, Integer>();
		List<TriadLink> triads = new Vector<TriadLink>();
		for (String s : pSet.keySet()) {
			count++;
			List<MRecord> contig_pairs = pSet.get(s);
			// remove the repeats
			if(paras.isRepMask())
			{
				List<MRecord> temp = new Vector<MRecord>(contig_pairs.size());
				try {
					for (MRecord m : contig_pairs) {
	//					String cName = m.gettName();
	//					if (!repeats.contains(cName))
	//						temp.add(m);
						String cName = m.gettName();
						if(!repeats.contains(cName))
							temp.add(m);
					}
				} catch (Exception e) {
					logger.debug(this.getClass().getName() + "\t" + count + e.getMessage() + "\t" + e.getClass().getName());
				}
				contig_pairs = null;
				contig_pairs = temp;
				temp = null;
			}
			// sorting the contig_pairs;
			Collections.sort(contig_pairs, new ByLocOrderComparator());
			// checking the similarity contigs 
			if(contig_pairs.size() >= 2)
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
				if(cpSize >= 3)
				{
					for(int i = 0; i <= cpSize - 3; i++)
					{
						MRecord m1 = contig_pairs.get(i);
						MRecord m2 = contig_pairs.get(i + 1);
						MRecord m3 = contig_pairs.get(i + 2);
						Contig pre = new Contig();
						pre.setID(m1.gettName());
//						pre.setLength(m1.gettLength());
						Contig mid = new Contig();
						mid.setID(m2.gettName());
//						mid.setLength(m2.gettLength());
						Contig lst = new Contig();
						lst.setID(m3.gettName());
//						lst.setLength(m3.gettLength());
						TriadLink tl = new TriadLink(pre, mid, lst);
						tl.setSupLinks(1);
						if(triads.contains(tl))
						{
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
						} else
						{
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
	
	// this method used to checking almost identical contigs linking by same pacbio reads
	private List<MRecord> validateContigPair(List<MRecord> data)
	{
		List<MRecord> values = new Vector<MRecord>(data.size());
		List<MRecord> simCNTs = new Vector<MRecord>(5);
		List<MRecord> delCNTs = new Vector<MRecord>(5);
//		int minOLLen = paras.getMinOLLen();// defined by parameter;
		int minOLLen = 300; // defined by myself;
		int start = 0;
		int end = 0;
		for(int i = 0; i < data.size(); i++)
		{
			MRecord m = data.get(i);
			if(i == 0)
			{
				start = m.getqStart();
				end = m.getqEnd();
				continue;
			}		
			int tStart = m.getqStart();
			int tEnd = m.getqEnd();
			int ol = 0;
			if(end > tEnd)
			{
				ol = tEnd - tStart;
			} else
			{
				ol =  end - tStart;
			}
			if(ol > minOLLen)
			{
				if(!simCNTs.contains(data.get(i - 1)))
					simCNTs.add(data.get(i - 1));
				if(!simCNTs.contains(m))
					simCNTs.add(m);
			} else
			{
				start = tStart;
				end = tEnd;
			}
		}
		if(simCNTs.size() != 0){
			double olweight = 0.6;
			double identweight = 0.4;
			long score = 0;
//			logger.debug("Similairty contigs:");
			for(int i = 0; i < simCNTs.size(); i++)
			{
				MRecord m = simCNTs.get(i);
//				logger.debug(m.gettName());
				int pS = m.getqStart();
				int pE = m.getqEnd();
				int ol = pE - pS;
				int ident = (int) Math.round(m.getIdentity() * 1000);
				if(i == 0)
				{
					score = Math.round(ol * olweight + ident * identweight);
					continue;
				}
				long temp = Math.round(ol * olweight + ident * identweight);
				if(score <= temp)
				{
					delCNTs.add(simCNTs.get(i - 1));
					score = temp;
				} else
				{
					delCNTs.add(m);
				}
			}
			data.removeAll(delCNTs);
		}
		return data;
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

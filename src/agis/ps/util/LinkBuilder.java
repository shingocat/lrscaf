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
import java.util.List;
import java.util.Map;
import java.util.Vector;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import agis.ps.M5Record;
import agis.ps.SimplePath;
import agis.ps.link.Contig;
import agis.ps.link.PBLink;
import agis.ps.link.PBLinkM5;

public class LinkBuilder {
	private static Logger logger = LoggerFactory.getLogger(LinkBuilder.class);
	private List<M5Record> m5s;
	private Parameter paras;

	public LinkBuilder() {
		// do nothing;
	}

	public LinkBuilder(List<M5Record> m5s, Parameter paras) {
		this.m5s = m5s;
		this.paras = paras;
	}

	public List<PBLinkM5> m5Record2Link() {
		if (m5s == null || m5s.size() == 0)
			throw new IllegalArgumentException("LinkBuilder: The M5Records could not be empty!");
		return this.m5Record2Link(m5s, paras);
	}

	// method to change valid m5record to link two contigs;
	public List<PBLinkM5> m5Record2Link(List<M5Record> m5Records, Parameter paras) {
		List<PBLinkM5> pbLinks = new Vector<PBLinkM5>();
		HashMap<String, List<M5Record>> pSet = new HashMap<String, List<M5Record>>();
		// get all valid M5Record, and storing in a hash map by the pacbio read
		// id;
		for (M5Record m5 : m5Records) {
			if(m5.gettName().equals("1025"))
				System.out.println(m5.toString());
			
			int minOLLen = paras.getMinOLLen();
			double minOLRatio = paras.getMinOLRatio();
			int maxOHLen = paras.getMaxOHLen();
			double maxOHRatio = paras.getMaxOHRatio();
			int maxEndLen = paras.getMaxEndLen();
			double maxEndRatio = paras.getMaxEndRatio();

			int pbLen = m5.getqLength();
			int pbStart = m5.getqStart();
			int pbEnd = m5.getqEnd() - 1;
			int contLen = m5.gettLength();
			int contStart = m5.gettStart();
			int contEnd = m5.gettEnd() - 1;
			Strand tStrand = m5.gettStrand();

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
			if (pSet.containsKey(m5.getqName())) {
				pSet.get(m5.getqName()).add(m5);
			} else {
				List<M5Record> records = new Vector<M5Record>();
				records.add(m5);
				pSet.put(m5.getqName(), records);
				records = null;
			}
		}

		logger.debug(this.getClass().getName() + "\t" + "Valid link: " + pSet.size());
		RepeatFinder rf = new RepeatFinder();
		List<String> repeats =  rf.findRepeat(pSet);
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
		for (String s : pSet.keySet()) {
			count++;
			List<M5Record> contig_pairs = pSet.get(s);

			List<M5Record> temp = new Vector<M5Record>(contig_pairs.size());
			// remove the repeats
			try {
				for (M5Record m : contig_pairs) {
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

			int cpSize = contig_pairs.size();

			// at least having two contigs under the same pacbio read;
			if (cpSize > 1) {
				// sorting the contig_pairs;
				Collections.sort(contig_pairs, new ByLocOrderComparator());
				// build only the successive link; A->B->C, it will build A->B
				// and B->C, omitted A->C
				for (int i = 0; i <= cpSize - 2; i++) {
					M5Record m1 = contig_pairs.get(i);
					M5Record m2 = contig_pairs.get(i + 1);

					if (m1.gettName().equalsIgnoreCase(m2.gettName())) {
						m1 = null;
						m2 = null;
						continue;
					}
					PBLinkM5 p = new PBLinkM5();
					p.setOrigin(m1);
					p.setTerminus(m2);
					p.setId(s);
					pbLinks.add(p);
					p = null;
					m1 = null;
					m2 = null;
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
		pSet = null;
		return pbLinks;
	}

}

class ByLocOrderComparator implements Comparator<Object> {
	@Override
	public int compare(Object cFirst, Object cSecond) {
		M5Record m1 = (M5Record) cFirst;
		M5Record m2 = (M5Record) cSecond;
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

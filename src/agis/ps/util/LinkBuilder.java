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
	
	public LinkBuilder()
	{
		// do nothing;
	}
	
	public LinkBuilder(List<M5Record> m5s, Parameter paras)
	{
		this.m5s = m5s;
		this.paras = paras;
	}
	
	public List<PBLinkM5> m5Record2Link()
	{
		if(m5s == null || m5s.size() == 0)
			throw new IllegalArgumentException("LinkBuilder: The M5Records could not be empty!");
		return this.m5Record2Link(m5s, paras);
	}
	
	// method to change valid m5record to link two contigs;
	public List<PBLinkM5> m5Record2Link(List<M5Record> m5Records, Parameter paras) {
		Integer minPBLen = paras.getMinPBLen();
		Integer minContLen = paras.getMinContLen();
		Integer minOLLen = paras.getMinOLLen();
		Double minOLRatio = paras.getMinOLRatio();
		Integer maxOHLen = paras.getMaxOHLen();
		Double maxOHRatio = paras.getMaxOHRatio();
		Integer maxEndLen = paras.getMaxEndLen();
		Double maxEndRatio = paras.getMaxEndRatio();
		Double identity = paras.getIdentity();
//		List<PBLink> pbLinks = new Vector<PBLink>();
		List<PBLinkM5> pbLinks = new Vector<PBLinkM5>();
		HashMap<String, List<M5Record>> pSet = new HashMap<String, List<M5Record>>();
		// get all valid M5Record, and storing in a hashmap by the pacbio read id;
		for (M5Record m5 : m5Records) {
			int pbLen = m5.getqLength();
			int pbStart = m5.getqStart();
			int pbEnd = m5.getqEnd() - 1;
			int contLen = m5.gettLength();
			int contStart = m5.gettStart();
			int contEnd = m5.gettEnd() - 1;
			
			// pacbio length less than specified value, next;
			if (pbLen < minPBLen.intValue())
				continue;
			// contig length less than specified value, next;
			if (contLen < minContLen.intValue())
				continue;
			// if the contig and pacbio identity less than specified value, next;
			if (m5.getIdentity() < identity)
				continue;
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
			// defined inner or not
			// if(pbStart < maxEndLen){
			// if(pbEnd > (pbLen - 1 - maxEndLen))
			// isInner = true;
			// } else if(pbStart > maxEndLen && pbStart < (pbLen - 1 -
			// maxEndLen))
			// {
			// if(pbEnd < (pbLen - 1 - maxEndLen))
			// isInner = true;
			// }
			if (pbStart >= maxEndLen && pbStart <= (pbLen - maxEndLen)) {
				if (pbEnd <= (pbLen - maxEndLen))
					isInner = true;
			}
			int ol_len = pbEnd - pbStart + 1;
			double ratio = (double) ol_len / contLen;
			int contLeftLen = contStart;
			int contRigthLen = contLen - contEnd - 1;
			int ohDefLen = (int) Math.round(contLen * maxOHRatio);
			if (ohDefLen < maxOHLen)
				maxOHLen = ohDefLen;
			if (isInner) {
				// the overlap length less than specified value, next;
				if (ol_len < minOLLen.intValue())
					continue;
				// if the overlap length enough for specified value, the ratio
				// is less than specified value, also next;
				if (ratio <= minOLRatio.doubleValue())
					continue;
				// for two end length of contig not allow larger than maxOHLen
				if (contLeftLen > maxOHLen || contRigthLen > maxOHLen)
					continue;
			} else {
				if (pbStart <= maxEndLen && pbEnd <= (pbLen - maxEndLen)) {
					// for left side, p1-p2 and p2-p3
					if (ol_len < minOLLen)
						continue;
					if (contRigthLen > maxOHLen)
						continue;
				} else if (pbStart <= maxEndLen && pbEnd >= (pbLen - maxEndLen)) { 
					// for p1-p2 and p3-p4 outer case!
					if (ol_len < minOLLen)
						continue;
				} else if (pbStart >= maxEndLen && pbEnd >= (pbLen - maxEndLen)) 
				{
					// for right side, p2-p3 and p3-p4
					if (ol_len < minOLLen)
						continue;
					if (contLeftLen > maxOHLen)
						continue;
				}
			}
			// put all valid records to a hashmap by key as the read id;
//			logger.debug(m5.getqName());
			if (pSet.containsKey(m5.getqName())) {
				pSet.get(m5.getqName()).add(m5);
			} else {
				List<M5Record> records = new Vector<M5Record>();
				records.add(m5);
				pSet.put(m5.getqName(), records);
			}
		}
		logger.debug("LinkBuilder: valid link " + pSet.size());
		// transform valid M5Record to Pacbio link
		for (String s : pSet.keySet()) {
			List<M5Record> contig_pairs = pSet.get(s);
			int cpSize = contig_pairs.size();
//			logger.debug("LinkBuilder: contig pairs size " + cpSize);
//			at least having two contigs under the same pacbio read;
			if (cpSize >= 2) {
				// sorting the contig_pairs;
				Collections.sort(contig_pairs, new ByLocOrderComparator());
				for (int i = 0; i < cpSize - 1; i++) {
					// remove the link which link itself
					PBLinkM5 pl = new PBLinkM5();
					pl.setId(s);
					pl.setOrigin(contig_pairs.get(i));
					pl.setTerminus(contig_pairs.get(i + 1));
					// remove the link which link itself
					if(pl.isSelfLink())
						continue;
					pbLinks.add(pl);
//					PBLink pl = new PBLink();
//					M5Record ar1 = (M5Record) contig_pairs.get(i);
//					M5Record ar2 = (M5Record) contig_pairs.get(i + 1);
//					Contig origin = new Contig();
//					Contig terminus = new Contig();
//					origin.setID(ar1.gettName());
//					origin.setLength(ar1.gettLength());
//					terminus.setID(ar2.gettName());
//					terminus.setLength(ar2.gettLength());
//					pl.setID(s);
//					pl.setOrigin(origin);
//					pl.setTerminus(terminus);
//					pl.setoStartLoc(ar1.getqStart());
//					pl.setoEndLoc(ar1.getqEnd());
//					pl.setoStrand(ar1.gettStrand());
//					pl.settStartLoc(ar2.getqStart());
//					pl.settEndLoc(ar2.getqEnd());
//					pl.settStrand(ar2.gettStrand());
//					pbLinks.add(pl);
				}
			}
		}

		// original code on Scaffolder;
		// List<SimplePath> simPaths = new ArrayList<SimplePath>();
		// if (simPaths != null) {
		// simPaths.clear();
		// }
		// for (String s : pSet.keySet()) {
		// String[] arrs = pSet.get(s).split(";");
		// if (arrs.length >= 2) {
		// for (int i = 0; i < arrs.length - 1; i++) {
		// SimplePath sp = new SimplePath();
		// String[] ar1 = arrs[i].split("==");
		// String[] ar2 = arrs[i + 1].split("==");
		// if (Integer.valueOf(ar1[1]) <= Integer.valueOf(ar2[1])) {
		// sp.setStart(ar1[0]);
		// sp.setEnd(ar2[0]);
		// sp.setLabel(ar1[1] + ":" + ar2[1]);
		// sp.setColor(Color.BLUE);
		// } else {
		// sp.setStart(ar2[0]);
		// sp.setEnd(ar1[0]);
		// sp.setLabel(ar2[1] + ":" + ar1[1]);
		// sp.setColor(Color.BLUE);
		// }
		// simPaths.add(sp);
		// }
		// }
		// }
		return pbLinks;
	}

	public List<M5Record> getM5s() {
		return m5s;
	}

	public void setM5s(List<M5Record> m5s) {
		this.m5s = m5s;
	}

	public Parameter getParas() {
		return paras;
	}

	public void setParas(Parameter paras) {
		this.paras = paras;
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

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

import agis.ps.M5Record;
import agis.ps.SimplePath;
import agis.ps.link.Contig;
import agis.ps.link.PBLink;

public class LinkBuilder {
	private static int MIN_CONTIG_LENGTH = 1000; // The contig length should be
													// larger than this value as
													// an valid contig;
	private static int MIN_OVERLAP_LENGTH = 2000; // The contig overlapped
													// PacBio should be larger
													// than this value and also
													// with the ratio larger
													// than 80;
	private static double MIN_OVERLAP_RATIO = 0.8; // The contig overlapped
													// region divied to whole
													// length of contigs;
	private static int MIN_PB_LENGTH = 10000; // the minimum length of pacbio;

	public static List<PBLink> m5Record2Link(List<M5Record> m5Records, Map<String, String> parameters) {
		List<PBLink> pbLinks = new Vector<PBLink>();
		HashMap<String, List<M5Record>> pSet = new HashMap<String, List<M5Record>>();
		for (M5Record m5 : m5Records) {
			if (m5.gettLength() < LinkBuilder.MIN_CONTIG_LENGTH)
				continue;
			int ol_len = m5.gettEnd() - m5.gettStart();
			double ratio = (double)ol_len / m5.gettLength();
			if (ol_len < LinkBuilder.MIN_OVERLAP_LENGTH || ratio <= LinkBuilder.MIN_OVERLAP_RATIO)
				continue;
			//String c = m5.gettName() + "==" + m5.getqStart() + ";";
			if (pSet.containsKey(m5.getqName())) {
				pSet.get(m5.getqName()).add(m5);
				//pSet.put(m5.getqName(), pSet.get(m5.getqName()) + c);
			} else {
				List<M5Record> records = new Vector<M5Record>();
				records.add(m5);
				pSet.put(m5.getqName(), records);
				//pSet.put(m5.getqName(), c);
			}
		}

		for (String s : pSet.keySet()) {
			List<M5Record> contig_pairs = pSet.get(s);
			int cpSize = contig_pairs.size();
			if (cpSize >= 2) {
				// sorting the contig_pairs;
				Collections.sort(contig_pairs, new ByLocOrderComparator());
				for (int i = 0; i < cpSize - 1; i++) {
					PBLink pl = new PBLink();
					M5Record ar1 = (M5Record)contig_pairs.get(i);
					M5Record ar2 = (M5Record)contig_pairs.get(i + 1);
					Contig origin = new Contig();
					Contig terminus = new Contig();
					origin.setID(ar1.gettName());
					origin.setLength(ar1.gettLength());
					terminus.setID(ar2.gettName());
					terminus.setLength(ar2.gettLength());
					pl.setID(s);
					pl.setOrigin(origin);
					pl.setTerminus(terminus);
					pl.setoStartLoc(ar1.getqStart());
					pl.setoEndLoc(ar1.getqEnd());
					pl.setoStrand(ar1.gettStrand());
					pl.settStartLoc(ar2.getqStart());
					pl.settEndLoc(ar2.getqEnd());
					pl.settStrand(ar2.gettStrand());
					pbLinks.add(pl);
				}
			}
		}

		// original code on Scaffolder;
//		List<SimplePath> simPaths = new ArrayList<SimplePath>();
//		if (simPaths != null) {
//			simPaths.clear();
//		}
//		for (String s : pSet.keySet()) {
//			String[] arrs = pSet.get(s).split(";");
//			if (arrs.length >= 2) {
//				for (int i = 0; i < arrs.length - 1; i++) {
//					SimplePath sp = new SimplePath();
//					String[] ar1 = arrs[i].split("==");
//					String[] ar2 = arrs[i + 1].split("==");
//					if (Integer.valueOf(ar1[1]) <= Integer.valueOf(ar2[1])) {
//						sp.setStart(ar1[0]);
//						sp.setEnd(ar2[0]);
//						sp.setLabel(ar1[1] + ":" + ar2[1]);
//						sp.setColor(Color.BLUE);
//					} else {
//						sp.setStart(ar2[0]);
//						sp.setEnd(ar1[0]);
//						sp.setLabel(ar2[1] + ":" + ar1[1]);
//						sp.setColor(Color.BLUE);
//					}
//					simPaths.add(sp);
//				}
//			}
//		}
		return pbLinks;
	}

}

class ByLocOrderComparator implements Comparator<Object> {
	@Override
	public int compare(Object cFirst, Object cSecond) {
		M5Record m1 = (M5Record)cFirst;
		M5Record m2 = (M5Record)cSecond;
		int diff = Integer.valueOf(m1.getqStart()) - Integer.valueOf(m2.getqStart());
		if(diff > 0)
			return 1;
		if(diff < 0)
			return -1;
		else
			return 0;
//		String[] ar1 = ((String) cFirst).split("==");
//		String[] ar2 = ((String) cSecond).split("==");
//		int diff = Integer.valueOf(ar1[1]) - Integer.valueOf(ar2[1]);
//		if (diff > 0)
//			return 1;
//		if (diff < 0)
//			return -1;
//		else
//			return 0;
	}

}

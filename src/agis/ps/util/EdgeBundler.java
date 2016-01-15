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

import agis.ps.Edge;
import agis.ps.link.PBLink;

public class EdgeBundler {
	private static int MIN_LINK_NUM = 3;
	static final Logger logger = LoggerFactory.getLogger(EdgeBundler.class);

	public static List<Edge> pbLinkBundling(List<PBLink> links, Map<String, String> paras) {
		List<Edge> edges = new Vector<Edge>();
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
			if (temp.get(s).size() >= EdgeBundler.MIN_LINK_NUM) {
				count += 1;
				// statistical analysis contig pairs distance
				List<PBLink> tlinks = temp.get(s);
				List<Integer> dists = new Vector<Integer>();
				for (PBLink pb : tlinks) {
					dists.add(pb.getDistance());
				}
				int mean = MathTool.mean(dists);
				int sd = MathTool.sd(dists);
				int upper = mean + 2 * sd;
				int low = mean - 2 * sd;
				// if the distance larger than mean + 2 * sd, then remove;
				try {
//					This will show iteraction error when remove element
//					for (PBLink pb : tlinks) {
//						if (pb.getDistance() > upper || pb.getDistance() < low)
//							tlinks.remove(pb);
//					}
					for(int i = 0; i < tlinks.size(); i++)
					{
						PBLink pb = tlinks.get(i);
						if (pb.getDistance() > upper || pb.getDistance() < low)
							tlinks.remove(pb);
					}
				} catch (Exception e) {
					logger.debug(String.valueOf(count));
					logger.debug(e.getMessage());
					logger.error(e.getMessage());
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
		
		//merge two opposite direction edge to only one edge base
		// not implemented now
		return edges;
	}
}

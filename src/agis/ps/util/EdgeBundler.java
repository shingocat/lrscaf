/*
*File: agis.ps.util.EdgeBundler.java
*User: mqin
*Email: mqin@ymail.com
*Date: 2016年1月13日
*/
package agis.ps.util;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import agis.ps.Edge;
import agis.ps.link.PBLink;

public class EdgeBundler {
	private static int MIN_LINK_NUM = 3;
	public static List<Edge> pbLinkBundling(List<PBLink> links, Map<String, String> paras)
	{
		List<Edge> edges = new ArrayList<Edge>();
		// storing all the same origin and terminus to a hash map;
		Map<String, List<PBLink>> temp = new HashMap<String, List<PBLink>>();
		for(PBLink pb : links)
		{
			String id = pb.getOrigin().getID() + ":" + pb.getTerminus().getID();
			if(temp.containsKey(id))
			{
				temp.get(id).add(pb);
			} else
			{
				List<PBLink> tLinks = new ArrayList<PBLink>();
				tLinks.add(pb);
				temp.put(id, tLinks);
			}
		}
		//build the edge for fitted for criterion;
		for(String s : temp.keySet())
		{
			if(temp.get(s).size() >= EdgeBundler.MIN_LINK_NUM)
			{
				Edge edge = new Edge();
				edge.setOrigin(temp.get(s).get(0).getOrigin());
				edge.setTerminus(temp.get(s).get(0).getTerminus());
				edge.setLinkNum(temp.get(s).size());
				edges.add(edge);
			}
		}
		return edges;
	}
}



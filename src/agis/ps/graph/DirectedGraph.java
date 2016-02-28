/*
*File: agis.ps.graph.DirectedGraph.java
*User: mqin
*Email: mqin@ymail.com
*Date: 2016年2月25日
*/
package agis.ps.graph;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.Vector;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import agis.ps.Edge;
import agis.ps.link.Contig;

public class DirectedGraph extends Graph implements Serializable {

	private static final long serialVersionUID = 1L;
	private static Logger logger = LoggerFactory.getLogger(DirectedGraph.class);
	// a map where the contig id as key, and the adjacent to vertices as value;
	private Map<String, List<Contig>> adjTos = new HashMap<String, List<Contig>>();
	// a map where the contig id as key, and the adjacent from vertices as
	// value;
	private Map<String, List<Contig>> adjFroms = new HashMap<String, List<Contig>>();

	public DirectedGraph(List<Edge> edges) {
		super(edges);
		// TODO Auto-generated constructor stub
		initAdjTos();
		initAdjFroms();
	}

	private void initAdjTos() {
		// initiated the vertex point to where
		if (adjTos == null)
			adjTos = Collections.synchronizedMap(new HashMap<String, List<Contig>>());
		adjTos.clear();
		for (int i = 0; i < getEdges().size(); i++) {
			Edge e = getEdges().get(i);
			String oId = e.getOrigin().getID();
			// for vertex point to other vertex
			if (adjTos.containsKey(oId)) {
				List<Contig> adjVetices = adjTos.get(oId);
				adjVetices.add(e.getTerminus());
				adjTos.put(oId, adjVetices);
			} else {
				List<Contig> adjVetices = new Vector<Contig>();
				adjVetices.add(e.getTerminus());
				adjTos.put(oId, adjVetices);
			}
		}
	}

	private void initAdjFroms() {
		// initiated the vertex point from where
		if (adjFroms == null)
			adjFroms = Collections.synchronizedMap(new HashMap<String, List<Contig>>());
		adjFroms.clear();
		for (int i = 0; i < getEdges().size(); i++) {
			Edge e = getEdges().get(i);
			String tId = e.getTerminus().getID();
			// for vertex point from other vertex
			if (adjFroms.containsKey(tId)) {
				List<Contig> adjVetices = adjFroms.get(tId);
				adjVetices.add(e.getOrigin());
				adjFroms.put(tId, adjVetices);
			} else {
				List<Contig> adjVetices = new Vector<Contig>();
				adjVetices.add(e.getOrigin());
				adjFroms.put(tId, adjVetices);
			}
		}
	}

	public int getVertexAdjVerticesNum(Contig cnt) {
		int num = 0;
		num = getAdjVertices(cnt).size();
		return num;
	}
	
	// get vertext by specified id;
	@Override
	public Contig getVertex(String id)
	{
		Contig cnt = null;
		List<Contig> cnts = getVertices();
		for(Contig c : cnts)
		{
			if(c.getID().equals(id))
			{
				cnt = c;
				break;
			}
		}
		return cnt;
	}

	// return the specified contig's adjacent contigs, including adjacent 
	// from other contigs or adjacent to other contigs;
	@Override
	public List<Contig> getAdjVertices(Contig cnt) {
		Set<Contig> adjs = new HashSet<Contig>();
		List<Contig> pTCnts = null;
		List<Contig> pFCnts = null;
		if(adjTos.containsKey(cnt.getID()))
			pTCnts = adjTos.get(cnt.getID());
		if(adjFroms.containsKey(cnt.getID()))
			pFCnts = adjFroms.get(cnt.getID());
		if(pTCnts != null)
			adjs.addAll(pTCnts);
		if(pFCnts != null)
			adjs.addAll(pFCnts);
		List<Contig> cnts = new Vector<Contig>();
		cnts = Arrays.asList(adjs.toArray(new Contig[adjs.size()]));
		return cnts;
	}
	
	// return the next vertex by the specified current and former contig;
	@Override
	public Contig getNextVertex(Contig current, Contig former)
	{
		Contig next = null;
		List<Contig> adjs = this.getAdjVertices(current);
		if(former == null)
		{
			if(adjs.size() == 1)
			{ // normal case, only one adjacent vertex
				next = adjs.get(0);
			} else if(adjs.size() == 2)
			{ // normal case, two adjacent vertex;
				// two temporary contigs variables  
				Contig tCnt1 = adjs.get(0);
				Contig tCnt2 = adjs.get(1);
				List<Edge> tEdg1 = getEdgesInfo(current, tCnt1);
				List<Edge> tEdg2 = getEdgesInfo(current, tCnt2);
				// according to the supported link number to decide 
				// which contig will be the next;
				int tEdgSL1 = 0;
				int tEdgSL2 = 0;
				if(!tEdg1.isEmpty())
				{
					for(Edge e : tEdg1)
						tEdgSL1 += tEdgSL1 + e.getLinkNum();
				}
				if(!tEdg2.isEmpty())
				{
					for(Edge e : tEdg2)
						tEdgSL2 += tEdgSL2 + e.getLinkNum();
				}
				if(tEdgSL1 > tEdgSL2)
					next = tCnt1;
				else
					next = tCnt2;
			} else
			{
				// abnormal case, larger than two vertices 
				// return null;
				next = null;
			}
		} else // former not null
		{
			if(adjs.size() == 1)
			{
				next = adjs.get(0);
				// checking whether the former vertex equal to next
				if(!next.equals(former))
					throw new IllegalArgumentException("DirectedGraph: The former vertex was not equal to next vertex when " +
							"the adjacent vertex of the current was only one!");
			} else if(adjs.size() == 2)
			{ // normal case, return the not selected vertex
				for(Contig c : adjs)
				{
					if(!c.equals(former))
						next = c;
				}
			} else
			{ // abnormal case, return the former
				next = former;
			}
		}
		return next;
	}
	
	// return edges info between the start and end contig, 
	@Override
	public List<Edge> getEdgesInfo(Contig start, Contig end)
	{	
		List<Edge> info = new Vector<Edge>();
		for(Edge e : getEdges())
		{
			Contig original = e.getOrigin();
			Contig terminus = e.getTerminus();
			if(original.equals(start) && terminus.equals(end))
				info.add(e);
			else if(original.equals(end) && terminus.equals(start))
				info.add(e);
		}
		return info;
	}
}

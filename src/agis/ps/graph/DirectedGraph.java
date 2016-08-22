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
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.Vector;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import agis.ps.link.Edge;
import agis.ps.seqs.Contig;
import agis.ps.util.MathTool;

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
		// initAdjTos();
		// initAdjFroms();
	}
	
	

	private void initAdjTos() {
		// initiated the vertex point to where
		adjTos = null;
		adjTos = Collections.synchronizedMap(new HashMap<String, List<Contig>>());
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
		adjFroms = null;
		adjFroms = Collections.synchronizedMap(new HashMap<String, List<Contig>>());
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
	public Contig getVertex(String id) {
		Contig cnt = null;
		List<Contig> cnts = getVertices();
		for (Contig c : cnts) {
			if (c.getID().equalsIgnoreCase(id)) {
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
		if (adjTos == null || adjTos.size() == 0)
			initAdjTos();
		if (adjTos.containsKey(cnt.getID()))
			pTCnts = adjTos.get(cnt.getID());
		if (adjFroms == null || adjFroms.size() == 0)
			initAdjFroms();
		if (adjFroms.containsKey(cnt.getID()))
			pFCnts = adjFroms.get(cnt.getID());
		if (pTCnts != null)
			adjs.addAll(pTCnts);
		if (pFCnts != null)
			adjs.addAll(pFCnts);
		List<Contig> cnts = new Vector<Contig>();
		for (Contig c : adjs) {
			if (!cnts.contains(c))
				cnts.add(c);
		}

		// cnts = Arrays.asList(adjs.toArray(new Contig[adjs.size()]));
		return cnts;
	}

	// return next vertex
	public List<Contig> getNextVertices(Contig current, Contig former) {
		List<Contig> nexts = null;
		List<Contig> adjs = adjTos.get(current.getID());
		if (adjs == null || adjs.isEmpty())
			return null;
		nexts = new Vector<Contig>(5);
		for (Contig c : adjs) {
			if (!c.equals(former))
				nexts.add(c);
		}
		return nexts;
	}
	
	// return the next vertex by the specified current and former contig;
	// including return divergence vertex; and former could not be null for triad requirement
	public Contig getNextVertex2(Contig current, Contig former)
	{
		Contig next = null;
		List<Contig> adjs = this.getAdjVertices(current);
		if (former == null) {
			throw new IllegalArgumentException(this.getClass().getName() + "\t" + "The former vertex could not be null!");
		}
		
		if (adjs.size() == 1) {
			next = adjs.get(0);
			// checking whether the former vertex equal to next
			if (!next.equals(former))
				throw new IllegalArgumentException(
						"DirectedGraph: The former vertex was not equal to next vertex when "
								+ "the adjacent vertex of the current was only one!");
		} else if (adjs.size() == 2) { 
			// normal case, return the not selected vertex
			for (Contig c : adjs) {
				if (!c.equals(former))
					next = c;
			}
		} else { 
			// abnormal case, return the former
			//				next = former;
			next = null;
		}
		return next;
	}
	
	// return the next vertex by the specified current and former contig;
	// excluding return divergence vertex;
	@Override
	public Contig getNextVertex(Contig current, Contig former) {
		Contig next = null;
		List<Contig> adjs = this.getAdjVertices(current);
		if (former == null) {
			if (adjs.size() == 1) { 
				// normal case, only one adjacent vertex
				next = adjs.get(0);
			} else if (adjs.size() == 2) { 
				// normal case, two adjacent vertex; two temporary contigs variables, return the more supported links vertex
				Contig tCnt1 = adjs.get(0);
				Contig tCnt2 = adjs.get(1);
				List<Edge> tEdg1 = getEdgesInfo(current, tCnt1);
				List<Edge> tEdg2 = getEdgesInfo(current, tCnt2);
				// according to the supported link number to decide which contig will be the next;
				int tEdgSL1 = tEdg1.get(0).getLinkNum();
				int tEdgSL2 = tEdg2.get(0).getLinkNum();
				
				if (tEdgSL1 > tEdgSL2)
					next = tCnt1;
				else
					next = tCnt2;
				tEdg1 = null;
				tEdg2 = null;
			} else {
				// abnormal case, larger than two vertices return null;
				next = null;
			}
		} else // former not null
		{
			if (adjs.size() == 1) {
				next = adjs.get(0);
				// checking whether the former vertex equal to next
				if (!next.equals(former))
					throw new IllegalArgumentException(
							"DirectedGraph: The former vertex was not equal to next vertex when "
									+ "the adjacent vertex of the current was only one!");
			} else if (adjs.size() == 2) { 
				// normal case, return the not selected vertex
				for (Contig c : adjs) {
					if (!c.equals(former))
						next = c;
				}
			} else { 
				// abnormal case, return the former
//				next = former;
				next = null;
			}
		}
		return next;
	}

	// return edges info between the start and end contig,
	@Override
	public List<Edge> getEdgesInfo(Contig start, Contig end) {
		List<Edge> info = new Vector<Edge>();
		for (Edge e : getEdges()) {
			Contig original = e.getOrigin();
			Contig terminus = e.getTerminus();
			if (original.equals(start) && terminus.equals(end))
				info.add(e);
			else if (original.equals(end) && terminus.equals(start))
				info.add(e);
		}
		return info;
	}

	@Override
	public void transitiveReducting() {
		long start = System.currentTimeMillis();
		if(adjTos == null || adjTos.size() == 0)
			initAdjTos();
		// it need to transitive reduction candidate contigs list;
		List<Contig> trCnds = new Vector<Contig>();
		// it store the replication contig in the same transitive reduction,
		// it do need to do transitive reduction again
		List<Contig> repCnts = new Vector<Contig>();
		for(String id : adjTos.keySet())
		{
			List<Contig> adjs = adjTos.get(id);
			int adjCount = adjs.size();
			// by far only considering three adjacent statement;
			if(adjCount >= 3)
			{
				trCnds.add(getVertex(id));
			}
		}
		// checking
		if(trCnds.size() == 0)
		{
			logger.debug(this.getClass().getName() + " The graph do not contain transitive reduction structure!");
			return;
		}
		// the depth for searching, the alternative path could only accept 5 node;
		int depth = 5;
		int count = 1;
		for(Contig cnt : trCnds)
		{
			if(repCnts.contains(cnt))
				continue;
			List<Contig> adjs = adjTos.get(cnt.getID());
			LinkedList<Contig> p = null;
			// finding transitive reduction path
outer:			for(Contig c : adjs)
			{
				Contig current = c;
				Contig previous = cnt;
				p = new LinkedList<Contig>();
				p.addLast(previous);
				p.addLast(current);
				List<Contig> nexts = this.getNextVertices(current, previous);
inner:				while(nexts != null)
				{
					// if the path more than 10 contigs and could not find the start contig;
					// break;
					if(count > depth)
					{
						count = 1;
						p = null;
						break inner;
					}
					if(nexts.size() >= 2)
					{
						if(!nexts.contains(cnt))
						{
							count = 1;
							p = null;
							break inner;
						} else
						{
							nexts = null;
							repCnts.add(current);
							current = null;
							previous = null;
							break outer;
						}
					} else if(nexts.size() == 1){
						if(nexts.contains(cnt))
						{
							nexts = null;
							current = null;
							previous = null;
							break outer;
						}
						previous = current;
						current = nexts.get(0);
						nexts = this.getNextVertices(current, previous);
						p.addLast(current);
					} else if(nexts.size() == 0)
					{
						count = 1;
						p = null;
						break inner;
					}
					count++;
				}
			}
			// transitive reduction
			if(p == null || p.isEmpty())
				continue;
			// checking which edge should be transitive reduction
			List<Edge> fEs = this.getEdgesInfo(cnt, p.get(1));
			List<Edge> aEs = this.getEdgesInfo(cnt, p.getLast());
			int fDist = fEs.get(0).getDistMean();
			int aDist = aEs.get(0).getDistMean();
			Contig end = null;
			boolean isForward = false;
			if(fDist > aDist)
			{
				Edge e = fEs.get(0);
				if(e.getOrigin().equals(cnt))
					end = e.getTerminus();
				else
					end = e.getOrigin();
				isForward = false;
			} else
			{
				Edge e = aEs.get(0);
				if(e.getOrigin().equals(cnt))
					end = e.getTerminus();
				else
					end = e.getOrigin();
				isForward = true;
			}
			
			List<Edge> trEs = this.getEdgesInfo(cnt, end);
			List<Edge> alEs = new Vector<Edge>();
			int alDist = 0;
			int alSd = 0;
			if(isForward)
			{
				List<Integer> alDists = new Vector<Integer>(p.size() * 2);
				List<Integer> alSds = new Vector<Integer>(p.size() * 2);
	
				for(int i = 0; i < p.size() - 1; i++)
				{
					Contig current = p.get(i);
					Contig next = p.get(i + 1);
					List<Edge> temp = this.getEdgesInfo(current, next);
					alEs.addAll(temp);
					
					alDists.add(temp.get(0).getDistMean());
					alSds.add(temp.get(0).getDistSd());
					
					if(i != 0)
						alDists.add(current.getLength());
				}
				
				alDist = MathTool.sum(alDists);
				alSd = MathTool.avgSd(alSds);
			} else
			{
				p.addLast(cnt);
				List<Integer> alDists = new Vector<Integer>(p.size() * 2);
				List<Integer> alSds = new Vector<Integer>(p.size() * 2);
	
				for(int i = p.size() - 1; i >= 1 ; i--)
				{
					Contig current = p.get(i);
					Contig next = p.get(i - 1);
					List<Edge> temp = this.getEdgesInfo(current, next);
					alEs.addAll(temp);
					
					alDists.add(temp.get(0).getDistMean());
					alSds.add(temp.get(0).getDistSd());
					
					if(i != p.size() - 1)
						alDists.add(current.getLength());
				}
			
				alDist = MathTool.sum(alDists);
				alSd = MathTool.avgSd(alSds);
			}


			int trDist = trEs.get(0).getDistMean();
			int trSd = trEs.get(0).getDistSd();
			int trSL = trEs.get(0).getLinkNum();
			
			int sd = trSd >= alSd ? trSd : alSd;
			int diff = trDist - alDist;
			int range = 10 * sd;
			if(diff >= -range && diff <= range)
				{
					for(Edge e: trEs)
					{
						this.edges.remove(e);
					}
					for(int i = 0; i < alEs.size(); i++)
					{
						Edge e = alEs.get(i);
						e.setLinkNum(e.getLinkNum() + trSL);
					}
				}
		}
		long end = System.currentTimeMillis();
		logger.info("Transitive Reducing, erase time: " + (end - start) + " ms");
		updateGraph();
	}

	@Override
	public void delErrorProneEdge(double ratio) {
		long start = System.currentTimeMillis();
		if(adjTos == null || adjTos.size() == 0)
			initAdjTos();
		List<Edge> rmEdges = new Vector<Edge>(100);
		for(String id : adjTos.keySet())
		{
			List<Contig> adjs = adjTos.get(id);
			int adjCount = adjs.size();
			// by far only considering three adjacent statement;
			if(adjCount >= 3)
			{
				Contig cnt = new Contig();
				cnt.setID(id);
				int [] sls = new int[adjCount];
				for(int i = 0; i < adjCount; i++)
				{
					Contig c = adjs.get(i);
					List<Edge> es = this.getEdgesInfo(cnt, c);
					sls[i] = es.get(0).getLinkNum();
				}
				Arrays.sort(sls);
				int max = sls[adjCount - 1];
				for(int i = 0; i < adjCount; i++)
				{
					Contig c = adjs.get(i);
					List<Edge> es = this.getEdgesInfo(cnt, c);
					int sl = es.get(0).getLinkNum();
					double r = (double)sl/max;
					if(r <= ratio)
					{
//						this.edges.removeAll(es);
						rmEdges.addAll(es);
					}
				}
			}
		}
		logger.info(this.getClass().getName() + "\tDelete error prone edges:" + rmEdges.size());
		this.edges.removeAll(rmEdges);
		long end = System.currentTimeMillis();
		logger.info("Error prone Edge deleting, erase time: " + (end - start) + " ms");
		updateGraph();
	}

	@Override
	public void linearMergin() {
		// TODO Auto-generated method stub

	}

	@Override
	public boolean removeEdge(Edge e) {
		// TODO Auto-generated method stub
		boolean isRemove = false;
		if (this.getEdges().remove(e))
			isRemove = true;
		return false;
	}
	

	public void updateGraph() {
		initAdjTos();
		initAdjFroms();
	}



	@Override
	public void delTips() {
		// TODO Auto-generated method stub
		
	}
}

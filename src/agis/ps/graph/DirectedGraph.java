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

import agis.ps.Edge;
import agis.ps.link.Contig;
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
		if (adjTos != null)
			adjTos = null;
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
		if (adjFroms != null)
			adjFroms = null;
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
		nexts = new Vector<Contig>();
		for (Contig c : adjs) {
			if (!c.equals(former))
				nexts.add(c);
		}
		return nexts;
	}

	// return the next vertex by the specified current and former contig;
	@Override
	public Contig getNextVertex(Contig current, Contig former) {
		Contig next = null;
		List<Contig> adjs = this.getAdjVertices(current);
		if (former == null) {
			if (adjs.size() == 1) { // normal case, only one adjacent vertex
				next = adjs.get(0);
			} else if (adjs.size() == 2) { // normal case, two adjacent vertex;
											// two temporary contigs variables
				Contig tCnt1 = adjs.get(0);
				Contig tCnt2 = adjs.get(1);
				List<Edge> tEdg1 = getEdgesInfo(current, tCnt1);
				List<Edge> tEdg2 = getEdgesInfo(current, tCnt2);
				// according to the supported link number to decide
				// which contig will be the next;
				int tEdgSL1 = 0;
				int tEdgSL2 = 0;
				if (!tEdg1.isEmpty()) {
					for (Edge e : tEdg1)
						tEdgSL1 += tEdgSL1 + e.getLinkNum();
				}
				if (!tEdg2.isEmpty()) {
					for (Edge e : tEdg2)
						tEdgSL2 += tEdgSL2 + e.getLinkNum();
				}
				if (tEdgSL1 > tEdgSL2)
					next = tCnt1;
				else
					next = tCnt2;
			} else {
				// abnormal case, larger than two vertices
				// return null;
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
			} else if (adjs.size() == 2) { // normal case, return the not
											// selected vertex
				for (Contig c : adjs) {
					if (!c.equals(former))
						next = c;
				}
			} else { // abnormal case, return the former
				next = former;
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
		// TODO Auto-generated method stub
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
			if(adjCount == 3)
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
		// the depth for searching, the alternative path could only accepty 5 node;
		int depth = 10;
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
					List<Integer> tDist = new Vector<Integer>(temp.size());
					List<Integer> tSd = new Vector<Integer>(temp.size());
					for(Edge e : temp)
					{
						tDist.add(e.getDistMean()); 
						tSd.add(e.getDistSd());
					}
					try{
						alDists.add(MathTool.mean(tDist));
						alSds.add(MathTool.avgSd(tSd));
					} catch(Exception e)
					{
						logger.debug(this.getClass().getName() + e.getMessage() + "\t" + e.getClass().getName());
						logger.info(this.getClass().getName() + e.getMessage() + "\t" + e.getClass().getName());
					}
					if(i != 0)
						alDists.add(current.getLength());
				}
				try{
					alDist = MathTool.mean(alDists);
					alSd = MathTool.avgSd(alSds);
				} catch(Exception e)
				{
					logger.debug(this.getClass().getName() + e.getMessage() + "\t" + e.getClass().getName());
					logger.info(this.getClass().getName() + e.getMessage() + "\t" + e.getClass().getName());
				}
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
					List<Integer> tDist = new Vector<Integer>(temp.size());
					List<Integer> tSd = new Vector<Integer>(temp.size());
					for(Edge e : temp)
					{
						tDist.add(e.getDistMean()); 
						tSd.add(e.getDistSd());
					}
					try{
						alDists.add(MathTool.mean(tDist));
						alSds.add(MathTool.avgSd(tSd));
					} catch(Exception e)
					{
						logger.debug(this.getClass().getName() + e.getMessage() + "\t" + e.getClass().getName());
						logger.info(this.getClass().getName() + e.getMessage() + "\t" + e.getClass().getName());
					}
					if(i != p.size() - 1)
						alDists.add(current.getLength());
				}
				try{
					alDist = MathTool.mean(alDists);
					alSd = MathTool.avgSd(alSds);
				} catch(Exception e)
				{
					logger.debug(this.getClass().getName() + e.getMessage() + "\t" + e.getClass().getName());
					logger.info(this.getClass().getName() + e.getMessage() + "\t" + e.getClass().getName());
				}
			}
			List<Integer> trDists = new Vector<Integer>(trEs.size());
			List<Integer> trSds = new Vector<Integer>(trEs.size());
			List<Integer> trSLs = new Vector<Integer>(trEs.size());
			for(Edge e : trEs)
			{
				trDists.add(e.getDistMean());
				trSds.add(e.getDistSd());
				trSLs.add(e.getLinkNum());
			}
			int trDist = 0;
			int trSd = 0;
			int trSL = 0;
			try
			{
				trDist = MathTool.mean(trDists);
				trSd = MathTool.avgSd(trSds);
				trSL = MathTool.mean(trSLs);
			} catch(Exception e)
			{
				logger.debug(this.getClass().getName() + e.getMessage() + "\t" + e.getClass().getName());
				logger.info(this.getClass().getName() + e.getMessage() + "\t" + e.getClass().getName());
			}
			
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
				
//			List<Edge> e1s = this.getEdgesInfo(p.getFirst(), p.getLast());
//			List<Integer> e1Dists = new Vector<Integer>(2);
//			List<Integer> e1Sds = new Vector<Integer>(2);
//			for(Edge e: e1s)
//			{
//				e1Dists.add(e.getDistMean());
//				e1Sds.add(e.getDistSd());
//			}
//			int e1Dist = 0;
//			int e1Sd = 0; // square_root(sum(square(sd)));
//			try {
//				e1Dist = MathTool.mean(e1Dists);
//				e1Sd = MathTool.avgSd(e1Sds);
//			} catch (Exception e) {
//				logger.debug(this.getClass().getName() + e.getMessage() + "\t" + e.getClass().getName());
//				logger.info(this.getClass().getName() + e.getMessage() + "\t" + e.getClass().getName());
//			}
//			
//			List<Integer> e2Dists = new Vector<Integer>(p.size());
//			List<Integer> e2Sds = new Vector<Integer>(p.size());
//			int e2Dist = 0;
//			int e2Sd = 0;
//			// bug, need to fix
//			for(int i = 0; i < p.size() - 1; i++)
//			{
//				Contig current = p.get(i);
//				Contig next = p.get(i + 1);
//				List<Edge> e2s = this.getEdgesInfo(current, next);
//				List<Integer> tDists = new Vector<Integer>(e2s.size());
//				for(Edge e : e2s)
//				{
////					e2Dists.add(e.getDistMean());
//					tDists.add(e.getDistMean());
//					e2Sds.add(e.getDistSd());
//				}
//				int mean = 0;
//				try{
//					mean = MathTool.mean(tDists);
//					e2Dists.add(mean);
//				} catch (Exception e)
//				{
//					logger.debug(this.getClass().getName() + e.getMessage() + "\t" + e.getClass().getName());
//					logger.info(this.getClass().getName() + e.getMessage() + "\t" + e.getClass().getName());
//				}
//				if(i != 0)
//					e2Dists.add(current.getLength());
//			}
//			try{
////				e2Dist = MathTool.mean(e2Dists);
//				e2Dist = MathTool.sum(e2Dists);
//				e2Sd = MathTool.avgSd(e2Sds);
//			} catch(Exception e)
//			{
//				logger.debug(this.getClass().getName() + e.getMessage() + "\t" + e.getClass().getName());
//				logger.info(this.getClass().getName() + e.getMessage() + "\t" + e.getClass().getName());
//			}
//			int sd = e1Sd >= e2Sd ? e1Sd : e2Sd;
//			int diff = e1Dist - e2Dist;
//			int range = 10 * sd;
//			if(diff >= -range && diff <= range)
//			{
//				for(Edge e: e1s)
//				{
//					this.edges.remove(e);
//				}
//				for(int i = 0; i < p.size() - 1; i++)
//				{
//					Contig current = p.get(i);
//					Contig next = p.get(i + 1);
//					for(Edge e : edges)
//					{
//						if(e.getOrigin().equals(current) && e.getTerminus().equals(next))
//							e.setLinkNum(e.getLinkNum() + 1);
//						else if(e.getOrigin().equals(next) && e.getTerminus().equals(current))
//							e.setLinkNum(e.getLinkNum() + 1);
//					}
//				}
//				
//			}
		}
		
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
}

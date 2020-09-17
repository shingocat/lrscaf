/*
*File: agis.ps.graph2.DirectedGraph.java
*User: mqin
*Email: mqin@ymail.com
*Date: 2016年7月5日
*/
package agis.ps.graph;

//import java.io.File;
//import java.io.IOException;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Vector;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import agis.ps.file.SimilarityCntReader;
import agis.ps.file.TriadLinkWriter;
//import agis.ps.link.CntFileEncapsulate;
import agis.ps.link.Edge;
//import agis.ps.seqs.Scaffold;
//import agis.ps.seqs.Contig;
import agis.ps.seqs.Sequence;
import agis.ps.util.MathTool;
import agis.ps.util.Parameter;
import agis.ps.util.Strand;

public class DirectedGraph extends Graph implements Serializable {

	private static final long serialVersionUID = 1L;
	private static Logger logger = LoggerFactory.getLogger(DirectedGraph.class);
	private static int TR_TIMES = 15; // for transitive reduction
	private int TIP_LENGTH = 1500; // 1000 bp for tip length;
	private Parameter paras = null;
//	private Map<String, List<Contig>> adjTos = Collections.synchronizedMap(new HashMap<String, List<Contig>>());
	private Map<String, List<Sequence>> adjTos = Collections.synchronizedMap(new HashMap<String, List<Sequence>>());
	// for store multiple in and multiple out contig vertex;
//	private List<Contig> mimos = new ArrayList<Contig>();
	private List<Sequence> mimos = new ArrayList<Sequence>();
	private TriadLinkWriter tlWriter = null;
//	private CntFileEncapsulate cntfile;
//	private Map<String, Integer> cntLens;
//	private Map<String, Contig> cnts;
	private Map<String, Sequence> seqs;

	public DirectedGraph(List<Edge> edges, Parameter paras) {
		super(edges);
		this.paras = paras; 
		this.TIP_LENGTH = paras.getTipLength();
		initAdjTos();
//		initCntIndexer();
		tlWriter = new TriadLinkWriter(paras);
		tlWriter.init(true);
	}
	
//	public DirectedGraph(List<Edge> edges, Parameter paras, Map<String, Contig> cnts) {
//		this(edges, paras);
//		this.cnts= cnts;
//	}
	
	public DirectedGraph(List<Edge> edges, Parameter paras, Map<String, Sequence> seqs) {
		this(edges, paras);
		this.seqs= seqs;
	}
	
	
//	public DirectedGraph(List<Edge> edges, Parameter paras, Map<String, Integer> cntLens)
//	{
//		this(edges, paras);
//		this.cntLens = cntLens;
//	}

	private void initAdjTos() {
		if (adjTos == null)
			adjTos = new HashMap<String, List<Sequence>>();
		adjTos.clear();
		Iterator<Map.Entry<String, Edge>> it = edges.entrySet().iterator();
		while (it.hasNext()) {
			Map.Entry<String, Edge> entry = it.next();
			Edge e = entry.getValue();
			Sequence origin = e.getOrigin();
			Sequence terminus = e.getTerminus();
			String id = origin.getId();
			if (adjTos.containsKey(id)) {
				List<Sequence> temp = adjTos.get(id);
				if (!temp.contains(terminus))
					temp.add(terminus);
				adjTos.replace(id, temp);
				if (temp.size() >= 3) {
					if (!mimos.contains(origin))
						mimos.add(origin);
				}
			} else {
				List<Sequence> temp = new Vector<Sequence>();
				temp.add(terminus);
				adjTos.put(id, temp);
			}

		}
	}

	public int getVertexAdjVerticesNum(Sequence cnt) {
		List<Sequence> adjCnts = this.getAdjVertices(cnt);
		if(adjCnts == null)
			return 0;
		else
			return adjCnts.size();
	}
	
	@Override
	public boolean isDivergenceVertex(Sequence cnt) {
		if(this.getVertexAdjVerticesNum(cnt) > 2)
			return true;
		else 
			return false;
	}

	// return next vertices by current and former vertex
	// but exclude former contig;
	@Override
	public List<Sequence> getNextVertices(Sequence current, Sequence former) {
		List<Sequence> adjs = this.adjTos.get(current.getId());
		List<Sequence> values = new Vector<Sequence>(adjs.size() - 1);
		Iterator<Sequence> it = adjs.iterator();
		while (it.hasNext()) {
			Sequence c = it.next();
			if (!c.equals(former))
				values.add(c);
		}
		if (values == null || values.isEmpty())
			return null;
		return values;
	}

	@Override
	public Sequence getVertex(String id) {
//		Contig cnt = null;
//		if(this.vertices.containsKey(id))
//			cnt = this.vertices.get(id);
		Sequence cnt = new Sequence();
		cnt.setId(id);
		if (this.vertices.contains(cnt))
			cnt = this.vertices.get(this.vertices.indexOf(cnt));
		return cnt;
	}

	@Override
	public void transitiveReducting() {
		long start = System.currentTimeMillis();
		// it store the replication contig in the same transitive reduction,
		// it do need to do transitive reduction again
		if (mimos.size() == 0) {
			logger.info(this.getClass().getName() + " The graph do not contain transitive reduction structure!");
			return;
		}
		// the depth for searching, the alternative path could only accept 5
		// node;
		int depth = 5;
//		int index = 0;
		for (Sequence origin : mimos) {
//			index++;
//			if(index == 232)
//				logger.debug("breakpoint");
			List<Sequence> cnts = this.getAdjVertices(origin);
			Iterator<Sequence> it = cnts.iterator();
			int indicator = cnts.size();
			try {
				while (it.hasNext()) {
					Sequence c = it.next();
//					if(c.getID().equalsIgnoreCase("jcf7180002720970"))
//						logger.debug("breakpoint");
//					if(!c.equals(origin))
//					{
						LinkedList<Sequence> path = new LinkedList<Sequence>();
						path.addLast(origin);
						// path.addLast(c);
						this.transitiveReducting(c, origin, origin, depth, path);
//					} else
//					{
//						List<Edge> cycEdges =  this.getEdgesInfo(origin, c);
//						this.removeEdges(cycEdges);
//					}
					int temp = cnts.size();
					if (indicator != temp) {
						it = this.getAdjVertices(origin).iterator();
						indicator = temp;
					}
					// else
					// break;
				}
			} catch (Exception e) {
				logger.debug("Error: ", e);
				logger.error(this.getClass().getName() + "\t" + e.getMessage());
//				logger.error("index " + index);
//				logger.error(origin.getID());
			}
		}
		long end = System.currentTimeMillis();
		tlWriter.close();
		logger.info("Transitive Reducing, elapsed time: " + (end - start) + " ms");
		updateGraph();
	}

	// return contig is the next contig, if the return contig equal to start
	// point,
	// the recurrence break or if the depth is equal to defined depth return;
	private boolean transitiveReducting(Sequence current, Sequence former, Sequence start, int depth,
			LinkedList<Sequence> path) {
//		logger.debug("current depth " + depth);
//		logger.debug("current node " + current.getID());
//		logger.debug("former node " + former.getID());
		path.addLast(current);
		List<Sequence> cnts = this.getNextVertices(current, former);
		if (cnts == null) {
			path.removeLast();
			return false;
		}
		if (cnts.contains(start)) {
			int size = path.size();
			List<Edge> fEs = this.getEdgesInfo(start, path.get(1));
			List<Edge> aEs = this.getEdgesInfo(start, path.getLast());
			if (fEs == null)
				return false;
			if (aEs == null)
				return false;
			int fDist = fEs.get(0).getDistMean();
			int aDist = aEs.get(0).getDistMean();
			Sequence end = null;
			LinkedList<Sequence> alPath = new LinkedList<Sequence>();
			if (fDist > aDist) {
				end = path.get(1);
				for (int i = size - 1; i >= 1; i--) {
					alPath.addLast(path.get(i));
				}
				alPath.addFirst(start);
			} else {
				end = path.getLast();
				alPath = path;
			}
			List<Edge> trEs = this.getEdgesInfo(start, end);
			List<Edge> alEs = new Vector<Edge>();
			long alDist = 0;
			int alSd = 0;
			int alSize = alPath.size();
			List<Integer> alDists = new Vector<Integer>(alSize * 2);
			List<Integer> alSds = new Vector<Integer>(alSize * 2);
			for (int i = 0; i < (alSize - 1); i++) {
				Sequence now = alPath.get(i);
				Sequence next = alPath.get(i + 1);
				List<Edge> temp = this.getEdgesInfo(now, next);
				alEs.addAll(temp);

				alDists.add(temp.get(0).getDistMean());
				alSds.add(temp.get(0).getDistSd());

				if (i > 0 && i < (alPath.size() - 1))
				{
//					alDists.add(now.getLength());
//					alDists.add(this.indexCntLength(now.getID()));
					String id = now.getId();
//					int dist = cntfile.getLengthByNewId(id);
//					int dist = cntLens.get(id);
					int dist = this.seqs.get(id).getLength();
					alDists.add(dist);
				}
			}
			alDist = MathTool.sum(alDists);
			alSd = MathTool.avgSd(alSds);
			
			// check whether the direction is conflict;
			// there are two cases;
			// 1: the conflict point is the end;
			// 2: the conflict point is the start;
			// first case:
			Strand trStrand = null;
			for(Edge e : trEs) {
				if(e.getOrigin().equals(start) && e.getTerminus().equals(end)) {
					trStrand = e.gettStrand();
					break;
				}
			}
			// considering the previous last point;
			Strand alStrand = null;
			Sequence preLst = alPath.get(alPath.size() - 2);
			List<Edge> lstEdges = this.getEdgesInfo(preLst, end);
			for(Edge e : lstEdges) {
				if(e.getOrigin().equals(preLst) && e.getTerminus().equals(end)) {
					alStrand = e.gettStrand();
					break;
				}
			}
			if(!trStrand.equals(alStrand)) {
				path.removeLast();
				return false;
			}
			// case 2:
			// considering the next first point;
			for(Edge e : trEs)
			{
				if(e.getOrigin().equals(start) && e.getTerminus().equals(end))
				{
					trStrand = e.getoStrand();
					break;
				}
			}
			Sequence nextFirst = alPath.get(1);
			List<Edge> nextEdges = this.getEdgesInfo(start, nextFirst);
			for(Edge e : nextEdges)
			{
				if(e.getOrigin().equals(start) && e.getTerminus().equals(nextFirst))
				{
					alStrand = e.getoStrand();
					break;
				}
			}
			if(!trStrand.equals(alStrand))
			{
				path.removeLast();
				return false;
			}

			int trDist = trEs.get(0).getDistMean();
			int trSd = trEs.get(0).getDistSd();
			int trSL = trEs.get(0).getLinkNum();

			int sd = trSd >= alSd ? trSd : alSd;
			long diff = trDist - alDist;
			int range = TR_TIMES * sd;
			
			// considering two standards for transitive reduction
			// case 1 : diff is in the confidence interval, 
			// case 2 : diff is less than the proportion of paths
			if (diff >= -range && diff <= range) {
				tlWriter.write4Edges(trEs);
				// remove tr edges
				this.removeEdges(trEs);
				// modify alternative path edges;
				for (Edge e : alEs) {
					e.setLinkNum(e.getLinkNum() + trSL);
				}
				path.removeLast();
				return true;
			} else {
				
				double trRatio = ((double)diff / trDist);
				double alRatio = ((double)diff / alDist);
				if(Math.abs(trRatio) <= 0.5 && Math.abs(alRatio) <= 0.5)
				{
					tlWriter.write4Edges(trEs);
					// remove tr edges
					this.removeEdges(trEs);
					// modify alternative path edges;
					for (Edge e : alEs) {
						e.setLinkNum(e.getLinkNum() + trSL);
					}
					path.removeLast();
					return true;
				} else {
					path.removeLast();
					return false;
				}
			}
		} else if (depth == 0) {
			path.removeLast();
			return false;
		} else {
			for (Sequence c : cnts) {
				transitiveReducting(c, current, start, depth - 1, path);
			}
			path.removeLast();
			return false;
		}
	}

	@Override
	public void linearMergin() {
		// TODO Auto-generated method stub

	}

	@Override
	public void delErrorProneEdge(double ratio) {
		long start = System.currentTimeMillis();
		List<Edge> rmEdges = new Vector<Edge>(100);
		try {
			Iterator<Sequence> it = mimos.iterator();
			while (it.hasNext()) {
				Sequence c = it.next();
				String id = c.getId();
				List<Sequence> adjs = adjTos.get(id);
				if (adjs == null)
					continue;
				int adjCount = adjs.size();
				Sequence cnt = new Sequence();
				cnt.setId(id);
				int[] sls = new int[adjCount];
				int[] lens = new int[adjCount + 1];
				for (int i = 0; i < adjCount; i++) {
					Sequence t = adjs.get(i);
					List<Edge> es = this.getEdgesInfo(cnt, t);
					sls[i] = es.get(0).getLinkNum();
//					lens[i] = this.indexCntLength(t.getID());
//					lens[i] = cntfile.getLengthByNewId(t.getID());
//					lens[i] = cntLens.get(t.getID());
					lens[i] = this.seqs.get(t.getId()).getLength();
				}
//				lens[adjCount] = this.indexCntLength(id);
//				lens[adjCount] = cntfile.getLengthByNewId(id);
//				lens[adjCount] = cntLens.get(id);
				lens[adjCount] = this.seqs.get(id).getLength();
				Arrays.sort(sls);
				Arrays.sort(lens);
				int max = sls[adjCount - 1];
				for (int i = 0; i < adjCount; i++) {
					Sequence t = adjs.get(i);
					List<Edge> es = this.getEdgesInfo(cnt, t);
					int sl = es.get(0).getLinkNum();
					double r = (double) sl / max;
					if (r <= ratio) {
						rmEdges.addAll(es);
					}
				}
				if(adjCount == 3)
				{
					int min = lens[0];
					int next = lens[1];
					if(min <= 500 && next >= 5000)
					{
						for (int i = 0; i < adjCount; i++) {
							Sequence t = adjs.get(i);
//							if(this.indexCntLength(t.getID()) <= 500){
//							if(cntfile.getLengthByNewId(t.getID()) <= 500){
//							if(cntLens.get(t.getID()) <= 500){
							if(this.seqs.get(t.getId()).getLength() <= 500){
								List<Edge> es = this.getEdgesInfo(cnt, t);
								int sl = es.get(0).getLinkNum();
								double r = (double) sl / max;
								if (r <= 0.3) {
									for(Edge e : es)
									{
										if(!rmEdges.contains(e))
											rmEdges.add(e);
									}
								}
							}
						}
					}
				}
			}
		} catch (Exception e) {
			logger.debug("Error: ", e);
			logger.error(this.getClass().getName() + "\t" + e.getMessage());
		}
		logger.info("Delete error prone edges:" + rmEdges.size());
		this.removeEdges(rmEdges);
		this.updateGraph();
		long end = System.currentTimeMillis();
		logger.info("Error prone Edge deleting, elapsed time: " + (end - start) + " ms");
	}

	// return the specified contig's adjacent contigs, including adjacent
	// from other contigs or adjacent to other contigs;
	@Override
	public List<Sequence> getAdjVertices(Sequence cnt) {
		List<Sequence> values = this.adjTos.get(cnt.getId());
		return values;
	}

	// return the next vertex by the specified current and former contig;
	// excluding return divergence vertex;
	@Override
	public Sequence getNextVertex(Sequence current, Sequence former) {
		Sequence next = null;
		List<Sequence> adjs = this.getAdjVertices(current);
		if (former == null) {
			if (adjs.size() == 1) {
				// normal case, only one adjacent vertex
				next = adjs.get(0);
			} else if (adjs.size() == 2) {
				// normal case, two adjacent vertex; two temporary contigs
				// variables, return the more supported links vertex
				Sequence tCnt1 = adjs.get(0);
				Sequence tCnt2 = adjs.get(1);
				List<Edge> tEdg1 = getEdgesInfo(current, tCnt1);
				List<Edge> tEdg2 = getEdgesInfo(current, tCnt2);
				// according to the supported link number to decide which contig
				// will be the next;
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
				for (Sequence c : adjs) {
					if (!c.equals(former))
						next = c;
				}
			} else {
				// abnormal case, return the former
				// next = former;
				next = null;
			}
		}
		return next;
	}

	// return the next vertex by the specified current and former contig;
	// including return divergence vertex; and former could not be null for
	// triad requirement
	public Sequence getNextVertex2(Sequence current, Sequence former) {
		Sequence next = null;
		List<Sequence> adjs = this.getAdjVertices(current);
		if (former == null) {
			throw new IllegalArgumentException(
					this.getClass().getName() + "\t" + "The former vertex could not be null!");
		}

		if (adjs.size() == 1) {
			next = adjs.get(0);
			// checking whether the former vertex equal to next
			if (!next.equals(former))
				throw new IllegalArgumentException("DirectedGraph: The former vertex was not equal to next vertex when "
						+ "the adjacent vertex of the current was only one!");
		} else if (adjs.size() == 2) {
			// normal case, return the not selected vertex
			for (Sequence c : adjs) {
				if (!c.equals(former))
					next = c;
			}
		} else {
			// abnormal case, return the former
			// next = former;
			next = null;
		}
		return next;
	}

	@Override
	public List<Edge> getEdgesInfo(Sequence start, Sequence end) {
		if(start == null || end == null)
			return null;
		List<Edge> info = new Vector<Edge>(2);
		String id1 = start.getId() + "->" + end.getId();
		String id2 = end.getId() + "->" + start.getId();
		if (this.edges.containsKey(id1))
			info.add(this.edges.get(id1));
		if (this.edges.containsKey(id2))
			info.add(this.edges.get(id2));
		if (info.size() == 0)
			return null;
		return info;
	}

	@Override
	public boolean removeEdge(Edge e) {
		boolean isRemove = false;
		Sequence origin = e.getOrigin();
		Sequence terminus = e.getTerminus();
		String id = origin.getId() + "->" + terminus.getId();
		if (this.edges.containsKey(id)) {
			if (this.edges.remove(id) != null) {
				List<Sequence> temp = this.adjTos.get(origin.getId());
				if (temp.contains(terminus))
					temp.remove(terminus);
				this.adjTos.replace(origin.getId(), temp);
				isRemove = true;
			}
		}
		return isRemove;
	}

	@Override
	public void removeEdges(List<Edge> edges) {
		Iterator<Edge> it = edges.iterator();
		while (it.hasNext()) {
			Edge e = it.next();
			Sequence origin = e.getOrigin();
			Sequence terminus = e.getTerminus();
			String id = origin.getId() + "->" + terminus.getId();
			if (this.edges.containsKey(id))
				this.edges.remove(id);
			List<Sequence> temp = this.adjTos.get(origin.getId());
			if (temp.contains(terminus))
				temp.remove(terminus);
			this.adjTos.replace(origin.getId(), temp);
		}
		// updateGraph();
	}

	public void updateGraph() {
		initAdjTos();
	}
	
	
	/**
	 * 
	 * the method interface for the call on deleting tips
	 */
	@Override
	public void delTips()
	{
		long start = System.currentTimeMillis();
		List<Edge> rmEdges = new Vector<Edge>(100);
		int totalDels = 0;
		Sequence c = null;
		try {
			Iterator<Sequence> it = mimos.iterator();
			while (it.hasNext()) {
				c = it.next();
				String id = c.getId();
//				if(id.equals("331"))
//					logger.debug("breakpoint");
				List<Sequence> adjs = adjTos.get(id);
				if (adjs == null)
					continue;
				int adjCount = adjs.size();
				// only considering divergence point
				// and considering the odd number diverence point
//				if((adjCount % 2) != 0)
//				{
					for(int i = 0; i < adjCount; i++)
					{
						Sequence next = adjs.get(i);
						LinkedList<Sequence> path = new LinkedList<Sequence>();
						path.addLast(c);
						this.delTips(next, c, 2, path, rmEdges);
					}
//				}
			}
		} catch (Exception e) {
			logger.debug("Error: ", e);
			logger.error(this.getClass().getName() + "\t" + e.getMessage());
		}
		if(rmEdges.size() != 0)
		{
			totalDels += rmEdges.size();
			this.removeEdges(rmEdges);
			this.updateGraph();
			rmEdges.clear();
		}
		logger.info("Delete tip edges:" + totalDels);
		long end = System.currentTimeMillis();
		logger.info("Tip edge deleting, elapsed time: " + (end - start) + " ms");
	}
	
	/**
	 * The inner method for deleting tips
	 * 
	 * @param current
	 * @param former
	 * @param depth
	 * @param path
	 * @param removes
	 */
	private boolean delTips(Sequence current, Sequence former, int depth, LinkedList<Sequence> path, List<Edge> removes)
	{
		path.addLast(current);
		List<Sequence> nexts = this.getNextVertices(current, former);
		if(nexts == null)
		{
			// check the tip length;
			int length = 0;
			int size = path.size() - 1; 
			List<Edge> all = new Vector<Edge>(size * 2);
			for(int i = 0; i < size; i++)
			{
				Sequence f = path.get(i);
				Sequence n = path.get(i + 1);
				List<Edge> es = this.getEdgesInfo(f, n);
				all.addAll(es);
				length += es.get(0).getDistMean();
//				length += cntfile.getLengthByNewId(n.getID());
//				length += cntLens.get(n.getID());
				length += this.seqs.get(n.getId()).getLength();
			}
			if(length <= TIP_LENGTH)
			{
				// reverse to check the path could be delete or not;
				// if meet the divergence point then break;
				Sequence c = path.get(size);
				Sequence f = path.get(size - 1);
				List<Edge> es = this.getEdgesInfo(c, f);
				for(Edge e : es)
					if(!removes.contains(e))
						removes.add(e);
				path.removeLast();
				return true;
			}
			path.removeLast();
			return false;
		}
		if(depth == 0)
		{
			if(nexts == null || nexts.size() == 0)
			{
				// check the tip length;
				int length = 0;
				int size = path.size() - 1; 
				List<Edge> all = new Vector<Edge>(size * 2);
				for(int i = 0; i < size; i++)
				{
					Sequence f = path.get(i);
					Sequence n = path.get(i + 1);
					List<Edge> es = this.getEdgesInfo(f, n);
					all.addAll(es);
					length += es.get(0).getDistMean();
//					length += cntfile.getLengthByNewId(n.getID());
//					length += cntLens.get(n.getID());
					length += this.seqs.get(n.getId()).getLength();
				}
				if(length <= TIP_LENGTH)
				{
					// reverse to check the path could be delete or not;
					// if meet the divergence point then break;
					Sequence c = path.get(size);
					Sequence f = path.get(size - 1);
					List<Edge> es = this.getEdgesInfo(c, f);
					for(Edge e : es)
						if(!removes.contains(e))
							removes.add(e);
					path.removeLast();
					return true;
				}
			}
			path.removeLast();
			return false;
		}
		// if the nexts contigs are large than 1;
		// if all the path is less than specified tip length;
		// then whole branch will be delete;
		if(nexts.size() > 1)
		{
			boolean isDel = true;
			for(Sequence c : nexts)
			{
				boolean isTrueDel = this.delTips(c, current, depth - 1, path, removes);
				if(!isTrueDel)
					isDel = false;
			}
			if(isDel)
			{
				int size = path.size() - 1;
				Sequence c = path.get(size);
				Sequence f = path.get(size - 1);
				List<Edge> es = this.getEdgesInfo(c, f);
				for(Edge e : es)
					if(!removes.contains(e))
						removes.add(e);
				path.removeLast();
				return true;
			}
			path.removeLast();
			return false;
		} else 
		{
			former = current;
			current = nexts.get(0);
			boolean isDel = this.delTips(current, former, depth - 1, path, removes);
			if(isDel)
			{
				int size = path.size() - 1;
				Sequence c = path.get(size);
				Sequence f = path.get(size - 1);
				List<Edge> es = this.getEdgesInfo(c, f);
				for(Edge e : es)
					if(!removes.contains(e))
						removes.add(e);
			}
			path.removeLast();
			return isDel;
		}
	}

	@Override
	public void delCyclerEdges() {
		List<Edge> rmEdges = new Vector<Edge>(100);
		int totalDels = 0;
		List<Edge> tempEdges = this.getEdges();
		for(Edge e : tempEdges)
		{
			if(e.getOrigin().equals(e.getTerminus()))
			{
				rmEdges.add(e);
			}
		}
		if(rmEdges.size() != 0)
		{
			totalDels += rmEdges.size();
			this.removeEdges(rmEdges);
			this.updateGraph();
			rmEdges.clear();
		}
	}

//	@Override
//	public void delSimCntEdges() {
//		List<LinkedList<String>> simcnts = this.getSimCnts();
//		List<Edge> removeEdges = new Vector<Edge>();
//		for(LinkedList<String> sims : simcnts)
//		{
//			int simsize = sims.size();
//			// only considering two elements case now
//			if(simsize == 2)
//			{
//				Contig s1 = new Contig();
//				s1.setID(sims.get(0));
//				Contig s2 = new Contig();
//				s2.setID(sims.get(1));
//				
//				List<Contig> s1adjs = this.getAdjVertices(s1);
//				List<Contig> s2adjs = this.getAdjVertices(s2);
//				for(Contig cs1 : s1adjs)
//				{
//					for(Contig cs2 : s2adjs)
//					{
//						if(cs1.equals(cs2))
//						{
//							if(this.isDivergenceVertex(cs1))
//							{
//								List<Edge> es1s = this.getEdgesInfo(s1, cs1);
//								List<Edge> es2s = this.getEdgesInfo(s2, cs2);
//								if(es1s.get(0).getLinkNum() > es2s.get(0).getLinkNum())
//								{
//									for(Edge e : es2s)
//									{
//										removeEdges.add(e);
//									}
//								} else if(es1s.get(0).getLinkNum() < es2s.get(0).getLinkNum())
//								{
//									for(Edge e: es1s)
//									{
//										removeEdges.add(e);
//									}
//								} else
//								{
//									// do nothing
//								}
//							}
//						}
//					}
//				}
//			}
//		}
//		this.removeEdges(removeEdges);
//	}
//	
//	private List<LinkedList<String>> getSimCnts()
//	{
//		List<LinkedList<String>> simcnts = new Vector<LinkedList<String>>();
//		SimilarityCntReader scr = new SimilarityCntReader(paras);
//		simcnts = scr.read();
//		return simcnts;
//	}
	
}

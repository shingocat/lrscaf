/*
*File: agis.ps.graph2.DirectedGraph.java
*User: mqin
*Email: mqin@ymail.com
*Date: 2016年7月5日
*/
package agis.ps.graph;

import java.io.File;
import java.io.IOException;
import java.io.Serializable;
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
import agis.ps.link.CntFileEncapsulate;
import agis.ps.link.Edge;
import agis.ps.seqs.Contig;
import agis.ps.util.MathTool;
import agis.ps.util.Parameter;
import agis.ps.util.Strand;

public class DirectedGraph extends Graph implements Serializable {

	private static final long serialVersionUID = 1L;
	private static Logger logger = LoggerFactory.getLogger(DirectedGraph.class);
	private static int TR_TIMES = 15; // for transitive reduction
	private int TIP_LENGTH = 1500; // 1000 bp for tip length;
	private Parameter paras = null;
	private Map<String, List<Contig>> adjTos = Collections.synchronizedMap(new HashMap<String, List<Contig>>());
	// for store multiple in and multiple out contig vertex;
	private List<Contig> mimos = new Vector<Contig>();
	private TriadLinkWriter tlWriter = null;
//	private CntFileEncapsulate cntfile;
	private Map<String, Integer> cntLens;

	public DirectedGraph(List<Edge> edges, Parameter paras) {
		super(edges);
		this.paras = paras; 
		this.TIP_LENGTH = paras.getTipLength();
		initAdjTos();
//		initCntIndexer();
		tlWriter = new TriadLinkWriter(paras);
		tlWriter.init();
	}
	
	
//	public DirectedGraph(List<Edge> edges, Parameter paras, CntFileEncapsulate cntfile)
//	{
//		this(edges, paras);
//		this.cntfile = cntfile;
//	}
	
	public DirectedGraph(List<Edge> edges, Parameter paras, Map<String, Integer> cntLens)
	{
		this(edges, paras);
		this.cntLens = cntLens;
	}

	private void initAdjTos() {
		if (adjTos == null)
			adjTos = Collections.synchronizedMap(new HashMap<String, List<Contig>>());
		if (adjTos != null)
			adjTos.clear();
		Collection<Map.Entry<String, Edge>> collections = this.edges.entrySet();
		Iterator<Map.Entry<String, Edge>> it = collections.iterator();
		while (it.hasNext()) {
			Map.Entry<String, Edge> entry = it.next();
			Edge e = entry.getValue();
			Contig origin = e.getOrigin();
			Contig terminus = e.getTerminus();
			String id = origin.getID();
			if (adjTos.containsKey(id)) {
				List<Contig> temp = adjTos.get(id);
				if (!temp.contains(terminus))
					temp.add(terminus);
				adjTos.replace(id, temp);
				if (temp.size() >= 3) {
					if (!mimos.contains(origin))
						mimos.add(origin);
				}
			} else {
				List<Contig> temp = new Vector<Contig>();
				temp.add(terminus);
				adjTos.put(id, temp);
			}

		}
	}

	public int getVertexAdjVerticesNum(Contig cnt) {
		List<Contig> adjCnts = this.getAdjVertices(cnt);
		if(adjCnts == null)
			return 0;
		else
			return adjCnts.size();
	}
	
	@Override
	public boolean isDivergenceVertex(Contig cnt)
	{
		if(this.getVertexAdjVerticesNum(cnt) > 2)
			return true;
		else 
			return false;
	}

	// return next vertices by current and former vertex
	// but exclude former contig;
	@Override
	public List<Contig> getNextVertices(Contig current, Contig former) {
		List<Contig> adjs = this.adjTos.get(current.getID());
		List<Contig> values = new Vector<Contig>(adjs.size() - 1);
		Iterator<Contig> it = adjs.iterator();
		while (it.hasNext()) {
			Contig c = it.next();
			if (!c.equals(former))
				values.add(c);
		}
		if (values == null || values.isEmpty())
			return null;
		return values;
	}

	@Override
	public Contig getVertex(String id) {
		Contig cnt = null;
		if (this.vertices.containsKey(id))
			cnt = this.vertices.get(id);
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
		for (Contig origin : mimos) {
			List<Contig> cnts = this.getAdjVertices(origin);
			Iterator<Contig> it = cnts.iterator();
			int indicator = cnts.size();
			try {
				while (it.hasNext()) {
					Contig c = it.next();
					LinkedList<Contig> path = new LinkedList<Contig>();
					path.addLast(origin);
					// path.addLast(c);
					this.transitiveReducting(c, origin, origin, depth, path);
					int temp = cnts.size();
					if (indicator != temp) {
						it = this.getAdjVertices(origin).iterator();
						indicator = temp;
					}
					// else
					// break;
				}
			} catch (Exception e) {
				logger.error(this.getClass().getName() + "\t" + e.getMessage());
			}
		}
		long end = System.currentTimeMillis();
//		tlWriter.close();
		logger.info("Transitive Reducing, erase time: " + (end - start) + " ms");
		updateGraph();
	}

	// return contig is the next contig, if the return contig equal to start
	// point,
	// the recurrence break or if the depth is equal to defined depth return;
	private boolean transitiveReducting(Contig current, Contig former, Contig start, int depth,
			LinkedList<Contig> path) {
		path.addLast(current);
		List<Contig> cnts = this.getNextVertices(current, former);
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
			Contig end = null;
			LinkedList<Contig> alPath = new LinkedList<Contig>();
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
			int alDist = 0;
			int alSd = 0;
			int alSize = alPath.size();
			List<Integer> alDists = new Vector<Integer>(alSize * 2);
			List<Integer> alSds = new Vector<Integer>(alSize * 2);
			for (int i = 0; i < (alSize - 1); i++) {
				Contig now = alPath.get(i);
				Contig next = alPath.get(i + 1);
				List<Edge> temp = this.getEdgesInfo(now, next);
				alEs.addAll(temp);

				alDists.add(temp.get(0).getDistMean());
				alSds.add(temp.get(0).getDistSd());

				if (i > 0 && i < (alPath.size() - 1))
				{
//					alDists.add(now.getLength());
//					alDists.add(this.indexCntLength(now.getID()));
					String id = now.getID();
//					int dist = cntfile.getLengthByNewId(id);
					int dist = cntLens.get(id);
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
			for(Edge e : trEs)
			{
				if(e.getOrigin().equals(start) && e.getTerminus().equals(end))
				{
					trStrand = e.gettStrand();
					break;
				}
			}
			// considering the previous last point;
			Strand alStrand = null;
			Contig preLst = alPath.get(alPath.size() - 2);
			List<Edge> lstEdges = this.getEdgesInfo(preLst, end);
			for(Edge e : lstEdges)
			{
				if(e.getOrigin().equals(preLst) && e.getTerminus().equals(end))
				{
					alStrand = e.gettStrand();
					break;
				}
			}
			if(!trStrand.equals(alStrand))
			{
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
			Contig nextFirst = alPath.get(1);
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
			int diff = trDist - alDist;
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
				if(Math.abs(trRatio) <= 0.15 && Math.abs(alRatio) <= 0.15)
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
			for (Contig c : cnts) {
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
			Iterator<Contig> it = mimos.iterator();
			while (it.hasNext()) {
				Contig c = it.next();
				String id = c.getID();
				List<Contig> adjs = adjTos.get(id);
				if (adjs == null)
					continue;
				int adjCount = adjs.size();
				Contig cnt = new Contig();
				cnt.setID(id);
				int[] sls = new int[adjCount];
				int[] lens = new int[adjCount + 1];
				for (int i = 0; i < adjCount; i++) {
					Contig t = adjs.get(i);
					List<Edge> es = this.getEdgesInfo(cnt, t);
					sls[i] = es.get(0).getLinkNum();
//					lens[i] = this.indexCntLength(t.getID());
//					lens[i] = cntfile.getLengthByNewId(t.getID());
					lens[i] = cntLens.get(t.getID());
				}
//				lens[adjCount] = this.indexCntLength(id);
//				lens[adjCount] = cntfile.getLengthByNewId(id);
				lens[adjCount] = cntLens.get(id);
				Arrays.sort(sls);
				Arrays.sort(lens);
				int max = sls[adjCount - 1];
				for (int i = 0; i < adjCount; i++) {
					Contig t = adjs.get(i);
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
							Contig t = adjs.get(i);
//							if(this.indexCntLength(t.getID()) <= 500){
//							if(cntfile.getLengthByNewId(t.getID()) <= 500){
							if(cntLens.get(t.getID()) <= 500){
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
			logger.error(this.getClass().getName() + "\t" + e.getMessage());
		}
		logger.info(this.getClass().getName() + "\tDelete error prone edges:" + rmEdges.size());
		this.removeEdges(rmEdges);
		this.updateGraph();
		long end = System.currentTimeMillis();
		logger.info("Error prone Edge deleting, erase time: " + (end - start) + " ms");
	}

	// return the specified contig's adjacent contigs, including adjacent
	// from other contigs or adjacent to other contigs;
	@Override
	public List<Contig> getAdjVertices(Contig cnt) {
		List<Contig> values = this.adjTos.get(cnt.getID());
		return values;
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
				// normal case, two adjacent vertex; two temporary contigs
				// variables, return the more supported links vertex
				Contig tCnt1 = adjs.get(0);
				Contig tCnt2 = adjs.get(1);
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
				for (Contig c : adjs) {
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
	public Contig getNextVertex2(Contig current, Contig former) {
		Contig next = null;
		List<Contig> adjs = this.getAdjVertices(current);
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
			for (Contig c : adjs) {
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
	public List<Edge> getEdgesInfo(Contig start, Contig end) {
		if(start == null || end == null)
			return null;
		List<Edge> info = new Vector<Edge>(2);
		String id1 = start.getID() + "->" + end.getID();
		String id2 = end.getID() + "->" + start.getID();
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
		Contig origin = e.getOrigin();
		Contig terminus = e.getTerminus();
		String id = origin.getID() + "->" + terminus.getID();
		if (this.edges.containsKey(id)) {
			if (this.edges.remove(id) != null) {
				List<Contig> temp = this.adjTos.get(origin.getID());
				if (temp.contains(terminus))
					temp.remove(terminus);
				this.adjTos.replace(origin.getID(), temp);
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
			Contig origin = e.getOrigin();
			Contig terminus = e.getTerminus();
			String id = origin.getID() + "->" + terminus.getID();
			if (this.edges.containsKey(id))
				this.edges.remove(id);
			List<Contig> temp = this.adjTos.get(origin.getID());
			if (temp.contains(terminus))
				temp.remove(terminus);
			this.adjTos.replace(origin.getID(), temp);
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
		Contig c = null;
		try {
			Iterator<Contig> it = mimos.iterator();
			while (it.hasNext()) {
				c = it.next();
				String id = c.getID();
//				if(id.equals("331"))
//					logger.debug("breakpoint");
				List<Contig> adjs = adjTos.get(id);
				if (adjs == null)
					continue;
				int adjCount = adjs.size();
				// only considering divergence point
				// and considering the odd number diverence point
				if((adjCount % 2) != 0)
				{
					for(int i = 0; i < adjCount; i++)
					{
						Contig next = adjs.get(i);
						LinkedList<Contig> path = new LinkedList<Contig>();
						path.addLast(c);
						this.delTips(next, c, 2, path, rmEdges);
					}
				}
			}
		} catch (Exception e) {
			logger.error(this.getClass().getName() + "\t" + e.getMessage());
		}
		if(rmEdges.size() != 0)
		{
			totalDels += rmEdges.size();
			this.removeEdges(rmEdges);
			this.updateGraph();
			rmEdges.clear();
		}
		logger.info(this.getClass().getName() + "\tDelete tip edges:" + totalDels);
		long end = System.currentTimeMillis();
		logger.info("Tip edge deleting, erase time: " + (end - start) + " ms");
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
	private boolean delTips(Contig current, Contig former, int depth, LinkedList<Contig> path, List<Edge> removes)
	{
		path.addLast(current);
		List<Contig> nexts = this.getNextVertices(current, former);
		if(nexts == null)
		{
			// check the tip length;
			int length = 0;
			int size = path.size() - 1; 
			List<Edge> all = new Vector<Edge>(size * 2);
			for(int i = 0; i < size; i++)
			{
				Contig f = path.get(i);
				Contig n = path.get(i + 1);
				List<Edge> es = this.getEdgesInfo(f, n);
				all.addAll(es);
				length += es.get(0).getDistMean();
//				length += cntfile.getLengthByNewId(n.getID());
				length += cntLens.get(n.getID());
			}
			if(length <= TIP_LENGTH)
			{
				// reverse to check the path could be delete or not;
				// if meet the divergence point then break;
				Contig c = path.get(size);
				Contig f = path.get(size - 1);
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
					Contig f = path.get(i);
					Contig n = path.get(i + 1);
					List<Edge> es = this.getEdgesInfo(f, n);
					all.addAll(es);
					length += es.get(0).getDistMean();
//					length += cntfile.getLengthByNewId(n.getID());
					length += cntLens.get(n.getID());
				}
				if(length <= TIP_LENGTH)
				{
					// reverse to check the path could be delete or not;
					// if meet the divergence point then break;
					Contig c = path.get(size);
					Contig f = path.get(size - 1);
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
			for(Contig c : nexts)
			{
				boolean isTrueDel = this.delTips(c, current, depth - 1, path, removes);
				if(!isTrueDel)
					isDel = false;
			}
			if(isDel)
			{
				int size = path.size() - 1;
				Contig c = path.get(size);
				Contig f = path.get(size - 1);
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
				Contig c = path.get(size);
				Contig f = path.get(size - 1);
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
	public void delSimCntEdges() {
		List<LinkedList<String>> simcnts = this.getSimCnts();
		List<Edge> removeEdges = new Vector<Edge>();
		for(LinkedList<String> sims : simcnts)
		{
			int simsize = sims.size();
			// only considering two elements case now
			if(simsize == 2)
			{
				Contig s1 = new Contig();
				s1.setID(sims.get(0));
				Contig s2 = new Contig();
				s2.setID(sims.get(1));
				
				List<Contig> s1adjs = this.getAdjVertices(s1);
				List<Contig> s2adjs = this.getAdjVertices(s2);
				for(Contig cs1 : s1adjs)
				{
					for(Contig cs2 : s2adjs)
					{
						if(cs1.equals(cs2))
						{
							if(this.isDivergenceVertex(cs1))
							{
								List<Edge> es1s = this.getEdgesInfo(s1, cs1);
								List<Edge> es2s = this.getEdgesInfo(s2, cs2);
								if(es1s.get(0).getLinkNum() > es2s.get(0).getLinkNum())
								{
									for(Edge e : es2s)
									{
										removeEdges.add(e);
									}
								} else if(es1s.get(0).getLinkNum() < es2s.get(0).getLinkNum())
								{
									for(Edge e: es1s)
									{
										removeEdges.add(e);
									}
								} else
								{
									// do nothing
								}
							}
						}
					}
				}
			}
		}
		this.removeEdges(removeEdges);
	}
	
	private List<LinkedList<String>> getSimCnts()
	{
		List<LinkedList<String>> simcnts = new Vector<LinkedList<String>>();
		SimilarityCntReader scr = new SimilarityCntReader(paras);
		simcnts = scr.read();
		return simcnts;
	}
	
}

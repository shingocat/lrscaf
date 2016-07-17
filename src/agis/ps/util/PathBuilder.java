/*
*File: agis.ps.util.PathBuilder.java
*User: mqin
*Email: mqin@ymail.com
*Date: 2016年1月14日
*/
package agis.ps.util;

import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Vector;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import agis.ps.file.DotGraphFileWriter;
import agis.ps.file.TriadLinkReader;
import agis.ps.graph2.DirectedGraph;
import agis.ps.graph2.Graph;
import agis.ps.link.ContInOut;
import agis.ps.link.Edge;
import agis.ps.link.TriadLink;
import agis.ps.link.TriadLinkComparator;
import agis.ps.path.Node;
import agis.ps.path.NodePath;
import agis.ps.seqs.Contig;

public class PathBuilder {
	public static Logger logger = LoggerFactory.getLogger(PathBuilder.class);
	private List<Edge> edges;
	private Parameter paras;
	private List<TriadLink> triads;
	private static int index = 0;

	public PathBuilder() {
		// do nothing;
	}

	public PathBuilder(List<Edge> edges, Parameter paras) {
		this.edges = edges;
		this.paras = paras;
	}

	public List<NodePath> buildPath() {
		if (edges == null || edges.size() == 0)
			throw new IllegalArgumentException("PathBuilder ： The Edges could not be empty!");
		return this.buildPath2(edges, paras);
	}
	
	private List<NodePath> buildPath2(List<Edge> edges, Parameter paras)
	{
		long start = System.currentTimeMillis();
		if (edges == null || edges.size() == 0)
			throw new IllegalArgumentException("PathBuilder: The Edges could not be empty!");
		List<NodePath> paths = new Vector<NodePath>();
		Graph diGraph = null;
		try {
			diGraph = new DirectedGraph(edges);
			// do transitive reduction
			diGraph.transitiveReducting();
			List<Edge> tempEdges = diGraph.getEdges();
			logger.info("Edges size after transitive reducing: " + tempEdges.size());
			String edgeFile = paras.getOutFolder() + System.getProperty("file.separator") + "edges_after_tr.info";
			DotGraphFileWriter.writeEdge(edgeFile, tempEdges);
			// delete error prone edge
			diGraph.delErrorProneEdge(paras.getRatio());
			tempEdges = diGraph.getEdges();
			logger.info("Edges size after error prone deleting: " + tempEdges.size());
			edgeFile = paras.getOutFolder() + System.getProperty("file.separator") + "edges_after_dep.info";
			DotGraphFileWriter.writeEdge(edgeFile, tempEdges);
			NodePath path = null;
//			TriadLinkReader tlr = new TriadLinkReader(paras);
//			List<TriadLink> triads = tlr.read();
			// travel the graph, random start
			// do not including the divergence end point in the path
			while (diGraph.isExistUnSelectedVertices()) {
//				index++;
//				if(index == 49)
//					logger.debug("index error " + index);
				Contig current = diGraph.getRandomVertex();
				// if the return conting is null and the
				// isExistUnSelectedVertices equal false then break;
				if (current == null)
					break;
				List<Contig> adjs = diGraph.getAdjVertices(current);
				if(adjs == null)
				{
					diGraph.setVertexAsSelected(current);
					continue;
				}
				int adjsSize = adjs.size();
				// random selected the frist 
				if(adjsSize == 0)
				{
					// orphan contig, only one element in path
					path = new NodePath();
					Node node = new Node();
					node.setCnt(current);
					node.setOrphan(true);
					path.push(node);
					paths.add(path);
					diGraph.setVertexAsSelected(current);
				} else if(adjsSize == 1)
				{
					// normal start point, always on the linear end point;
					path = new NodePath();
//					Contig next = diGraph.getNextVertex(current, null);
					Contig next = adjs.get(0);
					Contig startPoint = current;
					while(true)
					{   
						Node node = new Node();
						node.setCnt(current);
						node.setOrphan(false);
						diGraph.setVertexAsSelected(current);
						path.push(node);
						node = null;
						Contig previous = current;
						current = next;
						next = diGraph.getNextVertex(current, previous);
						if(next == null)
						{ // for the divergence point
							next = getTriadLinkNext(current, previous);
							// checking next is not null statement,
							// but the next will be not the adjacent vertex 
							if(next != null)
							{
								List<Contig> tAdjs = diGraph.getAdjVertices(current);
								if(!tAdjs.contains(next))
									next = null;
							}
							if(next == null){
								next = null;
								current = null;
								previous = null;
								break;
							} 
						} 
						if(next.equals(previous))
						{ // for the linear end point
							Node n = new Node();
							n.setCnt(current);
							n.setOrphan(false);
							diGraph.setVertexAsSelected(current);
							path.push(n);
							n = null;
							next = null;
							current = null;
							previous = null;
							break;
						}
					}
					paths.add(path);
				} else if(adjsSize == 2)
				{ // middle point;
					// normal start point, located in the linear path
					path = new NodePath();
					Node node = new Node();
					node.setCnt(current);
					path.push(node);
					diGraph.setVertexAsSelected(current);
					Contig startPoint = current;
					Contig c1 = adjs.get(0);
					Contig c2 = adjs.get(1);
					// c1--current--c2
					// directed by both direction;
					// for c1 direction; using c2 as previous point; checking whether c1 is valid point
					// all the element unshift into path;
					Contig previous = startPoint;
					current = c1;
					diGraph.setVertexAsSelected(current);
					Contig next = diGraph.getNextVertex(current, previous);
					while(true)
					{
						if(next != null)
						{
							if(next.equals(startPoint))
							{
								node = new Node();
								node.setCnt(current);
								path.unshift(node);
								diGraph.setVertexAsSelected(current);
								break;
							}
							if(next.equals(previous))
							{
								node = new Node();
								node.setCnt(current);
								path.unshift(node);
								diGraph.setVertexAsSelected(current);
								break;
							}
							node = new Node();
							node.setCnt(current);
							path.unshift(node);
							diGraph.setVertexAsSelected(current);
							previous = current;
							current = next;
							next = diGraph.getNextVertex(current, previous);
						} else
						{
							next = getTriadLinkNext(current, previous);
							// checking next is not null statement,
							// but the next will be not the adjacent vertex 
							if(next != null)
							{
								List<Contig> tAdjs = diGraph.getAdjVertices(current);
								if(!tAdjs.contains(next))
									next = null;
							}
							if(next == null)
							{
								diGraph.setVertexAsSelected(current);
								break;
							}
						}
					}
					// for c2 direction; using c1 as previous point; checking whether c2 is valid point;
					// all the valid element push into path;
					if(!diGraph.isExistUnSelectedVertices())
					{
						paths.add(path);
						break;
					}
					previous = startPoint;
					current = c2;
					next = diGraph.getNextVertex(current, previous);
					while(true)
					{
						if(next != null)
						{
							if(next.equals(startPoint))
							{
								node = new Node();
								node.setCnt(current);
								path.push(node);
								diGraph.setVertexAsSelected(current);
								break;
							}
							if(next.equals(previous))
							{
								node = new Node();
								node.setCnt(current);
								path.push(node);
								diGraph.setVertexAsSelected(current);
								break;
							}
							node = new Node();
							node.setCnt(current);
							path.push(node);
							diGraph.setVertexAsSelected(current);
							previous = current;
							current = next;
							next = diGraph.getNextVertex(current, previous);
						} else
						{
							next = getTriadLinkNext(current, previous);
							// checking next is not null statement,
							// but the next will be not the adjacent vertex 
							if(next != null)
							{
								List<Contig> tAdjs = diGraph.getAdjVertices(current);
								if(!tAdjs.contains(next))
									next = null;
							}
							if(next == null)
							{
								diGraph.setVertexAsSelected(current);
								break;
							}
						}
					}
					paths.add(path);
					
				} else if(adjsSize > 2)
				{ // divergence point
					diGraph.setVertexAsSelected(current);
					continue;
				}
			}
		} catch (Exception e) {
			logger.debug(this.getClass().getName() + "\t" + e.getMessage() + "\t" + e.getClass().getName());
			logger.error(this.getClass().getName() + "\t" + e.getMessage() + "\t" + e.getClass().getName());
		}
		// orientation contig in the paths;
		// define the first element in path is forward;
		try {
			for (NodePath np : paths) {
				int pathSize = np.getPathSize();
				if (pathSize == 1) {
					Node current = np.getElement(0);
					current.setStrand(Strand.FORWARD);
					current.setMeanDist2Next(0);
					current.setSdDist2Next(0);
					current.setSupportLinkNum(0);
				} else
				{
					for(int i = 0 ; i < pathSize - 1; i++)
					{
						Node current = np.getElement(i);
						Node next = np.getElement(i + 1);
						Contig cCnt = current.getCnt();
						Contig nCnt = next.getCnt();
						Edge e = null;
						List<Edge> es = diGraph.getEdgesInfo(current.getCnt(), next.getCnt());
						for(Edge t : es)
						{
							if(t.getOrigin().equals(cCnt) && t.getTerminus().equals(nCnt))
								e = t;
						}
						if(current.getStrand() == null)
							current.setStrand(e.getoStrand());
						if(next.getStrand() == null)
							next.setStrand(e.gettStrand());
						int meanSum = e.getDistMean();
						int sdSum = e.getDistSd();
						int slSum = e.getLinkNum();
						current.setMeanDist2Next(meanSum);
						current.setSdDist2Next(sdSum);
						current.setSupportLinkNum(slSum);
					}
				}
			}
		} catch (Exception e) {
			logger.debug(this.getClass().getName() + "\t" + e.getMessage() + "\t" + e.getClass().getName());
			logger.error(this.getClass().getName() + "\t" + e.getMessage() + "\t" + e.getClass().getName());
		}
		long end = System.currentTimeMillis();
		logger.info("Path Building, erase time: " + (end - start) + " ms");
		return paths;
	}

	@SuppressWarnings("unused")
	private List<NodePath> buildPath(List<Edge> edges, Parameter paras) {
		if (edges == null || edges.size() == 0)
			throw new IllegalArgumentException("PathBuilder: The Edges could not be empty!");
		List<NodePath> paths = new Vector<NodePath>();
		Graph diGraph = null;
		try {
			diGraph = new DirectedGraph(edges);
			// do transitive reduction
			diGraph.transitiveReducting();
			List<Edge> tempEdges = diGraph.getEdges();
			String edgeFile = paras.getOutFolder() + System.getProperty("file.separator") + "edges_after_tr.info";
			DotGraphFileWriter.writeEdge(edgeFile, tempEdges);
			// delete error prone edge
			diGraph.delErrorProneEdge(paras.getRatio());
			tempEdges = diGraph.getEdges();
			edgeFile = paras.getOutFolder() + System.getProperty("file.separator") + "edges_after_dep.info";
			DotGraphFileWriter.writeEdge(edgeFile, tempEdges);
			NodePath path = null;
			// travel the graph, random start
			// do not including the divergence end point in the path
			while (diGraph.isExistUnSelectedVertices()) {
				Contig current = diGraph.getRandomVertex();
				// Contig cnt = diGraph.getVertex("1413");
				// if the return conting is null and the
				// isExistUnSelectedVertices equal false then break;
				if (current == null)
					break;
				// checking the adjacent set of the specified contig;
				// if adjacent size == 0, it is orphan contig;
				// if adjacent size <= 2, it is normal start point;
				// if adjacent size >= 3, it is abnormal start point, next loop;
				List<Contig> adjSet = diGraph.getAdjVertices(current);
				int adjSetSize = adjSet.size();
				if (adjSetSize == 0) {
					// orphan contig, only one element in path
					path = new NodePath();
					Node node = new Node();
					node.setCnt(current);
					node.setOrphan(true);
					path.push(node);
					paths.add(path);
					diGraph.setVertexAsSelected(current);
				} else if (adjSetSize == 1) {
					// normal start point, always on the linear end point;
					path = new NodePath();
					Contig next = diGraph.getNextVertex(current, null);
					Contig startPoint = current;
					while(true)
					{   
						Node node = new Node();
						node.setCnt(current);
						node.setOrphan(false);
						diGraph.setVertexAsSelected(current);
						path.push(node);
						node = null;
						Contig previous = current;
						current = next;
						next = diGraph.getNextVertex(current, previous);
						if(next == null)
						{ // for the divergence point
							next = null;
							current = null;
							previous = null;
							break;
						} 
						if(next.equals(previous))
						{ // for the linear end point
							Node n = new Node();
							n.setCnt(current);
							n.setOrphan(false);
							diGraph.setVertexAsSelected(current);
							path.push(n);
							n = null;
							next = null;
							current = null;
							previous = null;
							break;
						}
					}
					paths.add(path);
//					while (next != null) {
//						Node node = new Node();
//						node.setCnt(cnt);
//						node.setOrphan(false);
//						diGraph.setVertexAsSelected(cnt);
//						path.push(node);
//						Contig temp = cnt;
//						cnt = next;
//						next = diGraph.getNextVertex(cnt, temp);
//						if (path.isNextExist(cnt, 0) && path.isNextExist(next, 0))
//						{
//							paths.add(path);
//							break;
//						}
//						if (next.equals(startPoint))
//							count = count + 1;
//					}
				} else if (adjSetSize == 2) {
					// normal start point, located in the linear path
					path = new NodePath();
					Node node = new Node();
					node.setCnt(current);
					path.push(node);
					diGraph.setVertexAsSelected(current);
					Contig startPoint = current;
					Contig c1 = adjSet.get(0);
					Contig c2 = adjSet.get(1);
					// c1--current--c2
					// directed by both direction;
					// for c1 direction; using c2 as previous point; checking whether c1 is valid point
					// all the element unshift into path;
					Contig previous = startPoint;
					current = c1;
					diGraph.setVertexAsSelected(current);
					Contig next = diGraph.getNextVertex(current, previous);
					while(true)
					{
						if(next != null)
						{
							if(next.equals(startPoint))
							{
								node = new Node();
								node.setCnt(current);
								path.unshift(node);
								diGraph.setVertexAsSelected(current);
								break;
							}
							if(next.equals(previous))
							{
								node = new Node();
								node.setCnt(current);
								path.unshift(node);
								diGraph.setVertexAsSelected(current);
								break;
							}
							node = new Node();
							node.setCnt(current);
							path.unshift(node);
							diGraph.setVertexAsSelected(current);
							previous = current;
							current = next;
							next = diGraph.getNextVertex(current, previous);
						} else
						{
							diGraph.setVertexAsSelected(current);
							break;
						}
					}
					// for c2 direction; using c1 as previous point; checking whether c2 is valid point;
					// all the valid element push into path;
					if(!diGraph.isExistUnSelectedVertices())
					{
						paths.add(path);
						break;
					}
					previous = startPoint;
					current = c2;
					next = diGraph.getNextVertex(current, previous);
					while(true)
					{
						if(next != null)
						{
							if(next.equals(startPoint))
							{
								node = new Node();
								node.setCnt(current);
								path.unshift(node);
								diGraph.setVertexAsSelected(current);
								break;
							}
							if(next.equals(previous))
							{
								node = new Node();
								node.setCnt(current);
								path.push(node);
								diGraph.setVertexAsSelected(current);
								break;
							}
							node = new Node();
							node.setCnt(current);
							path.push(node);
							diGraph.setVertexAsSelected(current);
							previous = current;
							current = next;
							next = diGraph.getNextVertex(current, previous);
						} else
						{
							diGraph.setVertexAsSelected(current);
							break;
						}
					}
					paths.add(path);
					
//					int count = 0;
//					boolean isReverse = false;
//					Contig next = diGraph.getNextVertex(current, null);
//					while (next != null) {
//						Node node = new Node();
//						node.setCnt(current);
//						node.setOrphan(false);
//						diGraph.setVertexAsSelected(current);
//						if (!path.isNextExist(current, 0)) {
//							if (isReverse)
//								path.unshift(node);
//							else
//								path.push(node);
//						}
//						Contig temp = current;
//						current = next;
//						next = diGraph.getNextVertex(current, temp);
//						if (next.equals(startPoint))
//							count += 1;
//						if (current.equals(startPoint))
//							isReverse = true;
//						if (current.equals(startPoint) && count == 2) {
//							paths.add(path);
//							break;
//						}
//					}
				} else if (adjSetSize >= 3) {
					diGraph.setVertexAsSelected(current);
					continue;
				}
			}
		} catch (Exception e) {
			logger.debug(this.getClass().getName() + "\t" + e.getMessage() + "\t" + e.getClass().getName());
			logger.error(this.getClass().getName() + "\t" + e.getMessage() + "\t" + e.getClass().getName());
		}
		// orientation contig in the paths;
		// define the first element in path is forward;
		Strand previousStrand = null;
		try {
			for (NodePath np : paths) {
				int pathSize = np.getPathSize();
				if (pathSize == 1) {
					Node current = np.getElement(0);
					current.setStrand(Strand.FORWARD);
					current.setMeanDist2Next(0);
					current.setSdDist2Next(0);
					current.setSupportLinkNum(0);
				} else
				{
					for(int i = 0 ; i < pathSize - 1; i++)
					{
						Node current = np.getElement(i);
						Node next = np.getElement(i + 1);
						Contig cCnt = current.getCnt();
						Contig nCnt = next.getCnt();
						Edge e = null;
						List<Edge> es = diGraph.getEdgesInfo(current.getCnt(), next.getCnt());
						for(Edge t : es)
						{
							if(t.getOrigin().equals(cCnt) && t.getTerminus().equals(nCnt))
								e = t;
						}
						if(current.getStrand() == null)
							current.setStrand(e.getoStrand());
						if(next.getStrand() == null)
							next.setStrand(e.gettStrand());
//						
//						
//						if (e.getOrigin().equals(current.getCnt())) {
//							current.setStrand(e.getoStrand());
//							next.setStrand(e.gettStrand());
//						} else {
//							// for the previous point
//							if (e.gettStrand().equals(Strand.FORWARD))
//								current.setStrand(Strand.REVERSE);
//							else
//								current.setStrand(Strand.FORWARD);
//							// for the following point
//							if (e.getoStrand().equals(Strand.FORWARD))
//								next.setStrand(Strand.REVERSE);
//							else
//								next.setStrand(Strand.FORWARD);
//						}
						// for the distance, sd and support links
						int meanSum = e.getDistMean();
						int sdSum = e.getDistSd();
						int slSum = e.getLinkNum();
						current.setMeanDist2Next(meanSum);
						current.setSdDist2Next(sdSum);
						current.setSupportLinkNum(slSum);
					}
				}
				
/*				else if (pathSize == 2) {
					Node current = np.getElement(0);
					Node next = np.getElement(1);
					List<Edge> eInfo = diGraph.getEdgesInfo(current.getCnt(), next.getCnt());
					Edge e = eInfo.get(0);
					if (e.getOrigin().equals(current.getCnt())) {
						current.setStrand(e.getoStrand());
						next.setStrand(e.gettStrand());
					} else {
						// for the previous point
						if (e.gettStrand().equals(Strand.FORWARD))
							current.setStrand(Strand.REVERSE);
						else
							current.setStrand(Strand.FORWARD);
						// for the following point
						if (e.getoStrand().equals(Strand.FORWARD))
							next.setStrand(Strand.REVERSE);
						else
							next.setStrand(Strand.FORWARD);
					}
					// for the distance, sd and support links
					int meanSum = e.getDistMean();
					int sdSum = e.getDistSd();
					int slSum = e.getLinkNum();
//					if (eInfo.size() == 1) {
//						meanSum = eInfo.get(0).getDistMean();
//						sdSum = eInfo.get(0).getDistSd();
//						slSum = eInfo.get(0).getLinkNum();
//					} else {
//						Edge e2 = eInfo.get(1);
//						meanSum = MathTool.mean(new Integer[] { e.getDistMean(), e2.getDistMean() });
//						sdSum = MathTool.mean(new Integer[] { e.getDistSd(), e2.getDistSd() });
//						slSum = e.getLinkNum() + e2.getLinkNum();
//					}
					current.setMeanDist2Next(meanSum);
					current.setSdDist2Next(sdSum);
					current.setSupportLinkNum(slSum);
				} else {
					for (int i = 0; i < pathSize - 1; i++) {
						Node previous = np.getElement(i);
						Node following = np.getElement(i + 1);
						List<Edge> eInfo = diGraph.getEdgesInfo(previous.getCnt(), following.getCnt());
						Edge e1 = eInfo.get(0);
						if (e1.getOrigin().equals(previous.getCnt())) {
							if (previous.getStrand() == null) {
								previous.setStrand(e1.getoStrand());
								following.setStrand(e1.gettStrand());
							} else {
								if (e1.getoStrand().equals(previous.getStrand())) {
									following.setStrand(e1.gettStrand());
								} else {
									// do some check, since the orientation is
									// wrong in path!
								}
							}
						} else {
							if (previous.getStrand() == null) {
								if (e1.gettStrand().equals(Strand.FORWARD))
									previous.setStrand(Strand.REVERSE);
								else
									previous.setStrand(Strand.FORWARD);
								if (e1.getoStrand().equals(Strand.FORWARD))
									following.setStrand(Strand.REVERSE);
								else
									following.setStrand(Strand.FORWARD);
							} else {
								if (!e1.gettStrand().equals(previous.getStrand())) {
									if (e1.getoStrand().equals(Strand.FORWARD))
										following.setStrand(Strand.REVERSE);
									else
										following.setStrand(Strand.FORWARD);
								} else {
									// do some check, since the orientation is
									// wrong in path!
									// it need to implemented!
								}
							}
						}

						// distance, sd and supported link
						int meanSum = 0;
						int sdSum = 0;
						int slSum = 0;
						if (eInfo.size() == 1) {
							meanSum = eInfo.get(0).getDistMean();
							sdSum = eInfo.get(0).getDistSd();
							slSum = eInfo.get(0).getLinkNum();
						} else {
							Edge e2 = eInfo.get(1);
							meanSum = MathTool.mean(new Integer[] { e1.getDistMean(), e2.getDistMean() });
							sdSum = MathTool.mean(new Integer[] { e1.getDistSd(), e2.getDistSd() });
							slSum = e1.getLinkNum() + e2.getLinkNum();
						}
						previous.setMeanDist2Next(meanSum);
						previous.setSdDist2Next(sdSum);
						previous.setSupportLinkNum(slSum);
					}
				}*/
			}
		} catch (Exception e) {
			logger.debug(this.getClass().getName() + "\t" + e.getMessage() + "\t" + e.getClass().getName());
			logger.error(this.getClass().getName() + "\t" + e.getMessage() + "\t" + e.getClass().getName());
		}
		return paths;
	}

	/*
	 * public List<Path> buildPath(List<Edge> edges, Parameter paras) { try {
	 * diGraph = new DiGraph(edges); Map<String, Integer> eStat =
	 * diGraph.getEdgesStatistics(); int lower =
	 * eStat.get("SUPPORT_LINKS_LOWER"); int upper =
	 * eStat.get("SUPPORT_LINKS_UPPER"); logger.debug("PathBuilder lower : " +
	 * lower); logger.debug("PathBuilder upper : " + upper); // original edges
	 * statistics for (String s : eStat.keySet()) { logger.debug(s + ":" +
	 * eStat.get(s)); } // remove outlier edges; diGraph.removeEdge(lower,
	 * upper); eStat = diGraph.getEdgesStatistics(); for (String s :
	 * eStat.keySet()) { logger.debug(s + ":" + eStat.get(s)); } // remove the
	 * edges by user specified value; lower = paras.getMinSupLinks(); upper =
	 * paras.getMaxSupLinks(); diGraph.removeEdge(lower, upper); eStat =
	 * diGraph.getEdgesStatistics(); for (String s : eStat.keySet()) {
	 * logger.debug(s + ":" + eStat.get(s)); } // pesudo edges;
	 * diGraph.addPesudoEdges();
	 * 
	 * // check each contig sorted indegree and outdegree List<ContInOut> values
	 * = diGraph.getCandVertices(); for (ContInOut c : values) {
	 * logger.debug(c.toString()); } diGraph = untangle(diGraph); // go through
	 * the graph // for random start; // String id =
	 * diGraph.getVertexByOrdering(); String id = diGraph.getOneRandomVertex();
	 * // DiGraph.selectedVertices.add(id); diGraph.addId2SletVerts(id);
	 * logger.debug("id: " + id); // String id = "1709"; String startId = id;
	 * boolean isReverse = false; List<Path> paths = new Vector<Path>(); Path
	 * path = new Path(); Strand strandStatus = Strand.FORWARD; // need to
	 * rethink implement for fast define valid edges;
	 * while(diGraph.getValidEdgeNums() >= 0) { List<Edge> pTEdges =
	 * diGraph.getPTValidEdgesById(id); if(pTEdges.size() == 0) { List<Edge>
	 * pFEdges = diGraph.getPFValidEdgesById(id); if(pFEdges.size() == 0) {
	 * diGraph.addId2SletVerts(id); id = diGraph.getOneRandomVertex(); continue;
	 * } else if(pFEdges.size() == 1) { Edge e = pFEdges.get(0); Contig o =
	 * e.getOrigin(); Contig t = e.getTerminus(); path.push(o); path.push(t);
	 * path.pushStrand(e.getoStrand()); path.pushStrand(e.gettStrand());
	 * 
	 * } else if(pFEdges.size() == 2) {
	 * 
	 * } else if(pFEdges.size() >2) {
	 * 
	 * }
	 * 
	 * } else if(pTEdges.size() == 1) {
	 * 
	 * } else if(pTEdges.size() == 2) {
	 * 
	 * } else if(pTEdges.size() > 2) {
	 * 
	 * } } // while (!diGraph.isEdgesEmpty()) { // List<Edge> pTEdges =
	 * diGraph.getEdgesBySpecifiedId(id); // // remove the reverse edge if in
	 * the path; // if (!path.isEmpty()) { // if (pTEdges.size() > 1) { // for
	 * (int i = 0; i < pTEdges.size(); i++) { // Edge e = pTEdges.get(i); // if
	 * (path.isExistReverseEdge(e)) { // pTEdges.remove(e); // continue; // } //
	 * if (path.isExistEdge(e)) { // pTEdges.remove(e); // continue; // } // }
	 * // } // } // if (pTEdges.isEmpty()) { // if (!path.isEmpty()) //
	 * paths.add(path); // path = new Path(); // // id =
	 * diGraph.getOneRandomVertex(); // id = diGraph.getVertexByOrdering(); //
	 * logger.debug("id: " + id); // if (id == null || id.length() == 0) //
	 * break; // DiGraph.selectedVertices.add(id); // if (isReverse) //
	 * isReverse = false; // } else { // // define the selected edge // Edge
	 * selectedE = null; // for (Edge e : pTEdges) { // if (selectedE == null) {
	 * // selectedE = e; // continue; // } else if (e.getLinkNum() >
	 * selectedE.getLinkNum() && e.getoStrand().equals(strandStatus)) { //
	 * selectedE = e; // continue; // } // } // // push or unshift vertex and
	 * Strand into Path // Contig origin = selectedE.getOrigin(); // Contig
	 * terminus = selectedE.getTerminus(); // // if (path.isEmpty()) { //
	 * path.push(origin); // path.pushStrand(selectedE.getoStrand()); //
	 * path.push(terminus); // path.pushStrand(selectedE.gettStrand()); // }
	 * else { // int oIndex = path.containVertex(origin); // int tIndex =
	 * path.containVertex(terminus); // if (oIndex == 0 && tIndex == -1) { // if
	 * (!isReverse) // isReverse = true; // path.unshift(terminus); //
	 * path.unshiftStrand(selectedE.gettStrand()); // } else if (oIndex ==
	 * path.getSize() - 1 && tIndex == -1) { // path.push(terminus); //
	 * path.pushStrand(selectedE.gettStrand()); // } // } // // storing the
	 * terminus Strand status; // strandStatus = selectedE.gettStrand(); // //
	 * remove edge in the digraph after push or unshift in the // // path //
	 * diGraph.removeEdge(origin.getID(), terminus.getID()); // id =
	 * terminus.getID(); // if (startId.equals(terminus.getID()) && isReverse) {
	 * // if (!path.isEmpty()) // paths.add(path); // path = new Path(); // //
	 * id = diGraph.getOneRandomVertex(); // id = diGraph.getVertexByOrdering();
	 * // logger.debug("id: " + id); // if (id == null || id.length() == 0) //
	 * break; // DiGraph.selectedVertices.add(id); // if (isReverse) //
	 * isReverse = false; // } // // if empty after remove edge, than put the
	 * path into paths // if (diGraph.isEdgesEmpty()) { // paths.add(path); //
	 * break; // } // } // } // logger.debug("Contain " + paths.size() +
	 * " Paths!"); int count = 0; for (Path p : paths) { logger.debug("Path " +
	 * count + ": " + p.toString()); count++; } return paths; } catch (Exception
	 * e) { logger.debug(e.getMessage()); logger.error(e.getMessage()); return
	 * null; } }
	 */

	// c1 for middle and c2 for two adjacent;
	public Contig getTriadLinkNext(Contig c1, Contig c2)
	{
		Contig next = null;
		if(triads == null)
		{
			TriadLinkReader tlr = new TriadLinkReader(paras);
			triads = tlr.read();
		}
		List<TriadLink> tls = new Vector<TriadLink>();
		for(TriadLink tl : triads)
		{
			Contig pre = tl.getPrevious();
			Contig mid = tl.getMiddle();
			Contig lst = tl.getLast();
			if(!c1.equals(mid))
				continue;
			if(c2.equals(pre) || c2.equals(lst))
			{
				tls.add(tl);
			}
		}
		if(tls.isEmpty())
			return null;
		TriadLinkComparator tlc = new TriadLinkComparator();
		Collections.sort(tls, tlc);
		TriadLink tl = tls.get(tls.size() - 1);
		if(tl.getPrevious().equals(c2))
			next = tl.getLast();
		else if(tl.getLast().equals(c2))
			next = tl.getPrevious();
		triads.remove(tl);
		return next;
	}

}

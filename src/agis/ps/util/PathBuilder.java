/*
*File: agis.ps.util.PathBuilder.java
*User: mqin
*Email: mqin@ymail.com
*Date: 2016年1月14日
*/
package agis.ps.util;

import java.util.Collection;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Vector;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import agis.ps.DiGraph;
import agis.ps.Edge;
import agis.ps.Path;
import agis.ps.link.ContInOut;
import agis.ps.link.Contig;

public class PathBuilder {
	public static Logger logger = LoggerFactory.getLogger(PathBuilder.class);
	private List<Edge> edges;
	private Parameter paras;
	private DiGraph diGraph;
	
	public PathBuilder()
	{
		// do nothing;
	}
	
	public PathBuilder(List<Edge> edges, Parameter paras)
	{
		this.edges = edges;
		this.paras = paras;
	}

	public List<Path> buildEulerPath(List<Edge> edges) {
		return null;
	}
	
	public List<Path> buildHamiltonPath()
	{
		if(edges == null || edges.size() == 0)
			throw new IllegalArgumentException("PathBuilder : The Edges could not be empty!");
		return this.buildHamiltonPath(edges);
	}
	
	public List<Path> buildHamiltonPath(List<Edge> edges) {
		try {
			DiGraph diGraph = new DiGraph(edges);
			// logger.debug("Vertices Num: " + diGraph.getVerNum());
			// logger.debug("Edges Num: " + diGraph.getEdgNum());
			// logger.debug("Indegree...");
			// Map<String, Integer> indegrees = diGraph.indegrees();
			// for(String s : indegrees.keySet())
			// {
			// logger.debug(s + ":" + indegrees.get(s));
			// }
			// logger.debug("Outdegree...");
			// Map<String, Integer> outdegrees = diGraph.outdegrees();
			// for(String s : outdegrees.keySet())
			// {
			// logger.debug(s + ":" + outdegrees.get(s));
			// }
			// logger.debug("smallest indegere...");
			// Map<String, Integer> sIns = diGraph.minIndegeres();
			// for(String s : sIns.keySet())
			// {
			// logger.debug(s + ":" + sIns.get(s));
			// }
			// statistics of edges info
			Map<String, Integer> eStat = diGraph.getEdgesStatistics();
			int lower = eStat.get("SUPPORT_LINKS_LOWER");
			int upper = eStat.get("SUPPORT_LINKS_UPPER");
			logger.debug("lower : " + lower);
			logger.debug("upper : " + upper);
			for (String s : eStat.keySet()) {
				logger.debug(s + ":" + eStat.get(s));
			}
			// simplest graph by remove lower and upper links number
			// edges former;
			// logger.debug("Before remove");
			// List<Edge> oE = diGraph.getEdges();
			// for (Edge e : oE) {
			// logger.debug("B: " + e.getOrigin().getID() + "->" +
			// e.getTerminus().getID() + " : " + e.getLinkNum());
			// }
			diGraph.removeEdge(lower, upper);
			// logger.debug("After process");
			// List<Edge> aE = diGraph.getEdges();
			// for (Edge e : aE) {
			// logger.debug("A: " + e.getOrigin().getID() + "->" +
			// e.getTerminus().getID() + " : " + e.getLinkNum());
			// }
			eStat = diGraph.getEdgesStatistics();
			for (String s : eStat.keySet()) {
				logger.debug(s + ":" + eStat.get(s));
			}
			// remove the support links less than 10
			diGraph.removeEdge(1);
			// logger.debug("After process less than 4");
			// aE = diGraph.getEdges();
			// for (Edge e : aE) {
			// logger.debug("A: " + e.getOrigin().getID() + "->" +
			// e.getTerminus().getID() + " : " + e.getLinkNum());
			// }
			eStat = diGraph.getEdgesStatistics();
			for (String s : eStat.keySet()) {
				logger.debug(s + ":" + eStat.get(s));
			}
			logger.debug(System.getProperty("user.dir"));
			DotGraphFileWriter.writeEdge(
					System.getProperty("user.dir") + System.getProperty("file.separator") + "removeLessThan10.txt",
					diGraph.getEdges());
			// check each contig sorted indegree and outdegree
			List<ContInOut> values = diGraph.getCandVertices();
			for (ContInOut c : values) {
				logger.debug(c.toString());
			}
			// go go go through the graph
			// for random start;
			String id = diGraph.getVertexByOrdering();
//			DiGraph.selectedVertices.add(id);
			diGraph.addId2SletVerts(id);
			logger.debug("id: " + id);
			// String id = "1709";
			String startId = id;
			boolean isReverse = false;
			List<Path> paths = new Vector<Path>();
			Path path = new Path();
			Strand strandStatus = Strand.FORWARD;
			while (!diGraph.isEdgesEmpty()) {
				List<Edge> pTEdges = diGraph.getEdgesBySpecifiedId(id);
				// remove the reverse edge if in the path;
				if (!path.isEmpty()) {
					if (pTEdges.size() > 1) {
						for (int i = 0; i < pTEdges.size(); i++) {
							Edge e = pTEdges.get(i);
							if (path.isExistReverseEdge(e)) {
								pTEdges.remove(e);
								continue;
							}
							if (path.isExistEdge(e)) {
								pTEdges.remove(e);
								continue;
							}
						}
					}
				}
				if (pTEdges.isEmpty()) {
					if (!path.isEmpty())
						paths.add(path);
					path = new Path();
					// id = diGraph.getOneRandomVertex();
					id = diGraph.getVertexByOrdering();
					logger.debug("id: " + id);
					if (id == null || id.length() == 0)
						break;
//					DiGraph.selectedVertices.add(id);
					diGraph.addId2SletVerts(id);
					if (isReverse)
						isReverse = false;
				} else {
					// define the selected edge
					Edge selectedE = null;
					for (Edge e : pTEdges) {
						if (selectedE == null) {
							selectedE = e;
							continue;
						} else if (e.getLinkNum() > selectedE.getLinkNum() && e.getoStrand().equals(strandStatus)) {
							selectedE = e;
							continue;
						}
					}
					// push or unshift vertex and Strand into Path
					Contig origin = selectedE.getOrigin();
					Contig terminus = selectedE.getTerminus();

					if (path.isEmpty()) {
						path.push(origin);
						path.pushStrand(selectedE.getoStrand());
						path.push(terminus);
						path.pushStrand(selectedE.gettStrand());
					} else {
						int oIndex = path.containVertex(origin);
						int tIndex = path.containVertex(terminus);
						if (oIndex == 0 && tIndex == -1) {
							if (!isReverse)
								isReverse = true;
							path.unshift(terminus);
							path.unshiftStrand(selectedE.gettStrand());
						} else if (oIndex == path.getSize() - 1 && tIndex == -1) {
							path.push(terminus);
							path.pushStrand(selectedE.gettStrand());
						}
					}
					// storing the terminus Strand status;
					strandStatus = selectedE.gettStrand();
					// remove edge in the digraph after push or unshift in the
					// path
					diGraph.removeEdge(origin.getID(), terminus.getID());
					id = terminus.getID();
					if (startId.equals(terminus.getID()) && isReverse) {
						if (!path.isEmpty())
							paths.add(path);
						path = new Path();
						// id = diGraph.getOneRandomVertex();
						id = diGraph.getVertexByOrdering();
						logger.debug("id: " + id);
						if (id == null || id.length() == 0)
							break;
//						DiGraph.selectedVertices.add(id);
						diGraph.addId2SletVerts(id);
						if (isReverse)
							isReverse = false;
					}
					// if empty after remove edge, than put the path into paths
					if (diGraph.isEdgesEmpty()) {
						paths.add(path);
						break;
					}
				}
			}

			logger.debug("Contain " + paths.size() + " Paths!");
			int count = 0;
			for (Path p : paths) {
				logger.debug("Path " + count + ": " + p.toString());
				count++;
			}
			DotGraphFileWriter.writePaths(
					System.getProperty("user.dir") + System.getProperty("file.separator") + "paths.txt", paths);
			return null;
		} catch (Exception e) {
			logger.debug(e.getMessage());
			logger.error(e.getMessage());
			return null;
		}
	}

	public List<Path> buildPath()
	{
		if(edges == null || edges.size() == 0)
			throw new IllegalArgumentException("PathBuilder ： The Edges could not be empty!");
		return this.buildPath(edges, paras);
	}
	
	public List<Path> buildPath(List<Edge> edges, Parameter paras) {
		try {
			diGraph = new DiGraph(edges);
			Map<String, Integer> eStat = diGraph.getEdgesStatistics();
			int lower = eStat.get("SUPPORT_LINKS_LOWER");
			int upper = eStat.get("SUPPORT_LINKS_UPPER");
			logger.debug("PathBuilder lower : " + lower);
			logger.debug("PathBuilder upper : " + upper);
			// original edges statistics
			for (String s : eStat.keySet()) {
				logger.debug(s + ":" + eStat.get(s));
			}
			// remove outlier edges;
			diGraph.removeEdge(lower, upper);
			eStat = diGraph.getEdgesStatistics();
			for (String s : eStat.keySet()) {
				logger.debug(s + ":" + eStat.get(s));
			}
			// remove the edges by user specified value;
			lower = paras.getMinSupLinks();
			upper = paras.getMaxSupLinks();
			diGraph.removeEdge(lower, upper);
			eStat = diGraph.getEdgesStatistics();
			for (String s : eStat.keySet()) {
				logger.debug(s + ":" + eStat.get(s));
			}
			// pesudo edges;
			diGraph.addPesudoEdges();
			
			// check each contig sorted indegree and outdegree
			List<ContInOut> values = diGraph.getCandVertices();
			for (ContInOut c : values) {
				logger.debug(c.toString());
			}
			diGraph = untangle(diGraph);
			// go through the graph
			// for random start;
			// String id = diGraph.getVertexByOrdering();
			String id = diGraph.getOneRandomVertex();
//			DiGraph.selectedVertices.add(id);
			diGraph.addId2SletVerts(id);
			logger.debug("id: " + id);
			// String id = "1709";
			String startId = id;
			boolean isReverse = false;
			List<Path> paths = new Vector<Path>();
			Path path = new Path();
			Strand strandStatus = Strand.FORWARD;
			// need to rethink implement for fast define valid edges;
			while(diGraph.getValidEdgeNums() >= 0)
			{
				List<Edge> pTEdges = diGraph.getPTValidEdgesById(id);
				if(pTEdges.size() == 0)
				{
					List<Edge> pFEdges = diGraph.getPFValidEdgesById(id);
					if(pFEdges.size() == 0)
					{
						diGraph.addId2SletVerts(id);
						id = diGraph.getOneRandomVertex();
						continue;
					} else if(pFEdges.size() ==  1)
					{
						Edge e = pFEdges.get(0);
						Contig o = e.getOrigin();
						Contig t = e.getTerminus();
						path.push(o);
						path.push(t);
						path.pushStrand(e.getoStrand());
						path.pushStrand(e.gettStrand());
						
					} else if(pFEdges.size() == 2)
					{
						
					} else if(pFEdges.size() >2)
					{
						
					}
						
				} else if(pTEdges.size() == 1)
				{
					
				} else if(pTEdges.size() == 2)
				{
					
				} else if(pTEdges.size() > 2)
				{
					
				}
			}
//			while (!diGraph.isEdgesEmpty()) {
//				List<Edge> pTEdges = diGraph.getEdgesBySpecifiedId(id);
//				// remove the reverse edge if in the path;
//				if (!path.isEmpty()) {
//					if (pTEdges.size() > 1) {
//						for (int i = 0; i < pTEdges.size(); i++) {
//							Edge e = pTEdges.get(i);
//							if (path.isExistReverseEdge(e)) {
//								pTEdges.remove(e);
//								continue;
//							}
//							if (path.isExistEdge(e)) {
//								pTEdges.remove(e);
//								continue;
//							}
//						}
//					}
//				}
//				if (pTEdges.isEmpty()) {
//					if (!path.isEmpty())
//						paths.add(path);
//					path = new Path();
//					// id = diGraph.getOneRandomVertex();
//					id = diGraph.getVertexByOrdering();
//					logger.debug("id: " + id);
//					if (id == null || id.length() == 0)
//						break;
//					DiGraph.selectedVertices.add(id);
//					if (isReverse)
//						isReverse = false;
//				} else {
//					// define the selected edge
//					Edge selectedE = null;
//					for (Edge e : pTEdges) {
//						if (selectedE == null) {
//							selectedE = e;
//							continue;
//						} else if (e.getLinkNum() > selectedE.getLinkNum() && e.getoStrand().equals(strandStatus)) {
//							selectedE = e;
//							continue;
//						}
//					}
//					// push or unshift vertex and Strand into Path
//					Contig origin = selectedE.getOrigin();
//					Contig terminus = selectedE.getTerminus();
//
//					if (path.isEmpty()) {
//						path.push(origin);
//						path.pushStrand(selectedE.getoStrand());
//						path.push(terminus);
//						path.pushStrand(selectedE.gettStrand());
//					} else {
//						int oIndex = path.containVertex(origin);
//						int tIndex = path.containVertex(terminus);
//						if (oIndex == 0 && tIndex == -1) {
//							if (!isReverse)
//								isReverse = true;
//							path.unshift(terminus);
//							path.unshiftStrand(selectedE.gettStrand());
//						} else if (oIndex == path.getSize() - 1 && tIndex == -1) {
//							path.push(terminus);
//							path.pushStrand(selectedE.gettStrand());
//						}
//					}
//					// storing the terminus Strand status;
//					strandStatus = selectedE.gettStrand();
//					// remove edge in the digraph after push or unshift in the
//					// path
//					diGraph.removeEdge(origin.getID(), terminus.getID());
//					id = terminus.getID();
//					if (startId.equals(terminus.getID()) && isReverse) {
//						if (!path.isEmpty())
//							paths.add(path);
//						path = new Path();
//						// id = diGraph.getOneRandomVertex();
//						id = diGraph.getVertexByOrdering();
//						logger.debug("id: " + id);
//						if (id == null || id.length() == 0)
//							break;
//						DiGraph.selectedVertices.add(id);
//						if (isReverse)
//							isReverse = false;
//					}
//					// if empty after remove edge, than put the path into paths
//					if (diGraph.isEdgesEmpty()) {
//						paths.add(path);
//						break;
//					}
//				}
//			}
//
			logger.debug("Contain " + paths.size() + " Paths!");
			int count = 0;
			for (Path p : paths) {
				logger.debug("Path " + count + ": " + p.toString());
				count++;
			}
			return paths;
		} catch (Exception e) {
			logger.debug(e.getMessage());
			logger.error(e.getMessage());
			return null;
		}
	}

	// untangling the digraph, remove the self-connected, parallel links, repeat
	// links;
	public static DiGraph untangle(DiGraph diGraph) {
		// remove self connected
		diGraph = removeSelfConnected(diGraph);
		return diGraph;
	}
	
	// remove the self connected edges;
	public static DiGraph removeSelfConnected(DiGraph diGraph) {
		// finding the self-connected edges;
		while (diGraph.hasSelfConnected()) {
			List<Edge> edges = diGraph.getSelfConnected();
			if (edges == null) {
				logger.debug("PathBuilder : There are no self connected edges in the DiGraph!");
				break;
			}
			for (int i = 0; i < edges.size(); i++) {
				Edge edge = edges.get(i);
				diGraph.removeEdge(edge.getOrigin().getID(), edge.getTerminus().getID());
			}
		}
		// finding the parallel edges

		return diGraph;
	}
	
	// remove tips 

	public DiGraph getDiGraph() {
		return diGraph;
	}

	public void setDiGraph(DiGraph diGraph) {
		this.diGraph = diGraph;
	}

	
}

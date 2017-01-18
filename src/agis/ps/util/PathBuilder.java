/*
*File: agis.ps.util.PathBuilder.java
*User: mqin
*Email: mqin@ymail.com
*Date: 2016年1月14日
*/
package agis.ps.util;

import java.util.Collections;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Vector;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import agis.ps.file.DotGraphFileWriter;
import agis.ps.file.TriadLinkReader;
import agis.ps.graph.DirectedGraph;
import agis.ps.graph.Graph;
import agis.ps.link.CntFileEncapsulate;
import agis.ps.link.Edge;
import agis.ps.link.TriadLink;
//import agis.ps.link.TriadLinkComparator;
import agis.ps.path.InternalNode2;
import agis.ps.path.InternalNode;
import agis.ps.path.InternalPath;
import agis.ps.path.InternalPathComparator;
import agis.ps.path.Node;
import agis.ps.path.NodePath;
import agis.ps.seqs.Contig;

public class PathBuilder {
	public static Logger logger = LoggerFactory.getLogger(PathBuilder.class);
	private static int MAXIMUM_INTERNAL_LENGTH = 10000; // 5000 bp for validating
	private static int INTERNAL_LENGTH = 0; // store the internal length;
	@SuppressWarnings("unused")
	private static int INDEX = 0;
	private List<Edge> edges;
	private Parameter paras;
	private List<TriadLink> triads;
	private Graph diGraph;
	private NodePath path;
//	private CntFileEncapsulate cntfile;
	private Map<String, Integer> cntLens;

	public PathBuilder(List<Edge> edges, Parameter paras, CntFileEncapsulate cntfile) {
		this.edges = edges;
		this.paras = paras;
//		this.cntfile = cntfile;
		TriadLinkReader tlr = new TriadLinkReader(paras);
		triads = tlr.read();
	}
	
	public PathBuilder(List<Edge> edges, Parameter paras, Map<String, Integer> cntLens)
	{
		this.edges = edges;
		this.paras = paras;
		this.cntLens = cntLens;
	}

	public List<NodePath> buildPath() {
		if (edges == null || edges.size() == 0)
			throw new IllegalArgumentException("PathBuilder ： The Edges could not be empty!");
		// return this.buildPath2(edges, paras);
//		return this.buildPathByCntFileEncapsulate2();
		return buildPath3();
	}
	
	private List<NodePath> buildPath3()
	{
		long start = System.currentTimeMillis();
		if (edges == null || edges.size() == 0)
			throw new IllegalArgumentException("PathBuilder: The Edges could not be empty!");
		List<NodePath> paths = new Vector<NodePath>();
		diGraph = null;
		try {
			logger.info("Original Edges size: " + edges.size());
			String edgeFile = paras.getOutFolder() + System.getProperty("file.separator") + "edges.info";
			DotGraphFileWriter.writeEdge(edgeFile, edges);
//			diGraph = new DirectedGraph(edges, paras, cntfile);
			diGraph = new DirectedGraph(edges, paras, cntLens);
			// do transitive reduction
			diGraph.transitiveReducting();
			List<Edge> tempEdges = diGraph.getEdges();
			logger.info("Edges size after transitive reducing: " + tempEdges.size());
			edgeFile = paras.getOutFolder() + System.getProperty("file.separator") + "edges_after_tr.info";
			DotGraphFileWriter.writeEdge(edgeFile, tempEdges);
			// delete error prone edge
			diGraph.delErrorProneEdge(paras.getRatio());
			tempEdges = diGraph.getEdges();
			logger.info("Edges size after error prone deleting: " + tempEdges.size());
			edgeFile = paras.getOutFolder() + System.getProperty("file.separator") + "edges_after_dep.info";
			DotGraphFileWriter.writeEdge(edgeFile, tempEdges);
			// delete tips
			diGraph.delTips();
			tempEdges = diGraph.getEdges();
			logger.info("Edges size after deleting tips: " + tempEdges.size());
			edgeFile = paras.getOutFolder() + System.getProperty("file.separator") + "edges_after_dt.info";
			DotGraphFileWriter.writeEdge(edgeFile, tempEdges);
			// delete similarity contigs error prone edges;
			// diGraph.delSimCntEdges();
			// tempEdges = diGraph.getEdges();
			// logger.info("Edges size after deleting similarity contigs error
			// prone edges : " + tempEdges.size());
			// edgeFile = paras.getOutFolder() +
			// System.getProperty("file.separator") + "edges_after_dsc.info";
			// DotGraphFileWriter.writeEdge(edgeFile, tempEdges);
			tempEdges = null;
			path = null;
			// travel the graph, random start
			// do not including the divergence end point in the path
			while (diGraph.isExistUnSelectedVertices()) {
//				INDEX++;
//				logger.debug(""+INDEX);
//				if(INDEX == 700)
//					logger.debug("breakpoint");
				Contig current = diGraph.getRandomVertex();
				if (current == null)
					break;
				List<Contig> adjs = diGraph.getAdjVertices(current);
				if (adjs == null) {
					if (this.isValidCnt(current)) {
						path = new NodePath();
						this.addNode(null, current, null, null, null, path, true);
						paths.add(path);
					} else {
						diGraph.setVertexAsSelected(current);
					}
					continue;
				}
				int adjsSize = adjs.size();
				if (adjsSize == 0) {
					// orphan contig, only one element in path
					if (this.isValidCnt(current)) {
						path = new NodePath();
						this.addNode(null, current, null, null, null, path, true);
						paths.add(path);
					} else {
						diGraph.setVertexAsSelected(current);
					}
				} else if (adjsSize == 1) {
					// normal start point, always on the linear end point;
					path = new NodePath();
					Contig next = adjs.get(0);
					this.travelgraph(null, current, next, path, true);
					paths.add(path);
				} else if (adjsSize == 2) {
					// normal start point, located in the linear path
					path = new NodePath();
					Contig c1 = adjs.get(0);
					Contig c2 = adjs.get(1);
					// adding the c2 into path first
					this.addNode(null, c2, null, null, null,  path, false);
					// travel backward 
					this.travelgraph(c2, current, c1, path, false);
					// travel forward
					int count = 0;
					while(path.hasLast())
					{
						if(count > 1)
							break;
						Node node = path.pop();
						diGraph.setVertextUnselected(node.getCnt());
						count++;
					}
					// checking whether is cicle case;
//					Contig last = null;
//					if(path.getPathSize() != 0)
//						last = path.getElement(path.getPathSize() - 1).getCnt();
//					
//
//					if(!c2.equals(last))
//					{
						this.travelgraph(c1, current, c2, path, true);
//					}
					// orientation conflict three point;
					if(path.getPathSize() == 1)
					{
						// for c1 
						path = new NodePath();
						this.addNode(null, c1, null, null, null,  path, false);
						paths.add(path);
						// for current
						path = new NodePath();
						this.addNode(null, current, null, null, null,  path, false);
						paths.add(path);
						// for c2
						path = new NodePath();
						this.addNode(null, c2, null, null, null,  path, false);
						paths.add(path);
					} else
					{
						paths.add(path);
					}
				} else if (adjsSize > 2) { 
					// divergence point, do not considering
					diGraph.setVertexAsSelected(current);
					continue;
				}
			}
		} catch (Exception e) {
			// logger.debug("INDEX = " + INDEX);
			logger.error(this.getClass().getName() + "\t" + e.getMessage() + "\t" + e.getClass().getName());
		}
		// orientation contig in the paths;
		// define the first element in path is forward;
		// require to check the orientation conflict!
		// cleaning and call gc
		this.diGraph = null;
		this.edges = null;
		this.triads = null;
		System.gc();
		long end = System.currentTimeMillis();
		logger.info("Path Building, erase time: " + (end - start) + " ms");
		return paths;
	}

//	private List<NodePath> buildPathByCntFileEncapsulate2() {
//		long start = System.currentTimeMillis();
//		if (edges == null || edges.size() == 0)
//			throw new IllegalArgumentException("PathBuilder: The Edges could not be empty!");
//		List<NodePath> paths = new Vector<NodePath>();
//		diGraph = null;
//		try {
//			logger.info("Original Edges size: " + edges.size());
//			String edgeFile = paras.getOutFolder() + System.getProperty("file.separator") + "edges.info";
//			DotGraphFileWriter.writeEdge(edgeFile, edges);
//			diGraph = new DirectedGraph(edges, paras, cntfile);
//			// do transitive reduction
//			diGraph.transitiveReducting();
//			List<Edge> tempEdges = diGraph.getEdges();
//			logger.info("Edges size after transitive reducing: " + tempEdges.size());
//			edgeFile = paras.getOutFolder() + System.getProperty("file.separator") + "edges_after_tr.info";
//			DotGraphFileWriter.writeEdge(edgeFile, tempEdges);
//			// delete error prone edge
//			diGraph.delErrorProneEdge(paras.getRatio());
//			tempEdges = diGraph.getEdges();
//			logger.info("Edges size after error prone deleting: " + tempEdges.size());
//			edgeFile = paras.getOutFolder() + System.getProperty("file.separator") + "edges_after_dep.info";
//			DotGraphFileWriter.writeEdge(edgeFile, tempEdges);
//			// delete tips
//			diGraph.delTips();
//			tempEdges = diGraph.getEdges();
//			logger.info("Edges size after deleting tips: " + tempEdges.size());
//			edgeFile = paras.getOutFolder() + System.getProperty("file.separator") + "edges_after_dt.info";
//			DotGraphFileWriter.writeEdge(edgeFile, tempEdges);
//			// delete similarity contigs error prone edges;
//			// diGraph.delSimCntEdges();
//			// tempEdges = diGraph.getEdges();
//			// logger.info("Edges size after deleting similarity contigs error
//			// prone edges : " + tempEdges.size());
//			// edgeFile = paras.getOutFolder() +
//			// System.getProperty("file.separator") + "edges_after_dsc.info";
//			// DotGraphFileWriter.writeEdge(edgeFile, tempEdges);
//			tempEdges = null;
//			path = null;
//			// travel the graph, random start
//			// do not including the divergence end point in the path
//			while (diGraph.isExistUnSelectedVertices()) {
////				INDEX++;
////				if(INDEX == 7)
////					logger.debug("breakpoint");
//				Contig current = diGraph.getRandomVertex();
//				if (current == null)
//					break;
//				List<Contig> adjs = diGraph.getAdjVertices(current);
//				if (adjs == null) {
//					if (this.isValidCnt(current)) {
//						path = new NodePath();
//						this.addNode2(null, current, null, null, null, path, true);
//						paths.add(path);
//					} else {
//						diGraph.setVertexAsSelected(current);
//					}
//					continue;
//				}
//				int adjsSize = adjs.size();
//				if (adjsSize == 0) {
//					// orphan contig, only one element in path
//					if (this.isValidCnt(current)) {
//						path = new NodePath();
//						this.addNode2(null, current, null, null, null, path, true);
//						paths.add(path);
//					} else {
//						diGraph.setVertexAsSelected(current);
//					}
//				} else if (adjsSize == 1) {
//					// normal start point, always on the linear end point;
//					path = new NodePath();
//					Contig next = adjs.get(0);
//					this.travelgraph(null, current, next, path, true);
//					paths.add(path);
//				} else if (adjsSize == 2) {
//					// normal start point, located in the linear path
//					path = new NodePath();
//					Contig c1 = adjs.get(0);
//					Contig c2 = adjs.get(1);
//					// adding the c2 into path first
//					this.addNode2(null, c2, null, null, null,  path, false);
//					// travel backward 
//					this.travelgraph(c2, current, c1, path, false);
//					// travel forward
//					int count = 0;
//					while(path.hasLast())
//					{
//						if(count > 1)
//							break;
//						Node node = path.pop();
//						diGraph.setVertextUnselected(node.getCnt());
//						count++;
//					}
//					// checking whether is cicle case;
////					Contig last = null;
////					if(path.getPathSize() != 0)
////						last = path.getElement(path.getPathSize() - 1).getCnt();
////					
////
////					if(!c2.equals(last))
////					{
//						this.travelgraph(c1, current, c2, path, true);
////					}
//					paths.add(path);
//				} else if (adjsSize > 2) { 
//					// divergence point, do not considering
//					diGraph.setVertexAsSelected(current);
//					continue;
//				}
//			}
//		} catch (Exception e) {
//			// logger.debug("INDEX = " + INDEX);
//			logger.error(this.getClass().getName() + "\t" + e.getMessage() + "\t" + e.getClass().getName());
//		}
//		// orientation contig in the paths;
//		// define the first element in path is forward;
//		// require to check the orientation conflict!
//		// cleaning and call gc
//		this.diGraph = null;
//		this.edges = null;
//		this.triads = null;
//		System.gc();
//		long end = System.currentTimeMillis();
//		logger.info("Path Building, erase time: " + (end - start) + " ms");
//		return paths;
//	}
	
	private void travelgraph(Contig previous, Contig current, Contig next, NodePath path, boolean isForward)
	{
		Contig startpoint = current;
		List<Edge> p2ces = null; // previous to current edges;
		Edge p2ce = null; // previous to current edge;
		List<Edge> c2nes = null; // current to next edges;
		Edge c2ne = null; // current to next edge;
		if(previous != null)
		{
			p2ces = diGraph.getEdgesInfo(previous, current);
			for(Edge e : p2ces)
			{
				if(e.getOrigin().equals(previous) && e.getTerminus().equals(current))
				{
					p2ce = e;
					break;
				}
			}
		}
		// get current to next edges;
		c2nes = diGraph.getEdgesInfo(current, next);
		for(Edge e : c2nes)
		{
			if(e.getOrigin().equals(current) && e.getTerminus().equals(next))
			{
				c2ne = e;
				break;
			}
		}
		// try to extend the path;
		while (true) {
			if (this.validateOrientation(previous, current, next, p2ce, c2ne)) {
				if (this.isValidCnt(current)) {
					this.addNode(previous, current, next, p2ce, c2ne, path, isForward);
				} else {
					break;
				}

			} else
			{
				next = null;
//				if (this.isValidCnt(current)) {
					this.addNode(previous, current, next, p2ce, c2ne, path, isForward);
//				} 
				break;
			}
			previous = current;
			current = next;
			if(current.getID().equals("1317"))
				logger.debug("breakpoint");
			next = diGraph.getNextVertex(current, previous);
			// get previous to current edges 
			p2ces = c2nes;
			p2ce = c2ne;
			// for the divergence point which could not determine
			// how-to get next
			if (next == null) {
				LinkedList<Contig> formers = new LinkedList<Contig>();
				int num = 0;
				// only considering two former contigs
				if(isForward)
				{
					while(true)
					{
						Node node = path.getElement(path.getPathSize() - (num + 2));
						if(node != null)
							formers.addLast(node.getCnt());
						else
							break;
						num++;
						if(num == 2)
							break;
					}
				} else
				{
					while(true)
					{
						Node node = path.getElement((num + 1));
						if(node != null)
							formers.addLast(node.getCnt());
						else
							break;
						num++;
						if(num == 2)
							break;
					}
				}
				// List<Contig> selectedCnts =
				// this.getTriadLinkNext3(current, previous);
				List<Contig> selectedCnts = this.getInternalPath(current, previous, formers);
				if (selectedCnts == null) {
					next = null;
					this.addNode(previous, current, next, p2ce, c2ne, path, isForward);
					break;
				}
				int size = selectedCnts.size();
				if (size == 1) {
					next = selectedCnts.get(0);
				} else
				{
					int index = 0;
					boolean isConflict = false;
					do{
						next = selectedCnts.get(index);
						// get current to next edges;
						c2nes = diGraph.getEdgesInfo(current, next);
						for(Edge e : c2nes)
						{
							if(e.getOrigin().equals(current) && e.getTerminus().equals(next))
							{
								c2ne = e;
								break;
							}
						}
						if (this.validateOrientation(previous, current, next, p2ce, c2ne)) {
							// because internal path will segment duplication;
							// so validate contig is travel again which is not divergence point
							// will not apply into this case;
//							if (this.isValidCnt(current)) {
							this.addNode(previous, current, next, p2ce, c2ne, path, isForward);
//							} else {
//								isConflict = true;
//								break;
//							}
						} else
						{
							next = null;
//							if (this.isValidCnt(current)) {
							this.addNode(previous, current, next, p2ce, c2ne, path, isForward);
//							} 
							isConflict = true;
							break;
						}
						p2ces = c2nes;
						p2ce = c2ne;
						previous = current;
						current = next;
						index++;
					} while(index <= (size - 2));
					// if exist orientation or other conflict 
					if(isConflict)
						break;
					next = selectedCnts.get(index);
				}
			} 
			// get current to next edges;
			c2nes = diGraph.getEdgesInfo(current, next);
			for(Edge e : c2nes)
			{
				if(e.getOrigin().equals(current) && e.getTerminus().equals(next))
				{
					c2ne = e;
					break;
				}
			}
			// check whether is reverse to path
			if (next.equals(previous)) { // for the linear end point
				next = null;
				this.addNode(previous, current, next, p2ce, c2ne, path, isForward);
				break;
			}
			// for the circle case 
			if (next.equals(startpoint))
			{
				next = null;
				if(this.isValidCnt(current))
					this.addNode(previous, current, next, p2ce, c2ne, path, isForward);
				break;
			}
		}
	}

	/**
	 * A method to build nodepath
	 * 
	 * @return
	 */
//	private List<NodePath> buildPathByCntFileEncapsulate() {
//		long start = System.currentTimeMillis();
//		if (edges == null || edges.size() == 0)
//			throw new IllegalArgumentException("PathBuilder: The Edges could not be empty!");
//		List<NodePath> paths = new Vector<NodePath>();
//		diGraph = null;
//		try {
//			diGraph = new DirectedGraph(edges, paras, cntfile);
//			// do transitive reduction
//			diGraph.transitiveReducting();
//			List<Edge> tempEdges = diGraph.getEdges();
//			logger.info("Edges size after transitive reducing: " + tempEdges.size());
//			String edgeFile = paras.getOutFolder() + System.getProperty("file.separator") + "edges_after_tr.info";
//			DotGraphFileWriter.writeEdge(edgeFile, tempEdges);
//			// delete error prone edge
//			diGraph.delErrorProneEdge(paras.getRatio());
//			tempEdges = diGraph.getEdges();
//			logger.info("Edges size after error prone deleting: " + tempEdges.size());
//			edgeFile = paras.getOutFolder() + System.getProperty("file.separator") + "edges_after_dep.info";
//			DotGraphFileWriter.writeEdge(edgeFile, tempEdges);
//			// delete tips
//			diGraph.delTips();
//			tempEdges = diGraph.getEdges();
//			logger.info("Edges size after deleting tips: " + tempEdges.size());
//			edgeFile = paras.getOutFolder() + System.getProperty("file.separator") + "edges_after_dt.info";
//			DotGraphFileWriter.writeEdge(edgeFile, tempEdges);
//			// delete similarity contigs error prone edges;
//			// diGraph.delSimCntEdges();
//			// tempEdges = diGraph.getEdges();
//			// logger.info("Edges size after deleting similarity contigs error
//			// prone edges : " + tempEdges.size());
//			// edgeFile = paras.getOutFolder() +
//			// System.getProperty("file.separator") + "edges_after_dsc.info";
//			// DotGraphFileWriter.writeEdge(edgeFile, tempEdges);
//			tempEdges = null;
//			path = null;
//			// TriadLinkReader tlr = new TriadLinkReader(paras);
//			// List<TriadLink> triads = tlr.read();
//			// travel the graph, random start
//			// do not including the divergence end point in the path
//			while (diGraph.isExistUnSelectedVertices()) {
//				INDEX++;
//				if (INDEX == 3)
//					logger.debug("breakpoint");
//				Contig current = diGraph.getRandomVertex();
//				// if(current.getID().equals("801"))
//				// logger.debug("breakpoint");
//				// if the return conting is null and the
//				// isExistUnSelectedVertices equal false then break;
//				if (current == null)
//					break;
//				List<Contig> adjs = diGraph.getAdjVertices(current);
//				if (adjs == null) {
//					diGraph.setVertexAsSelected(current);
//					continue;
//				}
//				int adjsSize = adjs.size();
//				// random selected the frist
//				if (adjsSize == 0) {
//					// orphan contig, only one element in path
//					if (this.isValidCnt(current)) {
//						path = new NodePath();
//						this.addNode(current, path, false);
//						paths.add(path);
//					} else {
//						diGraph.setVertexAsSelected(current);
//					}
//
//				} else if (adjsSize == 1) {
//					// normal start point, always on the linear end point;
//					path = new NodePath();
//					// Contig next = diGraph.getNextVertex(current, null);
//					Contig next = adjs.get(0);
//
//					while (true) {
//						if (this.isValidCnt(current)) {
//							this.addNode(current, path, false);
//						} else {
//							break;
//						}
//						Contig previous = current;
//						current = next;
//						if (current.getID().equals("2606"))
//							logger.debug("breakpoint");
//						next = diGraph.getNextVertex(current, previous);
//						// for the divergence point which could not determine
//						// how-to get next
//						if (next == null) {
//							LinkedList<Contig> formers = new LinkedList<Contig>();
//							if (path.getPathSize() > 2) {
//								formers.addLast(path.getElement(path.getPathSize() - 2).getCnt());
//								formers.addLast(path.getElement(path.getPathSize() - 3).getCnt());
//							} else if (path.getPathSize() == 2) {
//								formers.addLast(path.getElement(path.getPathSize() - 2).getCnt());
//							}
//							// List<Contig> selectedCnts =
//							// this.getTriadLinkNext3(current, previous);
//							List<Contig> selectedCnts = this.getInternalPath(current, previous, formers);
//							if (selectedCnts == null) {
//								this.addNode(current, path, false);
//								break;
//							}
//							int size = selectedCnts.size();
//							if (size == 1) {
//								next = selectedCnts.get(0);
//							} else if (size == 2) {
//								this.addNode(current, path, false);
//								previous = current;
//								current = selectedCnts.get(0);
//								next = selectedCnts.get(1);
//							} else {
//								this.addNode(current, path, false);
//								for (int i = 0; i < size - 2; i++) {
//									current = selectedCnts.get(i);
//									this.addNode(current, path, false);
//								}
//								previous = current;
//								current = selectedCnts.get(size - 2);
//								next = selectedCnts.get(size - 1);
//							}
//						}
//						if (next.equals(previous)) { // for the linear end point
//							this.addNode(current, path, false);
//							break;
//						}
//					}
//					paths.add(path);
//				} else if (adjsSize == 2) {
//					// middle point; normal start point, located in the linear
//					// path
//					path = new NodePath();
//					if (this.isValidCnt(current))
//						this.addNode(current, path, true);
//					else
//						continue;
//					Contig startPoint = current;
//					Contig c1 = adjs.get(0);
//					Contig c2 = adjs.get(1);
//					// c1--current--c2
//					// directed by both direction;
//					// for c1 direction; using c2 as previous point; checking
//					// whether c1 is valid point
//					// all the element unshift into path;
//					// checking the c1 and c2 whether travel
//					if (!this.isValidCnt(c1) && !this.isValidCnt(c2)) {
//						paths.add(path);
//						continue;
//					}
//					Contig previous = startPoint;
//					current = c1;
//					Contig next = diGraph.getNextVertex(current, previous);
//					while (true) {
//						if (next != null) {
//							if (next.equals(startPoint)) {
//								this.addNode(current, path, true);
//								break;
//							}
//							if (next.equals(previous)) {
//								this.addNode(current, path, true);
//								break;
//							}
//							// check the current contig whether selected befored
//							// if true break;
//							if (this.isValidCnt(current)) {
//								this.addNode(current, path, true);
//							} else {
//								break;
//							}
//							previous = current;
//							current = next;
//							if (current.getID().equals("2606"))
//								logger.debug("Breakpoint");
//							next = diGraph.getNextVertex(current, previous);
//						} else {
//							LinkedList<Contig> formers = new LinkedList<Contig>();
//							if (path.getPathSize() > 2) {
//								formers.addLast(path.getElement(1).getCnt());
//								formers.addLast(path.getElement(2).getCnt());
//							} else if (path.getPathSize() == 2) {
//								formers.addLast(path.getElement(1).getCnt());
//							} else if (path.getPathSize() == 1) {
//								formers.add(c2);
//							}
//							// to get the next point over the divergence point;
//							List<Contig> selectedCnts = this.getInternalPath(current, previous, formers);
//							if (selectedCnts == null) {
//								this.addNode(current, path, true);
//								break;
//							}
//							int size = selectedCnts.size();
//							if (size == 1) {
//								next = selectedCnts.get(0);
//							} else if (size == 2) {
//								this.addNode(current, path, true);
//								previous = current;
//								current = selectedCnts.get(0);
//								next = selectedCnts.get(1);
//							} else {
//								this.addNode(current, path, true);
//								for (int i = 0; i < size - 2; i++) {
//									current = selectedCnts.get(i);
//									this.addNode(current, path, true);
//								}
//								previous = current;
//								current = selectedCnts.get(size - 2);
//								next = selectedCnts.get(size - 1);
//							}
//						}
//					}
//					// for c2 direction; using c1 as previous point; checking
//					// whether c2 is valid point;
//					// all the valid element push into path;
//					if (!diGraph.isExistUnSelectedVertices()) {
//						paths.add(path);
//						break;
//					}
//					// if it is the loop, then the former current contig equal
//					// to
//					// next start contig c2, break;
//					if (current != null && c2 != null) {
//						if (current.equals(c2)) {
//							paths.add(path);
//							continue;
//						}
//					}
//					previous = startPoint;
//					current = c2;
//					next = diGraph.getNextVertex(current, previous);
//					while (true) {
//						if (next != null) {
//							if (next.equals(startPoint)) {
//								this.addNode(current, path, false);
//								break;
//							}
//							if (next.equals(previous)) {
//								this.addNode(current, path, false);
//								break;
//							}
//							// check the current contig whether selected before
//							if (this.isValidCnt(current))
//								this.addNode(current, path, false);
//							else
//								break;
//							previous = current;
//							current = next;
//							if (current.getID().equals("2606"))
//								logger.debug("breakpoint");
//							next = diGraph.getNextVertex(current, previous);
//						} else {
//							LinkedList<Contig> formers = new LinkedList<Contig>();
//							if (path.getPathSize() > 2) {
//								formers.addLast(path.getElement(path.getPathSize() - 2).getCnt());
//								formers.addLast(path.getElement(path.getPathSize() - 3).getCnt());
//							} else if (path.getPathSize() == 2) {
//								formers.addLast(path.getElement(path.getPathSize() - 2).getCnt());
//							}
//							List<Contig> selectedCnts = this.getInternalPath(current, previous, formers);
//							if (selectedCnts == null) {
//								this.addNode(current, path, false);
//								break;
//							}
//							int size = selectedCnts.size();
//							if (size == 1) {
//								next = selectedCnts.get(0);
//							} else if (size == 2) {
//								this.addNode(current, path, false);
//								previous = current;
//								current = selectedCnts.get(0);
//								next = selectedCnts.get(1);
//							} else {
//								this.addNode(current, path, false);
//								for (int i = 0; i < size - 2; i++) {
//									current = selectedCnts.get(i);
//									this.addNode(current, path, false);
//								}
//								previous = current;
//								current = selectedCnts.get(size - 2);
//								next = selectedCnts.get(size - 1);
//							}
//						}
//					}
//					paths.add(path);
//
//				} else if (adjsSize > 2) { // divergence point
//					diGraph.setVertexAsSelected(current);
//					continue;
//				}
//			}
//		} catch (Exception e) {
//			// logger.debug("INDEX = " + INDEX);
//			logger.error(this.getClass().getName() + "\t" + e.getMessage() + "\t" + e.getClass().getName());
//		}
//		// orientation contig in the paths;
//		// define the first element in path is forward;
//		// require to check the orientation conflict!
//		List<NodePath> tPaths = new Vector<NodePath>(paths.size());
//		try {
//			for (NodePath np : paths) {
//				int pathSize = np.getPathSize();
//				NodePath tNP = new NodePath();
//				Strand cStrand = null;
//				Strand nStrand = null;
//				if (pathSize == 1) {
//					Node current = np.getElement(0);
//					current.setStrand(Strand.FORWARD);
//					current.setMeanDist2Next(0);
//					current.setSdDist2Next(0);
//					current.setSupportLinkNum(0);
//					tNP.push(current);
//				} else {
//					Node current = null;
//					Node next = null;
//					for (int i = 0; i < pathSize - 1; i++) {
//						// INDEX++;
//						current = np.getElement(i);
//						next = np.getElement(i + 1);
//						Contig cCnt = current.getCnt();
//						Contig nCnt = next.getCnt();
//						Edge e = null;
//						List<Edge> es = diGraph.getEdgesInfo(cCnt, nCnt);
//						for (Edge t : es) {
//							if (t.getOrigin().equals(cCnt) && t.getTerminus().equals(nCnt)) {
//								e = t;
//								break;
//							}
//						}
//						// checking whether the edges links this two contigs are
//						// conflict;
//						if (cStrand == null)
//						// if (current.getStrand() == null)
//						{
//							// current.setStrand(e.getoStrand());
//							cStrand = e.getoStrand();
//						} else {
//							// conflict case;
//							// if(!current.getStrand().equals(e.getoStrand()))
//							// {
//							// tNP.push(current);
//							// tPaths.add(tNP);
//							// tNP = new NodePath();
//							// continue;
//							// }
//							if (!cStrand.equals(e.getoStrand())) {
//								Node cNode = new Node();
//								cNode.setCnt(cCnt);
//								cNode.setStrand(cStrand);
//								cNode.setOrphan(false);
//								tNP.push(cNode);
//								tPaths.add(tNP);
//								tNP = new NodePath();
//								cStrand = null;
//								nStrand = null;
//								i--;
//								continue;
//							}
//						}
//						// if (next.getStrand() == null)
//						// next.setStrand(e.gettStrand());
//						if (nStrand == null)
//							nStrand = e.gettStrand();
//						int meanSum = e.getDistMean();
//						int sdSum = e.getDistSd();
//						int slSum = e.getLinkNum();
//						Node cNode = new Node();
//						cNode.setCnt(cCnt);
//						cNode.setStrand(cStrand);
//						cNode.setMeanDist2Next(meanSum);
//						cNode.setSdDist2Next(sdSum);
//						cNode.setSupportLinkNum(slSum);
//						cNode.setOrphan(false);
//						tNP.push(cNode);
//						// current.setMeanDist2Next(meanSum);
//						// current.setSdDist2Next(sdSum);
//						// current.setSupportLinkNum(slSum);
//						// tNP.push(current);
//						// for the last node
//						if (i == (pathSize - 2)) {
//							Node nNode = new Node();
//							nNode.setCnt(nCnt);
//							nNode.setStrand(nStrand);
//							tNP.push(nNode);
//						}
//						cStrand = nStrand;
//						nStrand = null;
//					}
//				}
//				if (tNP != null && tNP.getPathSize() > 0)
//					tPaths.add(tNP);
//			}
//		} catch (Exception e) {
//			// logger.debug("INDEX\t" + INDEX);
//			logger.error(this.getClass().getName() + "\t" + e.getMessage() + "\t" + e.getClass().getName());
//		}
//		// cleaning and call gc
//		this.diGraph = null;
//		this.edges = null;
//		this.triads = null;
//		System.gc();
//		long end = System.currentTimeMillis();
//		logger.info("Path Building, erase time: " + (end - start) + " ms");
//		return tPaths;
//	}

	private List<Contig> getInternalPath(Contig internal, Contig external, List<Contig> formers) {
		// get adjacent contigs of internal excluding external;
		List<Contig> adjInternals = diGraph.getNextVertices(internal, external);
		if (adjInternals == null || adjInternals.isEmpty())
			return null;
		// build the list InternalNode
//		List<List<InternalNode>> inss = new Vector<List<InternalNode>>(); 
		List<List<InternalNode>> inss = new Vector<List<InternalNode>>();
		List<Edge> es = diGraph.getEdgesInfo(external, internal);
		Strand strand = null;
		for (Edge e : es) {
			if (e.getOrigin().equals(external) && e.getTerminus().equals(internal)) {
				strand = e.gettStrand();
				break;
			}
		}
		
		for (Contig c : adjInternals) {
			Contig convergence = internal;
			Map<String, Object> outcome = this.isBubbleAdjacent(internal, c, 4, adjInternals);
			boolean isBubble = (boolean) outcome.get("IsBubble");
			Contig tcon = (Contig) outcome.get("CNT");
			if(isBubble && tcon != null )
				convergence = tcon;
//			List<InternalNode> ins = new Vector<InternalNode>(30);
//			InternalNode in = new InternalNode();
//			in.setGrandfather(null);
//			in.setParent(null);
//			in.setChildren(internal);
//			in.setLeaf(false);
//			ins.add(in);
//			this.getInternalPath(c, internal, null, strand, isBubble, convergence, 6, ins);
//			inss.add(ins);
			List<InternalNode> ins = new Vector<InternalNode>(30);
			InternalNode in = new InternalNode();
			in.setParent(null);
			in.setChild(internal);
			in.setLeaf(false);
			this.getInternalPath(c, internal, null, strand, isBubble, convergence, 6, ins, in);
			inss.add(ins);
		}
		// build the path from internal node list;
		List<InternalPath> ips = buildInternalPathFromInternalNode(inss);
		if(ips == null)
			return null;
		this.computeInternalPathSupported(ips, external, formers);
		InternalPathComparator ic = new InternalPathComparator();
		Collections.sort(ips, ic);
		LinkedList<Contig> internalPath = null;
		if(ips == null || ips.isEmpty())
			return null;
		
		InternalPath ip1 = ips.get(0);
		int size = ips.size() - 1;
		int ip1score = 0;
		int ip2score = 0;
		InternalPath ip2 = null; // ip need to test
		while(true)
		{
			ip1score = ip1.getScore();
			if(ip1.getScore() <= 0)
			{
				ip1 = null;
				break;
			}
			for(int index = 1; index <= size; index++)
			{
				ip2 = ips.get(index);
				ip2score = ip2.getScore();
				if(ip1score > ip2score)
				{
					break;
				} else
				{
					ip1 = getBestInternalPath(ip1, ip2);
				}
			}
			if(ip1 != null)
				break;
//			if(isBestSupportedInternalPath(ip1, adjInternals))
//			{
//				break;
//			} else
//			{
//				ips.remove(ip1);
//				if(ips == null || ips.isEmpty())
//				{
//					ip1 = null;
//					break;
//				}
//				size = ips.size() - 1;
//				ip1 = ips.get(0);
//			}
		}
		if(ip1 != null)
		{
			internalPath = ip1.getPath();
			internalPath.remove(internal);
		}
		return internalPath;
	}
	
//	private boolean isBestSupportedInternalPath(InternalPath ip, List<Contig> adjs)
//	{
//		int size = ip.getPath().size();
//		if(size <= 1)
//			return true;
//		// only considering the uniques contigs
//		List<Contig> temp = new LinkedList<Contig>();
//		for(int i = size - 1; i >= 0; i--)
//		{
//			Contig c = ip.getCnt(i);
//			if(diGraph.isDivergenceVertex(c))
//				break;
//			temp.add(c);
//		}
//		int score = ip.getScore();
//		int value = 0;
//		for(Contig c : adjs)
//		{
//			if(ip.isContain(c))
//				continue;
//			for(TriadLink tl : triads)
//			{
//				if(tl.isContain(c))
//				{
//					Contig pre = tl.getPrevious();
//					Contig mid = tl.getMiddle();
//					Contig lst = tl.getLast();
//					if(temp.contains(pre))
//						value += tl.getSupLinks();
//					else if(temp.contains(lst))
//						value += tl.getSupLinks();
//					else if(mid != null && temp.contains(mid))
//						value += tl.getSupLinks();
//				}
//				if(value > score)
//					return false;
//			}
//		}
//		return true;
//	}
	
	/**
	 * @param previous : the previous contig
	 * @param current : the current contig need to check whether is bubble structure;
	 * @param depth : the depth need to check;
	 * @param adjacents : the adjacent contigs of internal contig;
	 * @return
	 */
	private Map<String, Object> isBubbleAdjacent(Contig previous, Contig current, int depth, List<Contig> adjcents)
	{
		boolean isBubble = false;
		Contig convergence = null;
		Map<String, Object> values = new HashMap<String, Object>();
		if(depth == 0)
		{
			isBubble = false;
			convergence = null;
			values.put("IsBubble", isBubble);
			values.put("CNT", convergence);
			return values;
		}
		List<Contig> nexts = diGraph.getNextVertices(current, previous);
		if(nexts == null || nexts.isEmpty())
		{
			isBubble = false;
			convergence = null;
			values.put("IsBubble", isBubble);
			values.put("CNT", convergence);
			return values;
		}
		if(nexts.size() > 1)
			convergence = current;
		boolean isExist = false;
		for(Contig c : nexts)
		{
			if(adjcents.contains(c))
			{
				isExist = true;
				isBubble = true;
				values.put("IsBubble", isBubble);
				values.put("CNT", convergence);
				break;
			}
		}
		if(isExist){
			return values;
		} else
		{
			for(Contig c : nexts)
			{
				Map<String, Object> outcome = isBubbleAdjacent(current, c, depth - 1, adjcents);
				if((boolean) outcome.get("IsBubble"))
				{
					isBubble = true;
//					convergence = (Contig) outcome.get("CNT");
					values.put("IsBubble", isBubble);
					if(outcome.get("CNT") == null)
						values.put("CNT", convergence);
					else 
						values.put("CNT", outcome.get("CNT"));
					break;
				}
				values.put("IsBubble", isBubble);
				values.put("CNT", convergence);
			}
		}
		return values;
	}
	
//	private Map<String, Object> isBubble(Contig current, Contig previous, List<Contig> adjs, int depth)
//	{
//		boolean isBubble = false;
//		Contig convergence = null;
//		Map<String, Object> values = new HashMap<String, Object>();
//		if(depth == 0)
//		{
//			isBubble = false;
//			convergence = null;
//			values.put("IsBubble", isBubble);
//			values.put("CNT", convergence);
//			return values;
//		}
//		List<Contig> nexts = diGraph.getNextVertices(current, previous);
//		if(nexts == null || nexts.isEmpty())
//		{
//			isBubble = false;
//			convergence = null;
//			values.put("IsBubble", isBubble);
//			values.put("CNT", convergence);
//			return values;
//		}
//		if(nexts.size() > 1)
//			convergence = current;
//		for(Contig c : nexts)
//		{
//			if(adjs.contains(c))
//			{
//				isBubble = true;
//				values.put("IsBubble", isBubble);
//				values.put("CNT", convergence);
//				break;
//			}
//			Map<String, Object> outcome = isBubble(c, current, adjs, depth - 1);
//			if((boolean) outcome.get("IsBubble"))
//			{
//				isBubble = true;
//				convergence = (Contig) outcome.get("CNT");
//				values.put("IsBubble", isBubble);
//				values.put("CNT", convergence);
//				break;
//			}
//		}
//		return values;
//	}
	
//	private boolean isInternalPathTraveled(List<Contig> ip)
//	{
//		if(ip == null || ip.isEmpty())
//			return false;
//		// do not considering the last point;
//		boolean isTraveled = true;
//		int size = ip.size() - 2;
//		for(int i = 0; i <= size; i++)
//		{
//			Contig c = ip.get(i);
//			if(!diGraph.isDivergenceVertex(c))
//			{
//				if(!diGraph.isVertexSelected(c))
//				{
//					isTraveled = false;
//					break;
//				}
//			}
//		}
//		return isTraveled;
//	}
	
	private InternalPath getBestInternalPath(InternalPath ip1, InternalPath ip2)
	{
		LinkedList<Contig> ip1cnts = ip1.getPath();
		Contig ip1c1 = ip1cnts.get(0);
		Contig ip1c2 = ip1cnts.get(1);
		LinkedList<Contig> ip2cnts = ip2.getPath();
		Contig ip2c1 = ip2cnts.get(0);
		Contig ip2c2 = ip2cnts.get(1);
		
		List<Edge> ip1es = diGraph.getEdgesInfo(ip1c1, ip1c2);
		int ip1esl = ip1es.get(0).getLinkNum();
		List<Edge> ip2es = diGraph.getEdgesInfo(ip2c1, ip2c2);
		int ip2esl = ip2es.get(0).getLinkNum();
		
		if(ip1esl > ip2esl)
			return ip1;
		else
			return ip2;
		
	}
	
	private void computeInternalPathSupported(List<InternalPath> ips, Contig external, List<Contig> formers)
	{
		if(ips == null || ips.isEmpty())
			return;
		if(triads == null)
		{
			TriadLinkReader tlr = new TriadLinkReader(paras);
			triads = tlr.read();
		}
		// external uniques 
		LinkedList<Contig> extUniqs = new LinkedList<Contig>();
		extUniqs.addLast(external);
		for(Contig c : formers)
		{
			extUniqs.addLast(c);
		}
		List<TriadLink> selectedTLs = new Vector<TriadLink>();
		// scoring each path;
		for (InternalPath ip : ips) {
			// logger.debug(ip.toString());
			for (TriadLink tl : triads) {
				Contig pre = tl.getPrevious();
				Contig mid = tl.getMiddle();
				Contig lst = tl.getLast();
//				if(pre.getID().equals("1317") || lst.getID().equals("1317"))
//					logger.debug("breakpoint");
				if(tl.isValid())
				{
//					Contig pre = tl.getPrevious();
//					Contig mid = tl.getMiddle();
//					Contig lst = tl.getLast();
//					if(pre.getID().equals("1317") || lst.getID().equals("1317"))
//						logger.debug("breakpoint");
					// there are two cases: 
					// 1)connect between unique contig and internal path;
					// 2)support internal path;
					if(extUniqs.contains(pre))
					{
						if(mid == null)
						{
							if(ip.isContain(lst))
							{
								ip.addScore(tl.getSupLinks());
								if(!selectedTLs.contains(tl))
									selectedTLs.add(tl);
							}
						} else
						{
							if(ip.isContain(mid))
							{
								if(ip.isContain(lst))
								{
									ip.addScore(tl.getSupLinks());
									if(!selectedTLs.contains(tl))
										selectedTLs.add(tl);
								}
							} else
							{
								if(ip.isContain(lst))
								{
									ip.addScore(tl.getSupLinks());
									if(!selectedTLs.contains(tl))
										selectedTLs.add(tl);
								}
							}
						}
					} else if(extUniqs.contains(mid))
					{
						if(ip.isContain(pre) || ip.isContain(lst))
						{
							ip.addScore(tl.getSupLinks());
							if(!selectedTLs.contains(tl))
								selectedTLs.add(tl);
						}
					} else if(extUniqs.contains(lst))
					{
						if(mid == null)
						{
							if(ip.isContain(pre))
							{
								ip.addScore(tl.getSupLinks());
								if(!selectedTLs.contains(tl))
									selectedTLs.add(tl);
							}
						} else
						{
							if(ip.isContain(mid))
							{
								if(ip.isContain(pre))
								{
									ip.addScore(tl.getSupLinks());
									if(!selectedTLs.contains(tl))
										selectedTLs.add(tl);
								}
							} else
							{
								if(ip.isContain(pre))
								{
									ip.addScore(tl.getSupLinks());
									if(!selectedTLs.contains(tl))
										selectedTLs.add(tl);
								}
							}
						}
					} else
					{
						//support internal path case
//						if(ip.isInternalPathContain(pre))
//						{
//							if(mid == null)
//							{
//								if(ip.isInternalPathContain(lst))
//								{
//									ip.addScore(tl.getSupLinks());
//									if(!selectedTLs.contains(tl))
//										selectedTLs.add(tl);
//								}
//							} else
//							{
//								if(ip.isInternalPathContain(mid))
//								{
//									ip.addScore(tl.getSupLinks());
//									if(!selectedTLs.contains(tl))
//										selectedTLs.add(tl);
//								}
//							}
//						} else if(ip.isInternalPathContain(mid))
//						{
//							if(ip.isInternalPathContain(lst))
//							{
//								ip.addScore(tl.getSupLinks());
//								if(!selectedTLs.contains(tl))
//									selectedTLs.add(tl);
//							}
//						} 
					}
				}
			}
		}
		
		for(TriadLink tl : selectedTLs)
		{
			tl.setValid(false);
//			logger.debug(tl.toString());
		}
	}

//	private boolean validateInternalPath2(InternalPath path) {
//		boolean isValid = true;
//		if (path == null || path.isEmpty())
//			return false;
//		List<Contig> tpath = path.getPath();
//		Contig current = null;
//		Contig next = null;
//		Strand cStrand = null;
//		int pSize = tpath.size();
//		for (int i = 0; i < pSize - 1; i++) {
//			current = tpath.get(i);
//			next = tpath.get(i + 1);
//			List<Edge> es = diGraph.getEdgesInfo(current, next);
//			Edge e = null;
//			for (Edge t : es) {
//				if (t.getOrigin().equals(current) && t.getTerminus().equals(next)) {
//					e = t;
//					break;
//				}
//			}
//			if (cStrand == null) {
//				cStrand = e.gettStrand();
//				continue;
//			} else {
//				if (!e.getoStrand().equals(cStrand)) {
//					isValid = false;
//					break;
//				} else {
//					cStrand = e.gettStrand();
//				}
//			}
//		}
//		return isValid;
//	}
	
	private List<InternalPath> buildInternalPathFromInternalNode(List<List<InternalNode>> inss) {
		if (inss == null || inss.isEmpty())
			return null;
		List<InternalPath> ips = new Vector<InternalPath>();
		for(int i = 0; i < inss.size(); i++)
		{
			List<InternalNode> ins = inss.get(i);
			Iterator<InternalNode> it = ins.iterator();
			while (it.hasNext()) {
				InternalNode in = it.next();
				if (in.isLeaf()) {
					InternalPath ip = new InternalPath();
					while (true) {
						Contig child = in.getChild();
						InternalNode parent = in.getParent();
						ip.addFirst(child);
						if (parent == null)
							break;
						in = parent;
//						in = findReverseInternalNode2(parent , ins);
					}
					//checking whether internal path is already travel;
					List<Contig> path = ip.getPath();
					boolean isTravel = true;
					for(Contig c : path)
					{
						if(!diGraph.isDivergenceVertex(c))
						{
							if(!diGraph.isVertexSelected(c))
							{
								isTravel = false;
								break;
							}
						}
					}
					if(!isTravel)
						ips.add(ip);
				}
			}
		}
		if(ips.isEmpty())
			return null;
		return ips;
	}

//	private List<InternalPath> buildInternalPathFromInternalNode(List<List<InternalNode>> inss) {
//		if (inss == null || inss.isEmpty())
//			return null;
//		List<InternalPath> ips = new Vector<InternalPath>();
//		for(int i = 0; i < inss.size(); i++)
//		{
//			List<InternalNode> ins = inss.get(i);
//			Iterator<InternalNode> it = ins.iterator();
//			while (it.hasNext()) {
//				InternalNode in = it.next();
//				if (in.isLeaf()) {
//					InternalPath ip = new InternalPath();
//					while (true) {
//						Contig child = in.getChildren();
//						Contig parent = in.getParent();
//						Contig grandfather = in.getGrandfather();
//						ip.addFirst(child);
//						if (parent == null)
//							break;
//						in = findReverseInternalNode(parent, grandfather, ins);
//					}
//					//checking whether internal path is already travel;
//					List<Contig> path = ip.getPath();
//					boolean isTravel = true;
//					for(Contig c : path)
//					{
//						if(!diGraph.isDivergenceVertex(c))
//						{
//							if(!diGraph.isVertexSelected(c))
//							{
//								isTravel = false;
//								break;
//							}
//						}
//					}
//					if(!isTravel)
//						ips.add(ip);
//				}
//			}
//		}
//		if(ips.isEmpty())
//			return null;
//		return ips;
//	}
//	
	
//	private InternalNode2 findReverseInternalNode2(InternalNode2 parent,  List<InternalNode2> ins) {
//		for (InternalNode2 in : ins) {
//			if(in.isLeaf())
//				continue;
//			if (in.equals(parent)) {
//				return in;
//			}
//		}
//		return null;
//	}

//	private InternalNode findReverseInternalNode(Contig parent, Contig grandfather, List<InternalNode> ins) {
//		InternalNode value = null;
//		for (InternalNode in : ins) {
//			if(in.isLeaf())
//				continue;
//			if (in.getChildren().equals(parent)) {
//				if (grandfather == null) {
//					if (in.getParent() == null) {
//						value = in;
//						break;
//					}
//				} else {
//					if (in.getParent() != null && in.getParent().equals(grandfather)) {
//						value = in;
//						break;
//					}
//				}
//			}
//		}
//		return value;
//	}
//	
	
	/**
	 * A method to get the unique contig aside to divergence contig;
	 * 
	 * @param child
	 *            - the current contig the check for unique contig
	 * @param father
	 *            - the former contig
	 * @param grandfather
	 *            - the former former contig
	 * @param fatherStrand
	 *            - former contig strand.
	 * @param depth
	 *            - the check depth
	 * @param ins
	 *            - the list to store the internal node;
	 *@param  innercnts
	 *			  - the list to store already travel contigs;          
	 *
	 */
	private void getInternalPath(Contig child, Contig father, Contig grandfather, Strand fatherStrand, 
			boolean isBubble, Contig convergence, int depth, List<InternalNode> ins, InternalNode parent) {
		// check orientation conflict
		Map<String, Object> checks = this.validateInternalPathOrientation(child, father, fatherStrand);
		if (!(boolean) checks.get("VALID"))
			return;
		Strand childStrand = (Strand) checks.get("STRAND");
		List<Contig> nextAdjs = diGraph.getNextVertices(child, father);
		InternalNode in = new InternalNode();
		in.setChild(child);
		in.setParent(parent);
		if (nextAdjs == null || nextAdjs.size() == 0) {
			in.setLeaf(true);
			ins.add(in);
			return;
		}
		if (depth == 0) {
			in.setLeaf(true);
			ins.add(in);
			return;
		}
		if(child.equals(convergence))
		{
			in.setLeaf(true);
			ins.add(in);
			return;
		}
		if(isBubble && father.equals(convergence))
		{
			in.setLeaf(true);
			ins.add(in);
			return;
		}

		// two cases: 1) the next point is divergence; 2) the next point is
		// unique;
		if (nextAdjs.size() > 1) {
			List<Edge> egs = diGraph.getEdgesInfo(child, father);
			int eLen = egs.get(0).getDistMean();
			INTERNAL_LENGTH += eLen;
			// INTERNAL_LENGTH += this.indexLen(current.getID());
//			int cLen = cntfile.getLengthByNewId(Integer.valueOf(child.getID()));
			int cLen = cntLens.get(child.getID());
			INTERNAL_LENGTH += cLen;
			if (INTERNAL_LENGTH <= MAXIMUM_INTERNAL_LENGTH) {
				boolean validAlls = false;
				for (Contig c : nextAdjs) {
					if(c.equals(convergence))
						continue;
					Map<String, Object> nValues = this.validateInternalPathOrientation(c, child, childStrand);
					if ((boolean) nValues.get("VALID")) {
						this.getInternalPath(c, child, father, childStrand, isBubble, convergence, depth - 1, ins, in);
						validAlls = true;
					} 
				}
				if (!validAlls)
					in.setLeaf(true);
				ins.add(in);
				INTERNAL_LENGTH -= cLen;
				INTERNAL_LENGTH -= eLen;
				return;
			} else {
				boolean validAlls = false;
				for (Contig c : nextAdjs) {
					if(c.equals(convergence))
						continue;
					Map<String, Object> nValues = this.validateInternalPathOrientation(c, child, childStrand);
					if ((boolean) nValues.get("VALID")) {
						InternalNode iin = new InternalNode();
						iin.setChild(c);
						iin.setParent(in);
						iin.setLeaf(true);
						ins.add(iin);
						validAlls = true;
					} 
				}
				if (!validAlls)
					in.setLeaf(true);
				ins.add(in);
				INTERNAL_LENGTH -= eLen;
				INTERNAL_LENGTH -= cLen;
				return;
			}
		} else {
			Contig c = nextAdjs.get(0);
			List<Edge> egs = diGraph.getEdgesInfo(child, father);
			int eLen = egs.get(0).getDistMean();
			INTERNAL_LENGTH += eLen;
			// INTERNAL_LENGTH += this.indexLen(current.getID());
//			int cLen = cntfile.getLengthByNewId(Integer.valueOf(child.getID()));
			int cLen = cntLens.get(child.getID());
			INTERNAL_LENGTH += cLen;
			// if the internal length of path already large than MAXIMUM_INTERNAL_LENGTH
			if (INTERNAL_LENGTH > MAXIMUM_INTERNAL_LENGTH) {
				in.setLeaf(true);
				ins.add(in);
				INTERNAL_LENGTH -= cLen;
				INTERNAL_LENGTH -= eLen;
				return;
			}
			if(c.equals(convergence))
			{
				in.setLeaf(true);
				ins.add(in);
				INTERNAL_LENGTH -= cLen;
				INTERNAL_LENGTH -= eLen;
				return;
			}
			Map<String, Object> nValues = this.validateInternalPathOrientation(c, child, childStrand);
			if ((boolean) nValues.get("VALID")) {
				this.getInternalPath(c, child, father, childStrand, isBubble, convergence, depth - 1, ins, in);
			} else {
				in.setLeaf(true);
			}
			ins.add(in);
			INTERNAL_LENGTH -= cLen;
			INTERNAL_LENGTH -= eLen;
		}
	}

	/**
	 * A method to get the unique contig aside to divergence contig;
	 * 
	 * @param child
	 *            - the current contig the check for unique contig
	 * @param father
	 *            - the former contig
	 * @param grandfather
	 *            - the former former contig
	 * @param fatherStrand
	 *            - former contig strand.
	 * @param depth
	 *            - the check depth
	 * @param ins
	 *            - the list to store the internal node;
	 *@param  innercnts
	 *			  - the list to store already travel contigs;          
	 *
	 */
//	private void getInternalPath(Contig child, Contig father, Contig grandfather, Strand fatherStrand, 
//			boolean isBubble, Contig convergence, int depth, List<InternalNode> ins) {
//		// check orientation conflict
//		Map<String, Object> checks = this.validateInternalPathOrientation(child, father, fatherStrand);
//		if (!(boolean) checks.get("VALID"))
//			return;
//		Strand childStrand = (Strand) checks.get("STRAND");
//		List<Contig> nextAdjs = diGraph.getNextVertices(child, father);
//		InternalNode in = new InternalNode();
//		in.setChildren(child);
//		in.setParent(father);
//		in.setGrandfather(grandfather);
//		if (nextAdjs == null || nextAdjs.size() == 0) {
//			in.setLeaf(true);
//			ins.add(in);
//			return;
//		}
//		if (depth == 0) {
//			in.setLeaf(true);
//			ins.add(in);
//			return;
//		}
//		if(child.equals(convergence))
//		{
//			in.setLeaf(true);
//			ins.add(in);
//			return;
//		}
//		if(isBubble && father.equals(convergence))
//		{
//			in.setLeaf(true);
//			ins.add(in);
//			return;
//		}
//
//		// two cases: 1) the next point is divergence; 2) the next point is
//		// unique;
//		if (nextAdjs.size() > 1) {
//			List<Edge> egs = diGraph.getEdgesInfo(child, father);
//			int eLen = egs.get(0).getDistMean();
//			INTERNAL_LENGTH += eLen;
//			// INTERNAL_LENGTH += this.indexLen(current.getID());
////			int cLen = cntfile.getLengthByNewId(Integer.valueOf(child.getID()));
//			int cLen = cntLens.get(child.getID());
//			INTERNAL_LENGTH += cLen;
//			if (INTERNAL_LENGTH <= MAXIMUM_INTERNAL_LENGTH) {
//				boolean validAlls = false;
//				for (Contig c : nextAdjs) {
//					if(c.equals(convergence))
//						continue;
//					Map<String, Object> nValues = this.validateInternalPathOrientation(c, child, childStrand);
//					if ((boolean) nValues.get("VALID")) {
//						this.getInternalPath(c, child, father, childStrand, isBubble, convergence, depth - 1, ins);
//						validAlls = true;
//					} 
//				}
//				if (!validAlls)
//					in.setLeaf(true);
//				ins.add(in);
//				INTERNAL_LENGTH -= cLen;
//				INTERNAL_LENGTH -= eLen;
//				return;
//			} else {
//				boolean validAlls = false;
//				for (Contig c : nextAdjs) {
//					if(c.equals(convergence))
//						continue;
//					Map<String, Object> nValues = this.validateInternalPathOrientation(c, child, childStrand);
//					if ((boolean) nValues.get("VALID")) {
//						InternalNode iin = new InternalNode();
//						iin.setChildren(c);
//						iin.setParent(child);
//						iin.setGrandfather(father);
//						iin.setLeaf(true);
//						ins.add(iin);
//						validAlls = true;
//					} 
//				}
//				if (!validAlls)
//					in.setLeaf(true);
//				ins.add(in);
//				INTERNAL_LENGTH -= eLen;
//				INTERNAL_LENGTH -= cLen;
//				return;
//			}
//		} else {
//			Contig c = nextAdjs.get(0);
//			List<Edge> egs = diGraph.getEdgesInfo(child, father);
//			int eLen = egs.get(0).getDistMean();
//			INTERNAL_LENGTH += eLen;
//			// INTERNAL_LENGTH += this.indexLen(current.getID());
////			int cLen = cntfile.getLengthByNewId(Integer.valueOf(child.getID()));
//			int cLen = cntLens.get(child.getID());
//			INTERNAL_LENGTH += cLen;
//			// if the internal length of path already large than MAXIMUM_INTERNAL_LENGTH
//			if (INTERNAL_LENGTH > MAXIMUM_INTERNAL_LENGTH) {
//				in.setLeaf(true);
//				ins.add(in);
//				INTERNAL_LENGTH -= cLen;
//				INTERNAL_LENGTH -= eLen;
//				return;
//			}
//			if(c.equals(convergence))
//			{
//				in.setLeaf(true);
//				ins.add(in);
//				INTERNAL_LENGTH -= cLen;
//				INTERNAL_LENGTH -= eLen;
//				return;
//			}
//			Map<String, Object> nValues = this.validateInternalPathOrientation(c, child, childStrand);
//			if ((boolean) nValues.get("VALID")) {
//				this.getInternalPath(c, child, father, childStrand, isBubble, convergence, depth - 1, ins);
//			} else {
//				in.setLeaf(true);
//			}
//			ins.add(in);
//			INTERNAL_LENGTH -= cLen;
//			INTERNAL_LENGTH -= eLen;
//		}
//	}

	private Map<String, Object> validateInternalPathOrientation(Contig child, Contig father, Strand fStrand) {

		List<Edge> es = diGraph.getEdgesInfo(father, child);
		if (es == null || es.isEmpty()) {
			Map<String, Object> values = new HashMap<String, Object>();
			values.put("VALID", false);
			values.put("STRAND", null);
			return values;
		}
		Strand oStrand = null;
		Strand tStrand = null;
		for (Edge e : es) {
			if (e.getOrigin().equals(father) && e.getTerminus().equals(child)) {
				oStrand = e.getoStrand();
				tStrand = e.gettStrand();
				break;
			}
		}
		if (oStrand.equals(fStrand)) {
			Map<String, Object> values = new HashMap<String, Object>();
			values.put("VALID", true);
			values.put("STRAND", tStrand);
			return values;
		} else {
			Map<String, Object> values = new HashMap<String, Object>();
			values.put("VALID", false);
			values.put("STRAND", null);
			return values;
		}
	}

	/**
	 * The fourth method implemented to get next point over the divergence
	 * point;
	 * 
	 * @param internal
	 *            - the divergence point
	 * @param external
	 *            - the former point adjacent to divergence point
	 * @param formers
	 *            - the formers contigs already in the path
	 * @return
	 */

//	private List<Contig> getTriadLinkNext(Contig internal, Contig external, List<Contig> formers) {
//		// contigs adjacents to internal but exclude external;
//		List<Contig> adjInternals = diGraph.getNextVertices(internal, external);
//		if (adjInternals.size() == 0)
//			return null;
//		// experience for yeast bubbles structure, if the adjInternals contig
//		// links path formers contigs then remove it;
//		List<Contig> tAdjInternals = new Vector<Contig>(adjInternals.size());
//		if (formers.size() > 0) {
//			for (Contig c : adjInternals) {
//				for (Contig i : formers) {
//					List<Edge> tes = diGraph.getEdgesInfo(c, i);
//					if (tes != null && tes.size() > 0) {
//						if (!tAdjInternals.contains(c))
//							tAdjInternals.add(c);
//					}
//				}
//			}
//		}
//
//		Iterator<Contig> it = adjInternals.iterator();
//		if (tAdjInternals != null && tAdjInternals.size() > 0) {
//			while (it.hasNext()) {
//				Contig c = it.next();
//				if (tAdjInternals.contains(c))
//					it.remove();
//			}
//		}
//		if (adjInternals.size() == 1)
//			return adjInternals;
//
//		// the unique contig list through path from external to internal;
//		List<Contig> uniques = new Vector<Contig>(5);
//		// loop the internal adjacent contigs
//		for (Contig c : adjInternals) {
//			List<Contig> tempUnique = new Vector<Contig>();
//			tempUnique.add(c);
//			this.getNextUniqueContigs(c, internal, 3, tempUnique);
//			for (Contig t : tempUnique) {
//				if (!uniques.contains(t))
//					uniques.add(t);
//			}
//		}
//		// remove unique contig already in path;
//		// remove unique contig is directly link former contigs
//		List<Contig> tempUniques = new Vector<Contig>(uniques.size());
//		for (Contig c : uniques) {
//			if (!path.isContain(c)) {
//				boolean isExist = false;
//				for (Contig in : formers) {
//					List<Edge> es = diGraph.getEdgesInfo(c, in);
//					if (es != null && !es.isEmpty()) {
//						isExist = true;
//						break;
//					}
//				}
//				if (!isExist) {
//					if (!c.equals(internal))
//						tempUniques.add(c);
//				}
//			}
//		}
//		uniques = null;
//		uniques = tempUniques;
//		tempUniques = null;
//
//		if (uniques.size() == 0)
//			return null;
//		// candidate triad link
//		List<TriadLink> canTls = new Vector<TriadLink>(10);
//		TriadLinkComparator tlc = new TriadLinkComparator();
//		// store requried deleting triadlinks;
//		Map<String, List<TriadLink>> delLinksMap = new HashMap<String, List<TriadLink>>();
//		for (Contig unique : uniques) {
//			List<TriadLink> tls = new Vector<TriadLink>();
//			int depth = 4;
//			this.getSupportedTraidLinks(external, internal, unique, internal, depth, tls);
//			// considering next external contigs
//			List<TriadLink> temp = null;
//			for (Contig c : formers) {
//				if (diGraph.isDivergenceVertex(c))
//					continue;
//				temp = this.findTriadLinks(c, internal, unique);
//				if (temp != null && temp.size() != 0) {
//					for (TriadLink t : temp) {
//						if (!tls.contains(t))
//							tls.add(t);
//					}
//				}
//			}
//
//			// if the inner contig is large than 5000 bp, considering it is
//			// unqiue;
//			// int inLen = this.indexLen(internal.getID());
//			int inLen = cntfile.getLengthByNewId(Integer.valueOf(internal.getID()));
//			if (inLen >= 5000) {
//				temp = this.findTriadLinks(internal, null, unique);
//				if (temp != null && !temp.isEmpty()) {
//					for (TriadLink t : temp) {
//						if (!tls.contains(t))
//							tls.add(t);
//					}
//				}
//			}
//
//			TriadLink tl = new TriadLink();
//			tl.setPrevious(external);
//			tl.setMiddle(internal);
//			tl.setLast(unique);
//
//			for (TriadLink t : tls) {
//				Contig middle = t.getMiddle();
//				Contig last = t.getLast();
//				// for special case when the middle is adjacent to internal
//				// point, it means
//				// this path is circle.
//				if (last != null && last.equals(unique)) {
//					if (middle != null && adjInternals.contains(middle))
//						continue;
//				}
//				if (middle != null) {
//					if (t.isContain(internal)) {
//						if (t.getMiddle().equals(internal)) {
//							if (t.isContain(unique))
//								tl.setSupLinks(tl.getSupLinks() + t.getSupLinks() + 2);
//							else
//								tl.setSupLinks(tl.getSupLinks() + t.getSupLinks());
//						} else {
//							if (t.isContain(external) && t.isContain(unique))
//								tl.setSupLinks(tl.getSupLinks() - t.getSupLinks());
//							else
//								tl.setSupLinks(tl.getSupLinks() + t.getSupLinks());
//						}
//					} else {
//						tl.setSupLinks(tl.getSupLinks() + t.getSupLinks());
//					}
//				} else {
//					tl.setSupLinks(tl.getSupLinks() + t.getSupLinks());
//				}
//			}
//			canTls.add(tl);
//			delLinksMap.put(unique.getID(), tls);
//		}
//
//		if (canTls.isEmpty())
//			return null;
//		Collections.sort(canTls, tlc);
//		// to get the best supported links;
//		// if there are two best links, it will check the revered direction
//		// and decide which one is the best else return null;
//		TriadLink tl = null;
//		for (TriadLink t : canTls) {
//			if (isBestSupportLink(t, uniques)) {
//				tl = t;
//				break;
//			}
//		}
//		// if (canTls.size() > 1) {
//		// // right now we only considered if there two best links;
//		// // the more complex case will implement if there are exist;
//		// TriadLink first = canTls.get(0);
//		// TriadLink second = canTls.get(1);
//		// if (first.getSupLinks() == second.getSupLinks())
//		// {
//		// // logic to check which one is the best;
//		// boolean p1 = isBestSupportLink(canTls.get(0), uniques);
//		// boolean p2 = isBestSupportLink(canTls.get(1), uniques);
//		// if(p1 && p2)
//		// { // both is the best case;
//		// return null;
//		// } else if(p1 && (!p2))
//		// { // the first one is the best;
//		// tl = first;
//		// } else if((!p1) && p2)
//		// { // the second one is the best;
//		// tl = second;
//		// } else
//		// {
//		// return null;
//		// }
//		// }else
//		// {
//		// tl = canTls.get(0);
//		// }
//		// } else {
//		// if (canTls.size() != 0)
//		// tl = canTls.get(0);
//		// else
//		// return null;
//		// }
//		if (tl == null || tl.getSupLinks() == 0)
//			return null;
//		// remove the supported triadlinks
//		// List<TriadLink> dels = delLinksMap.get(tl.getLast().getID());
//		// for (TriadLink t : dels) {
//		// triads.remove(t);
//		// }
//		// return the contig path;
//		LinkedList<Contig> path = new LinkedList<Contig>();
//		try {
//			// if (adjInternals.contains(tl.getLast())) {
//			// path.addLast(tl.getLast());
//			// } else {
//			// path = this.getInternalPath(external, internal, tl.getLast());
//			// }
//			path = this.getInternalPath(external, internal, tl.getLast());
//			if (path == null || path.size() == 0)
//				return null;
//		} catch (Exception e) {
//			logger.error(this.getClass().getName() + "\t" + e.getMessage());
//		}
//		return path;
//	}

	/**
	 * A method to check which is the best triadlink when there were two best
	 * support triadlink;
	 * 
	 * @param tl
	 * @param uniques
	 * @return
	 */
//	private boolean isBestSupportLink(TriadLink tl, List<Contig> uniques) {
//		// Contig pre = tl.getPrevious();
//		// Contig mid = tl.getMiddle();
//		Contig lst = tl.getLast();
//		int oSL = tl.getSupLinks();
//		// remove the lst contig in the uniques contigs list
//		if (uniques.contains(lst))
//			uniques.remove(lst);
//		// All contigs adjacents to mid;
//		for (Contig c : uniques) {
//			TriadLink temp = new TriadLink();
//			temp.setPrevious(lst);
//			// temp.setMiddle(mid);
//			temp.setLast(c);
//			for (TriadLink t : triads) {
//				// System.out.println(t.toString());
//				if (t.equals(temp))
//					temp.setSupLinks(temp.getSupLinks() + t.getSupLinks());
//				if (temp.getSupLinks() > oSL) {
//					uniques.add(lst);
//					return false;
//				}
//			}
//			// System.out.println(temp.toString());
//			// if(temp.getSupLinks() > oSL)
//			// {
//			// uniques.add(lst);
//			// return false;
//			// }
//		}
//		uniques.add(lst);
//		return true;
//	}

//	private void getSupportedTraidLinks(Contig external, Contig internal, Contig unique, Contig former, int depth,
//			List<TriadLink> links) {
//		List<Contig> adjs = diGraph.getAdjVertices(unique);
//		int size = adjs.size();
//		if (size == 1 || depth == 0) { // end point;
//			List<TriadLink> temp = findTriadLinks(external, internal, unique);
//			for (TriadLink t : temp) {
//				if (!links.contains(t))
//					links.add(t);
//			}
//			return;
//		} else if (size == 2) { // linear point;
//			// if the linear point, it should be have depth attribute;
//			// remove divergence point
//			List<Contig> adjsNext = new Vector<Contig>(adjs.size() - 1);
//			Iterator<Contig> it = adjs.iterator();
//			while (it.hasNext()) {
//				Contig c = it.next();
//				if (c.equals(internal))
//					continue;
//				if (c.equals(former))
//					continue;
//				if (diGraph.isDivergenceVertex(c))
//					continue;
//				adjsNext.add(c);
//			}
//			List<TriadLink> temp = findTriadLinks(external, internal, unique);
//			for (TriadLink t : temp) {
//				if (!links.contains(t))
//					links.add(t);
//			}
//			if (adjsNext.size() == 0)
//				return;
//			former = unique;
//			unique = adjsNext.get(0);
//			this.getSupportedTraidLinks(external, internal, unique, former, depth - 1, links);
//		} else if (size > 2) {// divergence point;
//			List<TriadLink> temp = findTriadLinks(external, internal, unique);
//			for (TriadLink t : temp) {
//				if (!links.contains(t))
//					links.add(t);
//			}
//			return;
//		}
//	}

	// finding triadlink by three contig;
//	private List<TriadLink> findTriadLinks(Contig external, Contig internal, Contig unique) {
//		List<TriadLink> links = new Vector<TriadLink>(10);
//		if (triads == null) {
//			TriadLinkReader tlr = new TriadLinkReader(paras);
//			triads = tlr.read();
//		}
//		// considering the next adjacent contig
//		Iterator<TriadLink> it = triads.iterator();
//		TriadLink t = null;
//		while (it.hasNext()) {
//			t = it.next();
//			if (t.isContain(external) && t.isContain(unique)) {
//				links.add(t);
//			}
//		}
//		return links;
//	}

//	private LinkedList<Contig> getInternalPath(Contig external, Contig internal, Contig unique) {
//		List<LinkedList<Contig>> paths = new Vector<LinkedList<Contig>>(5);
//		LinkedList<Contig> path = new LinkedList<Contig>();
//		List<Contig> nextAdjs = diGraph.getNextVertices(internal, external);
//		int depth = 5;
//
//		// checking whether the nextAdjs contains the unique
//		if (nextAdjs != null && nextAdjs.size() > 0) {
//			if (nextAdjs.contains(unique))
//				path.add(unique);
//		}
//		// if it is valid;
//		if (this.validateInternalPath(external, internal, path)) {
//			return path;
//		} else {
//			path = new LinkedList<Contig>();
//		}
//		// for more than one path could get to unique;
//		// get all the paths and then choose the shortest one;
//		for (Contig c : nextAdjs) {
//			boolean isExist = this.getInternalPath(internal, c, depth, unique, path);
//			if (isExist) {
//				paths.add(path);
//				path = new LinkedList<Contig>();
//			}
//		}
//		int sizes = paths.size();
//		// return the shortest one;
//		// if (sizes != 0)
//		// {
//		// int index = 0;
//		// int length = paths.get(0).size();
//		// int next = 0;
//		// for(int i = 1; i < sizes; i++)
//		// {
//		// next = paths.get(i).size();
//		// if(next < length)
//		// {
//		// index = i;
//		// length = next;
//		// }
//		// }
//		// return paths.get(index);
//		// } else
//		// {
//		// return null;
//		// }
//		// if there are more than one path could be arrive at unique contig
//		// then if the path orientation is not conflict, return it;
//
//		if (sizes != 0) {
//			path = new LinkedList<Contig>();
//			for (LinkedList<Contig> p : paths) {
//				if (this.validateInternalPath(external, internal, p)) {
//					path = p;
//					break;
//				}
//			}
//			if (path != null && !path.isEmpty())
//				return path;
//			else
//				return null;
//		} else {
//			return null;
//		}
//	}
//
//	private boolean getInternalPath(Contig previous, Contig current, int depth, Contig unique,
//			LinkedList<Contig> path) {
//		path.addLast(current);
//		List<Contig> nextAdjs = diGraph.getNextVertices(current, previous);
//		boolean isExist = false;
//		if (nextAdjs == null || nextAdjs.size() == 0) {
//			path.removeLast();
//			return false;
//		}
//		if (depth == 0) {
//			path.removeLast();
//			return false;
//		}
//		if (nextAdjs.contains(unique)) {
//			path.addLast(unique);
//			isExist = true;
//			return isExist;
//		}
//		if (nextAdjs.size() > 1) {
//			for (Contig c : nextAdjs) {
//				isExist = this.getInternalPath(current, c, depth - 1, unique, path);
//				if (isExist)
//					break;
//			}
//			if (!isExist)
//				path.removeLast();
//			return isExist;
//		}
//		Contig c = nextAdjs.get(0);
//		isExist = getInternalPath(current, c, depth - 1, unique, path);
//		if (!isExist)
//			path.removeLast();
//		return isExist;
//	}

	/**
	 * This method is used to check whether orientation of internal path is
	 * conflict;
	 * 
	 * @param path
	 * @return
	 */
//	private boolean validateInternalPath(Contig external, Contig internal, List<Contig> path) {
//		boolean isValid = true;
//		if (path == null || path.isEmpty())
//			return false;
//		List<Contig> tpath = new LinkedList<Contig>();
//		tpath.add(external);
//		tpath.add(internal);
//		tpath.addAll(path);
//		Contig current = null;
//		Contig next = null;
//		Strand cStrand = null;
//		int pSize = tpath.size();
//		for (int i = 0; i < pSize - 1; i++) {
//			current = tpath.get(i);
//			next = tpath.get(i + 1);
//			List<Edge> es = diGraph.getEdgesInfo(current, next);
//			Edge e = null;
//			for (Edge t : es) {
//				if (t.getOrigin().equals(current) && t.getTerminus().equals(next)) {
//					e = t;
//					break;
//				}
//			}
//			if (cStrand == null) {
//				cStrand = e.gettStrand();
//				continue;
//			} else {
//				if (!e.getoStrand().equals(cStrand)) {
//					isValid = false;
//					break;
//				} else {
//					cStrand = e.gettStrand();
//				}
//			}
//		}
//		return isValid;
//	}

	/**
	 * A method to get the unique contig aside to divergence contig;
	 * 
	 * @param current
	 *            - the current contig the check for unique contig
	 * @param previous
	 *            - the former contig
	 * @param depth
	 *            - the check depth
	 * @param uniques
	 *            - the list to store the uniques contigs;
	 */
//	private void getNextUniqueContigs(Contig current, Contig previous, int depth, List<Contig> uniques) {
//		// List<Contig> uniques = new Vector<Contig>(5);
//		List<Contig> nextAdjs = diGraph.getNextVertices(current, previous);
//		if (nextAdjs == null || nextAdjs.size() == 0) {
//			// INTERNAL_LENGTH = 0;
//			return;
//		}
//		if (depth == 0) {
//			// INTERNAL_LENGTH = 0;
//			return;
//		}
//		if (INTERNAL_LENGTH > MAXIMUM_INTERNAL_LENGTH) {
//			return;
//		}
//		// two cases: 1) the next point is divergence; 2) the next point is
//		// unique;
//		if (nextAdjs.size() > 1) {
//			List<Edge> egs = diGraph.getEdgesInfo(current, previous);
//			int eLen = egs.get(0).getDistMean();
//			INTERNAL_LENGTH += eLen;
//			// INTERNAL_LENGTH += this.indexLen(current.getID());
//			int cLen = cntfile.getLengthByNewId(Integer.valueOf(current.getID()));
//			INTERNAL_LENGTH += cLen;
//			if (INTERNAL_LENGTH <= MAXIMUM_INTERNAL_LENGTH) {
//				if (uniques != null && uniques.size() > 0)
//					uniques.remove(uniques.size() - 1);
//				for (Contig c : nextAdjs) {
//					if (!uniques.contains(c))
//						uniques.add(c);
//					this.getNextUniqueContigs(c, current, depth - 1, uniques);
//				}
//				INTERNAL_LENGTH -= cLen;
//				INTERNAL_LENGTH -= eLen;
//				return;
//			} else {
//				// current unique contig is divergence
//				if (nextAdjs != null && nextAdjs.size() > 1) {
//					if (cntfile.getLengthByNewId(Integer.valueOf(current.getID())) >= 5000) {
//						// large than 5k, then this contig could be as unique;
//					} else {
//						// if less than 5k, then remove it and used the next
//						// adjacent as unique;
//						if (uniques != null && uniques.size() > 0)
//							uniques.remove(uniques.size() - 1);
//						for (Contig c : nextAdjs)
//							uniques.add(c);
//					}
//				}
//				INTERNAL_LENGTH -= eLen;
//				INTERNAL_LENGTH -= cLen;
//				return;
//			}
//		} else {
//			Contig c = nextAdjs.get(0);
//			List<Edge> egs = diGraph.getEdgesInfo(current, previous);
//			int eLen = egs.get(0).getDistMean();
//			INTERNAL_LENGTH += eLen;
//			// INTERNAL_LENGTH += this.indexLen(current.getID());
//			int cLen = cntfile.getLengthByNewId(Integer.valueOf(current.getID()));
//			INTERNAL_LENGTH += cLen;
//			this.getNextUniqueContigs(c, current, depth - 1, uniques);
//			INTERNAL_LENGTH -= cLen;
//			INTERNAL_LENGTH -= eLen;
//
//		}
//	}

	/**
	 * A method to push cnt as node into path;
	 * 
	 * @param cnt
	 */
//	private void addNode(Contig cnt, NodePath path, boolean isForward) {
//		if (isForward) { // unshift contig on the forward direction;
//			Node node = new Node();
//			node.setCnt(cnt);
//			node.setOrphan(false);
//			path.unshift(node);
//			diGraph.setVertexAsSelected(cnt);
//		} else { // push contig on the backward direction;
//			Node node = new Node();
//			node.setCnt(cnt);
//			node.setOrphan(false);
//			path.push(node);
//			diGraph.setVertexAsSelected(cnt);
//		}
//	}

	/**
	 * A method to push cnt as node into path directly compute orientation and
	 * distance parameter
	 * 
	 * @param previous
	 * @param current
	 * @param next
	 * @param p2ce
	 * @param c2ne
	 * @param path
	 * @param isForward
	 * @param isFirstStart
	 *  
	 */
	private void addNode(Contig previous, Contig current, Contig next, Edge p2ce, Edge c2ne, NodePath path, 
			boolean isForward){
		// if it is forward direction, it will adding next contig into path and
		// their parameter of
		// current contig
		if (isForward) {
			if (next == null) {
				Node node = new Node();
				node.setCnt(current);
				node.setOrphan(false);
				// node.setMeanDist2Next(c2ne.getDistMean());
				// node.setSdDist2Next(c2ne.getDistSd());
				if(p2ce != null)
					node.setStrand(p2ce.gettStrand());
				else
					node.setStrand(Strand.FORWARD);
				path.push(node);
				diGraph.setVertexAsSelected(current);
			} else {
				Node node = new Node();
				node.setCnt(current);
				node.setOrphan(false);
				node.setMeanDist2Next(c2ne.getDistMean());
				node.setSdDist2Next(c2ne.getDistSd());
				node.setStrand(c2ne.getoStrand());
				path.push(node);
				diGraph.setVertexAsSelected(current);
			}
		} else {
			Node node = new Node();
			node.setCnt(current);
			node.setOrphan(false);
			if(p2ce != null)
			{
				node.setMeanDist2Next(p2ce.getDistMean());
				node.setSdDist2Next(p2ce.getDistSd());
				if(p2ce.gettStrand().equals(Strand.FORWARD))
					node.setStrand(Strand.REVERSE);
				else
					node.setStrand(Strand.FORWARD);
			}
			path.unshift(node);
			diGraph.setVertexAsSelected(current);
		}
	}

	/**
	 * The method used to check the next contig is not the divergence point and
	 * did not select before. If it is not the divergence point and do not
	 * seleect former, return true; else it is return false;
	 * 
	 * @return
	 */
	private boolean isValidCnt(Contig cnt) {
		if (diGraph.isDivergenceVertex(cnt)) {
			return true;
		} else {
			if (diGraph.isVertexSelected(cnt))
				return false;
			else
				return true;

		}
	}

	/**
	 * A method is used to check the orientation of three points
	 * 
	 * @param previous
	 * @param current
	 * @param next
	 * @return
	 */
	private boolean validateOrientation(Contig previous, Contig current, Contig next, Edge p2ce, Edge c2ne) {
		if (previous == null) {
			return true;
		} else {
			if (p2ce == null)
				return false;
//			Strand p2cPStrand = p2ce.getoStrand();
			Strand p2cCStrand = p2ce.gettStrand();
			Strand c2nCStrand = c2ne.getoStrand();
//			Strand c2nNStrand = c2ne.gettStrand();
			if (p2cCStrand.equals(c2nCStrand))
				return true;
			else
				return false;
		}
	}

}

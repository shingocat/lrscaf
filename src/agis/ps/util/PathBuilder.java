/*
*File: agis.ps.util.PathBuilder.java
*User: mqin
*Email: mqin@ymail.com
*Date: 2016年1月14日
*/
package agis.ps.util;

import java.io.File;
import java.util.Collections;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Vector;

import org.apache.lucene.analysis.Analyzer;
import org.apache.lucene.analysis.standard.StandardAnalyzer;
import org.apache.lucene.document.Document;
import org.apache.lucene.index.DirectoryReader;
import org.apache.lucene.index.IndexReader;
import org.apache.lucene.queryparser.classic.QueryParser;
import org.apache.lucene.search.IndexSearcher;
import org.apache.lucene.search.Query;
import org.apache.lucene.search.ScoreDoc;
import org.apache.lucene.search.TopDocs;
import org.apache.lucene.store.Directory;
import org.apache.lucene.store.SimpleFSDirectory;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import agis.ps.file.DotGraphFileWriter;
import agis.ps.file.TriadLinkReader;
import agis.ps.graph2.DirectedGraph;
import agis.ps.graph2.Graph;
import agis.ps.link.Edge;
import agis.ps.link.TriadLink;
import agis.ps.link.TriadLinkComparator;
import agis.ps.link2.CntFileEncapsulate;
import agis.ps.path.Node;
import agis.ps.path.NodePath;
import agis.ps.seqs.Contig;

public class PathBuilder {
	public static Logger logger = LoggerFactory.getLogger(PathBuilder.class);
	private static int MAXIMUM_INTERNAL_LENGTH = 5000; // 5000 bp for validating
														// segement duplication;
	private static int INTERNAL_LENGTH = 0; // store the internal length;
	private static int INDEX = 0;
	private List<Edge> edges;
	private Parameter paras;
	private List<TriadLink> triads;
	private Graph diGraph;
	private Directory directory;
	private IndexReader reader;
	private Analyzer analyzer;
	private QueryParser parser;
	private IndexSearcher searcher;
	private NodePath path;
	private CntFileEncapsulate cntfile;

	public PathBuilder() {
		// do nothing;
		this.initCntIndexer();
	}

	public PathBuilder(List<Edge> edges, Parameter paras) {
		this.edges = edges;
		this.paras = paras;
//		this.initCntIndexer();
	}
	
	public PathBuilder(List<Edge> edges, Parameter paras, CntFileEncapsulate cntfile)
	{
		this(edges, paras);
		this.cntfile = cntfile;
	}

	public List<NodePath> buildPath() {
		if (edges == null || edges.size() == 0)
			throw new IllegalArgumentException("PathBuilder ： The Edges could not be empty!");
//		return this.buildPath2(edges, paras);
		return this.buildPathByCntFileEncapsulate();
	}
	
	private List<NodePath> buildPathByCntFileEncapsulate()
	{
		long start = System.currentTimeMillis();
		if (edges == null || edges.size() == 0)
			throw new IllegalArgumentException("PathBuilder: The Edges could not be empty!");
		List<NodePath> paths = new Vector<NodePath>();
		diGraph = null;
		try {
			diGraph = new DirectedGraph(edges, paras, cntfile);
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
			// delete tips
			diGraph.delTips();
			tempEdges = diGraph.getEdges();
			logger.info("Edges size after deleting tips: " + tempEdges.size());
			edgeFile = paras.getOutFolder() + System.getProperty("file.separator") + "edges_after_dt.info";
			DotGraphFileWriter.writeEdge(edgeFile, tempEdges);
			tempEdges = null;
			path = null;
			// TriadLinkReader tlr = new TriadLinkReader(paras);
			// List<TriadLink> triads = tlr.read();
			// travel the graph, random start
			// do not including the divergence end point in the path
			while (diGraph.isExistUnSelectedVertices()) {
//				INDEX++;
//				if(INDEX == 33)
//					logger.debug("breakpoint");
				Contig current = diGraph.getRandomVertex();
				// if the return conting is null and the
				// isExistUnSelectedVertices equal false then break;
				if (current == null)
					break;
				List<Contig> adjs = diGraph.getAdjVertices(current);
				if (adjs == null) {
					diGraph.setVertexAsSelected(current);
					continue;
				}
				int adjsSize = adjs.size();
				// random selected the frist
				if (adjsSize == 0) {
					// orphan contig, only one element in path
					if(this.isValidCnt(current))
					{
						path = new NodePath();
						this.addNode(current, path, false);
						paths.add(path);
					} else
					{
						diGraph.setVertexAsSelected(current);
					}
					
				} else if (adjsSize == 1) {
					// normal start point, always on the linear end point;
					path = new NodePath();
					// Contig next = diGraph.getNextVertex(current, null);
					Contig next = adjs.get(0);

					while (true) {
						if(this.isValidCnt(current))
						{
							this.addNode(current, path, false);
						} else
						{
							break;
						}
						Contig previous = current;
						current = next;
						next = diGraph.getNextVertex(current, previous);
						// for the divergence point which could not determine how-to get next
						if (next == null) { 
							List<Contig> formers = new Vector<Contig>(3);
							if(path.getPathSize() > 2 )
							{
								formers.add(path.getElement(path.getPathSize() - 2).getCnt());
								formers.add(path.getElement(path.getPathSize() - 3).getCnt());
							} else if (path.getPathSize() == 2)
							{
								formers.add(path.getElement(path.getPathSize() - 2).getCnt());
							}
							// List<Contig> selectedCnts =
							// this.getTriadLinkNext3(current, previous);
							List<Contig> selectedCnts = this.getTriadLinkNext4(current, previous, formers);
							if (selectedCnts == null) {
								this.addNode(current, path, false);
								break;
							}
							int size = selectedCnts.size();
							if (size == 1) {
								next = selectedCnts.get(0);
							} else if (size == 2) {
								this.addNode(current, path, false);
								previous = current;
								current = selectedCnts.get(0);
								next = selectedCnts.get(1);
							} else {
								this.addNode(current, path, false);
								for (int i = 0; i < size - 2; i++) {
									current = selectedCnts.get(i);
									this.addNode(current, path, false);
								}
								previous = current;
								current = selectedCnts.get(size - 2);
								next = selectedCnts.get(size - 1);
							}
						}
						if (next.equals(previous)) { // for the linear end point
							this.addNode(current, path, false);
							break;
						}
					}
					paths.add(path);
				} else if (adjsSize == 2) {
					// middle point; normal start point, located in the linear
					// path
					path = new NodePath();
					if(this.isValidCnt(current))
						this.addNode(current, path, true);
					else
						continue;
					Contig startPoint = current;
					Contig c1 = adjs.get(0);
					Contig c2 = adjs.get(1);
					// c1--current--c2
					// directed by both direction;
					// for c1 direction; using c2 as previous point; checking
					// whether c1 is valid point
					// all the element unshift into path;
					// checking the c1 and c2 whether travel
					if(!this.isValidCnt(c1) && !this.isValidCnt(c2))
					{
						paths.add(path);
						continue;
					}
					Contig previous = startPoint;
					current = c1;
					Contig next = diGraph.getNextVertex(current, previous);
					while (true) {
						if (next != null) {
							if (next.equals(startPoint)) {
								this.addNode(current, path, true);
								break;
							}
							if (next.equals(previous)) {
								this.addNode(current, path, true);
								break;
							}
							// check the current contig whether selected befored
							// if true break;
							if(this.isValidCnt(current))
							{
								this.addNode(current, path, true);
							} else
							{
								break;
							}
							previous = current;
							current = next;
							next = diGraph.getNextVertex(current, previous);
						} else {
							List<Contig> formers = new Vector<Contig>(3);
							if(path.getPathSize() > 2 )
							{
								formers.add(path.getElement(1).getCnt());
								formers.add(path.getElement(2).getCnt());
							} else if (path.getPathSize() == 2)
							{
								formers.add(path.getElement(1).getCnt());
							}
							// to get the next point over the divergence point;
							List<Contig> selectedCnts = this.getTriadLinkNext4(current, previous, formers);
							if (selectedCnts == null) {
								this.addNode(current, path, true);
								break;
							}
							int size = selectedCnts.size();
							if (size == 1) {
								next = selectedCnts.get(0);
							} else if (size == 2) {
								this.addNode(current, path, true);
								previous = current;
								current = selectedCnts.get(0);
								next = selectedCnts.get(1);
							} else {
								this.addNode(current, path, true);
								for (int i = 0; i < size - 2; i++) {
									current = selectedCnts.get(i);
									this.addNode(current, path, true);
								}
								previous = current;
								current = selectedCnts.get(size - 2);
								next = selectedCnts.get(size - 1);
							}
						}
					}
					// for c2 direction; using c1 as previous point; checking
					// whether c2 is valid point;
					// all the valid element push into path;
					if (!diGraph.isExistUnSelectedVertices()) {
						paths.add(path);
						break;
					}
					// if it is the loop, then the former current contig equal to 
					// next start contig c2, break;
					if(current != null && c2 != null )
					{
						if(current.equals(c2))
						{
							paths.add(path);
							continue;
						}
					}
					previous = startPoint;
					current = c2;
					next = diGraph.getNextVertex(current, previous);
					while (true) {
						if (next != null) {
							if (next.equals(startPoint)) {
								this.addNode(current, path, false);
								break;
							}
							if (next.equals(previous)) {
								this.addNode(current, path, false);
								break;
							}
							// check the current contig whether selected before
							if(this.isValidCnt(current))
								this.addNode(current, path, false);
							else
								break;
							previous = current;
							current = next;
							next = diGraph.getNextVertex(current, previous);
						} else {
							List<Contig> formers = new Vector<Contig>(3);
							if(path.getPathSize() > 2 )
							{
								formers.add(path.getElement(path.getPathSize() - 2).getCnt());
								formers.add(path.getElement(path.getPathSize() - 3).getCnt());
							} else if (path.getPathSize() == 2)
							{
								formers.add(path.getElement(path.getPathSize() - 2).getCnt());
							}
							List<Contig> selectedCnts = this.getTriadLinkNext4(current, previous, formers);
							if (selectedCnts == null) {
								this.addNode(current, path, false);
								break;
							}
							int size = selectedCnts.size();
							if (size == 1) {
								next = selectedCnts.get(0);
							} else if (size == 2) {
								this.addNode(current, path, false);
								previous = current;
								current = selectedCnts.get(0);
								next = selectedCnts.get(1);
							} else {
								this.addNode(current, path, false);
								for (int i = 0; i < size - 2; i++) {
									current = selectedCnts.get(i);
									this.addNode(current, path, false);
								}
								previous = current;
								current = selectedCnts.get(size - 2);
								next = selectedCnts.get(size - 1);
							}
						}
					}
					paths.add(path);

				} else if (adjsSize > 2) { // divergence point
					diGraph.setVertexAsSelected(current);
					continue;
				}
			}
		} catch (Exception e) {
//			logger.debug("INDEX = " + INDEX);
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
				} else {
					for (int i = 0; i < pathSize - 1; i++) {
//						INDEX++;
//						if(INDEX == 66)
//							logger.debug("breakpoint");
						Node current = np.getElement(i);
						Node next = np.getElement(i + 1);
						Contig cCnt = current.getCnt();
						Contig nCnt = next.getCnt();
						Edge e = null;
						List<Edge> es = diGraph.getEdgesInfo(cCnt, nCnt);
						for (Edge t : es) {
							if (t.getOrigin().equals(cCnt) && t.getTerminus().equals(nCnt))
								e = t;
						}
						if (current.getStrand() == null)
							current.setStrand(e.getoStrand());
						if (next.getStrand() == null)
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
//			logger.debug("INDEX\t" + INDEX);
			logger.error(this.getClass().getName() + "\t" + e.getMessage() + "\t" + e.getClass().getName());
		}
		// cleaning and call gc
		this.diGraph = null;
		this.edges = null;
		this.triads = null;
		System.gc();
		long end = System.currentTimeMillis();
		logger.info("Path Building, erase time: " + (end - start) + " ms");
		return paths;
	}

	private List<NodePath> buildPath2(List<Edge> edges, Parameter paras) {
		long start = System.currentTimeMillis();
		if (edges == null || edges.size() == 0)
			throw new IllegalArgumentException("PathBuilder: The Edges could not be empty!");
		List<NodePath> paths = new Vector<NodePath>();
		diGraph = null;
		try {
			diGraph = new DirectedGraph(edges, paras);
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
			// delete tips
			diGraph.delTips();
			tempEdges = diGraph.getEdges();
			logger.info("Edges size after deleting tips: " + tempEdges.size());
			edgeFile = paras.getOutFolder() + System.getProperty("file.separator") + "edges_after_dt.info";
			DotGraphFileWriter.writeEdge(edgeFile, tempEdges);
			tempEdges = null;
			path = null;
			// TriadLinkReader tlr = new TriadLinkReader(paras);
			// List<TriadLink> triads = tlr.read();
			// travel the graph, random start
			// do not including the divergence end point in the path
			while (diGraph.isExistUnSelectedVertices()) {
//				INDEX++;
//				if(INDEX == 2)
//					logger.debug("breakpoint");
//				Contig current = new Contig();
//				current.setID("981");
				 Contig current = diGraph.getRandomVertex();
				// if the return conting is null and the
				// isExistUnSelectedVertices equal false then break;
				if (current == null)
					break;
				// if(current.getID().equals("1145"))
				// logger.debug("null exception");
				List<Contig> adjs = diGraph.getAdjVertices(current);
				if (adjs == null) {
					diGraph.setVertexAsSelected(current);
					continue;
				}
				int adjsSize = adjs.size();
				// random selected the frist
				if (adjsSize == 0) {
					// orphan contig, only one element in path
					path = new NodePath();
					Node node = new Node();
					node.setCnt(current);
					node.setOrphan(true);
					path.push(node);
					paths.add(path);
					diGraph.setVertexAsSelected(current);
				} else if (adjsSize == 1) {
					// normal start point, always on the linear end point;
					path = new NodePath();
					// Contig next = diGraph.getNextVertex(current, null);
					Contig next = adjs.get(0);
					// if the next and the next next contig is selected, 
					// then this point is not legal point
//					Contig nNext = diGraph.getNextVertex(next, current);
					if(diGraph.isDivergenceVertex(next) && diGraph.isVertexSelected(next))
					{
						Node node = new Node();
						node.setCnt(current);
						node.setOrphan(false);
						diGraph.setVertexAsSelected(current);
						path.push(node);
						paths.add(path);
						continue;
					}
					while (true) {
						Node node = new Node();
						node.setCnt(current);
						node.setOrphan(false);
						diGraph.setVertexAsSelected(current);
						path.push(node);
						node = null;
						Contig previous = current;
						current = next;
						next = diGraph.getNextVertex(current, previous);
						if(next != null){
							if(diGraph.isDivergenceVertex(next) && diGraph.isVertexSelected(next))
							{
								Node iNode = new Node();
								iNode.setCnt(current);
								iNode.setOrphan(false);
								diGraph.setVertexAsSelected(current);
								path.push(iNode);
//								paths.add(path);
								break;
							}
						}
						if (next == null) { // for the divergence point
							List<Contig> formers = new Vector<Contig>(3);
							if(path.getPathSize() > 2 )
							{
								formers.add(path.getElement(path.getPathSize() - 2).getCnt());
								formers.add(path.getElement(path.getPathSize() - 3).getCnt());
							} else if (path.getPathSize() == 2)
							{
								formers.add(path.getElement(path.getPathSize() - 2).getCnt());
							}
							// List<Contig> selectedCnts =
							// this.getTriadLinkNext3(current, previous);
							List<Contig> selectedCnts = this.getTriadLinkNext4(current, previous, formers);
							if (selectedCnts == null) {
								Node iNode = new Node();
								iNode.setCnt(current);
								iNode.setOrphan(false);
								diGraph.setVertexAsSelected(current);
								path.push(iNode);
								next = null;
								current = null;
								previous = null;
								break;
							}
							int size = selectedCnts.size();
							if (size == 1) {
								next = selectedCnts.get(0);
							} else if (size == 2) {
								Node iNode = new Node();
								iNode.setCnt(current);
								iNode.setOrphan(false);
								diGraph.setVertexAsSelected(current);
								path.push(iNode);
								previous = current;
								current = selectedCnts.get(0);
								next = selectedCnts.get(1);
							} else {
								Node iNode = new Node();
								iNode.setCnt(current);
								iNode.setOrphan(false);
								diGraph.setVertexAsSelected(current);
								path.push(iNode);
								for (int i = 0; i < size - 2; i++) {
									current = selectedCnts.get(i);
									iNode = new Node();
									iNode.setCnt(current);
									iNode.setOrphan(false);
									diGraph.setVertexAsSelected(current);
									path.push(iNode);
								}
								previous = current;
								current = selectedCnts.get(size - 2);
								next = selectedCnts.get(size - 1);
							}
							// next = getTriadLinkNext2(current, previous);
							// // checking next is not null statement,
							// // but the next will be not the adjacent vertex
							// if (next != null) {
							// List<Contig> tAdjs =
							// diGraph.getAdjVertices(current);
							// if (!tAdjs.contains(next))
							// next = null;
							// }
							// if (next == null) {
							// next = null;
							// current = null;
							// previous = null;
							// break;
							// }
						}
						if (next.equals(previous)) { // for the linear end point
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
				} else if (adjsSize == 2) {
					// middle point; normal start point, located in the linear
					// path
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
					// for c1 direction; using c2 as previous point; checking
					// whether c1 is valid point
					// all the element unshift into path;
					// checking the c1 and c2 whether travel
					if(diGraph.isVertexSelected(c1) && diGraph.isVertexSelected(c2))
					{
						paths.add(path);
						continue;
					}
					Contig previous = startPoint;
					current = c1;
					diGraph.setVertexAsSelected(current);
					Contig next = diGraph.getNextVertex(current, previous);
					while (true) {
						if (next != null) {
							if (next.equals(startPoint)) {
								node = new Node();
								node.setCnt(current);
								path.unshift(node);
								diGraph.setVertexAsSelected(current);
								break;
							}
							if (next.equals(previous)) {
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
//							if(current.getID().equals("997"))
//								logger.debug("breakpoint");
							next = diGraph.getNextVertex(current, previous);
						} else {
							List<Contig> formers = new Vector<Contig>(3);
							if(path.getPathSize() > 2 )
							{
								formers.add(path.getElement(1).getCnt());
								formers.add(path.getElement(2).getCnt());
							} else if (path.getPathSize() == 2)
							{
								formers.add(path.getElement(1).getCnt());
							}
							List<Contig> selectedCnts = this.getTriadLinkNext4(current, previous, formers);
							if (selectedCnts == null) {
								Node n = new Node();
								n.setCnt(current);
								n.setOrphan(false);
								diGraph.setVertexAsSelected(current);
								path.unshift(n);
								next = null;
								current = null;
								previous = null;
								break;
							}
							int size = selectedCnts.size();
							if (size == 1) {
								next = selectedCnts.get(0);
							} else if (size == 2) {
								Node iNode = new Node();
								iNode.setCnt(current);
								iNode.setOrphan(false);
								diGraph.setVertexAsSelected(current);
								path.unshift(iNode);
								previous = current;
								current = selectedCnts.get(0);
								next = selectedCnts.get(1);
							} else {
								// push right now current node
								Node iNode = new Node();
								iNode.setCnt(current);
								iNode.setOrphan(false);
								diGraph.setVertexAsSelected(current);
								path.unshift(iNode);
								for (int i = 0; i < size - 2; i++) {
									current = selectedCnts.get(i);
									iNode = new Node();
									iNode.setCnt(current);
									iNode.setOrphan(false);
									diGraph.setVertexAsSelected(current);
									path.unshift(iNode);
								}
								previous = current;
								current = selectedCnts.get(size - 2);
								next = selectedCnts.get(size - 1);
							}
							// List<Contig> unique =
							// this.getTriadLinkNext3(current, previous);
							// next = getTriadLinkNext2(current, previous);
							// // checking next is not null statement,
							// // but the next will be not the adjacent vertex
							// if (next != null) {
							// List<Contig> tAdjs =
							// diGraph.getAdjVertices(current);
							// if (!tAdjs.contains(next))
							// next = null;
							// }
							// if (next == null) {
							// diGraph.setVertexAsSelected(current);
							// break;
							// }
						}
					}
					// for c2 direction; using c1 as previous point; checking
					// whether c2 is valid point;
					// all the valid element push into path;
					if (!diGraph.isExistUnSelectedVertices()) {
						paths.add(path);
						break;
					}
					// if it is the loop, then the former current contig equal to 
					// next start contig c2, break;
					if(current != null && c2 != null )
					{
						if(current.equals(c2))
						{
							paths.add(path);
							continue;
						}
					}
					previous = startPoint;
					current = c2;
					next = diGraph.getNextVertex(current, previous);
					while (true) {
						if (next != null) {
							if (next.equals(startPoint)) {
								node = new Node();
								node.setCnt(current);
								path.push(node);
								diGraph.setVertexAsSelected(current);
								break;
							}
							if (next.equals(previous)) {
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
						} else {
							List<Contig> formers = new Vector<Contig>(3);
							if(path.getPathSize() > 2 )
							{
								formers.add(path.getElement(path.getPathSize() - 2).getCnt());
								formers.add(path.getElement(path.getPathSize() - 3).getCnt());
							} else if (path.getPathSize() == 2)
							{
								formers.add(path.getElement(path.getPathSize() - 2).getCnt());
							}
							List<Contig> selectedCnts = this.getTriadLinkNext4(current, previous, formers);
							if (selectedCnts == null) {
								Node n = new Node();
								n.setCnt(current);
								n.setOrphan(false);
								diGraph.setVertexAsSelected(current);
								path.push(n);
								next = null;
								current = null;
								previous = null;
								break;
							}
							int size = selectedCnts.size();
							if (size == 1) {
								next = selectedCnts.get(0);
							} else if (size == 2) {
								Node iNode = new Node();
								iNode.setCnt(current);
								iNode.setOrphan(false);
								diGraph.setVertexAsSelected(current);
								path.push(iNode);
								previous = current;
								current = selectedCnts.get(0);
								next = selectedCnts.get(1);
							} else {
								Node iNode = new Node();
								iNode.setCnt(current);
								iNode.setOrphan(false);
								diGraph.setVertexAsSelected(current);
								path.push(iNode);
								for (int i = 0; i < size - 2; i++) {
									current = selectedCnts.get(i);
									iNode = new Node();
									iNode.setCnt(current);
									iNode.setOrphan(false);
									diGraph.setVertexAsSelected(current);
									path.push(iNode);
								}
								previous = current;
								current = selectedCnts.get(size - 2);
								next = selectedCnts.get(size - 1);
							}
							// List<Contig> unique =
							// this.getTriadLinkNext3(current, previous);
							// next = getTriadLinkNext2(current, previous);
							// // checking next is not null statement,
							// // but the next will be not the adjacent vertex
							// if (next != null) {
							// List<Contig> tAdjs =
							// diGraph.getAdjVertices(current);
							// if (!tAdjs.contains(next))
							// next = null;
							// }
							// if (next == null) {
							// diGraph.setVertexAsSelected(current);
							// break;
							// }
						}
					}
					paths.add(path);

				} else if (adjsSize > 2) { // divergence point
					diGraph.setVertexAsSelected(current);
					continue;
				}
			}
		} catch (Exception e) {
//			logger.debug("INDEX = " + INDEX);
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
				} else {
					for (int i = 0; i < pathSize - 1; i++) {
//						INDEX++;
//						if(INDEX == 66)
//							logger.debug("breakpoint");
						Node current = np.getElement(i);
						Node next = np.getElement(i + 1);
						Contig cCnt = current.getCnt();
						Contig nCnt = next.getCnt();
						Edge e = null;
						List<Edge> es = diGraph.getEdgesInfo(cCnt, nCnt);
						for (Edge t : es) {
							if (t.getOrigin().equals(cCnt) && t.getTerminus().equals(nCnt))
								e = t;
						}
						if (current.getStrand() == null)
							current.setStrand(e.getoStrand());
						if (next.getStrand() == null)
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
//			logger.debug("INDEX\t" + INDEX);
			logger.error(this.getClass().getName() + "\t" + e.getMessage() + "\t" + e.getClass().getName());
		}
		// cleaning and call gc
		this.diGraph = null;
		this.edges = null;
		this.triads = null;
		System.gc();
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
			diGraph = new DirectedGraph(edges, paras);
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
					while (true) {
						Node node = new Node();
						node.setCnt(current);
						node.setOrphan(false);
						diGraph.setVertexAsSelected(current);
						path.push(node);
						node = null;
						Contig previous = current;
						current = next;
						next = diGraph.getNextVertex(current, previous);
						if (next == null) { // for the divergence point
							next = null;
							current = null;
							previous = null;
							break;
						}
						if (next.equals(previous)) { // for the linear end point
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
					// while (next != null) {
					// Node node = new Node();
					// node.setCnt(cnt);
					// node.setOrphan(false);
					// diGraph.setVertexAsSelected(cnt);
					// path.push(node);
					// Contig temp = cnt;
					// cnt = next;
					// next = diGraph.getNextVertex(cnt, temp);
					// if (path.isNextExist(cnt, 0) && path.isNextExist(next,
					// 0))
					// {
					// paths.add(path);
					// break;
					// }
					// if (next.equals(startPoint))
					// count = count + 1;
					// }
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
					// for c1 direction; using c2 as previous point; checking
					// whether c1 is valid point
					// all the element unshift into path;
					Contig previous = startPoint;
					current = c1;
					diGraph.setVertexAsSelected(current);
					Contig next = diGraph.getNextVertex(current, previous);
					while (true) {
						if (next != null) {
							if (next.equals(startPoint)) {
								node = new Node();
								node.setCnt(current);
								path.unshift(node);
								diGraph.setVertexAsSelected(current);
								break;
							}
							if (next.equals(previous)) {
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
						} else {
							diGraph.setVertexAsSelected(current);
							break;
						}
					}
					// for c2 direction; using c1 as previous point; checking
					// whether c2 is valid point;
					// all the valid element push into path;
					if (!diGraph.isExistUnSelectedVertices()) {
						paths.add(path);
						break;
					}
					previous = startPoint;
					current = c2;
					next = diGraph.getNextVertex(current, previous);
					while (true) {
						if (next != null) {
							if (next.equals(startPoint)) {
								node = new Node();
								node.setCnt(current);
								path.unshift(node);
								diGraph.setVertexAsSelected(current);
								break;
							}
							if (next.equals(previous)) {
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
						} else {
							diGraph.setVertexAsSelected(current);
							break;
						}
					}
					paths.add(path);

					// int count = 0;
					// boolean isReverse = false;
					// Contig next = diGraph.getNextVertex(current, null);
					// while (next != null) {
					// Node node = new Node();
					// node.setCnt(current);
					// node.setOrphan(false);
					// diGraph.setVertexAsSelected(current);
					// if (!path.isNextExist(current, 0)) {
					// if (isReverse)
					// path.unshift(node);
					// else
					// path.push(node);
					// }
					// Contig temp = current;
					// current = next;
					// next = diGraph.getNextVertex(current, temp);
					// if (next.equals(startPoint))
					// count += 1;
					// if (current.equals(startPoint))
					// isReverse = true;
					// if (current.equals(startPoint) && count == 2) {
					// paths.add(path);
					// break;
					// }
					// }
				} else if (adjSetSize >= 3) {
					diGraph.setVertexAsSelected(current);
					continue;
				}
			}
		} catch (Exception e) {
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
				} else {
					for (int i = 0; i < pathSize - 1; i++) {
						Node current = np.getElement(i);
						Node next = np.getElement(i + 1);
						Contig cCnt = current.getCnt();
						Contig nCnt = next.getCnt();
						Edge e = null;
						List<Edge> es = diGraph.getEdgesInfo(current.getCnt(), next.getCnt());
						for (Edge t : es) {
							if (t.getOrigin().equals(cCnt) && t.getTerminus().equals(nCnt))
								e = t;
						}
						if (current.getStrand() == null)
							current.setStrand(e.getoStrand());
						if (next.getStrand() == null)
							next.setStrand(e.gettStrand());
						//
						//
						// if (e.getOrigin().equals(current.getCnt())) {
						// current.setStrand(e.getoStrand());
						// next.setStrand(e.gettStrand());
						// } else {
						// // for the previous point
						// if (e.gettStrand().equals(Strand.FORWARD))
						// current.setStrand(Strand.REVERSE);
						// else
						// current.setStrand(Strand.FORWARD);
						// // for the following point
						// if (e.getoStrand().equals(Strand.FORWARD))
						// next.setStrand(Strand.REVERSE);
						// else
						// next.setStrand(Strand.FORWARD);
						// }
						// for the distance, sd and support links
						int meanSum = e.getDistMean();
						int sdSum = e.getDistSd();
						int slSum = e.getLinkNum();
						current.setMeanDist2Next(meanSum);
						current.setSdDist2Next(sdSum);
						current.setSupportLinkNum(slSum);
					}
				}

				/*
				 * else if (pathSize == 2) { Node current = np.getElement(0);
				 * Node next = np.getElement(1); List<Edge> eInfo =
				 * diGraph.getEdgesInfo(current.getCnt(), next.getCnt()); Edge e
				 * = eInfo.get(0); if (e.getOrigin().equals(current.getCnt())) {
				 * current.setStrand(e.getoStrand());
				 * next.setStrand(e.gettStrand()); } else { // for the previous
				 * point if (e.gettStrand().equals(Strand.FORWARD))
				 * current.setStrand(Strand.REVERSE); else
				 * current.setStrand(Strand.FORWARD); // for the following point
				 * if (e.getoStrand().equals(Strand.FORWARD))
				 * next.setStrand(Strand.REVERSE); else
				 * next.setStrand(Strand.FORWARD); } // for the distance, sd and
				 * support links int meanSum = e.getDistMean(); int sdSum =
				 * e.getDistSd(); int slSum = e.getLinkNum(); // if
				 * (eInfo.size() == 1) { // meanSum =
				 * eInfo.get(0).getDistMean(); // sdSum =
				 * eInfo.get(0).getDistSd(); // slSum =
				 * eInfo.get(0).getLinkNum(); // } else { // Edge e2 =
				 * eInfo.get(1); // meanSum = MathTool.mean(new Integer[] {
				 * e.getDistMean(), e2.getDistMean() }); // sdSum =
				 * MathTool.mean(new Integer[] { e.getDistSd(), e2.getDistSd()
				 * }); // slSum = e.getLinkNum() + e2.getLinkNum(); // }
				 * current.setMeanDist2Next(meanSum);
				 * current.setSdDist2Next(sdSum);
				 * current.setSupportLinkNum(slSum); } else { for (int i = 0; i
				 * < pathSize - 1; i++) { Node previous = np.getElement(i); Node
				 * following = np.getElement(i + 1); List<Edge> eInfo =
				 * diGraph.getEdgesInfo(previous.getCnt(), following.getCnt());
				 * Edge e1 = eInfo.get(0); if
				 * (e1.getOrigin().equals(previous.getCnt())) { if
				 * (previous.getStrand() == null) {
				 * previous.setStrand(e1.getoStrand());
				 * following.setStrand(e1.gettStrand()); } else { if
				 * (e1.getoStrand().equals(previous.getStrand())) {
				 * following.setStrand(e1.gettStrand()); } else { // do some
				 * check, since the orientation is // wrong in path! } } } else
				 * { if (previous.getStrand() == null) { if
				 * (e1.gettStrand().equals(Strand.FORWARD))
				 * previous.setStrand(Strand.REVERSE); else
				 * previous.setStrand(Strand.FORWARD); if
				 * (e1.getoStrand().equals(Strand.FORWARD))
				 * following.setStrand(Strand.REVERSE); else
				 * following.setStrand(Strand.FORWARD); } else { if
				 * (!e1.gettStrand().equals(previous.getStrand())) { if
				 * (e1.getoStrand().equals(Strand.FORWARD))
				 * following.setStrand(Strand.REVERSE); else
				 * following.setStrand(Strand.FORWARD); } else { // do some
				 * check, since the orientation is // wrong in path! // it need
				 * to implemented! } } }
				 * 
				 * // distance, sd and supported link int meanSum = 0; int sdSum
				 * = 0; int slSum = 0; if (eInfo.size() == 1) { meanSum =
				 * eInfo.get(0).getDistMean(); sdSum = eInfo.get(0).getDistSd();
				 * slSum = eInfo.get(0).getLinkNum(); } else { Edge e2 =
				 * eInfo.get(1); meanSum = MathTool.mean(new Integer[] {
				 * e1.getDistMean(), e2.getDistMean() }); sdSum =
				 * MathTool.mean(new Integer[] { e1.getDistSd(), e2.getDistSd()
				 * }); slSum = e1.getLinkNum() + e2.getLinkNum(); }
				 * previous.setMeanDist2Next(meanSum);
				 * previous.setSdDist2Next(sdSum);
				 * previous.setSupportLinkNum(slSum); } }
				 */
			}
		} catch (Exception e) {
			logger.error(this.getClass().getName() + "\t" + e.getMessage() + "\t" + e.getClass().getName());
		}
		return paths;
	}

	/**
	 * The fourth method implemented to get next point over the divergence point;
	 * @param internal - the divergence point
	 * @param external - the former point adjacent to divergence point
	 * @param formers - the formers contigs already in the path 
	 * @return
	 */

	private List<Contig> getTriadLinkNext4(Contig internal, Contig external, List<Contig> formers) {
		// contigs adjacents to internal but exclude external;
		List<Contig> adjInternals = diGraph.getNextVertices(internal, external);
		if (adjInternals.size() == 0)
			return null;
		// the unique contig list through path from external to internal;
		List<Contig> uniques = new Vector<Contig>(5);
		// loop the internal adjacent contigs
		for (Contig c : adjInternals) {
			List<Contig> tempUnique = new Vector<Contig>();
			tempUnique.add(c);
			this.getNextUniqueContigs(c, internal, 3, tempUnique);
			for(Contig t : tempUnique)
			{
				if(!uniques.contains(t))
					uniques.add(t);
			}
		}
		// remove unique contig already in path;
		// remove unique contig is directly link former contigs
		List<Contig> tempUniques = new Vector<Contig>(uniques.size());
		for(Contig c : uniques)
		{
			if(!path.isContain(c))
			{
				boolean isExist = false;
				for(Contig in : formers)
				{
					List<Edge> es = diGraph.getEdgesInfo(c, in);
					if(es != null && !es.isEmpty())
					{
						isExist = true;
						break;
					}
				}
				if(!isExist)
				{
					if(!c.equals(internal))
						tempUniques.add(c);
				}
			}
		}
		uniques = null;
		uniques = tempUniques;
		tempUniques = null;
		
		if(uniques.size() == 0)
			return null;
		// candidate triad link
		List<TriadLink> canTls = new Vector<TriadLink>(10);
		TriadLinkComparator tlc = new TriadLinkComparator();
		// store requried deleting triadlinks;
		Map<String, List<TriadLink>> delLinksMap = new HashMap<String, List<TriadLink>>();
		for (Contig unique : uniques) {
			List<TriadLink> tls = new Vector<TriadLink>();
			int depth = 4;
			this.getSupportedTraidLinks(external, internal, unique, internal, depth, tls);
			// considering next external contigs
			List<TriadLink> temp = null;
			for(Contig c : formers)
			{
				if(diGraph.isDivergenceVertex(c))
					continue;
				temp = this.findTriadLinks(c, internal, unique);
				if(temp != null && temp.size() != 0)
				{
					for(TriadLink t : temp)
					{
						if(!tls.contains(t))
							tls.add(t);
					}
				}
			}
			
			// if the inner contig is large than 5000 bp, considering it is unqiue;
//			int inLen = this.indexLen(internal.getID());
			int inLen = cntfile.getLengthByNewId(Integer.valueOf(internal.getID()));
			if(inLen >= 5000)
			{
				temp = this.findTriadLinks(internal, null, unique);
				if(temp != null && !temp.isEmpty())
				{
					for(TriadLink t : temp)
					{
						if(!tls.contains(t))
							tls.add(t);
					}
				}
			}

			TriadLink tl = new TriadLink();
			tl.setPrevious(external);
			tl.setMiddle(internal);
			tl.setLast(unique);

			for (TriadLink t : tls) {
				Contig middle = t.getMiddle();
				Contig last = t.getLast();
				// for special case when the middle is adjacent to internal point, it means
				// this path is circle.
				if(last != null && last.equals(unique))
				{
					if(middle != null && adjInternals.contains(middle))
						continue;
				}
				if(middle != null){
					if (t.isContain(internal)) {
						if (t.getMiddle().equals(internal)) {
							if(t.isContain(unique))
								tl.setSupLinks(tl.getSupLinks() + t.getSupLinks() + 2);
							else
								tl.setSupLinks(tl.getSupLinks() + t.getSupLinks());
						} else {
							if(t.isContain(external) && t.isContain(unique))
								tl.setSupLinks(tl.getSupLinks() - t.getSupLinks());
							else
								tl.setSupLinks(tl.getSupLinks() + t.getSupLinks());
						}
					} else {
						tl.setSupLinks(tl.getSupLinks() + t.getSupLinks());
					}
				} else
				{
					tl.setSupLinks(tl.getSupLinks() + t.getSupLinks());
				}
			}
			canTls.add(tl);
			delLinksMap.put(unique.getID(), tls);
		}

		if (canTls.isEmpty())
			return null;
		Collections.sort(canTls, tlc);
		// to get the best supported links;
		// if there are two best links, it will check the revered direction 
		// and decide which one is the best else return null;
		TriadLink tl = null;
		if (canTls.size() > 1) {
			// right now we only considered if there two best links;
			// the more complex case will implement if there are exist;
			TriadLink first = canTls.get(0);
			TriadLink second = canTls.get(1);
			if (first.getSupLinks() == second.getSupLinks())
			{
				// logic to check which one is the best;
				boolean p1 = isBestSupportLink(canTls.get(0), uniques);
				boolean p2 = isBestSupportLink(canTls.get(1), uniques);
				if(p1 && p2)
				{ // both is the best case;
					return null;
				} else if(p1 && (!p2))
				{ // the first one is the best;
					tl = first;
				} else if((!p1) && p2)
				{ // the second one is the best;
					tl = second;
				} else
				{
					return null;
				}
			}else
			{
				tl = canTls.get(0);
			}
		} else {
			if (canTls.size() != 0)
				tl = canTls.get(0);
			else
				return null;
		}
		if (tl.getSupLinks() == 0)
			return null;
		// remove the supported triadlinks
//		List<TriadLink> dels = delLinksMap.get(tl.getLast().getID());
//		for (TriadLink t : dels) {
//			triads.remove(t);
//		}
		// return the contig path;
		LinkedList<Contig> path = new LinkedList<Contig>();
		if (adjInternals.contains(tl.getLast())) {
			path.addLast(tl.getLast());
			return path;
		} else {
			return this.getInternalPath(external, internal, tl.getLast());
		}
	}
	
	
	/**
	 * A method to check which is the best triadlink when there were two
	 * best support triadlink; 
	 * @param tl
	 * @param uniques
	 * @return
	 */
	private boolean isBestSupportLink(TriadLink tl, List<Contig> uniques)
	{
//		Contig pre = tl.getPrevious();
//		Contig mid = tl.getMiddle();
		Contig lst = tl.getLast();
		int oSL = tl.getSupLinks();
		// remove the lst contig in the uniques contigs list
		if(uniques.contains(lst))
			uniques.remove(lst);
		// All contigs adjacents to mid;
		for(Contig c : uniques)
		{
			TriadLink temp = new TriadLink();
			temp.setPrevious(lst);
//			temp.setMiddle(mid);
			temp.setLast(c);
			for(TriadLink t : triads)
			{
//				System.out.println(t.toString());
				if(t.equals(temp))
					temp.setSupLinks(temp.getSupLinks() + t.getSupLinks());
			}
//			System.out.println(temp.toString());
			if(temp.getSupLinks() > oSL)
			{
				uniques.add(lst);
				return false;
			}
		}
		uniques.add(lst);
		return true;
	}

	/**
	 * The third method for determine which one is the next point to travel;
	 * @param internal
	 * @param external
	 * @return
	 */
	private List<Contig> getTriadLinkNext3(Contig internal, Contig external) {
		List<Contig> nexts = new Vector<Contig>();
		// if (triads == null) {
		// TriadLinkReader tlr = new TriadLinkReader(paras);
		// triads = tlr.read();
		// }
		// contigs adjacents to internal but exclude external;
		List<Contig> adjInternals = diGraph.getNextVertices(internal, external);
		if (adjInternals.size() == 0)
			return null;
		// the unique contig list through path from external to internal;
		List<Contig> uniques = new Vector<Contig>(5);
		// loop the internal adjacent contigs
		for (Contig c : adjInternals) {
			List<Contig> tempUnique = new Vector<Contig>();
			tempUnique.add(c);
			this.getNextUniqueContigs(c, internal, 3, tempUnique);
			uniques.addAll(tempUnique);
		}
		// based on unique contig to travel divergence node;

		// candidate triad link
		List<TriadLink> canTls = new Vector<TriadLink>(10);
		TriadLinkComparator tlc = new TriadLinkComparator();
		// store requried deleting triadlinks;
		Map<String, List<TriadLink>> delLinksMap = new HashMap<String, List<TriadLink>>();
		// for (Contig c : uniques) {
		// TriadLink tl = new TriadLink();
		// tl.setPrevious(external);
		// tl.setLast(c);
		// tl.setSupLinks(0);
		//
		// // the remainder triadlinks for next-next adjacents contig;
		// List<TriadLink> rmdTls = new Vector<TriadLink>(triads.size());
		//
		// // considering the next adjacent contig
		// Iterator<TriadLink> it = triads.iterator();
		// TriadLink t = null;
		// while (it.hasNext()) {
		// t = it.next();
		// // if(t.equals(tl))
		// if (t.isContain(external) && t.isContain(c)) {
		// if (t.isContain(internal)) {
		// // if the triad link contain middle contig,
		// // divide into two case:
		// // the middle contig is internal, ideal case;
		// // else the middle contig is not internal,
		// // then substract the supported links
		// Contig mid = t.getMiddle();
		// if (mid.equals(internal))
		// tl.setSupLinks(tl.getSupLinks() + t.getSupLinks());
		// else
		// tl.setSupLinks(tl.getSupLinks() - t.getSupLinks());
		// } else {
		// tl.setSupLinks(tl.getSupLinks() + t.getSupLinks());
		// }
		// if (delLinksMap.containsKey(c.getID())) {
		// List<TriadLink> tempLink = delLinksMap.get(c.getID());
		// if (!tempLink.contains(t)) {
		// tempLink.add(t);
		// delLinksMap.replace(c.getID(), tempLink);
		// }
		// } else {
		// List<TriadLink> tempLink = new Vector<TriadLink>(10);
		// tempLink.add(t);
		// delLinksMap.put(c.getID(), tempLink);
		// }
		// } else {
		// rmdTls.add(t);
		// }
		// }
		//
		// // considering the next next adjacent contig;
		// // only considering unique conditions;
		// List<Contig> canNextAdjs = diGraph.getNextVertices(c, internal);
		// if (canNextAdjs != null) {
		// // remove internal contig if canNextAdjs size large than 1
		// if (canNextAdjs.size() > 1) {
		// List<Contig> temp = new Vector<Contig>(canNextAdjs.size());
		// for (Contig ct : canNextAdjs) {
		// List<Contig> ctAdjs = diGraph.getNextVertices(ct, c);
		// if (ctAdjs == null || ctAdjs.size() == 1)
		// temp.add(ct);
		// }
		// canNextAdjs.clear();
		// for (Contig ct : temp) {
		// canNextAdjs.add(ct);
		// }
		// }
		// if (canNextAdjs.size() == 1) {
		// Contig cn = canNextAdjs.get(0);
		// // it = triads.iterator();
		// it = rmdTls.iterator();
		// while (it.hasNext()) {
		// t = it.next();
		// if (t.isContain(external) && t.isContain(cn)) {
		// if (t.isContain(internal)) {
		// // if the triad link contain middle contig,
		// // divide into two case:
		// // the middle contig is internal, ideal case;
		// // else the middle contig is not internal,
		// // then substract the supported links
		// Contig mid = t.getMiddle();
		// if (mid.equals(internal))
		// tl.setSupLinks(tl.getSupLinks() + t.getSupLinks());
		// else
		// tl.setSupLinks(tl.getSupLinks() - t.getSupLinks());
		// // it.remove();
		// } else {
		// tl.setSupLinks(tl.getSupLinks() + t.getSupLinks());
		// }
		// // break;
		// // it.remove();
		// if (delLinksMap.containsKey(c.getID())) {
		// List<TriadLink> tempLink = delLinksMap.get(c.getID());
		// if (!tempLink.contains(t)) {
		// tempLink.add(t);
		// delLinksMap.replace(c.getID(), tempLink);
		// }
		// } else {
		// List<TriadLink> tempLink = new Vector<TriadLink>(10);
		// tempLink.add(t);
		// delLinksMap.put(c.getID(), tempLink);
		// }
		// }
		// }
		// }
		// }
		// // adding to candidate triad links
		// // if(tl.getSupLinks() != 0)
		// canTls.add(tl);
		// }

		for (Contig unique : uniques) {
			List<TriadLink> tls = new Vector<TriadLink>();
			int depth = 4;
			this.getSupportedTraidLinks(external, internal, unique, internal, depth, tls);

			TriadLink tl = new TriadLink();
			tl.setPrevious(external);
			tl.setMiddle(internal);
			tl.setLast(unique);

			for (TriadLink t : tls) {
				if (t.isContain(internal)) {
					if (t.getMiddle().equals(internal)) {
						tl.setSupLinks(tl.getSupLinks() + t.getSupLinks());
					} else {
						tl.setSupLinks(tl.getSupLinks() - t.getSupLinks());
					}
				} else {
					tl.setSupLinks(tl.getSupLinks() + t.getSupLinks());
				}
			}
			canTls.add(tl);
			delLinksMap.put(unique.getID(), tls);
		}

		if (canTls.isEmpty())
			return null;
		Collections.sort(canTls, tlc);
		// if the first tl supported links equal to the second then return null;
		TriadLink tl = null;
		if (canTls.size() > 1) {
			if (canTls.get(0).getSupLinks() == canTls.get(1).getSupLinks())
				return null;
			tl = canTls.get(0);
		} else {
			if (canTls.size() != 0)
				tl = canTls.get(0);
			else
				return null;
		}
		if (tl.getSupLinks() == 0)
			return null;
		// remove the supported triadlinks
		List<TriadLink> dels = delLinksMap.get(tl.getLast().getID());
		for (TriadLink t : dels) {
			triads.remove(t);
		}
		// return the contig path;
		LinkedList<Contig> path = new LinkedList<Contig>();
		if (adjInternals.contains(tl.getLast())) {
			path.addLast(tl.getLast());
			return path;
		} else {
			return this.getInternalPath(external, internal, tl.getLast());
		}
	}

	private TriadLink getPesudoTriadLinks(Contig external, Contig internal, Contig unique) {
		TriadLink tl = new TriadLink();
		tl.setPrevious(external);
		tl.setMiddle(internal);
		tl.setLast(unique);

		int depth = 4;
		List<TriadLink> links = new Vector<TriadLink>(10);
		this.getSupportedTraidLinks(external, internal, unique, null, depth, links);
		for (TriadLink t : links) {
			if (t.isContain(internal)) {
				if (t.getMiddle().equals(internal)) {
					tl.setSupLinks(tl.getSupLinks() + t.getSupLinks());
				} else {
					tl.setSupLinks(tl.getSupLinks() - t.getSupLinks());
				}
			} else {
				tl.setSupLinks(tl.getSupLinks() + t.getSupLinks());
			}
		}
		return tl;
	}

	private void getSupportedTraidLinks(Contig external, Contig internal, Contig unique, Contig former, int depth,
			List<TriadLink> links) {
		List<Contig> adjs = diGraph.getAdjVertices(unique);
		int size = adjs.size();
		if (size == 1 || depth == 0) { // end point;
			List<TriadLink> temp = findTriadLinks(external, internal, unique);
			for (TriadLink t : temp) {
				if(!links.contains(t))
					links.add(t);
			}
			return;
		} else if (size == 2) { // linear point;
			// if the linear point, it should be have depth attribute;
			// remove divergence point
			List<Contig> adjsNext = new Vector<Contig>(adjs.size() - 1);
			Iterator<Contig> it = adjs.iterator();
			while (it.hasNext()) {
				Contig c = it.next();
				if (c.equals(internal))
					continue;
				if (c.equals(former))
					continue;
				if (diGraph.isDivergenceVertex(c))
					continue;
				adjsNext.add(c);
			}
			List<TriadLink> temp = findTriadLinks(external, internal, unique);
			for (TriadLink t : temp) {
				if(!links.contains(t))
					links.add(t);
			}
			if (adjsNext.size() == 0)
				return;
			former = unique;
			unique = adjsNext.get(0);
			this.getSupportedTraidLinks(external, internal, unique, former, depth - 1, links);
		} else if (size > 2) {// divergence point;
			List<TriadLink> temp = findTriadLinks(external, internal, unique);
			for (TriadLink t : temp) {
				if(!links.contains(t))
					links.add(t);
			}
			return;
		}
	}

	// finding triadlink by three contig;
	private List<TriadLink> findTriadLinks(Contig external, Contig internal, Contig unique) {
		List<TriadLink> links = new Vector<TriadLink>(10);
		if (triads == null) {
			TriadLinkReader tlr = new TriadLinkReader(paras);
			triads = tlr.read();
		}
		// considering the next adjacent contig
		Iterator<TriadLink> it = triads.iterator();
		TriadLink t = null;
		while (it.hasNext()) {
			t = it.next();
			if (t.isContain(external) && t.isContain(unique)) {
				links.add(t);
			}
		}
		return links;
	}

	private LinkedList<Contig> getInternalPath(Contig external, Contig internal, Contig unique) {
		LinkedList<Contig> path = new LinkedList<Contig>();
		List<Contig> nextAdjs = diGraph.getNextVertices(internal, external);
		int depth = 5;
		for (Contig c : nextAdjs) {
			this.getInternalPath(internal, c, depth, unique, path);
			if (path.contains(unique))
				break;
		}
		if (path.size() != 0)
			return path;
		else
			return null;
	}

	private void getInternalPath(Contig previous, Contig current, int depth, Contig unique, LinkedList<Contig> path) {
		path.addLast(current);
		List<Contig> nextAdjs = diGraph.getNextVertices(current, previous);
		if (nextAdjs == null || nextAdjs.size() == 0) {
			path.clear();
			return;
		}
		if (depth == 0) {
			path.clear();
			return;
		}
		if (nextAdjs.contains(unique)) {
			path.addLast(unique);
			return;
		}
		if (nextAdjs.size() > 1) {
			path.clear();
			return;
		}
		Contig c = nextAdjs.get(0);
		getInternalPath(current, c, depth - 1, unique, path);
	}

	/**
	 * A method to get the unique contig aside to divergence contig;
	 * @param current - the current contig the check for unique contig
	 * @param previous - the former contig 
	 * @param depth - the check depth
	 * @param uniques - the list to store the uniques contigs;
	 */
	private void getNextUniqueContigs(Contig current, Contig previous, int depth, List<Contig> uniques) {
		// List<Contig> uniques = new Vector<Contig>(5);
		List<Contig> nextAdjs = diGraph.getNextVertices(current, previous);
		if (nextAdjs == null || nextAdjs.size() == 0) {
			INTERNAL_LENGTH = 0;
			return;
		}
		if (depth == 0) {
			INTERNAL_LENGTH = 0;
			return;
		}
		// two cases: 1) the next point is divergence; 2) the next point is
		// unique;
		if (nextAdjs.size() > 1) {
			List<Edge> egs = diGraph.getEdgesInfo(current, previous);
			INTERNAL_LENGTH += egs.get(0).getDistMean();
//			INTERNAL_LENGTH += this.indexLen(current.getID());
			INTERNAL_LENGTH += cntfile.getLengthByNewId(Integer.valueOf(current.getID()));
			if (INTERNAL_LENGTH <= MAXIMUM_INTERNAL_LENGTH) {
				uniques.clear();
				for (Contig c : nextAdjs) {
					uniques.add(c);
				}
				INTERNAL_LENGTH = 0;
				return;
			} else {
				// do nothing, the first next point is unique;
				INTERNAL_LENGTH = 0;
				return;
			}
		} else {
			Contig c = nextAdjs.get(0);
			List<Edge> egs = diGraph.getEdgesInfo(current, previous);
			INTERNAL_LENGTH += egs.get(0).getDistMean();
//			INTERNAL_LENGTH += this.indexLen(current.getID());
			INTERNAL_LENGTH += cntfile.getLengthByNewId(Integer.valueOf(current.getID()));
			this.getNextUniqueContigs(c, current, depth - 1, uniques);
		}
	}

	private Contig getTriadLinkNext2(Contig internal, Contig external) {
		Contig next = null;
		if (triads == null) {
			TriadLinkReader tlr = new TriadLinkReader(paras);
			triads = tlr.read();
		}
		// the remainder triadlinks for next-next adjacents contig;
		List<TriadLink> rmdTls = new Vector<TriadLink>(triads.size());
		// candidate triad link
		List<TriadLink> canTls = new Vector<TriadLink>(10);
		// adjacent contig of internal contig
		List<Contig> adjs = diGraph.getAdjVertices(internal);
		// candidate adjacent contig exclude external;
		List<Contig> canAdjs = new Vector<Contig>(10);
		TriadLinkComparator tlc = new TriadLinkComparator();
		// store requried deleting triadlinks;
		Map<String, List<TriadLink>> delLinksMap = new HashMap<String, List<TriadLink>>();
		List<TriadLink> delLinks = new Vector<TriadLink>(10);
		// initiate canAdjs
		for (Contig c : adjs) {
			if (!c.equals(external))
				canAdjs.add(c);
		}
		// check whether exist triadlink to support
		for (Contig c : canAdjs) {
			TriadLink tl = new TriadLink();
			tl.setPrevious(external);
			tl.setMiddle(internal);
			tl.setLast(c);
			tl.setSupLinks(0);
			// considering the next adjacent contig
			Iterator<TriadLink> it = triads.iterator();
			while (it.hasNext()) {
				TriadLink t = it.next();
				// if(t.equals(tl))
				if (t.isContain(external) && t.isContain(c)) {
					if (t.isContain(internal)) {
						// if the triad link contain middle contig,
						// divide into two case:
						// the middle contig is internal, ideal case;
						// else the middle contig is not internal,
						// then substract the supported links
						Contig mid = t.getMiddle();
						if (mid.equals(internal))
							tl.setSupLinks(tl.getSupLinks() + t.getSupLinks());
						else
							tl.setSupLinks(tl.getSupLinks() - t.getSupLinks());
						// it.remove();
					} else {
						tl.setSupLinks(tl.getSupLinks() + t.getSupLinks());
					}
					// break;
					// it.remove();
					if (delLinksMap.containsKey(c.getID())) {
						List<TriadLink> temp = delLinksMap.get(c.getID());
						if (!temp.contains(t)) {
							temp.add(t);
							delLinksMap.replace(c.getID(), temp);
						}
					} else {
						List<TriadLink> temp = new Vector<TriadLink>(10);
						temp.add(t);
						delLinksMap.put(c.getID(), temp);
					}
				} else {
					rmdTls.add(t);
				}
			}

			// considering the next next adjacent contig;
			// only considering unique conditions;
			List<Contig> nextAdjs = diGraph.getAdjVertices(c);
			List<Contig> canNextAjds = new Vector<Contig>(3);
			for (Contig cn : nextAdjs) {
				if (!cn.equals(internal))
					canNextAjds.add(cn);
			}
			if (canNextAjds.size() == 1) {
				Contig cn = canNextAjds.get(0);
				// it = triads.iterator();
				it = rmdTls.iterator();
				while (it.hasNext()) {
					TriadLink t = it.next();
					if (t.isContain(external) && t.isContain(cn)) {
						if (t.isContain(internal)) {
							// if the triad link contain middle contig,
							// divide into two case:
							// the middle contig is internal, ideal case;
							// else the middle contig is not internal,
							// then substract the supported links
							Contig mid = t.getMiddle();
							if (mid.equals(internal))
								tl.setSupLinks(tl.getSupLinks() + t.getSupLinks());
							else
								tl.setSupLinks(tl.getSupLinks() - t.getSupLinks());
							// it.remove();
						} else {
							tl.setSupLinks(tl.getSupLinks() + t.getSupLinks());
						}
						// break;
						// it.remove();
						if (delLinksMap.containsKey(c.getID())) {
							List<TriadLink> temp = delLinksMap.get(c.getID());
							if (!temp.contains(t)) {
								temp.add(t);
								delLinksMap.replace(c.getID(), temp);
							}
						} else {
							List<TriadLink> temp = new Vector<TriadLink>(10);
							temp.add(t);
							delLinksMap.put(c.getID(), temp);
						}
					}
				}
			}
			// adding to candidate triad links
			// if(tl.getSupLinks() != 0)
			canTls.add(tl);
		}

		if (canTls.isEmpty())
			return null;
		Collections.sort(canTls, tlc);
		TriadLink tl = canTls.get(0);
		// remove the supported triadlinks
		List<TriadLink> dels = delLinksMap.get(tl.getLast().getID());
		for (TriadLink t : dels) {
			triads.remove(t);
		}
		return tl.getLast();
	}
	
	/**
	 * A method to push cnt as node into path;
	 * @param cnt
	 */
	private void addNode(Contig cnt, NodePath path, boolean isForward)
	{
		if(isForward)
		{ // unshift contig on the forward direction;
			Node node = new Node();
			node.setCnt(cnt);
			node.setOrphan(false);
			path.unshift(node);
			diGraph.setVertexAsSelected(cnt);
		} else
		{ // push contig on the backward direction;
			Node node = new Node();
			node.setCnt(cnt);
			node.setOrphan(false);
			path.push(node);
			diGraph.setVertexAsSelected(cnt);
		}
	}
	
	/**
	 * The method used to check the next contig is not
	 * the divergence point and did not select before.
	 * If it is not the divergence point and do not seleect
	 * former, return true; else it is return false;
	 * @return
	 */
	private boolean isValidCnt(Contig cnt)
	{
		if(diGraph.isDivergenceVertex(cnt)){
			return true;
		} else
		{
			if(diGraph.isVertexSelected(cnt))
				return false;
			else
				return true;
			
		}
	}

	// c1 for middle and c2 for two adjacent;
	public Contig getTriadLinkNext(Contig internal, Contig external) {
		Contig next = null;
		boolean isIdeal = true;
		if (triads == null) {
			TriadLinkReader tlr = new TriadLinkReader(paras);
			triads = tlr.read();
		}
		List<TriadLink> tls = new Vector<TriadLink>();
		List<Contig> adjs = diGraph.getAdjVertices(internal);
		List<Contig> tempAdjs = new Vector<Contig>(10);
		TriadLinkComparator tlc = new TriadLinkComparator();
		// remove the external contig in adjs;
		Iterator<Contig> itc = adjs.iterator();
		while (itc.hasNext()) {
			Contig c = itc.next();
			if (!c.equals(external))
				// itc.remove();
				tempAdjs.add(c);
		}
		// divide into two cases;
		// ideal case: the divergence contig is internal contig in triadlink and
		// the external
		// contig is beside to the internal contig;
		// non-ideal case: the triad link is missing the internal contig, but
		// connect to the
		// other external contig;
		// the ideal case:
		for (TriadLink tl : triads) {
			Contig pre = tl.getPrevious();
			Contig mid = tl.getMiddle();
			Contig lst = tl.getLast();
			if (!internal.equals(mid))
				continue;
			if (external.equals(pre) || external.equals(lst)) {
				tls.add(tl);
			}
		}

		if (!tls.isEmpty()) { // ideal case if the triad link is not empty
			// sort the collection from large to small based on supported link
			Collections.sort(tls, tlc);
			List<TriadLink> ctls = new Vector<TriadLink>(10);
			// merge the same path triad link
			for (TriadLink t1 : tls) {
				Contig pre = t1.getPrevious();
				Contig mid = t1.getMiddle();
				Contig lst = t1.getLast();

				// temp contig to check adjacent contigs;
				Contig temp = null;
				if (pre.equals(external)) {
					temp = lst;
				} else {
					temp = pre;
				}
				List<Contig> tAdjs = diGraph.getAdjVertices(temp);
				List<Contig> cAdjs = new Vector<Contig>(10);
				for (Contig c : tAdjs) {
					if (!c.equals(mid))
						cAdjs.add(c);
				}
				// check whether there are candidate path is in the triad link;
				for (TriadLink t2 : tls) {
					if (t1.equals(t2))
						continue;
					Contig p2 = t2.getPrevious();
					Contig m2 = t2.getMiddle();
					Contig l2 = t2.getLast();
					Contig c2 = null;
					if (p2.equals(external))
						c2 = l2;
					else
						c2 = p2;
					if (cAdjs.contains(c2)) {
						t1.setSupLinks(t1.getSupLinks() + t2.getSupLinks());
						triads.remove(t2);
						tls.remove(t2);
					}
				}
			}
			// sort again
			Collections.sort(tls, tlc);
			// using the maximum triad link as the best;
			TriadLink tl = tls.get(0);
			if (tl.getPrevious().equals(external))
				next = tl.getLast();
			else if (tl.getLast().equals(external))
				next = tl.getPrevious();
			triads.remove(tl);
		} else {
			// try non-ideal case, if ideal case does not exist;
			isIdeal = false;
			// if triads contain external have one of other contig;
			for (TriadLink tl : triads) {
				// if(external.equals(pre) || external.equals(lst))
				if (tl.isContain(external)) {
					for (Contig c : tempAdjs) {
						if (tl.isContain(c))
							tls.add(tl);
					}
				}
			}
			// if still null return;
			if (tls.isEmpty())
				return null;
			// not always the maximum triad link is best;
			// also considering the other external contig is the adjacent;
			TriadLink tl = null;
			OUTER: for (TriadLink t : tls) {
				for (Contig c : tempAdjs) {
					if (t.isContain(c)) {
						tl = t;
						break OUTER;
					}
				}
			}
			// if do not exist the best triad link, then the maximum is the best
			if (tl == null)
				tl = tls.get(0);
			// next = tl.getMiddle();
			Iterator<Contig> it = tempAdjs.iterator();
			while (it.hasNext()) {
				Contig c = it.next();
				if (tl.isContain(c))
					next = c;
			}
			triads.remove(tl);
		}
		return next;
	}

	private void initCntIndexer() {
		try {
			directory = new SimpleFSDirectory(
					new File(paras.getOutFolder() + System.getProperty("file.separator") + "cnt.index").toPath());
			reader = DirectoryReader.open(directory);
			searcher = new IndexSearcher(reader);
			analyzer = new StandardAnalyzer();
			parser = new QueryParser("id", analyzer);
		} catch (Exception e) {
			logger.error(this.getClass().getName() + "\t" + e.getMessage());
		}
	}

	private int indexLen(String id) {
		int len = 0;
		try {
			Query query = parser.parse(id);
			TopDocs tds = searcher.search(query, 10);
			for (ScoreDoc sd : tds.scoreDocs) {
				Document doc = searcher.doc(sd.doc);
				len = Integer.valueOf(doc.get("len"));
			}
		} catch (Exception e) {
			logger.error(this.getClass().getName() + "\t" + e.getMessage());
		}
		return len;
	}

}

/*
*File: agis.ps.Digraph.java
*User: mqin
*Email: mqin@ymail.com
*Date: 2016年1月14日
*/
package agis.ps;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.Set;
import java.util.Vector;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import agis.ps.link.ContInOut;
import agis.ps.link.Edge;
import agis.ps.seqs.Contig;
import agis.ps.util.MathTool;
import agis.ps.util.Strand;

//Directed Graph
public class DiGraph implements Serializable {
	private static final long serialVersionUID = 1L;
	private static Logger logger = LoggerFactory.getLogger(DiGraph.class);
	private Collection<String> sletVerts = Collections.synchronizedCollection(new HashSet<String>());
	private Collection<String> unSletVerts = Collections.synchronizedCollection(new HashSet<String>());
	private int verNum; // number of vertices in this directed graph;
	private int edgNum; // number of edges in this directed graph;
	// private HashSet<String> verticesId; // the id of vertices;
	private List<ContInOut> candiVertices; // the id of vertices for candidate
											// graph travel
	private Map<String, List<Edge>> pTMap; // outdegree, adjacency
														// list for
														// vertex v point to;
	private Map<String, List<Edge>> pFMap; // indegree, map of vertex
														// v point
														// from;
	private List<Edge> edges; // the edges list in this graph;
	private int validEdgeNums = 0;

	public DiGraph(List<Edge> edges) {
		if (edges == null || edges.size() == 0)
			throw new IllegalArgumentException("The parametes edges were empty for constructed Graph instance!");
		this.edges = edges;
		// getVerticesId();
		updateGraph();
	}

	// initiated UnSelected Vertices;
	private void initUnSletVerts() {
		// initiated the vertex id
		if (unSletVerts == null)
			unSletVerts = Collections.synchronizedCollection(new HashSet<String>());
		unSletVerts.clear();

		// initiated the above indegree and outdegree map;
		for (int i = 0; i < edges.size(); i++) {
			Edge e = edges.get(i);
			String oId = e.getOrigin().getID();
			String tId = e.getTerminus().getID();
			// add vertices id to hashset;
			unSletVerts.add(oId);
			unSletVerts.add(tId);
		}
	}
	
	// initiated valid edges numbers
	private void initEdgeNums()
	{
		if(validEdgeNums != 0)
			validEdgeNums = 0;
		for(Edge e : edges)
		{
			if(e.isValid())
				validEdgeNums += 1;
		}
	}

	// return all the vertices list based on indegree and outdegree increasing
	// order
	public List<ContInOut> getCandVertices() {
		if (candiVertices == null)
			candiVertices = Collections.synchronizedList(new LinkedList<ContInOut>());
		Map<String, Integer> indegrees = indegrees();
		Map<String, Integer> outdegrees = outdegrees();
		ContInOut cin = null;
		// build a map
		Map<String, Integer[]> values = new HashMap<String, Integer[]>();
		for (String s : indegrees.keySet()) {
			Integer[] value = new Integer[2];
			value[0] = indegrees.get(s);
			value[1] = 0;
			values.put(s, value);
		}
		for (String s : outdegrees.keySet()) {
			if (values.containsKey(s)) {
				Integer[] value = values.get(s);
				value[1] = outdegrees.get(s);
			} else {
				Integer[] value = new Integer[2];
				value[0] = 0;
				value[1] = outdegrees.get(s);
				values.put(s, value);
			}
		}
		// initiated the linkedlist
		for (String s : values.keySet()) {
			cin = new ContInOut();
			cin.setId(s);
			Integer[] value = values.get(s);
			cin.setIndegrees(value[0]);
			cin.setOutdegrees(value[1]);
			candiVertices.add(cin);
		}
		// sorted linked list by indegree and outdegree increasing order;
		Collections.sort(candiVertices, new Comparator<ContInOut>() {
			@Override
			public int compare(ContInOut o1, ContInOut o2) {
				Integer o1I = o1.getIndegrees();
				Integer o2I = o2.getIndegrees();
				Integer o1O = o1.getOutdegrees();
				Integer o2O = o2.getOutdegrees();

				int icomp = o1I.compareTo(o2I);
				if (icomp != 0) {
					return icomp;
				} else {
					return o1O.compareTo(o2O);
				}
			}
		});
		;
		return candiVertices;
	}

	// update the graph when initiated or modification
	private void updateGraph() {
		initUnSletVerts();
		initEdgeNums();
		this.getpFMap();
		this.getpTMap();
		verNum = pTMap.size();
		edgNum = edges.size();
	}

	// remove edges by criterion of support num;
	// this will remove all edges not larger than this values
	public void removeEdge(int supNum) {
		List<Edge> values = new Vector<Edge>(); // storing the valid values;
		for (int i = 0; i < edges.size(); i++) {
			if (edges.get(i).getLinkNum() >= supNum) {
				values.add(edges.get(i));
			}
		}
		this.setEdges(values);
		updateGraph();
	}

	// remove edges by criterion of support vertex id;
	public void removeEdge(String oId, String tTd) {
		List<Edge> values = new Vector<Edge>(); // storing the valid values;
		for (int i = 0; i < edges.size(); i++) {
			if ((!edges.get(i).getOrigin().getID().equals(oId)) | !(edges.get(i).getTerminus().getID().equals(tTd))) {
				values.add(edges.get(i));
			}
		}
		this.setEdges(values);
		updateGraph();
	}

	// remove edges by criterion of support link less or larger than the specify
	// value;
	public void removeEdge(int lower, int upper) {
		List<Edge> values = new Vector<Edge>(); // storing the valid values;
		for (int i = 0; i < edges.size(); i++) {
			Edge edge = edges.get(i);
			int linkNum = edge.getLinkNum();
			if (linkNum >= lower && linkNum <= upper) {
				values.add(edge);
			}
		}
		this.setEdges(values);
		updateGraph();
	}

	// return random vertex id;
	public String getOneRandomVertex(){
//		int random = new Random().nextInt(verNum);
//		String[] ss = pTMap.keySet().toArray(new String[pTMap.size()]);
//		String s = ss[random];
//		return s;
		if(unSletVerts.size() <= 0)
			throw new IllegalStateException("There are not more unselected vertices!");
		int random = new Random().nextInt(unSletVerts.size());
		String [] ss = unSletVerts.toArray(new String[unSletVerts.size()]);
		String s = ss[random];
		return s;
	}

	// return vertex by their indegree or outdegree order
	public String getVertexByOrdering() {
		String value = "";
		List<ContInOut> cins = getCandVertices();
		for (ContInOut c : cins) {
			if (sletVerts.contains(c.getId())) {
				continue;
			} else {
				value = c.getId();
				break;
			}
		}
		return value;
	}

	// add vertex id into sletVerts if used
	public boolean addId2SletVerts(String id) {
		boolean value = sletVerts.add(id);
		// delete vertex id from unSletVertx set if used
		if (value) {
			if (unSletVerts.contains(id))
				unSletVerts.remove(id);
		} else {
			if (sletVerts.contains(id))
				unSletVerts.remove(id);
		}
		return value;
	}

	public Map<String, Integer> minIndegeres() {
		Map<String, Integer> values = new LinkedHashMap<String, Integer>();
		Map<String, Integer> ins = indegrees();
		int min = Collections.min(ins.values());
		for (String s : ins.keySet()) {
			if (ins.get(s) == min)
				values.put(s, min);
		}
		return values;
	}

	public int indegree(String id) {
		if (pFMap.containsKey(id))
			return pFMap.get(id).size();
		else
			throw new IllegalArgumentException("ID of vertices in a Digraph was not exists!");
	}

	public Map<String, Integer> indegrees() {
		Map<String, Integer> values = new LinkedHashMap<String, Integer>();
		for (String s : pFMap.keySet()) {
			values.put(s, pFMap.get(s).size());
		}
		return values;
	}

	// return outdegree num of this vertex id
	public int outdegree(String id) {
		if (pTMap.containsKey(id))
			return pTMap.get(id).size();
		else
			throw new IllegalArgumentException("ID of vertices in a Digraph was not exists!");
	}

	// return outdegrees num of all vertex
	public Map<String, Integer> outdegrees() {
		Map<String, Integer> values = new LinkedHashMap<String, Integer>();
		for (String s : pTMap.keySet()) {
			values.put(s, pTMap.get(s).size());
		}
		return values;
	}

	public int getVerNum() {
		return pTMap.size();
	}

	public void setVerNum(int verNum) {
		this.verNum = verNum;
	}

	public int getEdgNum() {
		return edgNum;
	}

	public void setEdgNum(int edgNum) {
		this.edgNum = edgNum;
	}

	// original vertices list;
	// public HashSet<String> getVerticesId() {
	// // initiated the vertex id
	// if (verticesId == null)
	// verticesId = new HashSet<String>();
	// verticesId.clear();
	//
	// // initiated the above indegree and outdegree map;
	// for (int i = 0; i < edges.size(); i++) {
	// Edge e = edges.get(i);
	// String oId = e.getOrigin().getID();
	// String tId = e.getTerminus().getID();
	// // add vertices id to hashset;
	// verticesId.add(oId);
	// verticesId.add(tId);
	// }
	// return verticesId;
	// }
	// the new method for getting all vertices id, sum of selected and
	// unselected vertices
	public Collection<String> getVerticesId() {
		Collection<String> values = Collections.synchronizedCollection(new HashSet<String>());
		if (!unSletVerts.isEmpty())
			values.addAll(unSletVerts);
		if (!sletVerts.isEmpty())
			values.addAll(sletVerts);
		return values;
	}

//	public void setVerticesId(HashSet<String> verticesId) {
//		this.verticesId = verticesId;
//	}

	public Map<String, List<Edge>> getpTMap() {
		// initiated the vertex point to where
		if (pTMap == null)
			pTMap = Collections.synchronizedMap(new LinkedHashMap<String, List<Edge>>());
		pTMap.clear();
		// initiated the above indegree and outdegree map;
		for (int i = 0; i < edges.size(); i++) {
			Edge e = edges.get(i);
			String oId = e.getOrigin().getID();
			String tId = e.getTerminus().getID();
			// for vertex point to other vertex
			if (pTMap.containsKey(oId)) {
				List<Edge> adjVetices = pTMap.get(oId);
				adjVetices.add(e);
				pTMap.put(oId, adjVetices);
			} else {
				List<Edge> adjVetices = new ArrayList<Edge>();
				adjVetices.add(e);
				pTMap.put(oId, adjVetices);
			}
		}
		return pTMap;
	}

	public void setpTMap(LinkedHashMap<String, List<Edge>> pTMap) {
		this.pTMap = pTMap;
	}

	public Map<String, List<Edge>> getpFMap() {
		// initiated the vertex point from where
		if (pFMap == null)
			pFMap = Collections.synchronizedMap(new LinkedHashMap<String, List<Edge>>());
		pFMap.clear();
		// initiated the above indegree and outdegree map;
		for (int i = 0; i < edges.size(); i++) {
			Edge e = edges.get(i);
			String oId = e.getOrigin().getID();
			String tId = e.getTerminus().getID();
			// for vertex point from other vertex
			if (pFMap.containsKey(tId)) {
				List<Edge> adjVetices = pFMap.get(tId);
				adjVetices.add(e);
				pFMap.put(tId, adjVetices);
			} else {
				List<Edge> adjVetices = new ArrayList<Edge>();
				adjVetices.add(e);
				pFMap.put(tId, adjVetices);
			}
		}
		return pFMap;
	}

	public void setpFMap(LinkedHashMap<String, List<Edge>> pFMap) {
		this.pFMap = pFMap;
	}

	// return all adjacency vertices over each vertex;
	public Map<String, List<Edge>> getAdjs() {
		return pTMap;
	}

	// outdegree, return all adjacency vertices over the specific vertex id
	// point to;
	public LinkedHashMap<String, List<Edge>> getAdjsPT(String id) {
//		if (!verticesId.contains(id)) {
//			throw new IllegalArgumentException("ID of vertices in a Digraph was not exists!");
//		}
		if(!sletVerts.contains(id) && !unSletVerts.contains(id))
			throw new IllegalArgumentException("ID of vertices in a Digraph was not exists!");
		LinkedHashMap<String, List<Edge>> values = new LinkedHashMap<String, List<Edge>>();
		List<Edge> edges = getpTMap().get(id);
		values.put(id, edges);
		return values;
	}

	// indegere, return all adjacency vertices over the specific vertex id point
	// from other vertex;
	public LinkedHashMap<String, List<Edge>> getAdjsPF(String id) {
//		if (!verticesId.contains(id))
//			throw new IllegalArgumentException("Id of vertex in a digraph was not exists!");
		if(!sletVerts.contains(id) && !unSletVerts.contains(id))
			throw new IllegalArgumentException("Id of vertex in a digraph was not exists!");
		LinkedHashMap<String, List<Edge>> values = new LinkedHashMap<String, List<Edge>>();
		List<Edge> edges = getpFMap().get(id);
		values.put(id, edges);
		return values;
	}

	public void setAdjs(LinkedHashMap<String, List<Edge>> adjs) {
		this.pTMap = adjs;
	}

	// return edges info in this graph
	public LinkedHashMap<String, Integer> getEdgesStatistics() throws Exception {
		LinkedHashMap<String, Integer> values = new LinkedHashMap<String, Integer>();
		int mean = 0; // mean of support number of each edge;
		int sd = 0; // sd of support number of each edge;
		int edgesNum = 0; // the number of edges;
		int min = 0; // the minimum of support number of each edge;
		int max = 0; // the maximum of support number of each edge;
		int upper = 0; // the upper of support number of edges, 95% interval;
		int lower = 0; // the lower of support number of edges, 95% interval;
		ArrayList<Integer> supNums = new ArrayList<Integer>(); // the support
																// number of
																// each edge;
		for (int i = 0; i < edges.size(); i++) {
			supNums.add(edges.get(i).getLinkNum());
		}
		try {
			mean = MathTool.mean(supNums);
			sd = MathTool.sd(supNums);
		} catch (Exception e) {
			logger.debug("DiGraph: " + e.getMessage());
			logger.info("DiGraph: " + e.getMessage());
		}
		edgesNum = supNums.size();
		min = Collections.min(supNums);
		max = Collections.max(supNums);
		upper = mean + 2 * sd;
		lower = mean - 2 * sd;
		values.put("EDGES_COUNTS", edgesNum);
		values.put("SUPPORT_LINKS_MEAN", mean);
		values.put("SUPPORT_LINKS_SD", sd);
		values.put("SUPPORT_LINKS_MIN", min);
		values.put("SUPPORT_LINKS_MAX", max);
		values.put("SUPPORT_LINKS_UPPER", upper);
		values.put("SUPPORT_LINKS_LOWER", lower);
		return values;
	}

	public List<Edge> getEdges() {
		return edges;
	}

	public void setEdges(List<Edge> edges) {
		this.edges = edges;
	}

	public boolean isEdgesEmpty() {
		// if the edges in the digraph is empty, return true;
		if (this.getEdges().isEmpty())
			return true;
		return false;
	}

	// return true or false whether the specified id had edges
	public boolean hasEdges(String id) {
//		if (!verticesId.contains(id))
//			throw new IllegalArgumentException("The specified id do not exist in graph!");
		if (!sletVerts.contains(id) && !unSletVerts.contains(id))
			throw new IllegalArgumentException("The specified id do not exist in graph!");
		boolean value = false;
		if (getAdjsPT(id).size() > 0)
			return true;
		return value;
	}

	// return list of edges origin from the specified id
	public List<Edge> getEdgesBySpecifiedId(String id) {
//		if (!verticesId.contains(id))
//			throw new IllegalArgumentException("The specified id do not exist in graph!");
		if (!sletVerts.contains(id) && !unSletVerts.contains(id))
			throw new IllegalArgumentException("The specified id do not exist in graph!");
		List<Edge> values = new Vector<Edge>();
		values = getpTMap().get(id);
		if (values == null)
			return new Vector<Edge>();
		return values;
	}

	// checking the DiGraph has self connected edge or not;
	public boolean hasSelfConnected() {
		for (int i = 0; i < edges.size(); i++) {
			Contig origin = edges.get(i).getOrigin();
			Contig terminus = edges.get(i).getTerminus();
			if (origin.equals(terminus)) {
				return true;
			}
		}
		return false;
	}

	// return all the self connected edges
	public List<Edge> getSelfConnected() {
		List<Edge> edges = new Vector<Edge>();
		for (int i = 0; i < edges.size(); i++) {
			Contig origin = edges.get(i).getOrigin();
			Contig terminus = edges.get(i).getTerminus();
			if (origin.equals(terminus))
				edges.add(edges.get(i));
		}
		if (edges.size() == 0)
			return null;
		return edges;
	}

	// need to improve the implement!
	public boolean hasValidEdges() {
		for (int i = 0; i < edges.size(); i++) {
			if (edges.get(i).isValid())
				return true;
		}
		return false;
	}

	public List<Edge> getPTValidEdgesById(String id) {
//		if (!verticesId.contains(id))
//			throw new IllegalArgumentException("The specified id do not exist in graph!");
		if(!sletVerts.contains(id) && !unSletVerts.contains(id))
			throw new IllegalArgumentException("The specified id do not exist in graph!");
		List<Edge> values = new Vector<Edge>();
		values = getpTMap().get(id);
		// remove the unvalid edges
		List<Edge> reValues = new Vector<Edge>();
		for (Edge e : values) {
			if (e.isValid())
				reValues.add(e);
		}
		return reValues;
	}

	public List<Edge> getPFValidEdgesById(String id) {
//		if (!verticesId.contains(id))
//			throw new IllegalArgumentException("The specified id do not exist in graph!");
		if(!sletVerts.contains(id) && !unSletVerts.contains(id))
			throw new IllegalArgumentException("The specified id do not exist in graph!");
		List<Edge> values = new Vector<Edge>();
		values = getpFMap().get(id);
		// remove the unvalid edges;
		List<Edge> reValues = new Vector<Edge>();
		for (Edge e : values) {
			if (e.isValid())
				reValues.add(e);
		}
		return reValues;
	}

	public int getValidEdgeNums() {
		return validEdgeNums;
	}

	public void setValidEdgeNums(int validEdgeNums) {
		this.validEdgeNums = validEdgeNums;
	}
	
	public void setEdgeUnValid(Edge e)
	{
		if(!edges.contains(e))
			throw new IllegalArgumentException("DiGraph: The argument edge was not exists!");
		int index = edges.indexOf(e);
		edges.get(index).setValid(false);
		validEdgeNums = validEdgeNums - 1;
	}
	
	// adding pesudo edges for two contigs, when this digraph having only one direction edges;
	// it need to run after instance a DiGraph
	public List<Edge> addPesudoEdges()
	{
		List<Edge> values = Collections.synchronizedList(new Vector<Edge>());
		for(int i = 0; i < edges.size(); i++)
		{
			Edge e1 = edges.get(i); // the exist edge in the graph
			Edge e = new Edge(); // the pesudo edge, only turn over the origin and terminus contig
			e.setOrigin(e1.getTerminus());
			e.setTerminus(e1.getOrigin());
			if(!edges.contains(e))
			{
				e.setDistMean(e1.getDistMean());
				e.setDistSd(e1.getDistSd());
				e.setLinkNum(e1.getLinkNum());
				e.setOL(e1.isOL());
				e.setValid(e1.isValid());
				if(e1.getoStrand().equals(Strand.FORWARD) && e1.gettStrand().equals(Strand.REVERSE))
				{
					e.setoStrand(Strand.REVERSE);
					e.settStrand(Strand.REVERSE);
				} else if(e1.getoStrand().equals(Strand.FORWARD) && e1.gettStrand().equals(Strand.REVERSE))
				{
					e.setoStrand(Strand.FORWARD);
					e.settStrand(Strand.REVERSE);
				} else if(e1.getoStrand().equals(Strand.REVERSE) && e1.gettStrand().equals(Strand.FORWARD))
				{
					e.setoStrand(Strand.REVERSE);
					e.settStrand(Strand.FORWARD);
				} else if(e1.getoStrand().equals(Strand.REVERSE) && e1.gettStrand().equals(Strand.REVERSE))
				{
					e.setoStrand(Strand.FORWARD);
					e.settStrand(Strand.FORWARD);
				}
				e.setFake(true);
				values.add(e);
			} 
			values.add(e1);
		}
		edges = values;
		this.updateGraph();
		return values;
	}
}

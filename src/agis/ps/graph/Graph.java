/*
*File: agis.ps.graph.Graph.java
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
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Set;
import java.util.Vector;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import agis.ps.DiGraph;
import agis.ps.Edge;
import agis.ps.link.Contig;
import agis.ps.util.MathTool;

public class Graph implements Serializable {

	private static final long serialVersionUID = 1L;
	private static Logger logger = LoggerFactory.getLogger(Graph.class);
	private List<Edge> edges = new Vector<Edge>();
	private List<Contig> vertices = new Vector<Contig>();
	private List<Contig> selectedVertices = new Vector<Contig>();
	private List<Contig> unselectedVertices = new Vector<Contig>();

	public Graph(List<Edge> edges) {
		if (edges == null || edges.size() == 0)
			throw new IllegalArgumentException("Graph: The edges could not be null or empty when constructed graph!");
		if (edges != null)
			edges = new Vector<Edge>();
		this.edges = edges;
		unselectedVertices = this.getVertices();
	}

	public int getEdgeNum() {
		return this.getEdges().size();
	}

	public int getVertexNum() {
		return this.getVertices().size();
	}

	public List<Edge> getEdges() {
		return edges;
	}

	public List<Contig> getVertices() {
		if (vertices == null)
			vertices = new Vector<Contig>();
		for (Edge e : edges) {
			if (!vertices.contains(e.getOrigin()))
				vertices.add(e.getOrigin());
			if (!vertices.contains(e.getTerminus()))
				vertices.add(e.getTerminus());
		}
		return vertices;
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
			logger.debug("Graph: " + e.getMessage());
			logger.info("Graph: " + e.getMessage());
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

	// return true or false whether the graph has unselected edges
	public boolean hasUnselectedEdges() {
		boolean isTrue = false;

		return isTrue;
	}

	// checking whether there are existed unselected vertex
	public boolean isExistUnSelectedVertices() {
		boolean isExist = true;
		if (unselectedVertices.isEmpty())
			isExist = false;
		return isExist;
	}

	// if an vertex is selected, than set as selected status and add into
	// selectedVertices hashset;
	// and remove this vertex id from unSelectedVertices hashset;
	public void setVertexAsSelected(Contig cnt) {
		if (selectedVertices == null)
			selectedVertices = new Vector<Contig>();
		if (unselectedVertices.contains(cnt))
			unselectedVertices.remove(cnt);
		selectedVertices.add(cnt);
	}

	// return random vertex from graph
	public Contig getRandomVertex() {
		Contig cnt = null;
		if (this.isExistUnSelectedVertices())
			cnt = unselectedVertices.get(0);
		return cnt;
	}
}

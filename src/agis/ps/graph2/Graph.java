/*
*File: agis.ps.graph2.Graph.java
*User: mqin
*Email: mqin@ymail.com
*Date: 2016年7月5日
*/
package agis.ps.graph2;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Vector;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import agis.ps.graph.IUntangler;
import agis.ps.link.Edge;
import agis.ps.seqs.Contig;
import agis.ps.util.MathTool;

public abstract class Graph implements Serializable, IUntangler {

	private static final long serialVersionUID = 1L;
	private static Logger logger = LoggerFactory.getLogger(Graph.class);
	protected Map<String, Edge> edges = new HashMap<String, Edge>();
	protected Map<String, Contig> vertices = new HashMap<String, Contig>();
	protected Map<String, Contig> selectedVertices = new HashMap<String, Contig>();
	protected Map<String, Contig> unselectedVertices = new HashMap<String, Contig>();

	public Graph(List<Edge> edges) {
		if (edges == null || edges.size() == 0)
			throw new IllegalArgumentException(this.getClass().getName() + " , The edges" + " could not be empty!");
		if (this.edges == null)
			this.edges = new HashMap<String, Edge>();
		if (this.vertices == null)
			this.vertices = new HashMap<String, Contig>();
		if (this.selectedVertices == null)
			this.selectedVertices = new HashMap<String, Contig>();
		if (this.unselectedVertices == null)
			this.unselectedVertices = new HashMap<String, Contig>();
		this.edges.clear();
		this.vertices.clear();
		this.selectedVertices.clear();
		this.unselectedVertices.clear();
		for (Edge e : edges) {
			Contig origin = e.getOrigin();
			Contig terminus = e.getTerminus();
			String oId = origin.getID();
			String tId = terminus.getID();
			String id = oId + "->" + tId;
			this.edges.put(id, e);
			this.vertices.put(oId, origin);
			this.vertices.put(tId, terminus);
			this.unselectedVertices.put(oId, origin);
			this.unselectedVertices.put(tId, terminus);
			origin = null;
			terminus = null;
			oId = null;
			tId = null;
			id = null;
		}
	}

	public int getEdgeNum() {
		return this.edges.size();
	}

	public int getVertexNum() {
		return this.vertices.size();
	}

	public List<Edge> getEdges() {
		List<Edge> values = new Vector<Edge>(this.edges.size());
		Collection<Edge> temp = this.edges.values();
		Iterator<Edge> it = temp.iterator();
		while(it.hasNext())
		{
			values.add(it.next());
		}
		return values;
	}

	public List<Contig> getVertices() {
		return (List<Contig>) this.vertices.values();
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
		Iterator<Map.Entry<String, Edge>> it = this.edges.entrySet().iterator();
		while (it.hasNext()) {
			Map.Entry<String, Edge> entry = it.next();
			Edge e = entry.getValue();
			supNums.add(e.getLinkNum());
		}
		mean = MathTool.mean(supNums);
		sd = MathTool.sd(supNums);
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
		if (this.selectedVertices == null)
			this.selectedVertices = new HashMap<String, Contig>();
		if (this.unselectedVertices.containsKey(cnt.getID()))
			this.unselectedVertices.remove(cnt.getID());
		if (!this.selectedVertices.containsKey(cnt.getID()))
			this.selectedVertices.put(cnt.getID(), cnt);
	}

	// return random vertex from graph
	public Contig getRandomVertex() {
		Contig cnt = null;
		if (this.isExistUnSelectedVertices()) {
			Iterator<Contig> it = this.unselectedVertices.values().iterator();
			if (it.hasNext())
				cnt = it.next();
		}
		return cnt;
	}
	
	public boolean isVertexSelected(Contig c)
	{
		String id = c.getID();
		if(selectedVertices.containsKey(id))
			return true;
		else
			return false;
	}
	
	// return whether the contig is divergence contig;
	public abstract boolean isDivergenceVertex(Contig cnt);

	// return all the adjacent vertices to or from this contig;
	public abstract List<Contig> getAdjVertices(Contig cnt);

	// return the next vertex from the current vertex and former vertex;
	public abstract Contig getNextVertex(Contig current, Contig former);
	
	// return the next vertices from the current vertext and former vertext
	// but exclude the former vertex
	public abstract List<Contig> getNextVertices(Contig current, Contig former);

	// return the next vertex from the current vertex and former vertex;
	public abstract Contig getNextVertex2(Contig current, Contig former);

	// return the vertex by specified id;
	public abstract Contig getVertex(String id);

	// return edge info between these two contigs
	public abstract List<Edge> getEdgesInfo(Contig start, Contig end);

	// remove an edge
	public abstract boolean removeEdge(Edge e);
	
	// remove edges;
	public abstract void removeEdges(List<Edge> e);
}

/*
*File: agis.ps.graph2.Graph.java
*User: mqin
*Email: mqin@ymail.com
*Date: 2016年7月5日
*/
package agis.ps.graph;

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

//import org.slf4j.Logger;
//import org.slf4j.LoggerFactory;

import agis.ps.link.Edge;
//import agis.ps.seqs.Contig;
import agis.ps.seqs.Sequence;
import agis.ps.util.MathTool;

public abstract class Graph implements Serializable, IUntangler {

	private static final long serialVersionUID = 1L;
//	private static Logger logger = LoggerFactory.getLogger(Graph.class);
	protected Map<String, Edge> edges = null;
//	protected Map<String, Contig> vertices = new HashMap<String, Contig>();
//	protected Map<String, Contig> selectedVertices = new HashMap<String, Contig>();
//	protected Map<String, Contig> unselectedVertices = new HashMap<String, Contig>();
//	protected List<Contig> vertices = null;
//	protected List<Contig> selectedVertices = null;
//	protected List<Contig> unselectedVertices = null;
	protected List<Sequence> vertices = null;
	protected List<Sequence> selectedVertices = null;
	protected List<Sequence> unselectedVertices = null;

	public Graph(List<Edge> edges) {
		if (edges == null || edges.size() == 0)
			throw new IllegalArgumentException(this.getClass().getName() + " , The edges" + " could not be empty!");
		if (this.edges == null)
			this.edges = new HashMap<String, Edge>();
		if (this.vertices == null)
//			this.vertices = new HashMap<String, Contig>();
//			this.vertices = new ArrayList<Contig>();
			this.vertices = new ArrayList<Sequence>();
		if (this.selectedVertices == null)
//			this.selectedVertices = new HashMap<String, Contig>();
//			this.selectedVertices = new ArrayList<Contig>();
			this.selectedVertices = new ArrayList<Sequence>();
		if (this.unselectedVertices == null)
//			this.unselectedVertices = new HashMap<String, Contig>();
//			this.unselectedVertices = new ArrayList<Contig>();
			this.unselectedVertices = new ArrayList<Sequence>();
		this.edges.clear();
		this.vertices.clear();
		this.selectedVertices.clear();
		this.unselectedVertices.clear();
		for (Edge e : edges) {
			Sequence origin = e.getOrigin();
			Sequence terminus = e.getTerminus();
			String oId = origin.getId();
			String tId = terminus.getId();
			String id = oId + "->" + tId;
			this.edges.put(id, e);
			if(!this.vertices.contains(origin))
				this.vertices.add(origin);
			if(!this.unselectedVertices.contains(origin))
				this.unselectedVertices.add(origin);
			if(!this.vertices.contains(terminus))
				this.vertices.add(terminus);
			if(!this.unselectedVertices.contains(terminus))
				this.unselectedVertices.add(terminus);
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
//		List<Edge> values = new Vector<Edge>(this.edges.size());
//		Collection<Edge> temp = this.edges.values();
//		Iterator<Edge> it = temp.iterator();
//		while(it.hasNext())
//		{
//			values.add(it.next());
//		}
//		return values;
		return new ArrayList<Edge>(edges.values());
	}

	public List<Sequence> getVertices() {
		if(this.vertices == null)
			this.vertices = new ArrayList<Sequence>();
		this.vertices.clear();
		for(Map.Entry<String, Edge> entry : edges.entrySet()) {
			Edge e = entry.getValue();
			if(!this.vertices.contains(e.getOrigin()))
				this.vertices.add(e.getOrigin());
			if(!this.vertices.contains(e.getTerminus()))
				this.vertices.add(e.getTerminus());
		}
		return this.vertices;
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
	public void setVertexAsSelected(Sequence cnt) {
//		if (this.selectedVertices == null)
//			this.selectedVertices = new HashMap<String, Contig>();
//		if (this.unselectedVertices.containsKey(cnt.getID()))
//			this.unselectedVertices.remove(cnt.getID());
//		if (!this.selectedVertices.containsKey(cnt.getID()))
//			this.selectedVertices.put(cnt.getID(), cnt);
		if(this.selectedVertices == null)
			this.selectedVertices = new ArrayList<Sequence>();
		if(this.unselectedVertices.contains(cnt))
			this.unselectedVertices.remove(cnt);
		if(!this.selectedVertices.contains(cnt))
			this.selectedVertices.add(cnt);
	}
	
	// set vertext as unselect;
	public void setVertextUnselected(Sequence cnt)
	{
//		if (this.selectedVertices == null)
//			return;
//		if (this.selectedVertices.containsKey(cnt.getID()))
//			this.selectedVertices.remove(cnt.getID());
//		if (!this.unselectedVertices.containsKey(cnt.getID()))
//			this.unselectedVertices.put(cnt.getID(), cnt);
//		if(this.selectedVertices == null)
//			return;
		if(this.selectedVertices.contains(cnt))
			this.selectedVertices.remove(cnt);
		if(!this.unselectedVertices.contains(cnt))
			this.unselectedVertices.add(cnt);
	}

	// return random vertex from graph
	public Sequence getRandomVertex() {
		Sequence cnt = null;
		if (this.isExistUnSelectedVertices()) {
//			cnt = new ArrayList<Contig>(this.unselectedVertices.values()).get(0);
//			Iterator<Contig> it = this.unselectedVertices.values().iterator();
//			if (it.hasNext())
//				cnt = it.next();
			cnt = this.unselectedVertices.get(0);
		}
		return cnt;
	}
	
	public boolean isVertexSelected(Sequence c)
	{
//		String id = c.getID();
//		if(selectedVertices.containsKey(id))
//			return true;
//		else
//			return false;
		if(this.selectedVertices.contains(c))
			return true;
		else
			return false;
	}
	
	// return whether the contig is divergence contig;
	public abstract boolean isDivergenceVertex(Sequence cnt);

	// return all the adjacent vertices to or from this contig;
	public abstract List<Sequence> getAdjVertices(Sequence cnt);

	// return the next vertex from the current vertex and former vertex;
	public abstract Sequence getNextVertex(Sequence current, Sequence former);
	
	// return the next vertices from the current vertext and former vertext
	// but exclude the former vertex
	public abstract List<Sequence> getNextVertices(Sequence current, Sequence former);

	// return the next vertex from the current vertex and former vertex;
	public abstract Sequence getNextVertex2(Sequence current, Sequence former);

	// return the vertex by specified id;
	public abstract Sequence getVertex(String id);

	// return edge info between these two contigs
	public abstract List<Edge> getEdgesInfo(Sequence start, Sequence end);

	// remove an edge
	public abstract boolean removeEdge(Edge e);
	
	// remove edges;
	public abstract void removeEdges(List<Edge> e);

	public void delSimCntEdges() {
		// TODO Auto-generated method stub
		
	}

}

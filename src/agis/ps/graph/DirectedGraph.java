/*
*File: agis.ps.graph.DirectedGraph.java
*User: mqin
*Email: mqin@ymail.com
*Date: 2016年2月25日
*/
package agis.ps.graph;

import java.io.Serializable;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Vector;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import agis.ps.Edge;
import agis.ps.link.Contig;

public class DirectedGraph extends Graph implements Serializable{
	
	private static final long serialVersionUID = 1L;
	private static Logger logger = LoggerFactory.getLogger(DirectedGraph.class);
	private Map<String, List<Contig>> verPTMap = new HashMap<String, List<Contig>>();
	
	public DirectedGraph(List<Edge> edges) {
		super(edges);
		// TODO Auto-generated constructor stub
	}
	
	public int getVertexAdjVerticesNum(Contig cnt)
	{
		int num = 0;
		
		return num;
	}
	
	public List<Contig> getVertexAdjVertices(Contig cnt)
	{
		List<Contig> adjs = new Vector<Contig>();
		
		return adjs;
	}
}



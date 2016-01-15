/*
*File: agis.ps.Digraph.java
*User: mqin
*Email: mqin@ymail.com
*Date: 2016年1月14日
*/
package agis.ps;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.Set;

import agis.ps.util.MathTool;

//Directed Graph
public class DiGraph
{
	private int verNum; // number of vertices in this directed graph;
	private int edgNum; // number of edges in this directed graph;
	private HashSet<String> verticesId; // the id of vertices;
	private LinkedHashMap<String, List<Edge>> pTMap; // adjacency list for vertex v point to;
	private LinkedHashMap<String, List<Edge>> pFMap; // map of vertex v point from;
	private List<Edge> edges; // the edges list in this graph;
	
	public DiGraph(int verNum)
	{
		if(verNum < 0)
			throw new IllegalArgumentException("Number of vertices in a digraph must be nonnegative!");
		this.verNum = verNum;
	}
	
	public DiGraph(int verNum, int edgNum)
	{
		this(verNum);
		this.edgNum = edgNum;
	}
	
	public DiGraph(List<Edge> edges)
	{
		if(edges == null || edges.size() == 0)
			throw new IllegalArgumentException("The parametes is empty for constructed Graph instance!");
		this.edges = edges;
		// initiated the vertex point to where
		if(pTMap == null)
			pTMap = new LinkedHashMap<String, List<Edge>>();
		pTMap.clear();
		// initiated the vertex point from where
		if(pFMap == null)
			pFMap = new LinkedHashMap<String, List<Edge>>();
		pFMap.clear();
		// initiated the vertex id 
		if(verticesId == null)
			verticesId = new HashSet<String>();
		verticesId.clear();
		// initiated the above indegree and outdegree map;
		for(int i = 0; i < edges.size(); i++)
		{
			Edge e = edges.get(i);
			String oId = e.getOrigin().getID();
			String tId = e.getTerminus().getID();
			//add vertices id to hashset;
			verticesId.add(oId);
			verticesId.add(tId);
			//for vertex point to other vertex
			if(pTMap.containsKey(oId))
			{
				List<Edge> adjVetices = pTMap.get(oId);
				adjVetices.add(e);
				pTMap.put(oId, adjVetices);
			} else
			{
				List<Edge> adjVetices = new ArrayList<Edge>();
				adjVetices.add(e);
				pTMap.put(oId, adjVetices);
			}
			//for vertex point from other vertex
			if(pFMap.containsKey(tId))
			{
				List<Edge> adjVetices = pFMap.get(tId);
				adjVetices.add(e);
				pFMap.put(tId, adjVetices);
			} else
			{
				List<Edge> adjVetices = new ArrayList<Edge>();
				adjVetices.add(e);
				pFMap.put(tId, adjVetices);
			}
		}
		this.verNum = pTMap.size();
		// compute edges num in this graph
		int count = 0;
		for(String s : pTMap.keySet())
		{
			count += pTMap.get(s).size();
		}
		this.edgNum = count;
	}
	
	// return random vertex id;
	public String getOneRandomVertex()
	{
		int random = new Random().nextInt(verNum);
		String [] ss = pTMap.keySet().toArray(new String[pTMap.size()]);
		String s = ss[random];
		return s;
	}
	
	public Map<String, Integer> minIndegeres()
	{
		Map<String, Integer> values = new LinkedHashMap<String, Integer>();
		Map<String, Integer> ins = indegrees();
		int min = Collections.min(ins.values());
		for(String s : ins.keySet())
		{
			if(ins.get(s) == min)
				values.put(s, min);
		}
		return values;
	}
	
	public int indegree(String id)
	{
		if(pFMap.containsKey(id))
			return pFMap.get(id).size();
		else
			throw new IllegalArgumentException("ID of vertices in a Digraph was not exists!");
	}
	
	public Map<String, Integer> indegrees()
	{
		Map<String, Integer> values = new LinkedHashMap<String, Integer>();
		for(String s : pFMap.keySet())
		{
			values.put(s, pFMap.get(s).size());
		}
		return values;
	}
	// return outdegree num of this vertex id
	public int outdegree(String id)
	{
		if(pTMap.containsKey(id))
			return pTMap.get(id).size();
		else
			throw new IllegalArgumentException("ID of vertices in a Digraph was not exists!");
	}
	// return outdegrees num of all vertex
	public Map<String, Integer> outdegrees()
	{
		Map<String, Integer> values = new LinkedHashMap<String, Integer>();
		for(String s : pTMap.keySet())
		{
			values.put(s, pTMap.get(s).size());
		}
		return values;
	}

	public int getVerNum() {
		return verNum;
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
	
	public HashSet<String> getVerticesId() {
		return verticesId;
	}

	public void setVerticesId(HashSet<String> verticesId) {
		this.verticesId = verticesId;
	}

	public LinkedHashMap<String, List<Edge>> getpTMap() {
		return pTMap;
	}

	public void setpTMap(LinkedHashMap<String, List<Edge>> pTMap) {
		this.pTMap = pTMap;
	}

	public LinkedHashMap<String, List<Edge>> getpFMap() {
		return pFMap;
	}

	public void setpFMap(LinkedHashMap<String, List<Edge>> pFMap) {
		this.pFMap = pFMap;
	}

	//return all adjacency vertices over each vertex; 
	public LinkedHashMap<String, List<Edge>> getAdjs() {
		return pTMap;
	}
	//return all adjacency vertices over the specific vertex id;
	public LinkedHashMap<String, List<Edge>> getAdjs(String id)
	{
		if(!this.getVerticesId().contains(id))
		{
			throw new IllegalArgumentException("ID of vertices in a Digraph was not exists!");
		}
		LinkedHashMap<String, List<Edge>> values = new LinkedHashMap<String, List<Edge>>();
		List<Edge> edges = pFMap.get(id);
		values.put(id, edges);
		return values;
	}

	public void setAdjs(LinkedHashMap<String, List<Edge>> adjs) {
		this.pTMap = adjs;
	}
	
	public LinkedHashMap<String, Integer> getEdgesStatistics()
	{
		LinkedHashMap<String, Integer> values = new LinkedHashMap<String, Integer>();
		int mean = 0; // mean of support number of each edge;
		int sd = 0; // sd of support number of each edge;
		int edgesNum = 0; // the number of edges;
		int min = 0; // the minimum of support number of each edge;
		int max = 0; // the maximum of support number of each edge;
		ArrayList<Integer> supNums = new ArrayList<Integer>(); // the support number of each edge;
		for(int i = 0; i < edges.size(); i++)
		{
			supNums.add(edges.get(i).getLinkNum());
		}
		mean = MathTool.mean(supNums);
		sd = MathTool.sd(supNums);
		edgesNum = supNums.size();
		min = Collections.min(supNums);
		max = Collections.max(supNums);
		values.put("EDGES_COUNTS", edgesNum);
		values.put("SUPPORT_LINKS_MEAN", mean);
		values.put("SUPPORT_LINKS_SD", sd);
		values.put("SUPPORT_LINKS_MIN", min);
		values.put("SUPPORT_LINKS_MAX", max);
		return values;
	}
	
}

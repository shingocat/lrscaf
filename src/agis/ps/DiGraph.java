/*
*File: agis.ps.Digraph.java
*User: mqin
*Email: mqin@ymail.com
*Date: 2016年1月14日
*/
package agis.ps;

import java.util.ArrayList;
import java.util.Collections;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

import agis.ps.util.MathTool;

//Directed Graph
public class DiGraph
{
	private int verNum; // number of vertices in this directed graph;
	private int edgNum; // number of edges in this directed graph;
	private LinkedHashMap<String, List<Edge>> pTMap; // adjacency list for vertex v point to;
	private LinkedHashMap<String, List<Edge>> pFMap; // map of vertex v point from;
	
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
		// initiated the vertex point to where
		if(pTMap == null)
			pTMap = new LinkedHashMap<String, List<Edge>>();
		pTMap.clear();
		// initiated the vertex point from where
		if(pFMap == null)
			pFMap = new LinkedHashMap<String, List<Edge>>();
		pFMap.clear();
		// initiated the above indegree and outdegree map;
		for(int i = 0; i < edges.size(); i++)
		{
			Edge e = edges.get(i);
			String oId = e.getOrigin().getID();
			String tId = e.getTerminus().getID();
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
	
	public int outdegree(String id)
	{
		if(pTMap.containsKey(id))
			return pTMap.get(id).size();
		else
			throw new IllegalArgumentException("ID of vertices in a Digraph was not exists!");
	}
	
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

	public LinkedHashMap<String, List<Edge>> getAdjs() {
		return pTMap;
	}

	public void setAdjs(LinkedHashMap<String, List<Edge>> adjs) {
		this.pTMap = adjs;
	}
	
	
}

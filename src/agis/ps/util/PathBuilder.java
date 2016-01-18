/*
*File: agis.ps.util.PathBuilder.java
*User: mqin
*Email: mqin@ymail.com
*Date: 2016年1月14日
*/
package agis.ps.util;

import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import agis.ps.DiGraph;
import agis.ps.Edge;
import agis.ps.Path;
import agis.ps.Vertex;

public class PathBuilder {
	public static Logger logger = LoggerFactory.getLogger(PathBuilder.class);
	
	public static List<Path> buildEulerPath(List<Edge> edges)
	{
		return null;
	}
	
	public static List<Path> buildHamiltonPath(List<Edge> edges)
	{	
		try{
			DiGraph diGraph = new DiGraph(edges);
//			logger.debug("Vertices Num: " + diGraph.getVerNum());
//			logger.debug("Edges Num: " + diGraph.getEdgNum());
//			logger.debug("Indegree...");
//			Map<String, Integer> indegrees = diGraph.indegrees();
//			for(String s : indegrees.keySet())
//			{
//				logger.debug(s + ":" + indegrees.get(s));
//			}
//			logger.debug("Outdegree...");
//			Map<String, Integer> outdegrees = diGraph.outdegrees();
//			for(String s : outdegrees.keySet())
//			{
//				logger.debug(s + ":" + outdegrees.get(s));
//			}
//			logger.debug("smallest indegere...");
//			Map<String, Integer> sIns = diGraph.minIndegeres();
//			for(String s : sIns.keySet())
//			{
//				logger.debug(s + ":" + sIns.get(s));
//			}
			// statistics of edges info
			Map<String, Integer> eStat = diGraph.getEdgesStatistics();
			int lower = eStat.get("SUPPORT_LINKS_LOWER");
			int upper = eStat.get("SUPPORT_LINKS_UPPER");
			logger.debug("lower : " + lower);
			logger.debug("upper : " + upper);
			for(String s : eStat.keySet())
			{
				logger.debug(s + ":" + eStat.get(s));
			}
			// simplest graph by remove lower and upper links number 
			// edges former;
			logger.debug("Before remove");
			List<Edge> oE = diGraph.getEdges();
			for(Edge e : oE)
			{
				logger.debug("B: " + e.getOrigin().getID() + "->" + e.getTerminus().getID() + " : " + e.getLinkNum());
			}
			diGraph.removeEdge(lower, upper);
			logger.debug("After process");
			List<Edge> aE = diGraph.getEdges();
			for(Edge e : aE)
			{
				logger.debug("A: " + e.getOrigin().getID() + "->" + e.getTerminus().getID() + " : " + e.getLinkNum());
			}
			eStat = diGraph.getEdgesStatistics();
			for(String s : eStat.keySet())
			{
				logger.debug(s + ":" + eStat.get(s));
			}
			//remove the support links less than 4
			diGraph.removeEdge(4);
			logger.debug("After process less than 4");
			aE = diGraph.getEdges();
			for(Edge e : aE)
			{
				logger.debug("A: " + e.getOrigin().getID() + "->" + e.getTerminus().getID() + " : " + e.getLinkNum());
			}
			eStat = diGraph.getEdgesStatistics();
			for(String s : eStat.keySet())
			{
				logger.debug(s + ":" + eStat.get(s));
			}
			logger.debug(System.getProperty("user.dir"));
			DotGraphFileWriter.writeEdge(System.getProperty("user.dir") + System.getProperty("file.separator") + "removeLessThan10.txt", diGraph.getEdges());
			// go go go through the graph
			// for random start;
//			String id = diGraph.getOneRandomVertex();
//			logger.debug("id: " + id);
//			LinkedHashMap<String, List<Edge>> values = diGraph.getAdjs(id);
//			while(true)
//			{
//				if(diGraph.getAdjs(id).get(id).size() > 0)
//				{
//					
//				} else
//				{
//					break;
//				}
//			}
			return null;
		} catch(Exception e)
		{
			logger.debug(e.getMessage());
			logger.error(e.getMessage());
			return null;
		}
	}
	
	private void findNextVertex()
	{
		
	}
	
	private void findPreviousVertex()
	{
		
	}
}



/*
*File: agis.ps.Scaffolder2.java
*User: mqin
*Email: mqin@ymail.com
*Date: 2016年7月17日
*/
package agis.ps;

import java.util.List;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import agis.ps.file.DotGraphFileWriter;
import agis.ps.file.OutputFolderBuilder;
import agis.ps.file.ScaffoldWriter;
import agis.ps.link.Edge;
import agis.ps.link2.CntFileEncapsulate;
import agis.ps.path.NodePath;
import agis.ps.util.M4EdgeBundler;
import agis.ps.util.M5EdgeBundler;
import agis.ps.util.Parameter;
import agis.ps.util.PathBuilder;

/**
 * The second scaffolder by using java lucene to index
 * , change logic on handling links building, edges 
 * building and path building, also on build scaffolds; 
 * 
 * A new version by using encapsulate aligned and contig file;
 * 
 * @author mqin
 *
 */
public class Scaffolder2 {
	private static Logger logger = LoggerFactory.getLogger(Scaffolder.class);
	private Parameter paras;
	private List<Edge> edges;
	private CntFileEncapsulate cntfile;
	
	public Scaffolder2(Parameter paras)
	{
		this.paras = paras;
	}
	
	public void scaffolding()
	{
		try
		{
			// building output folder
			if(!buildOutputFolder())
				return;
			// building edges 
			this.buildEdges();
			if(edges == null || edges.size() == 0)
				return;
			this.writeEdgesInfo(edges, false);
			logger.info("Original Edges size: " + edges.size());
			PathBuilder pathBuilder = new PathBuilder(edges, paras, cntfile);
			List<NodePath> paths = pathBuilder.buildPath();
			this.writeNodePathInfo(paths);
			this.writeScaffolds(paras, paths);	
		} catch(Exception e)
		{
			logger.error(this.getClass().getName() + "\t" + e.getMessage());
		}
	}
	
	// building output folder 
	private boolean buildOutputFolder()
	{
		OutputFolderBuilder ofb = new OutputFolderBuilder(paras);
		return ofb.building();
	}
	
	// building edges
	private void buildEdges()
	{
		String type = paras.getType();
		if(type.equalsIgnoreCase("m5"))
		{
			// using m5 format to build edges;
			M5EdgeBundler bundler = new M5EdgeBundler(paras);
			this.edges = bundler.building();
			this.cntfile = bundler.getCntFile();
		} else if(type.equalsIgnoreCase("m4"))
		{
			// using m4 format to build edges;
			M4EdgeBundler bundler = new M4EdgeBundler(paras);
			this.edges = bundler.building();
		} else if(type.equalsIgnoreCase("sam") || type.equalsIgnoreCase("bam"))
		{
			// using sam format to build edges;
			// not implement right now;
		} else 
		{
			logger.info(this.getClass().getName() + "The aligned parameter do not set! only <m5>, <m4>, <sam> or <bam>");
			return;
		}
	}
	
	private void writeEdgesInfo(List<Edge> edges, boolean isPesudo) {
		// write the edges info into file.
		String edgeFile = paras.getOutFolder() + System.getProperty("file.separator") + "edges.info";
		DotGraphFileWriter.writeEdge(edgeFile, edges);
	}
	
	private void writeNodePathInfo(List<NodePath> paths)
	{
		String pathFile = paras.getOutFolder() + System.getProperty("file.separator") + "nodePaths.info";
		DotGraphFileWriter.writeNodePaths(pathFile, paths);
	}

	private void writeScaffolds(Parameter paras, List<NodePath> paths)
	{
		ScaffoldWriter sw = new ScaffoldWriter(paras, paths, cntfile);
		sw.writeByContigFileEncapsulate();
	}
	

}



/*
*File: agis.ps.Scaffolder2.java
*User: mqin
*Email: mqin@ymail.com
*Date: 2016年7月17日
*/
package agis.ps;

import java.util.List;
import java.util.Map;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import agis.ps.file.AlignmentFileReader;
import agis.ps.file.DotGraphFileWriter;
import agis.ps.file.M4Reader;
import agis.ps.file.M5Reader;
import agis.ps.file.OutputFolderBuilder;
import agis.ps.file.PBLinkWriter;
import agis.ps.file.ScaffoldWriter;
import agis.ps.file.TriadLinkWriter;
import agis.ps.link.Edge;
import agis.ps.link.MRecord;
import agis.ps.link.PBLink;
import agis.ps.link.TriadLink;
import agis.ps.path.NodePath;
import agis.ps.util.EdgeBundler;
import agis.ps.util.LinkBuilder;
import agis.ps.util.Parameter;
import agis.ps.util.PathBuilder;
import agis.ps.util.RepeatFinder;
//import sun.util.logging.resources.logging_zh_TW;

/**
 * The scaffolder
 * 
 * @author mqin
 *
 */
public class Scaffolder {
	private static Logger logger = LoggerFactory.getLogger(Scaffolder.class);
	private Parameter paras;
	private Map<String, List<MRecord>> records;
	private List<List<MRecord>> listRecords;
	private Map<String, Integer> cntCovs;
	private Map<String, Integer> cntLens;
	private List<String> repeats;
	private List<PBLink> links;
	private List<TriadLink> triads;
	private List<Edge> edges;
	private List<NodePath> paths;
	
	public Scaffolder(Parameter paras)
	{
		this.paras = paras;
	}
	
	public void scaffolding()
	{
		try
		{
			if(!buildOutputFolder())
				return;
			this.buildLinks();
			this.buildEdges();
			this.buildPaths();
			this.writeNodePathInfo();
			this.writeScaffolds();
		} catch(Exception e)
		{
			logger.error(this.getClass().getName() + "\t" + e.getMessage());
		}
	}
	
	private boolean buildOutputFolder()
	{
		OutputFolderBuilder ofb = new OutputFolderBuilder(paras);
		return ofb.building();
	}
	
	private void buildLinks()
	{
		String type = paras.getType();
		AlignmentFileReader reader = null;
		if(type.equalsIgnoreCase("m5"))
		{
			reader = new M5Reader(paras);
		} else if(type.equalsIgnoreCase("m4"))
		{
			reader = new M4Reader(paras);
			
		} else if(type.equalsIgnoreCase("sam") || type.equalsIgnoreCase("bam"))
		{
			logger.info(this.getClass().getName() + "Do not implement yet!");
			return;
		} else 
		{
			logger.info(this.getClass().getName() + "The aligned parameter do not set! only <m5>, <m4>, <sam> or <bam>");
			return;
		}
//		records = reader.read();
		listRecords = reader.read();
		cntCovs = reader.getCntCoverages();
		cntLens = reader.getCntLengths();
		RepeatFinder rf = new RepeatFinder(paras);
		repeats = rf.findRepeats(cntCovs);
		LinkBuilder linkBuilder = new LinkBuilder(paras);
//		links = linkBuilder.mRecords2Links(records, repeats);
		links = linkBuilder.mRecords2Links(listRecords, repeats);
		PBLinkWriter pblw = new PBLinkWriter(paras);
		pblw.write(links);
		triads = linkBuilder.getTriadLinks();
		TriadLinkWriter tlw = new TriadLinkWriter(paras);
		tlw.init();
		tlw.write(triads);
		tlw.close();
	}
	
	private void buildEdges()
	{
		long start = System.currentTimeMillis();
		EdgeBundler eb = new EdgeBundler(paras);
//		edges = eb.pbLink2Edges(links, cntLens);
		edges = eb.links2edges(links, cntLens);
		links = null;
		triads = null;
		long end = System.currentTimeMillis();
		logger.info("Building Edges, Erase Time: " + (end - start) + " ms");
	}
	
	private void buildPaths()
	{
		PathBuilder pb = new PathBuilder(edges, paras, cntLens);
		paths = pb.buildPath();
		cntLens = null;
		edges = null;
	}
	
	private void writeNodePathInfo()
	{
		String pathFile = paras.getOutFolder() + System.getProperty("file.separator") + "nodePaths.info";
		DotGraphFileWriter.writeNodePaths(pathFile, paths);
	}
	
	private void writeScaffolds()
	{
		ScaffoldWriter sw = new ScaffoldWriter(paras, paths);
		sw.write();
	}

}



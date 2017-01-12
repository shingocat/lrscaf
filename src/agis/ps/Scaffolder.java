/*
*File: agis.ps.Scaffolder2.java
*User: mqin
*Email: mqin@ymail.com
*Date: 2016年7月17日
*/
package agis.ps;

import java.util.List;
import java.util.Map;
import java.util.Vector;

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
import agis.ps.link.CntFileEncapsulate;
import agis.ps.link.Edge;
import agis.ps.link.MRecord;
import agis.ps.link.PBLink;
import agis.ps.link.PBLinkM;
import agis.ps.link.TriadLink;
import agis.ps.path.NodePath;
import agis.ps.util.EdgeBundler;
import agis.ps.util.LinkBuilder;
import agis.ps.util.M5EdgeBundler;
import agis.ps.util.Parameter;
import agis.ps.util.PathBuilder;
import agis.ps.util.RepeatFinder;

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
	private Map<String, Integer> cntCovs;
	private Map<String, Integer> cntLens;
	private List<String> repeats;
	private List<PBLink> links;
	private List<TriadLink> triads;
	private List<Edge> edges;
	private List<NodePath> paths;
	private CntFileEncapsulate cntfile;
	
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
			this.buildEdges();
			this.buildPaths2();
			this.writeNodePathInfo();
//			this.writeScaffolds(paras, paths);	
			this.writeScaffolds();
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
		AlignmentFileReader reader = null;
		if(type.equalsIgnoreCase("m5"))
		{
			// using m5 format to build edges;
//			M5EdgeBundler bundler = new M5EdgeBundler(paras);
//			this.edges = bundler.building();
//			this.cntfile = bundler.getCntFile();
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
		records = reader.read();
		cntCovs = reader.getCntCoverages();
		cntLens = reader.getCntLengths();
		RepeatFinder rf = new RepeatFinder(paras);
		repeats = rf.findRepeats(cntCovs);
		LinkBuilder linkBuilder = new LinkBuilder(paras);
		if(links == null)
			links = new Vector<PBLink>();
		links = linkBuilder.mRecords2Links(records, repeats);
		PBLinkWriter pblw = new PBLinkWriter(paras);
		pblw.write(links);
		triads = linkBuilder.getTriadLinks();
		TriadLinkWriter tlw = new TriadLinkWriter(paras);
		tlw.init();
		tlw.write(triads);
		tlw.close();
		if(edges == null)
			edges = new Vector<Edge>();
		EdgeBundler eb = new EdgeBundler(paras);
		edges = eb.pbLink2Edges(links, cntLens);
	}
	
	private void buildPaths2()
	{
		PathBuilder pb = new PathBuilder(edges, paras, cntLens);
		paths = pb.buildPath();
	}
	
	private void buildPaths()
	{
		PathBuilder pathBuilder = new PathBuilder(edges, paras, cntfile);
		paths = pathBuilder.buildPath();
	}
	
//	private void findRepeats()
//	{
//		RepeatFinder rf = new RepeatFinder(paras);
//		repeats = rf.findRepeat();
//	}
	
	private void writeNodePathInfo()
	{
		String pathFile = paras.getOutFolder() + System.getProperty("file.separator") + "nodePaths.info";
		DotGraphFileWriter.writeNodePaths(pathFile, paths);
	}

	private void writeScaffolds(Parameter paras, List<NodePath> paths)
	{
		ScaffoldWriter sw = new ScaffoldWriter(paras, paths, cntfile);
//		sw.writeByContigFileEncapsulate();
	}
	
	private void writeScaffolds()
	{
		ScaffoldWriter sw = new ScaffoldWriter(paras, paths);
		sw.write();
	}

}



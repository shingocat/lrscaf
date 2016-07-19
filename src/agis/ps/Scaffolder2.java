/*
*File: agis.ps.Scaffolder2.java
*User: mqin
*Email: mqin@ymail.com
*Date: 2016年7月17日
*/
package agis.ps;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.List;
import java.util.Map;
import java.util.Vector;

import org.apache.lucene.analysis.Analyzer;
import org.apache.lucene.analysis.standard.StandardAnalyzer;
import org.apache.lucene.document.Document;
import org.apache.lucene.document.StringField;
import org.apache.lucene.document.TextField;
import org.apache.lucene.document.Field.Store;
import org.apache.lucene.document.StoredField;
import org.apache.lucene.index.CorruptIndexException;
import org.apache.lucene.index.IndexWriter;
import org.apache.lucene.index.IndexWriterConfig;
import org.apache.lucene.index.IndexWriterConfig.OpenMode;
import org.apache.lucene.store.Directory;
import org.apache.lucene.store.LockObtainFailedException;
import org.apache.lucene.store.SimpleFSDirectory;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import agis.ps.file.ContigIndexer;
import agis.ps.file.ContigReader;
import agis.ps.file.DotGraphFileWriter;
import agis.ps.file.OutputFolderBuilder;
import agis.ps.file.ScaffoldWriter;
import agis.ps.link.Edge;
import agis.ps.link.M5Record;
import agis.ps.link.MRecord;
import agis.ps.path.NodePath;
import agis.ps.seqs.Contig;
import agis.ps.util.M4EdgeBundler;
import agis.ps.util.M5EdgeBundler;
import agis.ps.util.Parameter;
import agis.ps.util.PathBuilder;
import agis.ps.util.Strand;

public class Scaffolder2 {
	private static Logger logger = LoggerFactory.getLogger(Scaffolder.class);
	private Parameter paras;
	private List<Edge> edges;
	private Map<String, Contig> contigs;
	
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
			// reading contigs;
//			if(!this.readContigs(paras))
//				return;
			// index contigs
			if(!this.indexCnt())
				return;
			// building edges 
			edges = buildEdges();
			if(edges == null || edges.size() == 0)
				return;
			this.writeEdgesInfo(edges, false);
			logger.info("Edges size: " + edges.size());
			PathBuilder pathBuilder = new PathBuilder(edges, paras);
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
	
	private boolean indexCnt()
	{
		ContigIndexer ci = new ContigIndexer(paras);
		return ci.indexing();
	}
	
	private boolean readContigs(Parameter paras)
	{
		ContigReader cr = new ContigReader(paras);
		contigs = cr.read(); // original methods with filtering parameter;
		if(contigs == null)
			return false;
		else
			return true;
	}
	
	// building edges
	private List<Edge> buildEdges()
	{
		long start = System.currentTimeMillis();
		List<Edge> edges = null;
		String type = paras.getType();
		if(type.equalsIgnoreCase("m5"))
		{
			M5EdgeBundler bundler = new M5EdgeBundler(paras);
			edges = bundler.building();
		} else if(type.equalsIgnoreCase("m4"))
		{
			M4EdgeBundler bundler = new M4EdgeBundler(paras);
			edges = bundler.building();
		} else if(type.equalsIgnoreCase("sam") || type.equalsIgnoreCase("bam"))
		{
			
		} else 
		{
			logger.info(this.getClass().getName() + "The aligned parameter do not set! only <m5>, <m4>, <sam> or <bam>");
			return null;
		}
		long end = System.currentTimeMillis();
		logger.info("Building edges, erase time " + (end - start) / 1000 + " s");
		return edges;
	}
	
	private void writeEdgesInfo(List<Edge> edges, boolean isPesudo) {
		// write the edges info into file.
		String edgeFile = "";
		if(isPesudo)
			edgeFile = paras.getOutFolder() + System.getProperty("file.separator") + "edges_after_pesudo.info";
		else
			edgeFile = paras.getOutFolder() + System.getProperty("file.separator") + "edges.info";
		DotGraphFileWriter.writeEdge(edgeFile, edges);
	}
	
	private void writeNodePathInfo(List<NodePath> paths)
	{
		String pathFile = paras.getOutFolder() + System.getProperty("file.separator") + "nodePaths.info";
		DotGraphFileWriter.writeNodePaths(pathFile, paths);
	}

	private void writeScaffolds(Parameter paras, List<NodePath> paths)
	{
		ScaffoldWriter sw = new ScaffoldWriter(paras, paths);
		sw.write3();
	}
	

}



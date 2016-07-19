/*
*File: agis.ps.Scaffolder.java
*User: mqin
*Email: mqin@ymail.com
*Date: 2015年12月28日
*/
package agis.ps;

import java.io.IOException;
import java.net.MalformedURLException;
import java.io.File;
import java.io.FileNotFoundException;
import java.util.List;
import java.util.Map;
import java.util.HashMap;
import java.util.Iterator;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import agis.ps.file.ContigReader;
import agis.ps.file.DotGraphFileWriter;
import agis.ps.file.M4Reader;
import agis.ps.file.M5Reader;
import agis.ps.file.ScaffoldWriter;
import agis.ps.link.Edge;
import agis.ps.link.MRecord;
import agis.ps.link.PBLinkM;
import agis.ps.path.NodePath;
import agis.ps.seqs.Contig;
import agis.ps.util.EdgeBundler;
import agis.ps.util.LinkBuilder;
import agis.ps.util.Parameter;
import agis.ps.util.PathBuilder;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamInputResource;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;

public class Scaffolder {
	final static Logger logger = LoggerFactory.getLogger(Scaffolder.class);

	private Parameter paras;
	private String aFilePath;
	private String type;
	private Map<String, Contig> contigs;
	private SamReader samReader;
	private List<PBLinkM> pbLinks;
	private List<Edge> edges;

	public Scaffolder(Parameter paras) {
		this.paras = paras;
		this.aFilePath = paras.getAlgFile();
		this.type = paras.getType();
	}

	public void scaffolding() {
		long start = System.currentTimeMillis();
		logger.info("Starting....");
		try {
			// if could not build the output folder;
			// return;
			if(!buildOutputPath(paras))
				return;
			// if could not read the contigs file;
			// return;
			if(!(readContigs(paras)))
				return;
			List<MRecord> mRecords = null;
			if (type.equalsIgnoreCase("m5")) {
				mRecords = readM5Aligned(paras);
				// if there are some problem to return m5 records
				if(mRecords == null)
					return;
			} else if(type.equalsIgnoreCase("m4")){
				mRecords = readM4Alingned(paras);
				if(mRecords == null)
					return;

			} else if(type.equalsIgnoreCase("sam") || type.equalsIgnoreCase("bam")) {
				readSAMAligned(aFilePath);
			} else
			{
				logger.debug(this.getClass().getName() + "The aligned parameter do not set! only <m5>, <m4>, <sam> or <bam>");
				logger.debug(this.getClass().getName() + "The aligned parameter do not set! only <m5>, <m4>, <sam> or <bam>");
				return;
			}
			// links building
			LinkBuilder linkBuilder = new LinkBuilder(mRecords, paras);
			pbLinks = linkBuilder.mRecord2Link(); // original strictly building pb links
			linkBuilder = null;
			// edges building
			EdgeBundler edgeBundler = new EdgeBundler(pbLinks, paras, contigs);
			edges = edgeBundler.pbLinkM5Bundling(); // original strictly building edges;
			if(edges == null)
				return;
			this.writeEdgesInfo(edges, false);
			logger.info("Edges size: " + edges.size());
			PathBuilder pathBuilder = new PathBuilder(edges, paras);
			List<NodePath> paths = pathBuilder.buildPath();
			this.writeNodePathInfo(paths);
			writeScaffolds(paras,paths,contigs);
		} catch (NullPointerException e) {
			logger.debug(this.getClass().getName() + "\t" + e.getMessage() + "\t" + e.getClass().getName());
			logger.error(this.getClass().getName() + "\t" + e.getMessage() + "\t" + e.getClass().getName());
		} catch (MalformedURLException e) {
			logger.debug(this.getClass().getName() + "\t" + e.getMessage() + "\t" + e.getClass().getName());
			logger.error(this.getClass().getName() + "\t" + e.getMessage() + "\t" + e.getClass().getName());
		} catch (FileNotFoundException e) {
			logger.debug(this.getClass().getName() + "\t" + e.getMessage() + "\t" + e.getClass().getName());
			logger.error(this.getClass().getName() + "\t" + e.getMessage() + "\t" + e.getClass().getName());
		} catch (IOException e) {
			logger.debug(this.getClass().getName() + "\t" + e.getMessage() + "\t" + e.getClass().getName());
			logger.error(this.getClass().getName() + "\t" + e.getMessage() + "\t" + e.getClass().getName());
		} catch (Exception e) {
			logger.debug(this.getClass().getName() + "\t" + e.getMessage() + "\t" + e.getClass().getName());
			logger.error(this.getClass().getName() + "\t" + e.getMessage() + "\t" + e.getClass().getName());
		}
		logger.info("Ending....");
		long end = System.currentTimeMillis();
		logger.info("Scaffolding erase time: " + Double.valueOf((end - start)/1000) + " s");
	}
	
	public void writeScaffolds(Parameter paras, List<NodePath> paths, Map<String, Contig> cnts)
	{
		String filePath = paras.getOutFolder() + System.getProperty("file.separator") + "scaffolds.fasta";
		ScaffoldWriter sw = new ScaffoldWriter(paras, paths, cnts, filePath);
		sw.write2();
	}
	
	public void writeNodePathInfo(List<NodePath> paths)
	{
		String pathFile = paras.getOutFolder() + System.getProperty("file.separator") + "nodePaths.info";
		DotGraphFileWriter.writeNodePaths(pathFile, paths);
	}

	public void writeEdgesInfo(List<Edge> edges, boolean isPesudo) {
		// write the edges info into file.
		String edgeFile = "";
		if(isPesudo)
			edgeFile = paras.getOutFolder() + System.getProperty("file.separator") + "edges_after_pesudo.info";
		else
			edgeFile = paras.getOutFolder() + System.getProperty("file.separator") + "edges.info";
		DotGraphFileWriter.writeEdge(edgeFile, edges);
	}
	
	public boolean buildOutputPath(Parameter paras)
	{
		String path = paras.getOutFolder();
		long start = System.currentTimeMillis();
		boolean isValid = false;
		if (path == null || path.length() == 0) {
			logger.error(this.getClass().getName() + "The output path was not setted!");
			logger.debug(this.getClass().getName() + "The output path was not setted!");
			return isValid;
		}
		try {
			File output = new File(path);
			if (output.exists()) {
				logger.info(this.getClass().getName() + "The output folder was exist!");
				logger.debug(this.getClass().getName() + "The output folder was exist!");
				isValid = true;
			} else {
				isValid = output.mkdirs();
				if(isValid)
				{
					logger.info(this.getClass().getName() + "\t" + "Build output folder successfully!");
					logger.debug(this.getClass().getName() + "\t" + "Build output folder successfully!");
				} else
				{
					logger.info(this.getClass().getName() + "\t" + "Build output folder failed!");
					logger.debug(this.getClass().getName() + "\t" + "Build output folder failed!");
				}
			}
		} catch (SecurityException e) {
			logger.error(this.getClass().getName() + e.getMessage() + "\t" + e.getClass().getName());
			logger.debug(this.getClass().getName() + e.getMessage() + "\t" + e.getClass().getName());
		} catch (Exception e) {
			logger.error(this.getClass().getName() + e.getMessage() + "\t" + e.getClass().getName());
			logger.debug(this.getClass().getName() + e.getMessage() + "\t" + e.getClass().getName());
		}
		long end = System.currentTimeMillis();
		logger.info("Building output folder, erase time: " + (end - start) + " ms");
		return isValid;
		
	}
	
	public boolean readContigs(Parameter paras)
	{
		ContigReader cr = new ContigReader(paras);
		contigs = cr.read(); // original methods with filtering parameter;
		if(contigs == null)
			return false;
		else
			return true;
	}

	public List<MRecord> readM5Aligned(Parameter paras)
	{
		M5Reader reader = new M5Reader(paras);
		return reader.read(); // original method to read m5 record
//		return reader.readWithoutFiltering();
	}
	
	public List<MRecord> readM4Alingned(Parameter paras)
	{
		M4Reader reader = new M4Reader(paras);
		return reader.read();
	}

	public void readSAMAligned(String aFilePath) throws NullPointerException, MalformedURLException, IOException {
		if (aFilePath == null || aFilePath.length() == 0) {
			logger.debug("The aligned file was null or not setted!");
			logger.info("The aligned file was null or not setted!");
			return;
		}
		SamReaderFactory factory = SamReaderFactory.makeDefault()
				.enable(SamReaderFactory.Option.INCLUDE_SOURCE_IN_RECORDS,
						SamReaderFactory.Option.VALIDATE_CRC_CHECKSUMS)
				.validationStringency(ValidationStringency.LENIENT);
		SamInputResource resource = SamInputResource.of(new File(aFilePath)); 
		samReader = factory.open(resource);
		Iterator<SAMRecord> it = samReader.iterator();
		HashMap<String, String> pSet = new HashMap<String, String>();
		while (it.hasNext()) {
			SAMRecord r = it.next();
			String c = r.getReferenceName() + "==" + r.getAlignmentStart() + ";";
			if (pSet.containsKey(r.getReadName())) {
				pSet.put(r.getReadName(), pSet.get(r.getReadName()) + c);
			} else {
				pSet.put(r.getReadName(), c);
			}

		}
		for (String s : pSet.keySet()) {
			logger.debug(s + "\t" + pSet.get(s));
		}
	}
}

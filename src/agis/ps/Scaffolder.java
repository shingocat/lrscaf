/*
*File: agis.ps.Scaffolder.java
*User: mqin
*Email: mqin@ymail.com
*Date: 2015年12月28日
*/
package agis.ps;

import java.io.IOException;
import java.net.MalformedURLException;
import java.net.URL;
import java.io.File;
import java.io.FileNotFoundException;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;

import org.biojava.nbio.core.sequence.io.FastaReaderHelper;
import org.biojava.nbio.core.sequence.DNASequence;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import agis.ps.link.Contig;
import agis.ps.link.PBLink;
import agis.ps.link.PBLinkM5;
import agis.ps.path.NodePath;
import agis.ps.util.Color;
import agis.ps.util.ContigReader;
import agis.ps.util.DotGraphFileWriter;
import agis.ps.util.EdgeBundler;
import agis.ps.util.LinkBuilder;
import agis.ps.util.M5Reader;
import agis.ps.util.Parameter;
import agis.ps.util.PathBuilder;
import agis.ps.util.ScaffoldWriter;
import htsjdk.samtools.AlignmentBlock;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamInputResource;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;

public class Scaffolder {
	final static Logger logger = LoggerFactory.getLogger(Scaffolder.class);

	// private HashMap<String, Object> paras;
	private Parameter paras;
	private String cFilePath;
	private String aFilePath;
	private String gFilePath;
	private String outFolder;
	private String type;
	private boolean m5Header;
//	private LinkedHashMap<String, DNASequence> contigs;
	private Map<String, Contig> contigs;
	private SamReader samReader;
	private List<SimplePath> simPaths;
	private List<PBLinkM5> pbLinks;
	private List<Edge> edges;

	public Scaffolder(Parameter paras) {
		this.paras = paras;
		this.cFilePath = paras.getCntFile();
		this.aFilePath = paras.getAlgFile();
		this.type = paras.getType();
		this.outFolder = paras.getOutFolder();
	}

	// public Scaffolder(String cFilePath, String aFilePath)
	// {
	// this.cFilePath = cFilePath;
	// this.aFilePath = aFilePath;
	// }
	//
	// public Scaffolder(HashMap<String, Object> paras)
	// {
	// this.paras = paras;
	// this.cFilePath = (String)paras.get("CONTIG");
	// this.aFilePath = (String)paras.get("ALIGNED");
	// this.type = (String)paras.get("TYPE");
	//// this.gFilePath = (String)paras.get("DOTGRAPH");
	// this.outFolder = (String) paras.get("OUTPUT");
	// }

	public void scaffolding() {
		logger.info("Starting....");
		try {
			// if could not build the output folder;
			// return;
			if(!buildOutputPath(paras.getOutFolder()))
				return;
			// if could not read the contigs file;
			// return;
			if(!(readContigs(cFilePath, paras.getMinContLen())))
				return;
			List<M5Record> m5Records = null;
			if (type.equalsIgnoreCase("m")) {
//				m5Records = readM5Aligned(aFilePath, paras.getMinPBLen(), paras.getMinContLen());
				m5Records = readM5Aligned(paras);
				// if there are some problem to return m5 records
				if(m5Records == null)
					return;
//				pbLinks = LinkBuilder.m5Record2Link(m5Records, paras);
				LinkBuilder linkBuilder = new LinkBuilder(m5Records, paras);
				pbLinks = linkBuilder.m5Record2Link();
				linkBuilder = null;
			} else if(type.equalsIgnoreCase("sam") || type.equalsIgnoreCase("bam")) {
				readSAMAligned(aFilePath);
			} else
			{
				logger.debug(this.getClass().getName() + "The aligned parameter do not set! only <m>, <sam> or <bam>");
				logger.debug(this.getClass().getName() + "The aligned parameter do not set! only <m>, <sam> or <bam>");
				return;
			}
			EdgeBundler edgeBundler = new EdgeBundler(pbLinks, paras);
			edges = edgeBundler.pbLinkM5Bundling();
			if(edges == null)
				return;
//			edges = EdgeBundler.pbLinkM5Bundling(pbLinks, paras);
			this.writeEdgesInfo(edges, false);
			logger.debug("Edges size: " + edges.size());
			logger.debug("Pesudo edging: ");
			edges = null;
			edges = edgeBundler.pesudoEdging();
			this.writeEdgesInfo(edges, true);
			logger.debug("Edges size after pesudo edging :" + edges.size());
			edgeBundler = null;
			// List<Path> paths = PathBuilder.buildHamiltonPath(edges);
//			List<Path> paths = PathBuilder.buildPath(edges, paras);
			PathBuilder pathBuilder = new PathBuilder(edges, paras);
			List<NodePath> paths = pathBuilder.buildPath();
			this.writeNodePathInfo(paths);
			writeScaffolds(paths, contigs);
			//List<Path> paths = pathBuilder.buildPath();
			//this.writePathsInfo(paths);
			
			// listContigs();
			// listAligns();
			// if(gFilePath != null)
			// {
			// DotGraphFileWriter.writePBLink(gFilePath, pbLinks);
			// DotGraphFileWriter.writeEdge(gFilePath, edges);
			// DotGraphFileWriter dGFW = new DotGraphFileWriter(gFilePath,
			// simPaths);
			// dGFW.write();
			// }
		} catch (NullPointerException e) {
			logger.debug(this.getClass().getName() + e.getMessage() + "\t" + e.getClass().getName());
			logger.error(this.getClass().getName() + e.getMessage() + "\t" + e.getClass().getName());
		} catch (MalformedURLException e) {
			logger.debug(this.getClass().getName() + e.getMessage() + "\t" + e.getClass().getName());
			logger.error(this.getClass().getName() + e.getMessage() + "\t" + e.getClass().getName());
		} catch (FileNotFoundException e) {
			logger.debug(this.getClass().getName() + e.getMessage() + "\t" + e.getClass().getName());
			logger.error(this.getClass().getName() + e.getMessage() + "\t" + e.getClass().getName());
		} catch (IOException e) {
			logger.debug(this.getClass().getName() + e.getMessage() + "\t" + e.getClass().getName());
			logger.error(this.getClass().getName() + e.getMessage() + "\t" + e.getClass().getName());
		} catch (Exception e) {
			logger.debug(this.getClass().getName() + e.getMessage() + "\t" + e.getClass().getName());
			logger.error(this.getClass().getName() + e.getMessage() + "\t" + e.getClass().getName());
		}

		logger.info("Ending....");
	}
	
	public void writeScaffolds(List<NodePath> paths, Map<String, Contig> cnts)
	{
		String filePath = paras.getOutFolder() + System.getProperty("file.separator") + "scaffolds.fasta";
		ScaffoldWriter sw = new ScaffoldWriter(paths, cnts, filePath);
		sw.write();
	}
	
	public void writePathsInfo(List<Path> paths) {
		// write the paths info into file;
		String pathFile = paras.getOutFolder() + System.getProperty("file.separator") + "paths.info";
		DotGraphFileWriter.writePaths(pathFile, paths);
	}
	
	public void writeNodePathInfo(List<NodePath> paths)
	{
		String pathFile = paras.getOutFolder() + System.getProperty("file.separator") + "NodePaths.info";
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

	public boolean buildOutputPath(String path) {
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
				if (output.mkdirs()) {
					logger.info(this.getClass().getName() + "The output folder was created!");
					logger.debug(this.getClass().getName() + "The output folder was created!");
					isValid = true;
				}
			}
		} catch (SecurityException e) {
			logger.error(this.getClass().getName() + e.getMessage() + "\t" + e.getClass().getName());
			logger.debug(this.getClass().getName() + e.getMessage() + "\t" + e.getClass().getName());
		} catch (Exception e) {
			logger.error(this.getClass().getName() + e.getMessage() + "\t" + e.getClass().getName());
			logger.debug(this.getClass().getName() + e.getMessage() + "\t" + e.getClass().getName());
		}
		return isValid;
	}

	public boolean readContigs(String cFilePath, int cntLen) throws IOException {
//		if (cFilePath == null || cFilePath.length() == 0) {
//			logger.error("The contig file was null or not setted!");
//			logger.debug("The contig file was null or not setted!");
//			logger.info("The contig file was null or not setted!");
//			return;
//		}
//		File cFile = new File(cFilePath);
//		setContigs(FastaReaderHelper.readFastaDNASequence(cFile));
		ContigReader cr = new ContigReader(cFilePath, cntLen);
//		setContigs(cr.read());
		contigs = cr.read();
		if(contigs == null)
			return false;
		else
			return true;
	}
	
	public List<M5Record> readM5Aligned(Parameter paras)
	{
		M5Reader reader = new M5Reader(paras);
		return reader.read();
	}

	public List<M5Record> readM5Aligned(String aFilePath, int minPBLen, int minCNTLen) {
		if (aFilePath.isEmpty()) {
			logger.error(this.getClass().getName() + "The aligned file was null or not setted!");
			logger.debug(this.getClass().getName() + "The aligned file was null or not setted!");
			logger.info(this.getClass().getName() + "The aligned file was null or not setted!");
			return null;
		}
		M5Reader reader = new M5Reader(aFilePath, minPBLen, minCNTLen);
		List<M5Record> m5Records = reader.read();
		// HashMap<String, String> pSet = new HashMap<String, String>();
		//
		// for(M5Record m5 : m5Records)
		// {
		// String c = m5.gettName() + "==" + m5.getqStart() + ";";
		// if(pSet.containsKey(m5.getqName()))
		// {
		// pSet.put(m5.getqName(), pSet.get(m5.getqName()) + c);
		// } else
		// {
		// pSet.put(m5.getqName(), c);
		// }
		// }
		return m5Records;
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
		SamInputResource resource = SamInputResource.of(new File(aFilePath)); // .index(new
																				// URL("http://broadinstitute.org/my.bam.bai"));
		samReader = factory.open(resource);
		if (simPaths != null) {
			simPaths.clear();
		}
		simPaths = new ArrayList<SimplePath>();
		Iterator<SAMRecord> it = samReader.iterator();
		// HashMap<String, HashMap<String, Integer>> pSet = new HashMap<String,
		// HashMap<String, Integer>>();
		HashMap<String, String> pSet = new HashMap<String, String>();
		while (it.hasNext()) {
			SAMRecord r = it.next();
			// logger.debug("Read name:" + r.getReadName() + "\tReference name:"
			// + r.getReferenceName() +
			// "\tAligned start:" + r.getAlignmentStart() + "\tAligned end:" +
			// r.getAlignmentEnd() +
			// "\tFLAG:" + r.getFlags());
			String c = r.getReferenceName() + "==" + r.getAlignmentStart() + ";";
			// HashMap<String, Integer> cSet = new HashMap<String,Integer>();
			// cSet.put(r.getReferenceName(),
			// r.getReadPositionAtReferencePosition(r.getAlignmentStart()));
			if (pSet.containsKey(r.getReadName())) {
				pSet.put(r.getReadName(), pSet.get(r.getReadName()) + c);
			} else {
				pSet.put(r.getReadName(), c);
			}

		}
		for (String s : pSet.keySet()) {
			// logger.debug(s);
			logger.debug(s + "\t" + pSet.get(s));
		}
	}

	public void m5RecordToSimplePath(Map<String, String> pSet) {
		if (simPaths != null) {
			simPaths.clear();
		} else {
			simPaths = new ArrayList<SimplePath>();
		}
		for (String s : pSet.keySet()) {
			String[] arrs = pSet.get(s).split(";");
			if (arrs.length >= 2) {
				for (int i = 0; i < arrs.length - 1; i++) {
					SimplePath sp = new SimplePath();
					String[] ar1 = arrs[i].split("==");
					String[] ar2 = arrs[i + 1].split("==");
					if (Integer.valueOf(ar1[1]) <= Integer.valueOf(ar2[1])) {
						sp.setStart(ar1[0]);
						sp.setEnd(ar2[0]);
						sp.setLabel(ar1[1] + ":" + ar2[1]);
						sp.setColor(Color.BLUE);
					} else {
						sp.setStart(ar2[0]);
						sp.setEnd(ar1[0]);
						sp.setLabel(ar2[1] + ":" + ar1[1]);
						sp.setColor(Color.BLUE);
					}
					simPaths.add(sp);
				}
			}
		}
	}

	public void listContigs() {
		if (contigs == null || contigs.isEmpty()) {
			return;
		}
		for (String id : contigs.keySet()) {
			System.out.println("Id is:\t" + id);
		}
	}

	public void listAligns() {
		if (samReader == null)
			return;
		for (SAMRecord sr : samReader) {
			System.out.println("Read Name:\t" + sr.getReadName());
			System.out.println("Aligned Name:\t" + sr.getReferenceName());
		}
	}

	public String getcFilePath() {
		return cFilePath;
	}

	public void setcFilePath(String cFilePath) {
		this.cFilePath = cFilePath;
	}

	public String getaFilePath() {
		return aFilePath;
	}

	public void setaFilePath(String aFilePath) {
		this.aFilePath = aFilePath;
	}

//	public LinkedHashMap<String, DNASequence> getContigs() {
//		return contigs;
//	}
//
//	public void setContigs(LinkedHashMap<String, DNASequence> contigs) {
//		this.contigs = contigs;
//	}

}

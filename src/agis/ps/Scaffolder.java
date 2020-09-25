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
import agis.ps.file.SequenceReader;
import agis.ps.file.M4Reader;
import agis.ps.file.M5Reader;
import agis.ps.file.MMReader;
import agis.ps.file.NodePathWriter;
//import agis.ps.file.MisassembliesWriter;
import agis.ps.file.OutputFolderBuilder;
import agis.ps.file.RepeatWriter;
import agis.ps.file.LinkWriter;
import agis.ps.file.SamReader;
import agis.ps.file.ScaffoldWriter;
import agis.ps.file.TriadLinkWriter;
import agis.ps.link.Edge;
import agis.ps.link.ILink;
import agis.ps.link.MRecord;
import agis.ps.path.NodePath;
import agis.ps.seqs.Sequence;
import agis.ps.util.EdgeBundler;
import agis.ps.util.LinkBuilder;
//import agis.ps.util.MisassemblyChecker;
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
	private List<List<MRecord>> listRecords;
	private Map<String, Integer> seqCovs;
	private List<String> repeats;
//	private List<String> misassemblies;
	private List<ILink> links;
	private List<ILink> triads;
	private List<Edge> edges;
	private List<NodePath> paths;
//	private Map<String, Contig> cnts;
//	private Map<String, Scaffold> scafs;
	private Map<String, Sequence> seqs;
	// private List<Contig> unusedCnts; // contigs were not in building
	// scaffolds;

	public Scaffolder(Parameter paras) {
		this.paras = paras;
		// this.unusedCnts = new ArrayList<Contig>();
	}

	public void scaffolding() {
		try {
			if (!buildOutputFolder())
				return;
			this.readDraftFile();
			this.buildLinks();
			this.buildEdges();
			this.buildPaths();
			this.writeNodePathInfo();
			this.writeScaffolds();
//			this.writeSmallAndRepeatCnts();
			this.writeRepeatCnts();
//			this.writeMisassemblies();
		} catch (Exception e) {
			logger.debug("Error: ", e);
			logger.error(e.getMessage());
		}
	}

	private boolean buildOutputFolder() {
		OutputFolderBuilder ofb = new OutputFolderBuilder(paras);
		return ofb.building();
	}

	// read the contigs file;
	private void readDraftFile() {
		SequenceReader cr = new SequenceReader(paras);
//		cnts = cr.read();
		this.seqs = cr.read();
	}

	private void buildLinks() {
		String type = paras.getType();
		AlignmentFileReader reader = null;
		if (type.equalsIgnoreCase("m5")) {
			reader = new M5Reader(paras);
		} else if (type.equalsIgnoreCase("m4")) {
			reader = new M4Reader(paras);
		} else if (type.equalsIgnoreCase("sam")) {
			reader = new SamReader(paras);
		} else if (type.equalsIgnoreCase("mm")) {
			reader = new MMReader(paras);
		} else {
			logger.info("The aligned parameter should be not set! only <m5>, <m4>, <sam> or <bam>");
			return;
		}
//		listRecords = reader.read(cnts);
		listRecords = reader.read(this.seqs);
		seqCovs = reader.getCntCoverages();
		RepeatFinder rf = new RepeatFinder(paras);
		// original code 2020/9/16
//		repeats = rf.findRepeats(seqCovs, cnts);
		repeats = rf.findRepeats(seqCovs);
		// change sequence to repeat
		for(String id : repeats) {
			if(this.seqs.containsKey(id))
				this.seqs.get(id).setIsRepeat(true);
		}
//		misassemblies = MisassemblyChecker.findMisassemblies(paras, cnts);
		LinkBuilder linkBuilder = new LinkBuilder(paras, this.seqs);
		links = linkBuilder.mRecords2Links(listRecords, repeats);
//		links = linkBuilder.mRecords2Links(listRecords, repeats, misassemblies);
		LinkWriter pblw = new LinkWriter(paras);
		pblw.write(links);
		triads = linkBuilder.getTriadLinks();
		this.writeTriadLinks();
	}

	private void buildEdges() {
		EdgeBundler eb = new EdgeBundler(paras);
		// edges = eb.pbLink2Edges(links, cntLens);
		// edges = eb.links2edges(links, cntLens);
		edges = eb.links2edges(links);
		links = null;
		triads = null;
	}

	private void buildPaths() {
		// PathBuilder pb = new PathBuilder(edges, paras, cntLens);
		PathBuilder pb = new PathBuilder(edges, paras, this.seqs);
		paths = pb.buildPath();
		// cntLens = null;
		edges = null;
	}

	private void writeNodePathInfo() {
//		String pathFile = paras.getOutFolder() + System.getProperty("file.separator") + "nodePaths.info";
//		DotGraphFileWriter.writeNodePaths(pathFile, paths);
		NodePathWriter writer = new NodePathWriter(this.paras);
		writer.write(this.paths);
	}

	private void writeScaffolds() {
		ScaffoldWriter sw = new ScaffoldWriter(paras, paths, this.seqs);
		sw.write();
	}

	private void writeTriadLinks() {
		TriadLinkWriter tlw = new TriadLinkWriter(paras);
		tlw.init(false);
		tlw.write(triads);
		tlw.close();
	}
	
	private void writeRepeatCnts() {
		RepeatWriter rw = new RepeatWriter(paras, this.seqs);
		rw.write();
	}
	
//	private void writeMisassemblies()
//	{
//		MisassembliesWriter mw = new MisassembliesWriter(paras, cnts);
//		mw.write();
//	}
}

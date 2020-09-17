/*
*File: agis.ps.util.ContigReader.java
*User: mqin
*Email: mqin@ymail.com
*Date: 2016年2月28日
*/
package agis.ps.file;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
//import java.util.Vector;
//import java.util.regex.Matcher;
//import java.util.regex.Pattern;
import java.util.regex.PatternSyntaxException;

//import org.biojava.nbio.core.exceptions.CompoundNotFoundException;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

//import agis.ps.link.M5Record;
import agis.ps.seqs.Contig;
import agis.ps.seqs.Scaffold;
import agis.ps.seqs.Sequence;
import agis.ps.util.Parameter;
import agis.ps.util.SequenceUtils;

public class SequenceReader {

	public static Logger logger = LoggerFactory.getLogger(SequenceReader.class);
//	private Map<String, Contig> cnts = new HashMap<String, Contig>();
//	private Map<String, Scaffold> scafs = new HashMap<String, Scaffold>(); 
	private Map<String, Sequence> seqs = new HashMap<String, Sequence>();
	private String filePath;
	private List<Integer> lens;
	private Parameter paras;
	
	public SequenceReader(Parameter paras) {
		this.paras = paras;
		this.filePath = paras.getCntFile();
	}
	
	public Map<String, Sequence> read(){
		long start = System.currentTimeMillis();
		if (this.seqs == null)
			this.seqs = new HashMap<String, Sequence>();
		this.seqs.clear();
		if(lens == null)
			lens = new ArrayList<Integer>();
		lens.clear();
		BufferedReader br = null;
		try {
			File seqFile = new File(filePath);
			if (!seqFile.exists()) {
				logger.error("\t" + "The draft assembly file " + filePath + " do not exist!");
				return null;
			}

			br = Files.newBufferedReader(seqFile.toPath());
			String line = null;
			String id = null;
			StringBuilder sb = new StringBuilder();
			String temp;
			while (true) {
				line = br.readLine();
				if(line == null) {
					int length = sb.length();
					if(id != null && length >= 0) {
						id = SequenceUtils.formatId(id);
						Scaffold scaf = SequenceUtils.seq2scaf(id, sb.toString());
						if(length < paras.getMinContLen())
							scaf.setIsUsed(false);
						else
							scaf.setIsUsed(true);
						scaf.setIsRepeat(false);
						this.seqs.put(id, scaf);
						scaf = null;
						lens.add(length);
					}
					sb = null;
					temp = null;
					id = null;
					break;
				}
				line = line.trim();
				line = line.replaceAll(System.getProperty("line.separator"), "");
				if (line.startsWith(">")) {
					temp = id;
					id = line;
					int length = sb.length();
					if(temp != null && length >= 0) {
						temp = SequenceUtils.formatId(temp);
						Scaffold scaf = SequenceUtils.seq2scaf(temp, sb.toString());
						if(length < paras.getMinContLen())
							scaf.setIsUsed(false);
						else
							scaf.setIsUsed(true);
						scaf.setIsRepeat(false);
						this.seqs.put(temp, scaf);
						sb = null;
						scaf = null;
						temp = null;
						sb = new StringBuilder();
						lens.add(length);
					}
				} else {
					sb.append(line);
				}
			}

		} catch (ArrayIndexOutOfBoundsException e) {
			logger.debug("Error: ", e);
			logger.error(e.getMessage() + "\t" + e.getClass().getName());
		} catch (FileNotFoundException e) {
			logger.debug("Error: ", e);
			logger.error(e.getMessage() + "\t" + e.getClass().getName());
		} catch (IOException e) {
			logger.debug("Error: ", e);
			logger.error(e.getMessage() + "\t" + e.getClass().getName());
		} catch (PatternSyntaxException e) {
			logger.debug("Error: ", e);
			logger.error(e.getMessage() + "\t" + e.getClass().getName());
		} catch (Exception e) {
			logger.debug("Error: ", e);
			logger.error(e.getMessage() + "\t" + e.getClass().getName());
		} finally {
			try {
				if (br != null)
					br.close();
			} catch (IOException e) {
				logger.debug("Error: ", e);
				logger.error(e.getMessage() + "\t" + e.getClass().getName());
			}
		}
		long end = System.currentTimeMillis();
		logger.info("Reading contigs, elapsed times: " + (end - start) + " ms");
		this.writeCntsSummary();
		return this.seqs;
	}
	
//	public Map<String, Scaffold> readScaffolds(){
//		long start = System.currentTimeMillis();
//		if (this.scafs == null)
//			this.scafs = new HashMap<String, Scaffold>();
//		this.scafs.clear();
//		if(lens == null)
//			lens = new ArrayList<Integer>();
//		lens.clear();
//		BufferedReader br = null;
//		try {
//			File cntFile = new File(filePath);
//			if (!cntFile.exists()) {
//				logger.error("\t" + "The contig file " + filePath + " do not exist!");
//				return null;
//			}
//
//			br = Files.newBufferedReader(cntFile.toPath());
//			String line = null;
//			String id = null;
//			StringBuilder sb = new StringBuilder();
//			String temp;
//			while (true) {
//				line = br.readLine();
//				if(line == null) {
//					int length = sb.length();
//					if(id != null && length >= 0) {
//						id = SequenceUtils.formatId(id);
//						Scaffold scaf = SequenceUtils.seq2scaf(id, sb.toString());
//						if(length < paras.getMinContLen())
//							scaf.setIsUsed(false);
//						else
//							scaf.setIsUsed(true);
//						scaf.setIsRepeat(false);
//						this.scafs.put(id, scaf);
//						scaf = null;
//						lens.add(length);
//					}
//					sb = null;
//					temp = null;
//					id = null;
//					break;
//				}
//				line = line.trim();
//				line = line.replaceAll(System.getProperty("line.separator"), "");
//				if (line.startsWith(">")) {
//					temp = id;
//					id = line;
//					int length = sb.length();
//					if(temp != null && length >= 0) {
//						id = SequenceUtils.formatId(id);
//						Scaffold scaf = SequenceUtils.seq2scaf(id, sb.toString());
//						if(length < paras.getMinContLen())
//							scaf.setIsUsed(false);
//						else
//							scaf.setIsUsed(true);
//						scaf.setIsRepeat(false);
//						this.scafs.put(id, scaf);
//						sb = null;
//						scaf = null;
//						temp = null;
//						sb = new StringBuilder();
//						lens.add(length);
//					}
//				} else {
//					sb.append(line);
//				}
//			}
//
//		} catch (ArrayIndexOutOfBoundsException e) {
//			logger.debug("Error: ", e);
//			logger.error(e.getMessage() + "\t" + e.getClass().getName());
//		} catch (FileNotFoundException e) {
//			logger.debug("Error: ", e);
//			logger.error(e.getMessage() + "\t" + e.getClass().getName());
//		} catch (IOException e) {
//			logger.debug("Error: ", e);
//			logger.error(e.getMessage() + "\t" + e.getClass().getName());
//		} catch (PatternSyntaxException e) {
//			logger.debug("Error: ", e);
//			logger.error(e.getMessage() + "\t" + e.getClass().getName());
//		} catch (Exception e) {
//			logger.debug("Error: ", e);
//			logger.error(e.getMessage() + "\t" + e.getClass().getName());
//		} finally {
//			try {
//				if (br != null)
//					br.close();
//			} catch (IOException e) {
//				logger.debug("Error: ", e);
//				logger.error(e.getMessage() + "\t" + e.getClass().getName());
//			}
//		}
//		long end = System.currentTimeMillis();
//		logger.info("Reading contigs, erase times: " + (end - start) + " ms");
//		this.writeCntsSummary();
//		return this.scafs;
//	}
//	
//	// original method for reading contigs with filtering parameters
//	public Map<String, Contig> readContigs() {
//		long start = System.currentTimeMillis();
//		if (cnts == null)
//			cnts = new HashMap<String, Contig>();
//		cnts.clear();
//		if(lens == null)
//			lens = new ArrayList<Integer>();
//		lens.clear();
//		FileReader fr = null;
//		BufferedReader br = null;
//		try {
//			File cntFile = new File(filePath);
//			if (!cntFile.exists()) {
//				logger.error("\t" + "The contig file " + filePath + " do not exist!");
//				return null;
//			}
//
//			fr = new FileReader(cntFile);
//			br = new BufferedReader(fr);
//			String line = null;
//			String id = null;
//			StringBuilder sb = new StringBuilder();
//			String temp;
//			while (true) {
//				line = br.readLine();
//				if(line == null) {
//					int length = sb.length();
//					if(id != null && length >= 0) {
//						Contig cnt = new Contig();
//						cnt.setSeqs(sb.toString());
//						cnt.setLength(length);
//						if(length < paras.getMinContLen())
//							cnt.setIsUsed(false);
//						else
//							cnt.setIsUsed(true);
//						cnt.setIsRepeat(false);
//						id = id.replaceAll("^>", "");
//						id = id.split("\\s")[0];
//						id = id.trim();
//						cnt.setID(id);
//						cnts.put(id, cnt);
//						cnt = null;
//						lens.add(length);
//					}
//					sb = null;
//					temp = null;
//					id = null;
//					break;
//				}
//				line = line.trim();
//				line = line.replaceAll(System.getProperty("line.separator"), "");
//				if (line.startsWith(">")) {
//					temp = id;
//					id = line;
//					int length = sb.length();
//					if(temp != null && length >= 0)
//					{
//						temp = temp.replaceFirst("^>", "");
//						temp = temp.split("\\s")[0];
//						temp = temp.trim();
//						Contig cnt = new Contig();
//						cnt.setSeqs(sb.toString());
//						cnt.setLength(length);
//						cnt.setID(temp);
//						if(length < paras.getMinContLen())
//							cnt.setIsUsed(false);
//						else
//							cnt.setIsUsed(true);
//						cnt.setIsRepeat(false);
//						cnts.put(temp, cnt);
//						sb = null;
//						cnt = null;
//						temp = null;
//						sb = new StringBuilder();
//						lens.add(length);
//					}
//				} else {
//					sb.append(line);
//				}
//			}
//
//		} catch (ArrayIndexOutOfBoundsException e) {
//			logger.debug("Error: ", e);
//			logger.error(e.getMessage() + "\t" + e.getClass().getName());
//		} catch (FileNotFoundException e) {
//			logger.debug("Error: ", e);
//			logger.error(e.getMessage() + "\t" + e.getClass().getName());
//		} catch (IOException e) {
//			logger.debug("Error: ", e);
//			logger.error(e.getMessage() + "\t" + e.getClass().getName());
//		} catch (PatternSyntaxException e) {
//			logger.debug("Error: ", e);
//			logger.error(e.getMessage() + "\t" + e.getClass().getName());
//		} catch (Exception e) {
//			logger.debug("Error: ", e);
//			logger.error(e.getMessage() + "\t" + e.getClass().getName());
//		} finally {
//			try {
//				if (br != null)
//					br.close();
//			} catch (IOException e) {
//				logger.debug("Error: ", e);
//				logger.error(e.getMessage() + "\t" + e.getClass().getName());
//			}
//		}
//		long end = System.currentTimeMillis();
//		logger.info("Reading contigs, erase times: " + (end - start) + " ms");
//		this.writeCntsSummary();
//		return cnts;
//	}
	
	private void writeCntsSummary(){
		N50Writer writer = new N50Writer(paras, "draft", lens);
		writer.write();
	}
}

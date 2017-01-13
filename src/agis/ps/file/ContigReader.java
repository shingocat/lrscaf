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
import java.util.HashMap;
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
import agis.ps.util.Parameter;

public class ContigReader {

	public static Logger logger = LoggerFactory.getLogger(ContigReader.class);
	public Map<String, Contig> cnts = new HashMap<String, Contig>();
	private String filePath;
	
	public ContigReader(Parameter paras)
	{
		this.filePath = paras.getCntFile();
	}
	
	
	// original method for reading contigs with filtering parameters
	public Map<String, Contig> read() {
		long start = System.currentTimeMillis();
		if (cnts == null)
			cnts = new HashMap<String, Contig>();
		cnts.clear();
		FileReader fr = null;
		BufferedReader br = null;
		try {
			File cntFile = new File(filePath);
			if (!cntFile.exists()) {
				logger.debug(this.getClass().getName() + "\t" + "The contig file " + filePath + " do not exist!");
				logger.error(this.getClass().getName() + "\t" + "The contig file " + filePath + " do not exist!");
				return null;
			}

			fr = new FileReader(cntFile);
			br = new BufferedReader(fr);
			String line = null;
			String id = null;
			StringBuffer sb = new StringBuffer();
			String temp;
			while (true) {
				line = br.readLine();
				if(line == null)
				{
					int length = sb.length();
					if(id != null && length >= 0)
					{
						Contig cnt = new Contig();
						cnt.setSeqs(sb.toString());
						cnt.setLength(length);
						id = id.replaceAll("^>", "");
						id = id.split("\\s")[0];
						id = id.trim();
						cnt.setID(id);
						cnts.put(id, cnt);
						cnt = null;
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
					if(temp != null && length >= 0)
					{
						temp = temp.replaceFirst("^>", "");
						temp = temp.split("\\s")[0];
						temp = temp.trim();
						Contig cnt = new Contig();
						cnt.setSeqs(sb.toString());
						cnt.setLength(length);
						cnt.setID(temp);
						cnts.put(temp, cnt);
						sb = null;
						cnt = null;
						temp = null;
						sb = new StringBuffer();
					}
				} else {
					sb.append(line);
				}
			}

		} catch (ArrayIndexOutOfBoundsException e) {
			logger.error(this.getClass().getName() + "\t" + e.getMessage() + "\t" + e.getClass().getName());
		} catch (FileNotFoundException e) {
			logger.error(this.getClass().getName() + "\t" + e.getMessage() + "\t" + e.getClass().getName());
		} catch (IOException e) {
			logger.error(this.getClass().getName() + "\t" + e.getMessage() + "\t" + e.getClass().getName());
		} catch (PatternSyntaxException e) {
			logger.error(this.getClass().getName() + "\t" + e.getMessage() + "\t" + e.getClass().getName());
		} catch (Exception e) {
			logger.error(this.getClass().getName() + "\t" + e.getMessage() + "\t" + e.getClass().getName());
		} finally {
			try {
				if (br != null)
					br.close();
			} catch (IOException e) {
				logger.error(this.getClass().getName() + "\t" + e.getMessage() + "\t" + e.getClass().getName());
			}
		}
		long end = System.currentTimeMillis();
		logger.info("Reading contigs, erase times: " + (end - start) + " ms");
		return cnts;
	}
}

/*
*File: agis.ps.util.ContigReader.java
*User: mqin
*Email: mqin@ymail.com
*Date: 2016年2月28日
*/
package agis.ps.util;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.HashMap;
import java.util.Map;
import java.util.Vector;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.regex.PatternSyntaxException;

import org.biojava.nbio.core.exceptions.CompoundNotFoundException;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import agis.ps.M5Record;
import agis.ps.link.Contig;

public class ContigReader {

	public static Logger logger = LoggerFactory.getLogger(ContigReader.class);
	public Map<String, Contig> cnts = new HashMap<String, Contig>();
	private String filePath;

	public ContigReader(String filePath) {
		this.filePath = filePath;
	}

	public Map<String, Contig> read() {
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
					if(id != null && sb.length() != 0)
					{
						Contig cnt = new Contig(sb.toString());
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
					if(temp != null && sb.length() != 0)
					{
						temp = temp.replaceFirst("^>", "");
						temp = temp.split("\\s")[0];
						temp = temp.trim();
						if(temp.equalsIgnoreCase("1867"))
							System.out.println(sb.toString());
						Contig cnt = new Contig(sb.toString());
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

			// while((line = br.readLine()) != null)
			// {
			// line.replaceAll("\\s", "");
			// if(line.startsWith(">"))
			// {
			// if(sb.length() == 0)
			// {
			// id = line;
			// continue;
			// } else
			// {
			// id = id.replaceFirst("^>","");
			// id = id.split("\\s")[0];
			// Contig cnt = new Contig(sb.toString());
			// cnt.setID(id);
			// cnts.put(id, cnt);
			// sb = new StringBuffer();
			// id = line;
			// }
			// } else {
			// sb.append(line);
			// }
			// }

		} catch (ArrayIndexOutOfBoundsException e) {
			logger.debug(this.getClass().getName() + "\t" + e.getMessage() + "\t" + e.getClass().getName());
			logger.error(this.getClass().getName() + "\t" + e.getMessage() + "\t" + e.getClass().getName());
		} catch (FileNotFoundException e) {
			logger.debug(this.getClass().getName() + "\t" + e.getMessage() + "\t" + e.getClass().getName());
			logger.error(this.getClass().getName() + "\t" + e.getMessage() + "\t" + e.getClass().getName());
		} catch (IOException e) {
			logger.debug(this.getClass().getName() + "\t" + e.getMessage() + "\t" + e.getClass().getName());
			logger.error(this.getClass().getName() + "\t" + e.getMessage() + "\t" + e.getClass().getName());
		} catch (CompoundNotFoundException e) {
			logger.debug(this.getClass().getName() + "\t" + e.getMessage() + "\t" + e.getClass().getName());
			logger.error(this.getClass().getName() + "\t" + e.getMessage() + "\t" + e.getClass().getName());
		} catch (PatternSyntaxException e) {
			logger.debug(this.getClass().getName() + "\t" + e.getMessage() + "\t" + e.getClass().getName());
			logger.error(this.getClass().getName() + "\t" + e.getMessage() + "\t" + e.getClass().getName());
		} catch (Exception e) {
			logger.debug(this.getClass().getName() + "\t" + e.getMessage() + "\t" + e.getClass().getName());
			logger.error(this.getClass().getName() + "\t" + e.getMessage() + "\t" + e.getClass().getName());
		} finally {
			try {
				if (br != null)
					br.close();
			} catch (IOException e) {
				logger.debug(this.getClass().getName() + "\t" + e.getMessage() + "\t" + e.getClass().getName());
				logger.error(this.getClass().getName() + "\t" + e.getMessage() + "\t" + e.getClass().getName());
			}
		}
		return cnts;
	}
}

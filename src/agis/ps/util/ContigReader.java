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

import org.biojava.nbio.core.exceptions.CompoundNotFoundException;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import agis.ps.M5Record;
import agis.ps.link.Contig;

public class ContigReader {
	
	public static Logger logger = LoggerFactory.getLogger(ContigReader.class);
	public Map<String, Contig> cnts = new HashMap<String, Contig>();
	private String filePath;
	
	public ContigReader(String filePath)
	{
		this.filePath = filePath;
	}
	
	public Map<String, Contig> read()
	{
		if(cnts == null)
			cnts = new HashMap<String, Contig>();
		cnts.clear();
		FileReader fr = null;
		BufferedReader br = null;
		try
		{
			File cntFile = new File(filePath);
			if(!cntFile.exists())
			{
				logger.debug("The contig file, " + filePath + ", do not exist!");
				logger.error("The contig file, " + filePath + ", do not exist!");
				return null;
			}
			
			fr = new FileReader(cntFile);
			br = new BufferedReader(fr);
			String line = null;
			String id = null;
			StringBuilder sb = new StringBuilder();
			while((line = br.readLine()) != null)
			{
				line.replaceAll("\\s", "");
				if(line.startsWith(">"))
				{
					if(sb.length() == 0)
					{	
						id = line;
						continue;
					} else
					{
						id = id.replaceFirst("^>","");
						id = id.split("\\s")[0];
						Contig cnt = new Contig(sb.toString());
						cnt.setID(id);
						cnts.put(id, cnt);
						sb = new StringBuilder();
						id = line;
					}
				} else {
					sb.append(line);
				}
			}
			br.close();
			fr.close();
		} catch(ArrayIndexOutOfBoundsException e)
		{
			logger.debug(e.getMessage());
			logger.error(e.getMessage());
		} catch(FileNotFoundException e)
		{
			logger.debug(e.getMessage());
			logger.error(e.getMessage());
		} catch(IOException e)
		{
			logger.debug(e.getMessage());
			logger.error(e.getMessage());
		} catch (CompoundNotFoundException e) {
			logger.debug(e.getMessage());
			logger.error(e.getMessage());
		} finally
		{
			if(br != null)
			{
				try {
					br.close();
				} catch (IOException e) {
					logger.debug(e.getMessage());
					logger.error(e.getMessage());
				}
			} 
			if(fr != null)
			{
				try {
					fr.close();
				} catch (IOException e) {
					logger.debug(e.getMessage());
					logger.error(e.getMessage());
				}
			}
		}
		return cnts;
	}
}



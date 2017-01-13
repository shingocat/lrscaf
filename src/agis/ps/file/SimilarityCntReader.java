/*
*File: agis.ps.file.SimilarityCntReader.java
*User: mqin
*Email: mqin@ymail.com
*Date: 2016年12月12日
*/
package agis.ps.file;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.LinkedList;
import java.util.List;
import java.util.Vector;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

//import agis.ps.link.TriadLink;
//import agis.ps.seqs.Contig;
import agis.ps.util.Parameter;

public class SimilarityCntReader {
	public static Logger logger = LoggerFactory.getLogger(SimilarityCntReader.class);
	public List<LinkedList<String>> simcnts = null;
	public Parameter paras = null;
	
	public SimilarityCntReader (Parameter paras)
	{
		this.paras = paras;
	}
	
	public List<LinkedList<String>> read()
	{
		simcnts = new Vector<LinkedList<String>>(15);
		String outFolder = paras.getOutFolder();
		String fileName = outFolder + System.getProperty("file.separator") + "simcnts.info";
		File file = null;
		FileReader fr = null;
		BufferedReader br = null;
		try
		{
			file = new File(fileName);
			if(!file.exists())
				throw new IOException("The triad link info file do not exist!");
			fr = new FileReader(file);
			br = new BufferedReader(fr);
			String line = "";
			while((line = br.readLine()) != null)
			{
				String [] strs = line.split("\t");
				LinkedList<String> sims = new LinkedList<String>();
				for(String s : strs)
				{
					sims.add(s);
				}
				simcnts.add(sims);
			}
			
		} catch(IOException e)
		{
			logger.error(this.getClass().getName() + "\t" + e.getMessage() + "\t" + e.getClass().getName());
		} finally
		{
			try
			{
				if(br != null)
					br.close();
			}  catch(IOException e)
			{
				logger.error(this.getClass().getName() + "\t" + e.getMessage() + "\t" + e.getClass().getName());
			} 
		}
		return simcnts;
	}
}



/*
*File: agis.ps.file.TriadLinkReader.java
*User: mqin
*Email: mqin@ymail.com
*Date: 2016年5月13日
*/
package agis.ps.file;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.List;
import java.util.Vector;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import agis.ps.link.TriadLink;
import agis.ps.seqs.Contig;
import agis.ps.util.Parameter;

public class TriadLinkReader {
	public static Logger logger = LoggerFactory.getLogger(TriadLinkReader.class);
	public List<TriadLink> triads = null;
	public Parameter paras = null;
	
	public TriadLinkReader (Parameter paras)
	{
		this.paras = paras;
	}
	
	public List<TriadLink> read()
	{
		triads = new Vector<TriadLink>(100);
		String outFolder = paras.getOutFolder();
		String fileName = outFolder + System.getProperty("file.separator") + "triadlinks.info";
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
			int index = 0;
			while((line = br.readLine()) != null)
			{
//				if(index == 386)
//					logger.debug(index + " " + line);
				String [] strs = line.split(",");
				Contig pre = new Contig();
				pre.setID(strs[0]);
				Contig lst = new Contig();
				lst.setID(strs[2]);
				Contig mid = null;
				TriadLink tl = null;
				if(!strs[1].equalsIgnoreCase("-"))
				{
					mid = new Contig();
					mid.setID(strs[1]);
//					tl = new TriadLink(pre, mid, lst);
//					tl.setSupLinks(Integer.valueOf(strs[3]));
				}
				tl = new TriadLink(pre, mid, lst);
				tl.setSupLinks(Integer.valueOf(strs[3]));
				triads.add(tl);
				index++;
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
		return triads;
	}
}



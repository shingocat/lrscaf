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

import agis.ps.link.Contig;
import agis.ps.link.TriadLink;
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
		String fileName = outFolder + System.getProperty("file.separator") + "triad_link.info";
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
				String [] strs = line.split(",");
				String [] pres = strs[0].split("\\(length=");
				Contig pre = new Contig();
				pre.setID(pres[0]);
				pre.setLength(Integer.valueOf(pres[1].substring(0, pres[1].length() - 1)));
				String [] mids = strs[1].split("\\(length=");
				Contig mid = new Contig();
				mid.setID(mids[0]);
				mid.setLength(Integer.valueOf(mids[1].substring(0, mids[1].length() - 1)));
				String [] lsts = strs[2].split("\\(length=");
				Contig lst = new Contig();
				lst.setID(lsts[0]);
				lst.setLength(Integer.valueOf(lsts[1].substring(0, lsts[1].length() - 1)));
				String [] sls = strs[3].split("=");
				int sl = Integer.valueOf(sls[1]);
				TriadLink tl = new TriadLink(pre, mid, lst);
				tl.setSupLinks(sl);
				triads.add(tl);
			}
			
		} catch(IOException e)
		{
			logger.debug(this.getClass().getName() + "\t" + e.getMessage() + "\t" + e.getClass().getName());
			logger.error(this.getClass().getName() + "\t" + e.getMessage() + "\t" + e.getClass().getName());
		} finally
		{
			try
			{
				if(br != null)
					br.close();
			}  catch(IOException e)
			{
				logger.debug(this.getClass().getName() + "\t" + e.getMessage() + "\t" + e.getClass().getName());
				logger.error(this.getClass().getName() + "\t" + e.getMessage() + "\t" + e.getClass().getName());
			} 
		}
		return triads;
	}
}



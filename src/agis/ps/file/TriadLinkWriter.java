/*
*File: agis.ps.file.TriadLinkWriter.java
*User: mqin
*Email: mqin@ymail.com
*Date: 2016年5月13日
*/
package agis.ps.file;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Arrays;
import java.util.List;
import java.util.Map;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import agis.ps.link.TriadLink;
import agis.ps.seqs.Contig;
import agis.ps.util.Parameter;

public class TriadLinkWriter {
	public static Logger logger = LoggerFactory.getLogger(TriadLinkWriter.class);
	private Parameter paras;
	private List<TriadLink> triads;
	
	public TriadLinkWriter(Parameter paras, List<TriadLink> triads)
	{
		this.paras = paras;
		this.triads = triads;
	}

	public void write()
	{
		String outFolder = paras.getOutFolder();
		String fileName = outFolder + System.getProperty("file.separator") + "triad_link.info";
		File file = null; 
		FileWriter fw = null;
		BufferedWriter bw = null;
		try
		{
			file = new File(fileName);
			if (file.exists()) {
				logger.debug("The output file of scaffold is exist! It will not be overwrited!");
				logger.info("The output file of scaffold is exist! It will not be overwrited!");
				return;
			}
			if(!file.createNewFile())
			{
				logger.debug("ScaffoldWriter: The output file of scaffolds could not create!");
				logger.info("ScaffoldWriter: The output file of scaffolds could not create!");
				return;
			}
			fw = new FileWriter(file);
			bw = new BufferedWriter(fw);
			for(TriadLink tl : triads)
			{
				Contig pre = tl.getPrevious();
				Contig mid = tl.getMiddle();
				Contig lst = tl.getLast();
				bw.write(pre.getID() + "(length=" + pre.getLength() + ")," + 
						mid.getID() + "(length=" + mid.getLength() + ")," + 
						lst.getID() + "(length=" + lst.getLength() + ")," +
						"supLinks=" + tl.getSupLinks());
				bw.newLine();
			}
			bw.flush();
		} catch(IOException e)
		{
			logger.debug(this.getClass().getName() + "\t" + e.getMessage() + "\t" + e.getClass().getName());
			logger.error(this.getClass().getName() + "\t" + e.getMessage() + "\t" + e.getClass().getName());
		} finally
		{
			try{
				if(bw != null)
					bw.close();
			} catch(IOException e)
			{
				logger.debug(this.getClass().getName() + "\t" + e.getMessage() + "\t" + e.getClass().getName());
				logger.error(this.getClass().getName() + "\t" + e.getMessage() + "\t" + e.getClass().getName());
			}
		}
	}
}



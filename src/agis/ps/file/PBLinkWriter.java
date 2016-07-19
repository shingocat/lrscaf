/*
*File: agis.ps.file.PBLinkWriter.java
*User: mqin
*Email: mqin@ymail.com
*Date: 2016年7月19日
*/
package agis.ps.file;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Iterator;
import java.util.List;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import agis.ps.link.MRecord;
import agis.ps.link.PBLinkM;
import agis.ps.util.Parameter;

public class PBLinkWriter {
	
	private static Logger logger = LoggerFactory.getLogger(PBLinkWriter.class);
	private Parameter paras = null;
	private File file = null;
	private FileWriter fw = null;
	private BufferedWriter bw = null;
	
	public PBLinkWriter(Parameter paras)
	{
		this.paras = paras;
	}
	
	public void init()
	{
		try{
			file = new File(paras.getOutFolder() + System.getProperty("file.separator") + "links.info");
			fw = new FileWriter(file, true);
			bw = new BufferedWriter(fw);
		} catch(IOException e)
		{
			logger.error(this.getClass().getName() + "\t" + e.getMessage());
		} catch(Exception e)
		{
			logger.error(this.getClass().getName() + "\t" + e.getMessage());
		}
	}
	
	public void write(List<PBLinkM> links)
	{
		if(links != null && links.size() > 0)
		{
			Iterator<PBLinkM> it = links.iterator();
			while(it.hasNext())
			{
				PBLinkM link = it.next();
				MRecord origin = link.getOrigin();
				MRecord terminus = link.getTerminus();
				String id = link.getId();
				int dist = link.getDistance();
				String line = origin.gettName() + "\t" + origin.gettStrand() + "\t" + terminus.gettName() + "\t" +
				terminus.gettStrand() + "\t" + dist + "\t" + id + "\t" + origin.getqStart() + "\t" + origin.getqEnd()
				+ "\t" + terminus.getqStart() + "\t" + terminus.getqEnd();
				try {
					bw.write(line);
					bw.newLine();
				} catch (IOException e) {
					logger.error(this.getClass().getName() + "\t" + e.getMessage());
				}
			}
		}
	}
	
	public void close()
	{
		try{
			if(bw != null)
				bw.close();
		} catch(Exception e)
		{
			logger.error(this.getClass().getName() + "\t" + e.getMessage());
		}
	}

}



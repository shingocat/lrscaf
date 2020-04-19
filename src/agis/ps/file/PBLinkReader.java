/*
*File: agis.ps.file.PBLinkReader.java
*User: mqin
*Email: mqin@ymail.com
*Date: 2016年7月19日
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

import agis.ps.link.PBLink;
import agis.ps.util.Parameter;

public class PBLinkReader {
	private static Logger logger = LoggerFactory.getLogger(PBLinkReader.class);
	private Parameter paras = null;
	private File file = null;
	private FileReader fr = null;
	private BufferedReader br = null;
	
	public PBLinkReader(Parameter paras)
	{
		this.paras = paras;
	}
	
	public List<PBLink> read()
	{
		List<PBLink> links = new Vector<PBLink>(1000);
		try
		{
			String path = paras.getOutFolder();
			file = new File(path + System.getProperty("file.separator") + "links.info");
			fr = new FileReader(file);
			br = new BufferedReader(fr);
			String line = "";
			while((line = br.readLine()) != null)
			{
				line = line.trim();
				String arrs[] = line.split("\t");
				PBLink link = new PBLink();
				link.setOrigin(arrs[0]);
				link.setOStrand(arrs[1]);
				link.setTerminus(arrs[2]);
				link.setTStrand(arrs[3]);
				link.setDist(Integer.valueOf(arrs[4]));
				link.setPbId(arrs[5]);
				link.setOPStart(Integer.valueOf(arrs[6]));
				link.setOPEnd(Integer.valueOf(arrs[7]));
				link.setTPStart(Integer.valueOf(arrs[8]));
				link.setTPEnd(Integer.valueOf(arrs[9]));
				links.add(link);
			}
			br.close();
		} catch(IOException e)
		{
			logger.debug("Error: ", e);
			logger.error(this.getClass().getName() + "\t" + e.getMessage());
		} catch(Exception e)
		{
			logger.debug("Error: ", e);
			logger.error(this.getClass().getName() + "\t" + e.getMessage());
		} finally{
			try{
				
			} catch(Exception e)
			{
				logger.debug("Error: ", e);
				logger.error(this.getClass().getName() + "\t" + e.getMessage());
			}
		}		
		return links;
	}

}



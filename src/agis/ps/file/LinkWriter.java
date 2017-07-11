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

import agis.ps.link.Link;
//import agis.ps.link.MRecord;
import agis.ps.link.PBLink;
//import agis.ps.link.PBLinkM;
import agis.ps.util.Parameter;

public class LinkWriter {
	
	private static Logger logger = LoggerFactory.getLogger(LinkWriter.class);
	private Parameter paras = null;
	private File file = null;
	private FileWriter fw = null;
	private BufferedWriter bw = null;
	
	public LinkWriter(Parameter paras)
	{
		this.paras = paras;
	}
	
//	public void init()
//	{
//		try{
//			file = new File(paras.getOutFolder() + System.getProperty("file.separator") + "links.info");
//			fw = new FileWriter(file, true);
//			bw = new BufferedWriter(fw);
//		} catch(IOException e)
//		{
//			logger.error(this.getClass().getName() + "\t" + e.getMessage());
//		} catch(Exception e)
//		{
//			logger.error(this.getClass().getName() + "\t" + e.getMessage());
//		}
//	}
//	
//	public void write(List<PBLinkM> links)
//	{
//		if(links != null && links.size() > 0)
//		{
//			Iterator<PBLinkM> it = links.iterator();
//			while(it.hasNext())
//			{
//				PBLinkM link = it.next();
//				MRecord origin = link.getOrigin();
//				MRecord terminus = link.getTerminus();
//				String id = link.getId();
//				int dist = link.getDistance();
//				String line = origin.gettName() + "\t" + origin.gettStrand() + "\t" + terminus.gettName() + "\t" +
//				terminus.gettStrand() + "\t" + dist + "\t" + id + "\t" + origin.getqStart() + "\t" + origin.getqEnd()
//				+ "\t" + terminus.getqStart() + "\t" + terminus.getqEnd();
//				try {
//					bw.write(line);
//					bw.newLine();
//				} catch (IOException e) {
//					logger.error(this.getClass().getName() + "\t" + e.getMessage());
//				}
//			}
//		}
//	}
	
//	public void write(List<PBLink> links)
//	{
//		try{
//			file = new File(paras.getOutFolder() + System.getProperty("file.separator") + "links.info");
//			fw = new FileWriter(file, true);
//			bw = new BufferedWriter(fw);
//			if(links != null && links.size() > 0)
//			{
//				Iterator<PBLink> it = links.iterator();
//				while(it.hasNext())
//				{
//					PBLink link = it.next();
//					String line = link.getOrigin() + "\t" + link.getOStrand() + "\t" + link.getTerminus() + "\t" +
//					link.getTStrand() + "\t" + link.getDist() + "\t" + link.getPbId() + "\t" + link.getOPStart() + "\t" 
//					+ link.getOPEnd() + "\t" + link.getTPStart() + "\t" + link.getTPEnd();
//					bw.write(line);
//					bw.newLine();
//				}
//				bw.flush();
//				bw.close();
//			}
//		} catch(IOException e)
//		{
//			logger.error(this.getClass().getName() + "\t" + e.getMessage());
//		} catch(Exception e)
//		{
//			logger.error(this.getClass().getName() + "\t" + e.getMessage());
//		} finally
//		{
//			try{
//				if(bw != null)
//					bw.close();
//			} catch(Exception e)
//			{
//				logger.error(this.getClass().getName() + "\t" + e.getMessage());
//			}
//		}
//		
//	}
	
	public void write(List<Link> links)
	{
		try{
			file = new File(paras.getOutFolder() + System.getProperty("file.separator") + "links.info");
			fw = new FileWriter(file, true);
			bw = new BufferedWriter(fw);
			if(links != null && links.size() > 0)
			{
				Iterator<Link> it = links.iterator();
				while(it.hasNext())
				{
					Link link = it.next();
					String line = link.getOriginal().getID() + "\t" + link.getoStrand() + "\t" + link.getTerminus().getID() + "\t" +
					link.gettStrand() + "\t" + link.getDist() + "\t" + link.getLrId() + "\t" + link.getoStart() + "\t" 
					+ link.getoEnd() + "\t" + link.gettStart() + "\t" + link.gettEnd();
					bw.write(line);
					bw.newLine();
				}
				bw.flush();
				bw.close();
			}
		} catch(IOException e)
		{
			logger.error(this.getClass().getName() + "\t" + e.getMessage());
		} catch(Exception e)
		{
			logger.error(this.getClass().getName() + "\t" + e.getMessage());
		} finally
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
	
//	public void close()
//	{
//		try{
//			if(bw != null)
//				bw.close();
//		} catch(Exception e)
//		{
//			logger.error(this.getClass().getName() + "\t" + e.getMessage());
//		}
//	}

}



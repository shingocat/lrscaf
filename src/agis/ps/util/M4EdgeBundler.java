/*
*File: agis.ps.util.M4EdgeBundler.java
*User: mqin
*Email: mqin@ymail.com
*Date: 2016年7月17日
*/
package agis.ps.util;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.List;
import java.util.Vector;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import agis.ps.link.Edge;
import agis.ps.link.MRecord;

public class M4EdgeBundler {
	private static Logger logger = LoggerFactory.getLogger(M4EdgeBundler.class);
	private Parameter paras = null;
	
	public M4EdgeBundler(Parameter  paras)
	{
		this.paras = paras;
	}
	
	public List<Edge> building()
	{
		File file = null;
		FileReader fr = null;
		BufferedReader br = null;
		String alnFile = paras.getAlgFile();
		List<Edge> edges = new Vector<Edge>(200);
		int minPBLen = paras.getMinPBLen();
		int minCNTLen = paras.getMinContLen();
		double identity = paras.getIdentity();
		try{
			file = new File(alnFile);
			if(!file.exists())
			{
				logger.error(this.getClass().getName() + "The m5 file do not exist!");
				return null;
			}
			fr = new FileReader(file);
			br = new BufferedReader(fr);
			String line = null;
			String [] arrs = null;
			List<MRecord> records = new Vector<MRecord>(10);
			String id = null; // pacbio long read id
			while(true)
			{
				line = br.readLine();
				if(line != null){
					line.trim();
					line = line.replaceAll(System.getProperty("line.separator"), "");
					arrs = line.split("\\s+");
					if(arrs[0].equalsIgnoreCase("qName") && arrs[1].equalsIgnoreCase("qLength"))
						continue;
					if(id == null)
					{
						id = arrs[0];
						// if less than minimum pacbio length
						if(Integer.valueOf(arrs[1]) < minPBLen)
							continue;
						// if less tan minimum contig length
						if(Integer.valueOf(arrs[6]) < minCNTLen)
							continue;
						// if the identity less than specified value;
						double sum = Double.valueOf(arrs[11]) + Double.valueOf(arrs[12]) + Double.valueOf(arrs[13]) + Double.valueOf(arrs[14]);
						double value = Double.valueOf(arrs[11]) / sum;
						if(value < identity)
							continue;
					} else 
					{
						if(id.equals(arrs[0]))
						{
							
						} else
						{
							
						}
					}
				} else{
					// if line == null;
				}
			}
			
		} catch(IOException e)
		{
			logger.error(this.getClass().getName() + "\t" + e.getMessage());
		} catch(Exception e)
		{
			logger.error(this.getClass().getName() + "\t" + e.getMessage());
		}
		return edges;
	}
}



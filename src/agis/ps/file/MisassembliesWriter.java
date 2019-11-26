/** 
** Usage: TODO
** Author: mqin
** Email: mqin@outlook.com
** Date: 2017年8月9日
*/
package agis.ps.file;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Map;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import agis.ps.seqs.Contig;
import agis.ps.util.Parameter;

public class MisassembliesWriter {
	private static Logger logger = LoggerFactory.getLogger(MisassembliesWriter.class);
	private Parameter paras;
	private Map<String, Contig> cnts;
	
	public MisassembliesWriter(Parameter paras, Map<String, Contig> cnts)
	{
		this.paras = paras;
		this.cnts = cnts;
	}
	
	public void write()
	{
		File misassembly = null;
		FileWriter fwMisassembly = null;
		BufferedWriter bwMisassembly = null;
		try{
			misassembly = new File(paras.getOutFolder() + System.getProperty("file.separator") + "misassembly.contigs");
			if(misassembly.exists()) {
				logger.info("The output file " + misassembly.getCanonicalPath() + " existed. It will overwrite.");
			} else {
				if(!misassembly.createNewFile()) {
					logger.error("The output file " + misassembly.getCanonicalPath() + " could not create.");
					return;
				}
			}
			fwMisassembly = new FileWriter(misassembly, false);
			bwMisassembly = new BufferedWriter(fwMisassembly);
			for(Map.Entry<String, Contig> entry : cnts.entrySet())
			{
				Contig c = entry.getValue();
				String id = entry.getKey();
				if(c.getIsMisassembly())
				{
					bwMisassembly.write(">" + id);
					bwMisassembly.newLine();
					bwMisassembly.write(c.getForwardSeqs());
					bwMisassembly.newLine();
					continue;
				}
			}
			bwMisassembly.flush();
			bwMisassembly.close();
		} catch(IOException e)
		{
			logger.error(e.getMessage());
		} catch(Exception e)
		{
			logger.error(e.getMessage());
		} finally{
			try{
				if(bwMisassembly != null)
					bwMisassembly.close();
			} catch(IOException e)
			{
				logger.error(e.getMessage());
			}
			
		}
	}
}

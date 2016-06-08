/*
*File: agis.ps.util.ContigCoverageWriter.java
*User: mqin
*Email: mqin@ymail.com
*Date: 2016年5月12日
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

import agis.ps.seqs.Contig;
import agis.ps.util.Parameter;

public class ContigCoverageWriter {
	public static Logger logger = LoggerFactory.getLogger(ContigCoverageWriter.class);
	private Parameter paras;
	private Map<String, List<String>> args;
	
	public ContigCoverageWriter(Parameter paras, Map<String, List<String>> args)
	{
		this.paras = paras;
		this.args = args;
	}
	
	public void write()
	{
		String outFolder = paras.getOutFolder();
		String fileName = outFolder + System.getProperty("file.separator") + "contig_coverage.info";
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
			for(String s : args.keySet())
			{
				bw.write(s + "\t" + args.get(s).size() + "\t" + Arrays.toString(args.get(s).toArray()));
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



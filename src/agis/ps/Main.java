/*
*File: agis.ps.Scaffolder.java
*User: mqin
*Email: mqin@ymail.com
*Date: 2015Äê12ÔÂ28ÈÕ
*The main function interface for scaffolding
*/
package agis.ps;

import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.log4j.PropertyConfigurator;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.MissingArgumentException;
import org.apache.commons.cli.ParseException;
import org.apache.commons.cli.CommandLine;

import java.util.HashMap;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

public class Main {
	
	final static Logger logger = LoggerFactory.getLogger(Main.class);
	
	public static void main(String [] args)
	{
		//System.setProperty("log4j.configuration", "log4j.properties");
		try{
			Options ops = new Options();
			ops.addOption("c", "contig", true, "The file of contig in fasta format!");
			ops.addOption("a", "align", true, "The file in sam|bam format after aligned by using blasr!");
			ops.addOption("t", "type", true, "The aligned file type, s for sam or bam, m for m5!");
			ops.addOption("g", "dotgraph", true, "The output file path for dot graph file!");
			ops.addOption("h", "help", false, "Print help info!");
			
			CommandLineParser parser = new DefaultParser();
			CommandLine cl = parser.parse(ops, args);
			
			if(cl.hasOption("h") || !cl.hasOption("a") || !cl.hasOption("c"))
			{
				HelpFormatter f = new HelpFormatter();
				f.printHelp("Usage: java agisps.jar -c <Contig File> -a <Aligned File>", ops);
				System.exit(0);
			} else
			{
				String cFile = cl.getOptionValue("c");
				String aFile = cl.getOptionValue("a");
				String type = cl.getOptionValue("type");
				
				HashMap<String, Object> paras = new HashMap<String, Object>();
				paras.put("CONTIG", cFile);
				paras.put("ALIGNED", aFile);
				
				if(type != null)
				{
					paras.put("TYPE", type);
				} else
				{
					paras.put("TYPE", "s");
				}
				
				String dotPath = cl.getOptionValue("g");
				paras.put("DOTGRAPH", dotPath);
				
				Scaffolder scaffolder = new Scaffolder(paras);
				scaffolder.scaffolding();
			}
		} catch(MissingArgumentException e)
		{
			logger.error(e.getMessage());
			logger.debug(e.getMessage());
			logger.info(e.getMessage());
		}catch(ParseException e)
		{
			logger.error(e.getMessage());
			logger.debug(e.getMessage());
			logger.info(e.getMessage());
		} 
		
	}
}



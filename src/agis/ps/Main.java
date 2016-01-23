/*
*File: agis.ps.Scaffolder.java
*User: mqin
*Email: mqin@ymail.com
*Date: 2015年12月28日
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

import java.io.PrintWriter;
import java.util.HashMap;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import agis.ps.util.Parameter;
import agis.ps.util.XMLParser;

public class Main {

	final static Logger logger = LoggerFactory.getLogger(Main.class);

	public static void main(String[] args) {
		// System.setProperty("log4j.configuration", "log4j.properties");
		try {
			Options ops = new Options();
			ops.addOption("x", "xml", true, "The XML file for specified parameters on scaffolding!");
			ops.addOption("c", "contig", true, "The file of contig in fasta format!");
			ops.addOption("a", "align", true, "The file in sam|bam format after aligned by using blasr!");
			ops.addOption("t", "type", true, "The aligned file type, s for sam or bam, m for m5!");
//			ops.addOption("mh", "m5header", false, "The indicator of m5 file whether has header!");
//			ops.addOption("g", "dotgraph", true, "The output file path for dot graph file!");
			ops.addOption("o", "output", true, "The output file path!");
			ops.addOption("h", "help", false, "Print help info!");

			CommandLineParser parser = new DefaultParser();
			CommandLine cl = parser.parse(ops, args);
			HelpFormatter f = new HelpFormatter();
			Parameter paras = null;
			
			if (cl.hasOption("h")) {
				f.printHelp("Usage: java -jar agisps.jar -c <Contig File> -a <Aligned File> or"
						+ "\t java -jar agisps.jar -x <Parameters File>\n", ops);
				System.exit(0);
			} else {
				// for parse the xml configure, all the other parameters set by command line will dismissed
				if (cl.hasOption("x")) {
					logger.info("parse the xml configure, all the other parameters set by command line will dismissed");
					String xmlFile = cl.getOptionValue("x");
					XMLParser xmlParser = new XMLParser();
					paras = xmlParser.parseXML(xmlFile);
					Scaffolder scaffolder = new Scaffolder(paras);
					scaffolder.scaffolding();
				} else {
					// checking the aligned setting valid or not
					if(!cl.hasOption("a"))
					{
						logger.debug("The argument '-a' setting were not valid!");
						logger.error("The argument '-a' setting were not valid!");
						f.printHelp("Usage: java -jar agisps.jar -c <Contig File> -a <Aligned File> or"
								+ "\t java -jar agisps.jar -x <Parameters File>\n", ops);
						System.exit(0);
					}
					// checking the contig setting or not
					if(!cl.hasOption("c"))
					{
						logger.debug("The argument '-c' setting were not valid!");
						logger.error("The argument '-c' setting were not valid!");
						f.printHelp("Usage: java -jar agisps.jar -c <Contig File> -a <Aligned File> or"
								+ "\t java -jar agisps.jar -x <Parameters File>\n", ops);
						System.exit(0);
					}
					
					String cFile = cl.getOptionValue("c");
					String aFile = cl.getOptionValue("a");
					String type = cl.getOptionValue("t");

//					HashMap<String, Object> paras = new HashMap<String, Object>();
//					paras.put("CONTIG", cFile);
//					paras.put("ALIGNED", aFile);
					paras.setCntFile(cFile);
					paras.setAlgFile(aFile);

					if (type != null) {
						if (!(type.matches("(?i:[ms])"))) {
							logger.info("The type parameter is not valid, setted " + type + ", only accepted s or m!");
							logger.debug("The type parameter is not valid, setted " + type + ", only accepted s or m!");
							System.exit(0);
						}
						paras.setType(type);
//						paras.put("TYPE", type);
//						if (type.equalsIgnoreCase("m")) {
//							if (cl.hasOption("mh"))
//								paras.put("M5HEADER", true);
//							else
//								paras.put("M5HEADER", false);
//						}
					} else {
						type = "s";
						paras.setType(type);
//						paras.put("TYPE", "s");
					}

//					String dotPath = cl.getOptionValue("g");
//					paras.put("DOTGRAPH", dotPath);
					String outPath = cl.getOptionValue("o");
//					paras.put("OUTPUT", outPath);
					paras.setOutFolder(outPath);

					Scaffolder scaffolder = new Scaffolder(paras);
					scaffolder.scaffolding();
				}
			}
		} catch (MissingArgumentException e) {
			logger.error(e.getMessage());
			logger.debug(e.getMessage());
			logger.info(e.getMessage());
		} catch (ParseException e) {
			logger.error(e.getMessage());
			logger.debug(e.getMessage());
			logger.info(e.getMessage());
		} catch (Exception e) {
			logger.error(e.getMessage());
			logger.debug(e.getMessage());
			logger.info(e.getMessage());
		}
	}
}

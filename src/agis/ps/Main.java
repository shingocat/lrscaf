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
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.MissingArgumentException;
import org.apache.commons.cli.CommandLine;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import agis.ps.file.XMLParser;
import agis.ps.util.Parameter;

public class Main {

	final static Logger logger = LoggerFactory.getLogger(Main.class);

	public static void main(String[] args) {
		// System.setProperty("log4j.configuration", "log4j.properties");
		long start = System.currentTimeMillis();
		logger.info("Launching...");
		try {
			Options opts = Main.initOpts();
			CommandLineParser parser = new DefaultParser();
			CommandLine cl = parser.parse(opts, args);
			HelpFormatter f = new HelpFormatter();
			Parameter paras = null;
			
			if (cl.hasOption("h")) {
				f.printHelp("java -jar agisps.jar -c <Contig File> -a <Aligned File> -o <Output Folder> or"
						+ " java -jar agisps.jar -x <Parameters File>\n", opts);
				System.exit(0);
			} else {
				// for parse the xml configure, all the other parameters set by command line will dismissed
				if (cl.hasOption("x")) {
					logger.info("Parsing the xml configure, all the other parameters set by command line will dismissed");
					String xmlFile = cl.getOptionValue("x");
					XMLParser xmlParser = new XMLParser();
					paras = xmlParser.parseXML(xmlFile);
				} else {				
					paras = Main.parsering(cl, f, opts);
				}
				// The body for executing scaffolding process
				Scaffolder scaffolder = new Scaffolder(paras);
				scaffolder.scaffolding();
			}
		} catch (MissingArgumentException e) {
			logger.error(Main.class.getName() + "\t" + e.getMessage());
		} catch (ParseException e) {
			logger.error(Main.class.getName() + "\t"+ e.getMessage());
		} catch (Exception e) {
			logger.error(Main.class.getName() + "\t" + e.getMessage());
		}
		logger.info("Ending...");
		long end = System.currentTimeMillis();
		logger.info("Scaffolding erase time: " + (end - start)/1000 + " s.");
	}
	
	private static Options initOpts()
	{
		Options opts = new Options();
		// xml file
		opts.addOption("x", "xml", true, "The parameter XML file! All command-line parameters would be omitted if set.");
		// contig
		opts.addOption("c", "contig", true, "The file of pre-assembled contigs in fasta format.");
		// pacbio reads
		opts.addOption("p", "pacbio", true, "The PacBio Long Reads file for gap filling if requried.");
		// m5 format
		opts.addOption("m5", "m5", true, "The aligned file in -m 5 format of blasr.");
		// m4 format
		opts.addOption("m4", "m4", true, "The aligned file in -m 4 format of blasr.");
		// minimap format
		opts.addOption("mm", "m4", true, "The aligned file in minimap output format.");
		// sam format
		opts.addOption("sam", "sam", true, "The aligned file in sam format.");
		// output folder
		opts.addOption("o", "output", true, "The scaffolding output folder.");
		// help 
		opts.addOption("h", "help", false, "The help infomation!");
		// minimum contig length
		opts.addOption("micl", "miniCntLen", true, "The minimum contig's length for scaffolding! Default:<3000>.");
		// minimum pacbio long read length
		opts.addOption("mipl", "miniPBLen", true, "The minimum PacBio long read's length for scaffolding! Default:<5000>.");
		// identity
		opts.addOption("i", "identity", true, "The identity threshold for blasr alignment! Default: <0.8>.");
		// minimum overlap length
		opts.addOption("mioll", "miniOLLen", true, "The minimum overlap length threshold for blasr alignment! Default: <2400>.");
		// minimum overlap ratio
		opts.addOption("miolr", "miniOLRatio", true, "The minimum overlap ratio threshold for blasr alignment! Default: <0.8>.");
		// maximum overhang length;
		opts.addOption("maohl", "maOHLen", true, "The maximum overhang length threshold for blasr alignment! Default: <300>.");
		// maximum overhang ratio;
		opts.addOption("maohr", "maOHRatio", true, "The maximum overhang ratio threshold for blasr alignment! Default: <0.1>.");
		// maximum ending length
		opts.addOption("mael", "maELen", true, "The maximum ending length of PacBio's Long Read! Default: <300>.");
		// maximum ending ratio
		opts.addOption("maer", "maERatio", true, "The maximum ending ratio of PacBio's Long Read! Default: <0.1>.");
		// minimum support links;
		opts.addOption("misl", "miSLN", true, "The minimum support links number! Default: <1>.");
		// use overlap link
		opts.addOption("doll", "doll", false, "The indicator for using overlap relationship of contig! Default: <f>, discard overlap relationship.");
		// deleting error prone edges ratio;
		opts.addOption("r", "ratio", true, "The ratio for deleting error prone edges! Default: <0.2>.");
		// masking repeat
		opts.addOption("mr", "mr", false, "The indicator for masking repeats! Default: <f>.");
		// gap filling
		opts.addOption("gf", "gf", false, "The indicator for gap filling! Default: <f>.");
		// tip length
		opts.addOption("tl", "tiplength", false, "The maximum tip length.");
		// iqr time
		opts.addOption("iqrt", "iqrtime", false, "The IQR time for defined repeat outlier.");
		// only for minimap format , default 8;
		opts.addOption("mmcm", "mmcm", true, "The filter parameter only for last column format of minimap, default:8.");
		return opts;
	}
	
	private static Parameter parsering(CommandLine cl, HelpFormatter f, Options opts) throws MissingArgumentException
	{
		Parameter paras = new Parameter();
		// checking the aligned setting valid or not
		if((!cl.hasOption("m5")) && (!cl.hasOption("m4")))
		{
			logger.error("The aligned file could not be null!");
			f.printHelp("Usage: java -jar agisps.jar -c <Contig File> -a <Aligned File> or"
					+ "\t java -jar agisps.jar -x <Parameters File>\n", opts);
			throw new MissingArgumentException("The aligned file could not be null!");
		}
		// checking the contig setting or not
		if(!cl.hasOption("c"))
		{
			logger.error("The argument '-c' setting could not be null!");
			f.printHelp("Usage: java -jar agisps.jar -c <Contig File> -a <Aligned File> or"
					+ "\t java -jar agisps.jar -x <Parameters File>\n", opts);
			throw new MissingArgumentException("The aligned file could not be null!");
		}
		// checking the out
		if(!cl.hasOption("o"))
		{
			logger.error("The argument '-o' setting could not be null!");
			f.printHelp("Usage: java -jar agisps.jar -c <Contig File> -a <Aligned File> or"
					+ "\t java -jar agisps.jar -x <Parameters File>\n", opts);
			throw new MissingArgumentException("The aligned file could not be null!");
		}
		// parsering m4
		if(cl.hasOption("m5"))
		{
			paras.setAlgFile(cl.getOptionValue("m5"));
			paras.setType("m5");
		} else if(cl.hasOption("m4"))
		{
			paras.setAlgFile(cl.getOptionValue("m4"));
			paras.setType("m4");
		} else if(cl.hasOption("sam") || cl.hasOption("bam"))
		{
			paras.setAlgFile(cl.getOptionValue("sam"));
			paras.setType("sam");
		}
		// parsering contig;
		if(cl.hasOption("c"))
		{
			paras.setCntFile(cl.getOptionValue("c"));;
		}
		// parsering pb reads;
		if(cl.hasOption("p"))
		{
			paras.setPbFile(cl.getOptionValue("p"));
		}
		// parsering output folder;
		if(cl.hasOption("o"))
		{
			paras.setOutFolder(cl.getOptionValue("o"));
		}
		// parsering minimum contig length
		if(cl.hasOption("micl"))
		{
			paras.setMinContLen(Integer.valueOf(cl.getOptionValue("micl")));
		}
		// parsering minimum pacbio lng read
		if(cl.hasOption("mipl"))
		{
			paras.setMinPBLen(Integer.valueOf(cl.getOptionValue("mipl")));
		}
		// parsering identity
		if(cl.hasOption("i"))
		{
			paras.setIdentity(Double.valueOf(cl.getOptionValue("i")));
		} 
		// parsering minimum overlap length
		if(cl.hasOption("mioll"))
		{
			paras.setMinOLLen(Integer.valueOf(cl.getOptionValue("mioll")));
		}
		// parsering minimum overlap ratio
		if(cl.hasOption("miolr"))
		{
			paras.setMinOLRatio(Double.valueOf(cl.getOptionValue("miolr")));
		}
		// parsering maximum overhang length
		if(cl.hasOption("maohl"))
		{
			paras.setMaxOHLen(Integer.valueOf(cl.getOptionValue("maohl")));
		}
		// parsering maximum overhang ratio
		if(cl.hasOption("maohr"))
		{
			paras.setMaxOHRatio(Double.valueOf(cl.getOptionValue("maohr")));
		}
		// parsering maximum ending length
		if(cl.hasOption("mael"))
		{
			paras.setMaxEndLen(Integer.valueOf(cl.getOptionValue("mael")));
		}
		// parsering maximum ending ratio
		if(cl.hasOption("maer"))
		{
			paras.setMaxEndRatio(Double.valueOf(cl.getOptionValue("maer")));
		}
		// parsering minimum support links
		if(cl.hasOption("misl"))
		{
			paras.setMinSupLinks(Integer.valueOf(cl.getOptionValue("misl")));
		}
		// parsering use overlap link
		if(cl.hasOption("uoll"))
		{
			
			paras.setUseOLLink(true);
		}
		// parsering delete error prone edges ratio
		if(cl.hasOption("r"))
		{
			paras.setRatio(Double.valueOf(cl.getOptionValue("r")));
		}
		// parsering masking repeat
		if(cl.hasOption("mr"))
		{
			paras.setRepMask(true);
		}
		// parsering gap filling
		if(cl.hasOption("gf"))
		{
			paras.setGapFilling(true);
		}
		// parsing tip length
		if(cl.hasOption("tl"))
		{
			paras.setTipLength(Integer.valueOf(cl.getOptionValue("tl")));
		}
		// parsing iqr time
		if(cl.hasOption("iqrt"))
		{
			paras.setIqrTime(Double.valueOf(cl.getOptionValue("iqrt")));
		}
		return paras;
	}
}

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
//import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.MissingArgumentException;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.OptionGroup;
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
		try {
			Options opts = Main.initOpts();
			CommandLineParser parser = new DefaultParser();
			CommandLine cl = parser.parse(opts, args);
//			HelpFormatter f = new HelpFormatter();
			Parameter paras = null;
			
			if (cl.hasOption("h")) {
//				f.printHelp("java -jar agisps.jar -c <Contig File> -a <Aligned File> -o <Output Folder> or"
//						+ " java -jar agisps.jar -x <Parameters File>\n", opts, true);
				logger.info(printUsageInfo());
				System.exit(0);
			} else {
				// for parse the xml configure, all the other parameters set by command line will dismissed
				if (cl.hasOption("x")) {
					logger.info("Parsing the xml configure, all the other parameters set by command line will dismissed!");
					String xmlFile = cl.getOptionValue("x");
					XMLParser xmlParser = new XMLParser();
					paras = xmlParser.parseXML(xmlFile);
				} else {				
					paras = Main.parsering(cl,  opts);
				}
				// The body for executing scaffolding process
				logger.info("Launching...");
				Scaffolder scaffolder = new Scaffolder(paras);
				scaffolder.scaffolding();
				logger.info("Ending...");
				long end = System.currentTimeMillis();
				logger.info("Scaffolding erase time: " + (end - start)/1000 + " s.");
			}
		} catch (MissingArgumentException e) {
			logger.error(Main.class.getName() + "\t" + e.getMessage());
			logger.info(printUsageInfo());
		} catch (ParseException e) {
			logger.error(Main.class.getName() + "\t"+ e.getMessage());
			logger.info(printUsageInfo());
		} catch (Exception e) {
			logger.error(Main.class.getName() + "\t" + e.getMessage());
			logger.info(printUsageInfo());
		}
	}
	
	private static Options initOpts()
	{
		Options opts = new Options();
		// parameters sepcified by xml
		//OptionGroup xmlog = new OptionGroup();
		// xml file
		Option xml = new Option("x", "xml", true, "The parameter XML file! All command-line parameters would be omitted if set.");
		//xmlog.addOption(xml);
		//opts.addOptionGroup(xmlog);
		opts.addOption(xml);
		
		// cml options
		//OptionGroup inputog = new OptionGroup();
		// contig
		Option cnt = new Option("c", "contig", true, "The file of pre-assembled contigs in fasta format.");
		opts.addOption(cnt);
		// pacbio reads
//		opts.addOption("p", "pacbio", true, "The PacBio Long Reads file for gap filling if requried.");
		// m5 format
//		Option m5 = new Option("m5", "m5", true, "The aligned file in -m 5 format of blasr.");
//		opts.addOption("m5", "m5", true, "The aligned file in -m 5 format of blasr.");
//		inputog.addOption(m5);
		// m4 format
//		opts.addOption("m4", "m4", true, "The aligned file in -m 4 format of blasr.");
//		Option m4 = new Option("m4", "m4", true, "The aligned file in -m 4 format of blasr.");
//		inputog.addOption(m4);
		// minimap format
//		opts.addOption("mm", "mm", true, "The aligned file in minimap output format.");
//		Option mm = new Option("mm", "mm", true, "The aligned file in minimap output format.");
//		inputog.addOption(mm);
		// sam format
//		opts.addOption("sam", "sam", true, "The aligned file in sam format.");
//		Option sam = new Option("sam", "sam", true, "The aligned file in sam format.");
//		inputog.addOption(sam);
		Option alnFile = new Option("a", "alignedFile", true, "Required. The aligned file by using Minimap or BLASR mapper!");
		//inputog.addOption(alnFile);
		opts.addOption(alnFile);
		Option type = new Option("t", "type", true, "Requried. The aligned file format, supported: mm for Minimap, m4, m5 or sam for BLASR.");
		//inputog.addOption(type);
		//opts.addOptionGroup(inputog);
		opts.addOption(type);
		// output folder
		opts.addOption("o", "output", true, "The scaffolding output folder.");
		// help 
		opts.addOption("h", "help", false, "The help infomation!");
		// minimum contig length
		opts.addOption("micl", "miniCntLen", true, "The minimum contig's length for scaffolding! Default:<3000>.");
		// minimum pacbio long read length
//		opts.addOption("mipl", "miniPBLen", true, "The minimum PacBio long read's length for scaffolding! Default:<5000>.");
		// identity
		opts.addOption("i", "identity", true, "The identity threshold for blasr alignment! Default: <0.8>.");
		// minimum overlap length
		opts.addOption("mioll", "miniOLLen", true, "The minimum overlap length threshold for blasr alignment! Default: <2400>.");
		// minimum overlap ratio
		opts.addOption("miolr", "miniOLRatio", true, "The minimum overlap ratio threshold for blasr alignment! Default: <0.8>.");
		// maximum overhang length;
		opts.addOption("mxohl", "maxOHLen", true, "The maximum overhang length threshold for blasr alignment! Default: <300>.");
		// maximum overhang ratio;
		opts.addOption("mxohr", "maxOHRatio", true, "The maximum overhang ratio threshold for blasr alignment! Default: <0.1>.");
		// maximum ending length
		opts.addOption("mxel", "maxEndLen", true, "The maximum ending length of PacBio's Long Read! Default: <300>.");
		// maximum ending ratio
		opts.addOption("mxer", "maxEndRatio", true, "The maximum ending ratio of PacBio's Long Read! Default: <0.1>.");
		// minimum support links;
		opts.addOption("misl", "miniSupLinks", true, "The minimum support links number! Default: <1>.");
		// use overlap link
//		opts.addOption("doll", "doll", false, "The indicator for using overlap relationship of contig! Default: <f>, discard overlap relationship.");
		// deleting error prone edges ratio;
		opts.addOption("r", "ratio", true, "The ratio for deleting error prone edges! Default: <0.2>.");
		// masking repeat
//		opts.addOption("mr", "mr", false, "The indicator for masking repeats! Default: <f>.");
		// gap filling
//		opts.addOption("gf", "gf", false, "The indicator for gap filling! Default: <f>.");
		// tip length
//		opts.addOption("tl", "tiplength", false, "The maximum tip length.");
		// iqr time
		opts.addOption("iqrt", "iqrtime", false, "The IQR time for defined repeat outlier.");
		// only for minimap format , default 8;
		opts.addOption("mmcm", "mmcm", true, "The filter parameter only for last column format of minimap, default:8.");
		return opts;
	}
	
	private static Parameter parsering(CommandLine cl, Options opts) throws MissingArgumentException
	{
		Parameter paras = new Parameter();
		// checking the aligned setting valid or not
//		if((!cl.hasOption("m5")) && (!cl.hasOption("m4")))
//		{
//			//logger.error("The aligned file could not be null!");
////			f.printHelp("Usage: java -jar agisps.jar -c <Contig File> -a <Aligned File> or"
////					+ "\t java -jar agisps.jar -x <Parameters File>\n", opts);
//			throw new MissingArgumentException("The aligned file could not be null!");
//		}
		if(!cl.hasOption("a"))
		{
			throw new MissingArgumentException("The aligned file could not be null!");
		}
		if(!cl.hasOption("t"))
		{
			throw new MissingArgumentException("The type of aligned file could not be null!");
		}
		// checking the contig setting or not
		if(!cl.hasOption("c"))
		{
			//logger.error("The argument '-c' setting could not be null!");
//			f.printHelp("Usage: java -jar agisps.jar -c <Contig File> -a <Aligned File> or"
//					+ "\t java -jar agisps.jar -x <Parameters File>\n", opts);
			throw new MissingArgumentException("The draft assemblied genome could not be null!");
		}
		// checking the out
		if(!cl.hasOption("o"))
		{
			//logger.error("The argument '-o' setting could not be null!");
//			f.printHelp("Usage: java -jar agisps.jar -c <Contig File> -a <Aligned File> or"
//					+ "\t java -jar agisps.jar -x <Parameters File>\n", opts);
			throw new MissingArgumentException("The ouput folder could not be null!");
		}
		// parsing file type
//		if(cl.hasOption("m5"))
//		{
//			paras.setAlgFile(cl.getOptionValue("m5"));
//			paras.setType("m5");
//		} else if(cl.hasOption("m4"))
//		{
//			paras.setAlgFile(cl.getOptionValue("m4"));
//			paras.setType("m4");
//		} else if(cl.hasOption("sam") || cl.hasOption("bam"))
//		{
//			paras.setAlgFile(cl.getOptionValue("sam"));
//			paras.setType("sam");
//		} else if(cl.hasOption("mm"))
//		{
//			paras.setAlgFile(cl.getOptionValue("mm"));
//			paras.setType("mm");
//		}
		if(cl.hasOption("a"))
		{
			paras.setAlgFile(cl.getOptionValue("a"));
		}
		if(cl.hasOption("t"))
		{
			paras.setType(cl.getOptionValue("t"));
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
			int value = Integer.valueOf(cl.getOptionValue("micl"));
			if(value < 200)
			{
				paras.setMinContLen(200);
				logger.info("The minimum contig's length should be large than 200 bp!");
				logger.info("Mandatorily set to 200 bp!");
			} else
			{
				paras.setMinContLen(value);
			}
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
	
	private static String printUsageInfo()
	{
		StringBuilder sb = new StringBuilder();
		sb.append("\nLRScaf is a scaffolder by using Thired Generation Sequencing data to "
				+ "scaffold any draft assemblied genomes. it supports command line argument or XML configure file!\n");
		sb.append("Usage:\nXML: java -jar LRScaf.jar -x <XML_File>\n");
		sb.append("CML: java -jar LRScaf.jar -c <Draft_Assemblies_File> -a <Alinged_File> -t <m5|m4|mm|sam> [Other_Options]\n");
		sb.append("Argument Details:\n");
		sb.append("XML configure file:\n");
		sb.append("-x\t--xml\t<arg>\tThe parameter XML file! All command-line parameters would be omitted if set.\n");
		sb.append("Command line options:\n");
		sb.append("-c\t--contig\t<arg>\tRequired. The file of pre-assembled contigs in fasta format.\n");
		sb.append("## Input file:\n");
		sb.append("-a\t--alignedFile\tRequired. The aligned file by using Minimap or BLASR mapper!\n");
		sb.append("-t\t--type\tRequried. The aligned file format, supported: 1) mm for Minimap; 2)m4, m5 or sam for BLASR.\n");
		sb.append("## Output folder:\n");
		sb.append("-o\t--output\t<arg>\tRequired. The scaffolding output folder.\n");
		sb.append("## Other options:\n");
		sb.append("-i\t--identity\t<arg>\tThe identity threshold for blasr alignment! Default: <0.8>.\n");
		sb.append("-r\t--ratio\t<arg>\tThe ratio for deleting error prone edges! Default: <0.2>.\n");
		sb.append("-misl\t--miniSupLinks\t<arg>\tThe minimum support links number! Default: <1>.\n");
		sb.append("-iqrt\t--iqrtime\tThe IQR time for defined repeat outlier.Default: <1.5>\n");
		sb.append("-micl\t--miniCntLen\t<arg>\tThe minimum contig's length (must be larget than 200 bp) for scaffolding! Default:<200>.\n");
		sb.append("-mioll\t--miniOLLen\t<arg>\tThe minimum overlap length threshold for alignment! Default: <160>.\n");
		sb.append("-miolr\t--miniOLRatio\t<arg>\tThe minimum overlap ratio threshold for alignment! Default: <0.8>.\n");
		sb.append("-mxel\t--maxEndLen\t<arg>\tThe maximum ending length of Long Read! Default: <300>.\n");
		sb.append("-mxer\t--maxEndRatio\t<arg>\tThe maximum ending ratio of Long Read! Default: <0.1>.\n");
		sb.append("-mxohl\t--maxOHLen\t<arg>\tThe maximum overhang length threshold for alignment! Default: <300>.\n");
		sb.append("-mxohr\t--maxOHRatio\t<arg>\tThe maximum overhang ratio threshold for alignment! Default: <0.1>.\n");
		sb.append("-mmcm\t--mmcm\t<arg>\tThe filter parameter only for last column format of Minimap, default:<8>.\n");
 		sb.append("-h\t--help\tThe help infomation!\n");
		sb.append("Please report issues at <github_project_website>!");
		return sb.toString();
	}
}

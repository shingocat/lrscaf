/*
*File: agis.ps.Scaffolder.java
*User: mqin
*Email: mqin@ymail.com
*Date: 2015年12月28日
*/
package agis.ps;

import java.io.IOException;
import java.net.MalformedURLException;
import java.net.URL;
import java.io.File;
import java.io.FileNotFoundException;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;

import org.biojava.nbio.core.sequence.io.FastaReaderHelper;
import org.biojava.nbio.core.sequence.DNASequence;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import agis.ps.link.PBLink;
import agis.ps.util.Color;
import agis.ps.util.DotGraphFileWriter;
import agis.ps.util.EdgeBundler;
import agis.ps.util.LinkBuilder;
import agis.ps.util.M5Reader;
import agis.ps.util.Parameter;
import agis.ps.util.PathBuilder;
import htsjdk.samtools.AlignmentBlock;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamInputResource;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;

public class Scaffolder {
	final static Logger logger = LoggerFactory.getLogger(Scaffolder.class);
	
	private HashMap<String, Object> paras;
	private String cFilePath;
	private String aFilePath;
	private String gFilePath;
	private String outFolder;
	private String type;
	private boolean m5Header;
	private LinkedHashMap<String, DNASequence> contigs;
	private SamReader samReader;
	private List<SimplePath> simPaths;
	private List<PBLink> pbLinks;
	private List<Edge> edges;
	
	public Scaffolder(Parameter paras)
	{
		this.cFilePath = paras.getCntFile();
		this.aFilePath = paras.getAlgFile();
		this.type = paras.getType();
		this.outFolder = paras.getOutFolder();
	}
	
	public Scaffolder(String cFilePath, String aFilePath)
	{
		this.cFilePath = cFilePath;
		this.aFilePath = aFilePath;
	}
	
	public Scaffolder(HashMap<String, Object> paras)
	{
		this.paras = paras;
		this.cFilePath = (String)paras.get("CONTIG");
		this.aFilePath = (String)paras.get("ALIGNED");
		this.type = (String)paras.get("TYPE");
//		this.gFilePath = (String)paras.get("DOTGRAPH");
		this.outFolder = (String) paras.get("OUTPUT");
	}
	
	public void scaffolding()
	{
		System.out.println("Starting....");
		try
		{
			readContigs(cFilePath);
			List<M5Record> m5Records = null;
			if(type.equalsIgnoreCase("m"))
			{
				m5Records = readM5Aligned(aFilePath);
				pbLinks = LinkBuilder.m5Record2Link(m5Records, null);
			} else
			{
				readSAMAligned(aFilePath);
			}
			edges = EdgeBundler.pbLinkBundling(pbLinks, null);
			logger.debug("edges size: " + edges.size());
			PathBuilder.buildHamiltonPath(edges);
			//listContigs();
			//listAligns();
//			if(gFilePath != null)
//			{
//				DotGraphFileWriter.writePBLink(gFilePath, pbLinks);
//				DotGraphFileWriter.writeEdge(gFilePath, edges);
//				DotGraphFileWriter dGFW = new DotGraphFileWriter(gFilePath, simPaths);
//				dGFW.write();
//			}
		} catch(NullPointerException e)
		{
			logger.debug(e.getMessage());
		} catch(MalformedURLException e)
		{
			logger.debug(e.getMessage());
		} catch(FileNotFoundException e)
		{
			logger.debug(e.getMessage());
		} catch(IOException e)
		{
			logger.debug(e.getMessage());
		} 
		
		System.out.println("Ending....");
	}
	
	public void readContigs(String cFilePath) throws IOException
	{
		if(cFilePath == null || cFilePath.length() == 0)
		{
			logger.error("The contig file was null or not setted!");
			logger.debug("The contig file was null or not setted!");
			logger.info("The contig file was null or not setted!");
			return;
		}
		File cFile = new File(cFilePath);
		setContigs(FastaReaderHelper.readFastaDNASequence(cFile));
	}
	
	public List<M5Record> readM5Aligned(String aFilePath)
	{
		if(aFilePath.isEmpty())
		{	
			logger.error("The aligned file was null or not setted!");
			logger.debug("The aligned file was null or not setted!");
			logger.info("The aligned file was null or not setted!");
			return null;
		}
		M5Reader reader = new M5Reader(aFilePath);	
		List<M5Record> m5Records = reader.read();
		HashMap<String, String> pSet = new HashMap<String, String>();
		
		for(M5Record m5 : m5Records)
		{
			String c = m5.gettName() + "==" + m5.getqStart() + ";";
			if(pSet.containsKey(m5.getqName()))
			{
				pSet.put(m5.getqName(), pSet.get(m5.getqName()) + c);
			} else
			{
				pSet.put(m5.getqName(), c);
			}
		}

//		for(String s : pSet.keySet())
//		{
//			//logger.debug(s);
//			logger.debug(s + "\t" + pSet.get(s));
//		}
		//this.m5RecordToSimplePath(pSet);
		return m5Records;
	}
	
	public void readSAMAligned(String aFilePath) throws NullPointerException, MalformedURLException,IOException
	{
		if(aFilePath == null || aFilePath.length() == 0)
		{
			logger.debug("The aligned file was null or not setted!");
			logger.info("The aligned file was null or not setted!");
			return;
		}
		SamReaderFactory factory = SamReaderFactory.makeDefault().enable(SamReaderFactory.Option.INCLUDE_SOURCE_IN_RECORDS, 
				SamReaderFactory.Option.VALIDATE_CRC_CHECKSUMS).validationStringency(ValidationStringency.LENIENT);
		SamInputResource resource = SamInputResource.of(new File(aFilePath)); //.index(new URL("http://broadinstitute.org/my.bam.bai"));
		samReader = factory.open(resource);
		if(simPaths != null)
		{
			simPaths.clear();
		}
		simPaths = new ArrayList<SimplePath>();
		Iterator<SAMRecord> it = samReader.iterator();
//		HashMap<String, HashMap<String, Integer>> pSet = new HashMap<String, HashMap<String, Integer>>();
		HashMap<String, String> pSet = new HashMap<String, String>();
		while(it.hasNext())
		{
			SAMRecord r = it.next();
//			logger.debug("Read name:" + r.getReadName() + "\tReference name:" + r.getReferenceName() +
//					"\tAligned start:" + r.getAlignmentStart() + "\tAligned end:" + r.getAlignmentEnd() +
//					"\tFLAG:" + r.getFlags());
			String c =  r.getReferenceName() + "==" +  r.getAlignmentStart() + ";";
//			HashMap<String, Integer> cSet = new HashMap<String,Integer>();
//			cSet.put(r.getReferenceName(), r.getReadPositionAtReferencePosition(r.getAlignmentStart()));
			if(pSet.containsKey(r.getReadName()))
			{
				pSet.put(r.getReadName(), pSet.get(r.getReadName()) + c);
			} else
			{
				pSet.put(r.getReadName(), c);
			}
			
		}
		for(String s : pSet.keySet())
		{
			//logger.debug(s);
			logger.debug(s + "\t" + pSet.get(s));
		}
	}
	
	public void m5RecordToSimplePath(Map<String, String> pSet)
	{
		if(simPaths != null)
		{
			simPaths.clear();
		} else
		{
			simPaths = new ArrayList<SimplePath>();
		}
		for(String s : pSet.keySet())
		{
			String [] arrs = pSet.get(s).split(";");
			if(arrs.length >= 2)
			{
				for(int i = 0; i < arrs.length - 1; i++)
				{
					SimplePath sp = new SimplePath();
					String [] ar1 = arrs[i].split("==");
					String [] ar2 = arrs[i + 1].split("==");
					if(Integer.valueOf(ar1[1]) <= Integer.valueOf(ar2[1]))
					{
						sp.setStart(ar1[0]);
						sp.setEnd(ar2[0]);
						sp.setLabel(ar1[1] + ":" + ar2[1]);
						sp.setColor(Color.BLUE);
					} else
					{
						sp.setStart(ar2[0]);
						sp.setEnd(ar1[0]);
						sp.setLabel(ar2[1] + ":" + ar1[1]);
						sp.setColor(Color.BLUE);
					}
					simPaths.add(sp);
				}
			}
		}
	}
	
	public void listContigs()
	{
		if(contigs == null || contigs.isEmpty())
		{
			return;
		}
		for(String id : contigs.keySet())
		{
			System.out.println("Id is:\t" + id);
		}
	}
	
	public void listAligns()
	{
		if(samReader == null)
			return;
		for(SAMRecord sr : samReader)
		{
			System.out.println("Read Name:\t" + sr.getReadName());
			System.out.println("Aligned Name:\t" + sr.getReferenceName());
		}
	}

	public String getcFilePath() {
		return cFilePath;
	}

	public void setcFilePath(String cFilePath) {
		this.cFilePath = cFilePath;
	}

	public String getaFilePath() {
		return aFilePath;
	}

	public void setaFilePath(String aFilePath) {
		this.aFilePath = aFilePath;
	}

	public LinkedHashMap<String, DNASequence> getContigs() {
		return contigs;
	}

	public void setContigs(LinkedHashMap<String, DNASequence> contigs) {
		this.contigs = contigs;
	}
	
	

}



/*
*File: agis.ps.util.RepeatFinder.java
*User: mqin
*Email: mqin@ymail.com
*Date: 2016年5月12日
*/
package agis.ps.util;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Vector;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import agis.ps.file.ContigCoverageWriter;
import agis.ps.link.M5Record;
import agis.ps.link.MRecord;

public class RepeatFinder {
	
	public static Logger logger = LoggerFactory.getLogger(RepeatFinder.class);
	private Parameter paras;
	private List<String> repeats;
	private Map<String, List<MRecord>> cntMaps; 
	private int minCNTLen = 0;
	private int minPBLen = 0;
	private double identity = 0.0d;
	
	public RepeatFinder(Parameter paras)
	{
		this.paras = paras;
	}
	
	public List<String> findRepeat()
	{
		long start = System.currentTimeMillis();
		if(cntMaps == null)
			cntMaps = new HashMap<String, List<MRecord>>();
		cntMaps.clear();
		File file = null;
		FileReader fr = null;
		BufferedReader br = null;
		String alnFile = paras.getAlgFile();
		try {
			file = new File(alnFile);
			if (!file.exists()) {
				logger.error(this.getClass().getName() + "\t" + "The m5 file do not exist!");
				 return null;
			}
			fr = new FileReader(file);
			br = new BufferedReader(fr);
			String line = null;
			String[] arrs = null;
			while(true)
			{
				line = br.readLine();
				if(line != null)
				{
					arrs = line.split("\\s+");
					if (arrs[0].equalsIgnoreCase("qName") && arrs[1].equalsIgnoreCase("qLength"))
						continue;
					MRecord m = MRecordValidator.validate(arrs, paras);
					if(m != null)
					{
						String pbId = m.getqName();
						if(cntMaps.containsKey(pbId))
						{
							List<MRecord> ms = cntMaps.get(pbId);
							if(!ms.contains(m))
							{
								ms.add(m);
								cntMaps.replace(pbId, ms);
							}
						} else
						{
							List<MRecord> ms = new Vector<MRecord>(10);
							ms.add(m);
							cntMaps.put(pbId, ms);
						}
					}
				} else
				{
					break;
				}
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
				if(br != null)
					br.close();
			} catch(Exception e)
			{
				logger.error(this.getClass().getName() + "\t" + e.getMessage());
			}
		}
		repeats = this.findRepeat(cntMaps);
		long end = System.currentTimeMillis();
		logger.info("Finding repeat, erase time: " + (end - start) / 1000 + " s");
		return repeats;
	}
	
	public List<String> findRepeat(Map<String, List<MRecord>> args)
	{
		List<String> repeats = new Vector<String>(30);
		// contig id and pacbio id in List<String> 
		Map<String, List<String>> cntCounts = new HashMap<String, List<String>>();
		for(String s : args.keySet())
		{
			List<MRecord> cnts = args.get(s);
			for(MRecord m : cnts)
			{
				if(cntCounts.containsKey(m.gettName()))
				{
					List<String> temp = cntCounts.get(m.gettName());
					int tLen = m.gettLength();
					int tStart = m.gettStart();
					int tEnd = m.gettEnd();
					int tLeftLen = tStart;
					int tRightLen = tLen - tEnd;
					int maxOHLen = paras.getMaxOHLen();
					int defOHLen = (int) (tLen * paras.getMaxOHRatio());
					if(maxOHLen <= defOHLen)
						defOHLen = maxOHLen;
					if(tLeftLen > defOHLen )
						continue;
					if(tRightLen > defOHLen)
						continue;
					if(!temp.contains(s))
						temp.add(s);
				} else
				{
					Vector<String> temp = new Vector<String>(30);
					int tLen = m.gettLength();
					int tStart = m.gettStart();
					int tEnd = m.gettEnd();
					int tLeftLen = tStart;
					int tRightLen = tLen - tEnd;
					int maxOHLen = paras.getMaxOHLen();
					int defOHLen = (int) (tLen * paras.getMaxOHRatio());
					if(maxOHLen <= defOHLen)
						defOHLen = maxOHLen;
					if(tLeftLen > defOHLen )
						continue;
					if(tRightLen > defOHLen)
						continue;
					temp.add(s);
					cntCounts.put(m.gettName(), temp);
					temp = null;
				}
			}
		}
		List<Integer> covs = new Vector<Integer>();
		for(String id : cntCounts.keySet())
		{
//			logger.debug("id = " +  id + ", cov = " + cntCounts.get(id).size());
			covs.add(cntCounts.get(id).size());
		}
		int mean =  MathTool.mean(covs);
		int sd = MathTool.sd(covs);
		int median = MathTool.median(covs);
		int upper = median + 2 * sd;
		logger.info("Mean cov = " + mean);
		logger.info("Median cov = " + median);
		logger.info("S.D. = " + sd);
		logger.info("Mean Range 95%: [" + (mean - 2 * sd) + " : " + (mean + 2 * sd) + "]");
		logger.info("Mean Range 99%: [" + (mean - 3 * sd) + " : " + (mean + 3 * sd) + "]");
		logger.info("Median Range 95%: [" + (median - 2 * sd) + " : " + (median + 2 * sd) + "]");
		logger.info("Median Range 99%: [" + (median - 3 * sd) + " : " + (median + 3 * sd) + "]");
		logger.info("Pesudo repeat contigs");
//		Map<String, List<String>> repeat = new HashMap<String, List<String>>();
		for(String id : cntCounts.keySet())
		{
			if(cntCounts.get(id).size() > upper)
			{
//				repeat.put(id, cntCounts.get(id));
				repeats.add(id);
				logger.info(id + " might be repeat!");
			}
		}
		logger.info("Repeat count: " + repeats.size());
		ContigCoverageWriter ccw = new ContigCoverageWriter(paras, cntCounts);
		ccw.write();;
		return repeats;
	}
}



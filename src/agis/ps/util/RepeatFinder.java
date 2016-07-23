/*
*File: agis.ps.util.RepeatFinder.java
*User: mqin
*Email: mqin@ymail.com
*Date: 2016年5月12日
*/
package agis.ps.util;

import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Vector;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import agis.ps.file.ContigCoverageWriter;
import agis.ps.link.MRecord;

public class RepeatFinder {
	
	public static Logger logger = LoggerFactory.getLogger(RepeatFinder.class);
	private Parameter paras;
	
	public RepeatFinder(Parameter paras)
	{
		this.paras = paras;
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
		logger.debug("Mean cov = " + mean);
		logger.debug("Median cov = " + median);
		logger.debug("S.D. = " + sd);
		logger.debug("Mean Range 95%: [" + (mean - 2 * sd) + " : " + (mean + 2 * sd) + "]");
		logger.debug("Mean Range 99%: [" + (mean - 3 * sd) + " : " + (mean + 3 * sd) + "]");
		logger.debug("Median Range 95%: [" + (median - 2 * sd) + " : " + (median + 2 * sd) + "]");
		logger.debug("Median Range 99%: [" + (median - 3 * sd) + " : " + (median + 3 * sd) + "]");
		logger.debug("Pesudo repeat contigs");
//		Map<String, List<String>> repeat = new HashMap<String, List<String>>();
		for(String id : cntCounts.keySet())
		{
			if(cntCounts.get(id).size() > upper)
			{
//				repeat.put(id, cntCounts.get(id));
				repeats.add(id);
				logger.debug(id + " might be repeat!");
			}
		}
		logger.debug("Repeat count: " + repeats.size());
		ContigCoverageWriter ccw = new ContigCoverageWriter(paras, cntCounts);
		ccw.write();;
		return repeats;
	}
}



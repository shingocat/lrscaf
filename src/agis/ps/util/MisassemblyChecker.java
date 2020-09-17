/** 
** Usage: TODO
** Author: mqin
** Email: mqin@outlook.com
** Date: 2017年7月13日
*/
package agis.ps.util;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.Map;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import agis.ps.file.SequenceCoverageWriter;
import agis.ps.file.ContigMisassemblyWriter;
import agis.ps.link.MRecord;
//import agis.ps.seqs.Contig;
import agis.ps.seqs.Sequence;

public class MisassemblyChecker {
	
	private static Logger logger = LoggerFactory.getLogger(MisassemblyChecker.class);
	private Parameter paras;
//	private Map<String, Contig> cnts;
	private Map<String, Sequence> seqs;
	
	public MisassemblyChecker(Parameter paras)
	{
		this.paras = paras;
	}
	
	public static boolean checking2(Parameter paras, Map<String, Sequence> seqs, List<MRecord> maRs)
	{
		boolean isTrue = false;
		int maxOHLen = paras.getMaxOHLen();
		double maxOHRatio = paras.getMaxOHRatio();
		if(maRs == null || maRs.isEmpty())
			return isTrue;
		List<MRecord> tempRs = new ArrayList<MRecord>();
		int size = maRs.size();
		boolean isExist = false;
		double threshold = 0.2;
		for(int i = 0; i < size; i++)
		{ // define overlap region is repeat;
			MRecord target = maRs.get(i);
			for(int j = 0; j < size; j++)
			{
				if(j != i)
				{
					MRecord query = maRs.get(j);
					int tLRStart = target.getqStart();
					int tLREnd = target.getqEnd();
					int qLRStart = query.getqStart();
					int qLREnd = query.getqEnd();
					if(qLRStart < tLREnd && qLREnd > tLRStart)
					{
						int olRegion = 0;
						int [] positions = {tLRStart, tLREnd, qLRStart, qLREnd};
						Arrays.sort(positions);
						olRegion = positions[2] - positions[1];
						double tOLRatio = (double) olRegion / (tLREnd - tLRStart);
						double qOLRatio = (double) olRegion / (qLREnd - qLRStart);
						if(tOLRatio > threshold || qOLRatio > threshold)
						{
							isExist = true;
							break;
						}
					}
				}
			}
			if(!isExist)
			{
				tempRs.add(target);
				isExist = false;
			}
		}
		maRs = tempRs;
		size = maRs.size();
		if(size == 1)
			return isTrue;
		for(int i = 0; i < size; i++)
		{
			MRecord target = maRs.get(i);
			Strand tCntStrand = target.gettStrand();
			int lrLen = target.getqLength();
			int tLRStart = target.getqStart();
			int tLREnd = target.getqEnd();
			int tCntLen = target.gettLength();
			int tCntStart = target.gettStart();
			int tCntEnd = target.gettEnd();
			int defOHLen = (int) (lrLen * maxOHRatio);
			if(defOHLen > maxOHLen)
				defOHLen = maxOHLen;
			int tLRLeftLen = tLRStart;
			int tLRRightLen = lrLen - tLREnd;
			int tCntLeftLen = tCntStart;
			int tCntRightLen = tCntLen - tCntEnd;
			if(tCntStrand.equals(Strand.REVERSE))
			{
				int temp = tCntLeftLen;
				tCntLeftLen = tCntRightLen;
				tCntRightLen = temp;
			}
			if(tLRLeftLen > defOHLen && tLRRightLen < defOHLen)
			{
				if(tCntLeftLen > defOHLen)
				{
					// change to store all regions, and then sort them and compute the supported long reads 
					for(int j = i + 1; j < size;j++)
					{
						MRecord query = maRs.get(j);
						Strand qCntStrand = query.gettStrand();
						int qLRStart = query.getqStart();
						int qLREnd = query.getqEnd();
						int qCntLen = query.gettLength();
						int qCntStart = query.gettStart();
						int qCntEnd = query.gettEnd();
						int qLRLeftLen = qLRStart;
						int qLRRightLen = lrLen - qLREnd;
						int qCntLeftLen = qCntStart;
						int qCntRightLen = qCntLen - qCntEnd;
						if(qCntStrand.equals(Strand.REVERSE))
						{
							int temp = qCntLeftLen;
							qCntLeftLen = qCntRightLen;
							qCntRightLen = temp;
						}
						int internal = tLRStart - qLREnd;
						if(qLRRightLen > defOHLen)
						{
							if(tCntLeftLen > internal && qCntLeftLen > internal)
							{
								// target 
								MisassemblyRegion mr = new MisassemblyRegion();
								mr.setStart(tCntStart);
								mr.setEnd(tCntEnd);
								mr.setSupportLRs(1);
								Sequence cnt = seqs.get(target.gettName());
								List<MisassemblyRegion> mrs = cnt.getMisassemblies();
								if(!mrs.contains(mr))
									mrs.add(mr);
								// query
								mr = new MisassemblyRegion();
								mr.setStart(qCntStart);
								mr.setEnd(qCntEnd);
								mr.setSupportLRs(1);
								cnt = seqs.get(query.gettName());
								mrs = cnt.getMisassemblies();
								if(!mrs.contains(mr))
									mrs.add(mr);
							}
						}
					}
				}
			} else if(tLRLeftLen < defOHLen && tLRRightLen > defOHLen)
			{
				if(tCntRightLen > defOHLen)
				{
					// change to store all regions, and then sort them and compute the supported long reads 
					for(int j = i + 1; j < size;j++)
					{
						MRecord query = maRs.get(j);
						Strand qCntStrand = query.gettStrand();
						int qLRStart = query.getqStart();
						int qLREnd = query.getqEnd();
						int qCntLen = query.gettLength();
						int qCntStart = query.gettStart();
						int qCntEnd = query.gettEnd();
						int qLRLeftLen = qLRStart;
						int qLRRightLen = lrLen - qLREnd;
						int qCntLeftLen = qCntStart;
						int qCntRightLen = qCntLen - qCntEnd;
						if(qCntStrand.equals(Strand.REVERSE))
						{
							int temp = qCntLeftLen;
							qCntLeftLen = qCntRightLen;
							qCntRightLen = temp;
						}
						int internal = qLRStart - tLREnd;
						if(qLRLeftLen > defOHLen)
						{
							if(tCntRightLen > internal && qCntLeftLen > internal)
							{
								// target 
								MisassemblyRegion mr = new MisassemblyRegion();
								mr.setStart(tCntStart);
								mr.setEnd(tCntEnd);
								mr.setSupportLRs(1);
								Sequence cnt = seqs.get(target.gettName());
								List<MisassemblyRegion> mrs = cnt.getMisassemblies();
								if(!mrs.contains(mr))
									mrs.add(mr);
								// query
								mr = new MisassemblyRegion();
								mr.setStart(qCntStart);
								mr.setEnd(qCntEnd);
								mr.setSupportLRs(1);
								cnt = seqs.get(query.gettName());
								mrs = cnt.getMisassemblies();
								if(!mrs.contains(mr))
									mrs.add(mr);
							}
						}
					}
				}
			} else if(tLRLeftLen > defOHLen && tLRRightLen > defOHLen)
			{
				if(tCntLeftLen > defOHLen && tCntRightLen < defOHLen)
				{ // for right side ending
					for(int j = i + 1; j < size;j++)
					{
						MRecord query = maRs.get(j);
						Strand qCntStrand = query.gettStrand();
						int qLRStart = query.getqStart();
						int qLREnd = query.getqEnd();
						int qCntLen = query.gettLength();
						int qCntStart = query.gettStart();
						int qCntEnd = query.gettEnd();
						int qLRLeftLen = qLRStart;
						int qLRRightLen = lrLen - qLREnd;
						int qCntLeftLen = qCntStart;
						int qCntRightLen = qCntLen - qCntEnd;
						if(qCntStrand.equals(Strand.REVERSE))
						{
							int temp = qCntLeftLen;
							qCntLeftLen = qCntRightLen;
							qCntRightLen = temp;
						}
						int internal = tLRStart - qLREnd;
						if(qLRRightLen > defOHLen)
						{
							if(tCntLeftLen > internal && qCntRightLen > internal)
							{
								// target 
								MisassemblyRegion mr = new MisassemblyRegion();
								mr.setStart(tCntStart);
								mr.setEnd(tCntEnd);
								mr.setSupportLRs(1);
								Sequence seq = seqs.get(target.gettName());
								List<MisassemblyRegion> mrs = seq.getMisassemblies();
								if(!mrs.contains(mr))
									mrs.add(mr);
								// query
								mr = new MisassemblyRegion();
								mr.setStart(qCntStart);
								mr.setEnd(qCntEnd);
								mr.setSupportLRs(1);
								seq = seqs.get(query.gettName());
								mrs = seq.getMisassemblies();
								if(!mrs.contains(mr))
									mrs.add(mr);
							}
						}
					}
				} else if(tCntLeftLen < defOHLen && tCntRightLen > defOHLen)
				{// for left side ending
					for(int j = i + 1; j < size;j++)
					{
						MRecord query = maRs.get(j);
						Strand qCntStrand = query.gettStrand();
						int qLRStart = query.getqStart();
						int qLREnd = query.getqEnd();
						int qCntLen = query.gettLength();
						int qCntStart = query.gettStart();
						int qCntEnd = query.gettEnd();
						int qLRLeftLen = qLRStart;
						int qLRRightLen = lrLen - qLREnd;
						int qCntLeftLen = qCntStart;
						int qCntRightLen = qCntLen - qCntEnd;
						if(qCntStrand.equals(Strand.REVERSE))
						{
							int temp = qCntLeftLen;
							qCntLeftLen = qCntRightLen;
							qCntRightLen = temp;
						}
						int internal = qLRStart - tLREnd;
						if(qLRLeftLen > defOHLen)
						{
							if(tCntRightLen > internal && qCntLeftLen > internal)
							{
								// target 
								MisassemblyRegion mr = new MisassemblyRegion();
								mr.setStart(tCntStart);
								mr.setEnd(tCntEnd);
								mr.setSupportLRs(1);
								Sequence seq = seqs.get(target.gettName());
								List<MisassemblyRegion> mrs = seq.getMisassemblies();
								if(!mrs.contains(mr))
									mrs.add(mr);
								// query
								mr = new MisassemblyRegion();
								mr.setStart(qCntStart);
								mr.setEnd(qCntEnd);
								mr.setSupportLRs(1);
								seq = seqs.get(query.gettName());
								mrs = seq.getMisassemblies();
								if(!mrs.contains(mr))
									mrs.add(mr);
							}
						}
					}
				}
			}
		}
		return isTrue;
	}
	
	public static boolean checking(Parameter paras, Map<String, Sequence> seqs, MRecord record)
	{
		boolean isTrue = false;
		int miniCntLen = paras.getMinContLen();
		String tName = record.gettName();
		int tLength = record.gettLength();
		if(tLength < miniCntLen)
		{
			return isTrue;
		} else
		{
			Strand tStrand = record.gettStrand();
			int lrLen = record.getqLength();
			int lrStart = record.getqStart();
			int lrEnd = record.getqEnd();
			if((lrEnd - lrStart) < 1000)
				return isTrue;
			int cntStart = record.gettStart();
			int cntEnd = record.gettEnd();
			int maxEndLen = paras.getMaxEndLen();
			double maxEndRatio = paras.getMaxEndRatio();
			int defEndLen = (int) (lrLen * maxEndRatio);
			int maxOHLen = paras.getMaxOHLen();
			double maxOHRatio = paras.getMaxOHRatio() * 2;
			int defOHLen = (int) (lrLen * maxOHRatio);
			if(defOHLen > maxOHLen)
				defOHLen = maxOHLen;
			int lrLeftLen = lrStart;
			int lrRightLen = lrLen - lrEnd;
			int cntLeftLen = cntStart;
			int cntRightLen = tLength - cntEnd;
			Sequence contig = seqs.get(tName);
			if(tStrand.equals(Strand.REVERSE))
			{
				int temp = cntLeftLen;
				cntLeftLen = cntRightLen;
				cntRightLen = temp;
			}
			if(lrLeftLen > defOHLen && lrRightLen < defOHLen)
			{
				if(cntLeftLen > defOHLen)
				{
					// change to store all regions, and then sort them and compute the supported long reads 
					MisassemblyRegion mr = new MisassemblyRegion();
					mr.setStart(cntStart);
					mr.setEnd(cntEnd);
					mr.setSupportLRs(1);
					contig.addMisassemblyRegion(mr);
					logger.info(record.getqName());
					isTrue = true;
				}
			} else if(lrLeftLen < defOHLen && lrRightLen > defOHLen)
			{
				if(cntRightLen > defOHLen)
				{
					// change to store all regions, and then sort them and compute the supported long reads 
					MisassemblyRegion mr = new MisassemblyRegion();
					mr.setStart(cntStart);
					mr.setEnd(cntEnd);
					mr.setSupportLRs(1);
					contig.addMisassemblyRegion(mr);
					logger.info(record.getqName());
					isTrue = true;
				}
			} else if(lrLeftLen > defOHLen && lrRightLen > defOHLen)
			{
				if(cntLeftLen > defOHLen && cntRightLen < defOHLen)
				{ // for right side ending
					MisassemblyRegion mr = new MisassemblyRegion();
					mr.setStart(cntStart);
					mr.setEnd(cntEnd);
					mr.setSupportLRs(1);
					contig.addMisassemblyRegion(mr);
					logger.info(record.getqName());
					isTrue = true;
				} else if(cntLeftLen < defOHLen && cntRightLen > defOHLen)
				{// for left side ending
					MisassemblyRegion mr = new MisassemblyRegion();
					mr.setStart(cntStart);
					mr.setEnd(cntEnd);
					mr.setSupportLRs(1);
					contig.addMisassemblyRegion(mr);
					logger.info(record.getqName());
					isTrue = true;
				}
			}
			return isTrue;
		}
	}

	@SuppressWarnings("unchecked")
	public static List<String> findMisassemblies(Parameter paras, Map<String, Sequence> seqs) {
		// TODO Auto-generated method stub
		long start = System.currentTimeMillis();
		List<String> misassemblies = null;
		double iqrTime = paras.getIqrTime();
		// sort the misassebmly regions;
		for(Map.Entry<String, Sequence> entry : seqs.entrySet())
		{
			Sequence cnt = entry.getValue();
			List<MisassemblyRegion> regions = cnt.getMisassemblies();
			if(!regions.isEmpty())
			{
//				Collections.sort(regions, new MisassemblyRegionSorter());
				Collections.sort(regions, new Comparator<Object>(){
					@Override
					public int compare(Object o1, Object o2) {
						MisassemblyRegion mr1 = (MisassemblyRegion) o1;
						MisassemblyRegion mr2 = (MisassemblyRegion) o2;
						if(mr1.getStart() > mr2.getStart())
							return 1;
						else if(mr1.getStart() == mr2.getStart())
							return 0;
						else 
							return -1;
					}
					
				});
				MisassemblyRegion target = regions.get(0);
				ArrayList<MisassemblyRegion> tempMRs = new ArrayList<MisassemblyRegion>();
				int size = regions.size();
				for(int i = 1; i < size; i++)
				{
					MisassemblyRegion query = regions.get(i);
					if(target.isOverlap(query))
					{
						target.updateRegion(query);
					} else
					{
						tempMRs.add(target);
						target = query;
					}
				}
				if(!tempMRs.contains(target))
					tempMRs.add(target);
//				cnt.setMisassemblyRegions(tempMRs);
				cnt.setMisassemblies(tempMRs);
			}
		}
		List<Integer> covs = new ArrayList<Integer>();
		for(Map.Entry<String, Sequence> entry : seqs.entrySet())
		{
			List<MisassemblyRegion> regions = entry.getValue().getMisassemblies();
			for(MisassemblyRegion mr : regions)
			{
				covs.add(mr.getSupportLRs());
			}
		}
		if (misassemblies == null)
			misassemblies = new ArrayList<String>();
		if (!covs.isEmpty()) {
			Map<String, Double> values = MathTool.summary(covs);
			double firstQ = values.get("FIRSTQ");
			double thirdQ = values.get("THIRDQ");
			double iqr = thirdQ - firstQ;
			double upper = iqr * iqrTime + values.get("THIRDQ");
			double median = values.get("MEDIAN");
			logger.info("Finding Misassemblies:");
			logger.info("MIN: " + values.get("MIN"));
			logger.info("First Quartile: " + firstQ);
			logger.info("Median = " + median);
			logger.info("Third Quartile: " + thirdQ);
			logger.info("MAX: " + values.get("MAX"));
			logger.info("Interquartile Range: " + iqr);
			logger.info(iqrTime + "'s IQR " + ", Outlier Threshold: " + upper);
			// logger.info("Pesudo repeat contigs");
			// Map<String, List<String>> repeat = new HashMap<String,
			for (Map.Entry<String, Sequence> entry : seqs.entrySet()) {
				Sequence cnt = entry.getValue();
				List<MisassemblyRegion> regions = cnt.getMisassemblies();
				for(MisassemblyRegion mr : regions)
				{
					if (mr.getSupportLRs() > upper) {
						if(!misassemblies.contains(entry.getKey()))
							misassemblies.add(entry.getKey());
						cnt.setIsMisassembly(true);
					}
				}
			}
			ContigMisassemblyWriter cmw = new ContigMisassemblyWriter(paras);
			cmw.write(seqs);
		}
		logger.info("Misassemblies count: " + misassemblies.size());
		long end = System.currentTimeMillis();
		logger.info("Finding misassemblies, erase time: " + (end - start) + " ms");
		return misassemblies;
	}
	
}


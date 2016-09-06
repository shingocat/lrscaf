/*
*File: agis.ps.util.MRecordValidator.java
*User: mqin
*Email: mqin@ymail.com
*Date: 2016年7月28日
*/
package agis.ps.util;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import agis.ps.link.M5Record;
import agis.ps.link.MRecord;

public class MRecordValidator {
	private static Logger logger = LoggerFactory.getLogger(MRecordValidator.class);

	public static MRecord validate(String [] arrs, Parameter paras) {
		M5Record m5 = null;
		try{ 
			int minPBLen = paras.getMinPBLen();
			int minCNTLen = paras.getMinContLen();
			double identity = paras.getIdentity();
			// if less than minimum pacbio length
			if (Integer.valueOf(arrs[1]) < minPBLen)
				return null;
			// if less tan minimum contig length
			if (Integer.valueOf(arrs[6]) < minCNTLen)
				return null;
			// if the identity less than specified value;
			double sum = Double.valueOf(arrs[11]) + Double.valueOf(arrs[12]) + Double.valueOf(arrs[13])
					+ Double.valueOf(arrs[14]);
			double value = Double.valueOf(arrs[11]) / sum;
			if (value < identity)
				return null;
			m5 = new M5Record();
			m5.setqName(arrs[0]);
			m5.setqLength(Integer.valueOf(arrs[1]));
			m5.setqStart(Integer.valueOf(arrs[2]));
			m5.setqEnd(Integer.valueOf(arrs[3]));
			m5.setqStrand(arrs[4].equals("+") ? Strand.FORWARD : Strand.REVERSE);
			m5.settName(arrs[5]);
			m5.settLength(Integer.valueOf(arrs[6]));
			m5.settStart(Integer.valueOf(arrs[7]));
			m5.settEnd(Integer.valueOf(arrs[8]));
			m5.settStrand(arrs[9].equals("+") ? Strand.FORWARD : Strand.REVERSE);
			m5.setScore(Integer.valueOf(arrs[10]));
			m5.setNumMatch(Integer.valueOf(arrs[11]));
			m5.setNumMismatch(Integer.valueOf(arrs[12]));
			m5.setNumIns(Integer.valueOf(arrs[13]));
			m5.setNumDel(Integer.valueOf(arrs[14]));
			m5.setMapQV(Integer.valueOf(arrs[15]));
			// m5.setqAlignedSeq(arrs[16]);
			// m5.setMatchPattern(arrs[17]);
			// m5.settAlignedSeq(arrs[18]);
			// setting identity
			m5.setIdentity(value);
		}catch(Exception e)
		{
			logger.error(MRecordValidator.class.getName() + "\t" + e.getMessage());
		}
		return m5;
	}
	
	public static MRecord validate4Repeats(String [] arrs, Parameter paras)
	{
		M5Record m5 = null;
		try{ 
			int minPBLen = paras.getMinPBLen();
			int minCNTLen = paras.getMinContLen();
			double identity = paras.getIdentity();
			int maxOHLen = paras.getMaxOHLen();
			double maxOHRatio = paras.getMaxOHRatio();

			
			int tLen = Integer.valueOf(arrs[6]);
			int tStart = Integer.valueOf(arrs[7]);
			int tEnd = Integer.valueOf(arrs[8]);
			int tLeftLen = tStart;
			int tRightLen = tLen - tEnd;
			int defOHLen = (int) (tLen * maxOHRatio);
			if(maxOHLen <= defOHLen)
				defOHLen = maxOHLen;
			if(tLeftLen > defOHLen )
				return null;
			if(tRightLen > defOHLen)
				return null;
			// for finding repeats, do not considering filter standard;
//			// if less than minimum pacbio length
//			if (Integer.valueOf(arrs[1]) < minPBLen)
//				return null;
//			// if less tan minimum contig length
//			if (Integer.valueOf(arrs[6]) < minCNTLen)
//				return null;
//			// if the identity less than specified value;
			double sum = Double.valueOf(arrs[11]) + Double.valueOf(arrs[12]) + Double.valueOf(arrs[13])
					+ Double.valueOf(arrs[14]);
			double value = Double.valueOf(arrs[11]) / sum;
//			if (value < identity)
//				return null;
			m5 = new M5Record();
			m5.setqName(arrs[0]);
			m5.setqLength(Integer.valueOf(arrs[1]));
			m5.setqStart(Integer.valueOf(arrs[2]));
			m5.setqEnd(Integer.valueOf(arrs[3]));
			m5.setqStrand(arrs[4].equals("+") ? Strand.FORWARD : Strand.REVERSE);
			m5.settName(arrs[5]);
			m5.settLength(tLen);
			m5.settStart(tStart);
			m5.settEnd(tEnd);
			m5.settStrand(arrs[9].equals("+") ? Strand.FORWARD : Strand.REVERSE);
			m5.setScore(Integer.valueOf(arrs[10]));
			m5.setNumMatch(Integer.valueOf(arrs[11]));
			m5.setNumMismatch(Integer.valueOf(arrs[12]));
			m5.setNumIns(Integer.valueOf(arrs[13]));
			m5.setNumDel(Integer.valueOf(arrs[14]));
			m5.setMapQV(Integer.valueOf(arrs[15]));
			// m5.setqAlignedSeq(arrs[16]);
			// m5.setMatchPattern(arrs[17]);
			// m5.settAlignedSeq(arrs[18]);
			// setting identity
			m5.setIdentity(value);
		}catch(Exception e)
		{
			logger.error(MRecordValidator.class.getName() + "\t" + e.getMessage());
		}
		return m5;
	}

}

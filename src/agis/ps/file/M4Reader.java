/*
*File: agis.ps.file.M4Reader.java
*User: mqin
*Email: mqin@ymail.com
*Date: 2016年5月30日
*/
package agis.ps.file;

//import org.slf4j.Logger;
//import org.slf4j.LoggerFactory;

import agis.ps.link.MRecord;
import agis.ps.util.Parameter;
import agis.ps.util.Strand;

public class M4Reader extends AlignmentFileReader{
//	private final static Logger logger = LoggerFactory.getLogger(M4Reader.class);
	
	public M4Reader(Parameter paras) {
		super(paras);
	}
	
	@Override
	protected MRecord initMRecord(String [] arrs)
	{
		MRecord record = new MRecord();
		record.setqName(arrs[0]);
		record.settName(arrs[1]);
//		record.setScore(Integer.valueOf(arrs[2]));
		double identity = Double.valueOf(arrs[3]) / 100;
		identity = Double.valueOf(String.format("%.4f", identity));
		record.setIdentity(identity);
		record.setqStrand(Strand.FORWARD);
		record.setqStart(arrs[5]);
		record.setqEnd(arrs[6]);
		record.setqLength(arrs[7]);
		if(arrs[8].equals("0"))
		{
			record.settStrand(Strand.FORWARD);
			record.settStart(arrs[9]);
			record.settEnd(arrs[10]);
			record.settLength(arrs[11]);
		} else
		{
			record.settStrand(Strand.REVERSE);
			int tStart = Integer.valueOf(arrs[9]);
			int tEnd = Integer.valueOf(arrs[10]);
			int tLen = Integer.valueOf(arrs[11]);
			record.settStart(tLen - tEnd);
			record.settEnd(tLen - tStart);
			record.settLength(tLen);
		}
//		record.setMapQV(Integer.valueOf(arrs[12]));
		return record;
	}
}



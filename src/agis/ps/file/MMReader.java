/*
*File: agis.ps.file.MMReader.java
*User: mqin
*Email: mqin@ymail.com
*Date: 2017年2月22日
*/
package agis.ps.file;

import agis.ps.link.MRecord;
import agis.ps.util.Parameter;
import agis.ps.util.Strand;

public class MMReader extends AlignmentFileReader {

	public MMReader(Parameter paras) {
		super(paras);
	}
	
	@Override
	protected MRecord initMRecord(String [] arrs)
	{
		MRecord record = new MRecord();
		record.setqName(arrs[0]);
		record.setqLength(arrs[1]);
		record.setqStart(arrs[2]);
		record.setqEnd(arrs[3]);
		record.settName(arrs[5]);
		record.settLength(arrs[6]);
		int match = Integer.valueOf(arrs[9]);
		int allMatch = Integer.valueOf(arrs[10]);
		double identity = Double.valueOf(((double)match) /allMatch);
		identity = Double.valueOf(String.format("%.4f", identity));
		record.setIdentity(identity);
		if(arrs[4].equals("+"))
			record.settStrand(Strand.FORWARD);
		else
			record.settStrand(Strand.REVERSE);
		record.settStart(arrs[7]);
		record.settEnd(arrs[8]);
		int cms = Integer.valueOf(arrs[12].split("cm:i:")[1]);
		if(cms <= paras.getMmcm())
			return null;
		return record;
	}
}



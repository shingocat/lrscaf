/*
*File: agis.ps.file.SamReader.java
*User: mqin
*Email: mqin@ymail.com
*Date: 2017年2月22日
*/
package agis.ps.file;

import java.util.regex.Matcher;
import java.util.regex.Pattern;

import agis.ps.link.MRecord;
import agis.ps.util.Parameter;
import agis.ps.util.Strand;

public class SamReader extends AlignmentFileReader {

	private String matchRegex = "\\d+M";
	private String insRegex = "\\d+I";
	private String delRegex = "\\d+D";
	
	public SamReader(Parameter paras) {
		super(paras);
	}
	
	@Override
	protected MRecord initMRecord(String [] arrs)
	{
		MRecord record = new MRecord();
		String qName = arrs[0];
		int flag = Integer.valueOf(arrs[1]);
		String tName = arrs[2];
		int tStart = Integer.valueOf(arrs[3]) - 1;
//		int mapQV = Integer.valueOf(arrs[4]);
		String cigar = arrs[5];
//		String rnext = arrs[6];
//		String pnext = arrs[7];
		int tOverlapLeng = Integer.valueOf(arrs[8]);
//		String seq = arrs[9];
//		String qual = arrs[10];
//		String rg = arrs[11];
//		String as = arrs[12]; //score;
		int xs = Integer.valueOf(arrs[13].split(":i:")[1]) - 1; // query start position
		int xe = Integer.valueOf(arrs[14].split(":i:")[1]) - 1; // query end position
//		int ys = Integer.valueOf(arrs[15].split(":i:")[1]); // query length start position;
		int ye = Integer.valueOf(arrs[16].split(":i:")[1]); // query length;
//		int zm = Integer.valueOf(arrs[17].split(":i:")[1]);
//		int xl = Integer.valueOf(arrs[18].split(":i:")[1]);
//		int xt = Integer.valueOf(arrs[19].split(":i:")[1]);
		int nm = Integer.valueOf(arrs[20].split(":i:")[1]);
//		int fi = Integer.valueOf(arrs[21].split(":i:")[1]);
//		int xq = Integer.valueOf(arrs[22].split(":i:")[1]);
		// set record
		record.setqName(qName);
		record.setqLength(ye);
		record.setqStart(xs);
		record.setqEnd(xe);
		record.setqStrand(Strand.FORWARD);
		record.settName(tName);
		record.settLength(0);
		record.settStart(tStart);
		record.settEnd(tStart + tOverlapLeng);
		record.settStrand((flag & 16) == 16 ? Strand.REVERSE : Strand.FORWARD);
//		m.setScore(Integer.valueOf(arrs[10]));
		record.setIdentity(this.getIdentity(cigar, nm));
		return record;
	}
	
	// compute identity from cigar seq
	private double getIdentity(String cigar, int sumDisMatchs)
	{
		double ident = 0.0;
		Pattern matPat = Pattern.compile(matchRegex, Pattern.CASE_INSENSITIVE);
		Pattern insPat = Pattern.compile(insRegex, Pattern.CASE_INSENSITIVE);
		Pattern delPat = Pattern.compile(delRegex, Pattern.CASE_INSENSITIVE);
		Integer matchs = 0;
		Integer ins = 0;
		Integer dels = 0;
		Matcher matMat = matPat.matcher(cigar);
		while(matMat.find())
		{
			String temp = matMat.group();
			temp = temp.replace("M", "");
			matchs += Integer.valueOf(temp);
		}
		Matcher insMat = insPat.matcher(cigar);
		while(insMat.find())
		{
			String temp = insMat.group();
			temp = temp.replace("I", "");
			ins += Integer.valueOf(temp);
		}
		Matcher delMat = delPat.matcher(cigar);
		while(delMat.find())
		{
			String temp = delMat.group();
			temp = temp.replace("D", "");
			dels += Integer.valueOf(temp);
		}
		
		int misMatchs = sumDisMatchs - ins - dels;
		ident = (double)(matchs - misMatchs) / (matchs + ins + dels);
		ident = Double.valueOf(String.format("%.4f", ident));
		return ident;
	}
}



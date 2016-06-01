/*
*File: agis.ps.M5Record.java
*User: mqin
*Email: mqin@ymail.com
*Date: 2015年12月28日
*/
package agis.ps.link;

import java.io.Serializable;
import java.text.DecimalFormat;

import agis.ps.util.Strand;

public class M5Record extends MRecord implements Serializable {
	
	private static final long serialVersionUID = 1L;	
	private Integer numMatch;
	private Integer numMismatch;
	private Integer numIns;
	private Integer numDel;
	private String qAlignedSeq;
	private String matchPattern;
	private String tAlignedSeq;
	
	public M5Record()
	{
		// do not thing
	}

	public Integer getNumMatch() {
		return numMatch;
	}

	public void setNumMatch(Integer numMatch) {
		this.numMatch = numMatch;
	}

	public Integer getNumMismatch() {
		return numMismatch;
	}

	public void setNumMismatch(Integer numMismatch) {
		this.numMismatch = numMismatch;
	}

	public Integer getNumIns() {
		return numIns;
	}

	public void setNumIns(Integer numIns) {
		this.numIns = numIns;
	}

	public Integer getNumDel() {
		return numDel;
	}

	public void setNumDel(Integer numDel) {
		this.numDel = numDel;
	}

	public String getqAlignedSeq() {
		return qAlignedSeq;
	}

	public void setqAlignedSeq(String qAlignedSeq) {
		this.qAlignedSeq = qAlignedSeq;
	}

	public String getMatchPattern() {
		return matchPattern;
	}

	public void setMatchPattern(String matchPattern) {
		this.matchPattern = matchPattern;
	}

	public String gettAlignedSeq() {
		return tAlignedSeq;
	}

	public void settAlignedSeq(String tAlignedSeq) {
		this.tAlignedSeq = tAlignedSeq;
	}

	public Double getIdentity()
	{
		int sum = this.getNumMatch() + this.getNumMismatch() + this.getNumIns() + this.getNumDel();
		double value = (double)this.getNumMatch() / sum;
//		DecimalFormat df = new DecimalFormat("0.00");
//		value = Double.valueOf(df.format(value));
		return value;
	}


	@Override
	public String toString()
	{
		return super.toString() + " M5Record [numMismatch=" + numMismatch
				+ ", numIns=" + numIns + ", numDel=" + numDel + " ]";
	}
}



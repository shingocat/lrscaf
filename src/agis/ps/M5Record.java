/*
*File: agis.ps.M5Record.java
*User: mqin
*Email: mqin@ymail.com
*Date: 2015年12月28日
*/
package agis.ps;

import java.io.Serializable;
import java.text.DecimalFormat;

import agis.ps.util.Strand;

public class M5Record implements Serializable {
	
	private static final long serialVersionUID = 1L;
	
	private String qName;
	private Integer qLength;
	private Integer qStart;
	private Integer qEnd;
	private Strand qStrand;
	private String tName;
	private Integer tLength;
	private Integer tStart;
	private Integer tEnd;
	private Strand tStrand;
	private Integer score;
	private Integer numMatch;
	private Integer numMismatch;
	private Integer numIns;
	private Integer numDel;
	private Integer mapQV;
	private String qAlignedSeq;
	private String matchPattern;
	private String tAlignedSeq;
	
	public M5Record()
	{
		// do not thing
	}

	public String getqName() {
		return qName;
	}

	public void setqName(String qName) {
		this.qName = qName;
	}

	public Integer getqLength() {
		return qLength;
	}

	public void setqLength(Integer qLength) {
		this.qLength = qLength;
	}

	public Integer getqStart() {
		return qStart;
	}

	public void setqStart(Integer qStart) {
		this.qStart = qStart;
	}

	public Integer getqEnd() {
		return qEnd;
	}

	public void setqEnd(Integer qEnd) {
		this.qEnd = qEnd;
	}

	public Strand getqStrand() {
		return qStrand;
	}

	public void setqStrand(Strand qStrand) {
		this.qStrand = qStrand;
	}

	public String gettName() {
		return tName;
	}

	public void settName(String tName) {
		this.tName = tName;
	}

	public Integer gettLength() {
		return tLength;
	}

	public void settLength(Integer tLength) {
		this.tLength = tLength;
	}

	public Integer gettStart() {
		return tStart;
	}

	public void settStart(Integer tStart) {
		this.tStart = tStart;
	}

	public Integer gettEnd() {
		return tEnd;
	}

	public void settEnd(Integer tEnd) {
		this.tEnd = tEnd;
	}

	public Strand gettStrand() {
		return tStrand;
	}

	public void settStrand(Strand tStrand) {
		this.tStrand = tStrand;
	}

	public Integer getScore() {
		return score;
	}

	public void setScore(Integer score) {
		this.score = score;
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

	public Integer getMapQV() {
		return mapQV;
	}

	public void setMapQV(Integer mapQV) {
		this.mapQV = mapQV;
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

	/*@Override
	public String toString() {
		return "M5Record [qName=" + qName + ", qLength=" + qLength + ", qStart=" + qStart + ", qEnd=" + qEnd
				+ ", qStrand=" + qStrand + ", tName=" + tName + ", tLength=" + tLength + ", tStart=" + tStart
				+ ", tEnd=" + tEnd + ", tStrand=" + tStrand + ", score=" + score + ", numMismatch=" + numMismatch
				+ ", numIns=" + numIns + ", numDel=" + numDel + ", mapQV=" + mapQV + ", qAlignedSeq=" + qAlignedSeq
				+ ", matchPattern=" + matchPattern + ", tAlignedSeq=" + tAlignedSeq + "]";
	}*/
	
	@Override
	public String toString()
	{
		return "M5Record [qName=" + qName + ", qLength=" + qLength + ", qStart=" + qStart + ", qEnd=" + qEnd
				+ ", qStrand=" + qStrand + ", tName=" + tName + ", tLength=" + tLength + ", tStart=" + tStart
				+ ", tEnd=" + tEnd + ", tStrand=" + tStrand + ", score=" + score + ", numMismatch=" + numMismatch
				+ ", numIns=" + numIns + ", numDel=" + numDel + ", mapQV=" + mapQV + "]";
	}
	
}



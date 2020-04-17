/*
*File: agis.ps.link.MRecord.java
*User: mqin
*Email: mqin@ymail.com
*Date: 2016年5月30日
*/
package agis.ps.link;

import java.io.Serializable;

import agis.ps.util.Strand;

public class MRecord implements Serializable{

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
//	private Integer score;
	private Double identity;
//	private Integer mapQV;
	
	public String getqName() {
		return qName;
	}
	public void setqName(String qName) {
		this.qName = qName;
	}
	public Integer getqLength() {
		return qLength;
	}
	public void setqLength(int qLength)
	{
		this.qLength = qLength;
	}
	public void setqLength(String qLength)
	{
		this.qLength = Integer.valueOf(qLength);
	}
	public Integer getqStart() {
		return qStart;
	}
	public void setqStart(int qStart)
	{
		this.qStart = qStart;
	}
	public void setqStart(String qStart) {
		this.qStart = Integer.valueOf(qStart);
	}
	public Integer getqEnd() {
		return qEnd;
	}
	public void setqEnd(int qEnd)
	{
		this.qEnd = qEnd;
	}
	public void setqEnd(String qEnd) {
		this.qEnd = Integer.valueOf(qEnd);
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
	public void settLength(String tLength) {
		this.tLength = Integer.valueOf(tLength);
	}
	public void settLength(Integer tLength)
	{
		this.tLength = tLength;
	}
	public Integer gettStart() {
		return tStart;
	}
	public void settStart(String tStart) {
		this.tStart = Integer.valueOf(tStart);
	}
	public void settStart(Integer tStart)
	{
		this.tStart = tStart;
	}
	public Integer gettEnd() {
		return tEnd;
	}
	public void settEnd(String tEnd) {
		this.tEnd = Integer.valueOf(tEnd);
	}
	public void settEnd(Integer tEnd)
	{
		this.tEnd = tEnd;
	}
	public Strand gettStrand() {
		return tStrand;
	}
	public void settStrand(Strand tStrand) {
		this.tStrand = tStrand;
	}
//	public Integer getScore() {
//		return score;
//	}
//	public void setScore(Integer score) {
//		this.score = score;
//	}
	public Double getIdentity() {
		return identity;
	}
	public void setIdentity(Double identity) {
		this.identity = identity;
	}
//	public Integer getMapQV() {
//		return mapQV;
//	}
//	public void setMapQV(Integer mapQV) {
//		this.mapQV = mapQV;
//	}
	
	
//	@Override
//	public String toString() {
//		return "MRecord [qName=" + qName + ", qLength=" + qLength + ", qStart=" + qStart + ", qEnd=" + qEnd
//				+ ", qStrand=" + qStrand + ", tName=" + tName + ", tLength=" + tLength + ", tStart=" + tStart
//				+ ", tEnd=" + tEnd + ", tStrand=" + tStrand + ", score=" + score + ", identity=" + identity + ", mapQV="
//				+ mapQV + "]";
//	}
	
	
	
	@Override
	public int hashCode() {
		final int prime = 31;
		int result = 1;
		result = prime * result + ((qName == null) ? 0 : qName.hashCode());
		result = prime * result + ((tName == null) ? 0 : tName.hashCode());
		return result;
	}
	
	@Override
	public String toString() {
		return "MRecord [qName=" + qName + ", qLength=" + qLength + ", qStart=" + qStart + ", qEnd=" + qEnd
				+ ", qStrand=" + qStrand + ", tName=" + tName + ", tLength=" + tLength + ", tStart=" + tStart
				+ ", tEnd=" + tEnd + ", tStrand=" + tStrand + ", identity=" + identity + "]";
	}
	@Override
	public boolean equals(Object obj) {
		if (this == obj)
			return true;
		if (obj == null)
			return false;
		if (getClass() != obj.getClass())
			return false;
		MRecord other = (MRecord) obj;
		if (qName == null) {
			if (other.qName != null)
				return false;
		} else if (!qName.equals(other.qName))
			return false;
		if (tName == null) {
			if (other.tName != null)
				return false;
		} else if (!tName.equals(other.tName))
			return false;
		return true;
	}
	
	
}



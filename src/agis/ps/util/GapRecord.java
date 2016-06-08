/*
*File: agis.ps.util.GapRecord.java
*User: mqin
*Email: mqin@ymail.com
*Date: 2016年6月2日
*/
package agis.ps.util;

import java.io.Serializable;
import java.util.List;
import java.util.Vector;

import agis.ps.seqs.PBGapSeq;

public class GapRecord implements Serializable {

	private static final long serialVersionUID = 1L;
	private String start;
	private String end;
	private List<PBGapSeq> seqs = new Vector<PBGapSeq>(10);
	private int nums;
	
	public String getStart() {
		return start;
	}
	public void setStart(String start) {
		this.start = start;
	}
	public String getEnd() {
		return end;
	}
	public void setEnd(String end) {
		this.end = end;
	}
	public List<PBGapSeq> getSeqs() {
		return seqs;
	}
	public void setSeqs(List<PBGapSeq> seqs) {
		this.seqs = seqs;
	}
	public int getNums() {
		return seqs.size();
	}
//	public void setNums(int nums) {
//		this.nums = nums;
//	}
	public void addSeq(PBGapSeq seq)
	{
		if(seqs == null)
			seqs = new Vector<PBGapSeq>();
		seqs.add(seq);
	}
	
	@Override
	public int hashCode() {
		final int prime = 31;
		int result = 1;
		result = prime * result + ((end == null) ? 0 : end.hashCode());
		result = prime * result + ((start == null) ? 0 : start.hashCode());
		return result;
	}
	
	@Override
	public boolean equals(Object obj) {
		if (this == obj)
			return true;
		if (obj == null)
			return false;
		if (getClass() != obj.getClass())
			return false;
		GapRecord other = (GapRecord) obj;
		if(start.equalsIgnoreCase(other.getStart()) && end.equalsIgnoreCase(other.getEnd()))
		{
			return true;
		} else if(start.equalsIgnoreCase(other.getEnd()) && end.equalsIgnoreCase(other.getStart()))
			return true;
		else 
			return false;
//		if (end == null) {
//			if (other.end != null)
//				return false;
//		} else if (!end.equals(other.end))
//			return false;
//		if (start == null) {
//			if (other.start != null)
//				return false;
//		} else if (!start.equals(other.start))
//			return false;
//		return true;
	}
	
	@Override
	public String toString() {
		return "GapRecord [start=" + start + ", end=" + end + ", nums=" + nums + "]";
	}
}



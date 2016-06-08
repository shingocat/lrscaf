/*
*File: agis.ps.seqs.PBGapSeq.java
*User: mqin
*Email: mqin@ymail.com
*Date: 2016年6月2日
*/
package agis.ps.seqs;

import java.io.Serializable;

import agis.ps.util.Strand;

public class PBGapSeq implements Serializable {

	private static final long serialVersionUID = 1L;
	private String id;
	private int start;
	private int end;
	private Strand strand;
	
	public String getId() {
		return id;
	}
	public void setId(String id) {
		this.id = id;
	}
	public int getStart() {
		return start;
	}
	public void setStart(int start) {
		this.start = start;
	}
	public int getEnd() {
		return end;
	}
	public void setEnd(int end) {
		if(end < start)
			throw new IllegalArgumentException(this.getClass().getName() + "\t:The end position could not be less than start!");
		this.end = end;
	}
	public Strand getStrand() {
		return strand;
	}
	public void setStrand(Strand strand) {
		this.strand = strand;
	}
	@Override
	public String toString() {
		return "PBGapSeq [id=" + id + ", start=" + start + ", end=" + end + ", strand=" + strand + "]";
	}
}



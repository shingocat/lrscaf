package agis.ps.seqs;


/**
* @author Mao Qin
* @version 2020年9月16日 上午10:18:58
* @Email mqin@outlook.com
* @Description class usage.
* @Copyright All Right Reserved 2020.
*/
public class SeqRegion {
	private Integer index;
	private Integer start;
	private Integer end;
	private Integer length;
	private SequenceType seqType = SequenceType.NORMAL;
	
	public Integer getLength() {
		if(this.end >= this.start)
			return this.end - this.start;
		else 
			return 0;
	}
	
	public SequenceType getSeqType() {
		return seqType;
	}
	public void setSeqType(SequenceType seqType) {
		this.seqType = seqType;
	}
	public int getIndex() {
		return index;
	}
	public void setIndex(int index) {
		this.index = index;
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
		this.end = end;
	}
	@Override
	public String toString() {
		return "SeqRegion [index=" + index + ", start=" + start + ", end=" + end + ", length= " + this.getLength() + ", seqType=" + seqType + "]";
	}
}

/** 
** Usage: TODO
** Author: mqin
** Email: mqin@outlook.com
** Date: 2017年7月14日
*/
package agis.ps.util;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import agis.ps.seqs.Contig;

public class MisassemblyRegion {
	private static Logger logger = LoggerFactory.getLogger(MisassemblyRegion.class);
//	private Contig contig; // belong to which contig;
	private int start; // misassembled region start, always forward strand position;
	private int end; // misassembled region end, as start position;
	private int supportLRs; // the number of long reads supported this region which is misassembled;
	
	public MisassemblyRegion()
	{
		
	}
	
//	public Contig getContig() {
//		return contig;
//	}
//
//	public void setContig(Contig contig) {
//		this.contig = contig;
//	}

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

	public int getSupportLRs() {
		return supportLRs;
	}

	public void setSupportLRs(int supportLRs) {
		this.supportLRs = supportLRs;
	}
	
	public void updateRegion(int start, int end)
	{
		if(start < this.start)
			this.start = start;
		if(end > this.end)
			this.end = end;
		this.supportLRs = this.supportLRs + 1;
	}
	
	public void updateRegion(MisassemblyRegion mr)
	{
		int start = mr.getStart();
		int end = mr.getEnd();
		this.updateRegion(start, end);
		
	}
	
	public boolean isOverlap(int start, int end)
	{
		boolean isTrue = false;
		if(start <= this.start && end >= this.start)
		{
			isTrue = true;
		} else if(start > this.start && start <= this.end)
		{
			isTrue = true;
		}
		return isTrue;
	}
	
	public boolean isOverlap(MisassemblyRegion mr)
	{
		int start = mr.getStart();
		int end = mr.getEnd();
		return this.isOverlap(start, end);
	}

	@Override
	public String toString() {
		return "MisassemblyRegion [start=" + start + 
				", end=" + end + ", supportLRs=" + supportLRs + "]";
	}

	@Override
	public int hashCode() {
		final int prime = 31;
		int result = 1;
		result = prime * result + end;
		result = prime * result + start;
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
		MisassemblyRegion other = (MisassemblyRegion) obj;
		if (end != other.end)
			return false;
		if (start != other.start)
			return false;
		return true;
	}
	
}

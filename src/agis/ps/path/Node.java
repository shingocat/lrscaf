/*
*File: agis.ps.path.Node.java
*User: mqin
*Email: mqin@ymail.com
*Date: 2016年2月26日
*/
package agis.ps.path;

import java.io.Serializable;

import agis.ps.link.Contig;
import agis.ps.util.Strand;

public class Node implements Serializable {
	
	private static final long serialVersionUID = 1L;
	// the contig 
	private Contig cnt;
	// the contig strand in the path;
	// depend on the first start point
	private Strand strand;
	// the mean distance to the next vertex;
	// mean of two opposite edges if exist;
	private int meanDist2Next;
	// the sd distance to the next vertext;
	// mean of two opposite edges if exist;
	private int sdDist2Next;
	// the support link number;
	// sum of two opposite edges if exist;
	private int supportLinkNum;
	// indicated whether orphan node
	boolean isOrphan;
	
	public Contig getCnt() {
		return cnt;
	}
	public void setCnt(Contig cnt) {
		this.cnt = cnt;
	}
	public Strand getStrand() {
		return strand;
	}
	public void setStrand(Strand strand) {
		this.strand = strand;
	}
	public int getMeanDist2Next() {
		return meanDist2Next;
	}
	public void setMeanDist2Next(int meanDist2Next) {
		this.meanDist2Next = meanDist2Next;
	}
	public int getSdDist2Next() {
		return sdDist2Next;
	}
	public void setSdDist2Next(int sdDist2Next) {
		this.sdDist2Next = sdDist2Next;
	}
	public int getSupportLinkNum() {
		return supportLinkNum;
	}
	public void setSupportLinkNum(int supportLinkNum) {
		this.supportLinkNum = supportLinkNum;
	}
	public boolean isOrphan() {
		return isOrphan;
	}
	public void setOrphan(boolean isOrphan) {
		this.isOrphan = isOrphan;
	}
	
	@Override
	public String toString() {
		return "Node [cnt=" + cnt + ", strand=" + strand + ", meanDist2Next=" + meanDist2Next + ", sdDist2Next="
				+ sdDist2Next + ", supportLinkNum=" + supportLinkNum + ", isOrphan=" + isOrphan + "]";
	}
}



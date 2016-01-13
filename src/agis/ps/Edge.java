/*
*File: agis.ps.Edge.java
*User: mqin
*Email: mqin@ymail.com
*Date: 2015年12月28日
*/
package agis.ps;

import agis.ps.link.Contig;

public class Edge {
	private Contig origin;
	private Contig terminus;
	private int linkNum;

	public Contig getOrigin() {
		return origin;
	}

	public void setOrigin(Contig origin) {
		this.origin = origin;
	}

	public Contig getTerminus() {
		return terminus;
	}

	public void setTerminus(Contig terminus) {
		this.terminus = terminus;
	}

	public int getLinkNum() {
		return linkNum;
	}

	public void setLinkNum(int linkNum) {
		this.linkNum = linkNum;
	}

}

/*
*File: agis.ps.Edge.java
*User: mqin
*Email: mqin@ymail.com
*Date: 2015年12月28日
*/
package agis.ps;

import java.io.Serializable;

import agis.ps.link.Contig;
import agis.ps.util.Strand;

public class Edge implements Serializable{
	private Contig origin;
	private Contig terminus;
	private Strand oStrand;
	private Strand tStrand;
	private int linkNum;
	private int distMean;
	private int distSd;

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

	public Strand getoStrand() {
		return oStrand;
	}

	public void setoStrand(Strand oStrand) {
		this.oStrand = oStrand;
	}

	public Strand gettStrand() {
		return tStrand;
	}

	public void settStrand(Strand tStrnad) {
		this.tStrand = tStrnad;
	}

	public int getDistMean() {
		return distMean;
	}

	public void setDistMean(int distMean) {
		this.distMean = distMean;
	}

	public int getDistSd() {
		return distSd;
	}

	public void setDistSd(int distSd) {
		this.distSd = distSd;
	}

	@Override
	public String toString() {
		return "Edge [origin=" + origin.getID() + ", terminus=" + terminus.getID() + ", oStrand=" + oStrand + ", tStrand=" + tStrand
				+ ", linkNum=" + linkNum + ", distMean=" + distMean + ", distSd=" + distSd + "]";
	}
	
	

}

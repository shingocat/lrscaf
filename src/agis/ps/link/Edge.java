/*
*File: agis.ps.Edge.java
*User: mqin
*Email: mqin@ymail.com
*Date: 2015年12月28日
*/
package agis.ps.link;

import java.io.Serializable;

import agis.ps.seqs.Contig;
import agis.ps.util.Strand;

public class Edge implements Serializable{
	private static final long serialVersionUID = 1L;
	private Contig origin;
	private Contig terminus;
	private Strand oStrand;
	private Strand tStrand;
	private int linkNum;
	private int distMean;
	private int distSd;
	private boolean isOL; // overlap or not;
	private boolean isValid; // setting this edges is valid for build path;
	private boolean isFake; // setting this edge is fake or pesudo by program to travel graph!
	
	public boolean isOL() {
		return isOL;
	}

	public void setOL(boolean isOL) {
		this.isOL = isOL;
	}

	public boolean isValid() {
		return isValid;
	}

	public void setValid(boolean isValid) {
		this.isValid = isValid;
	}

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
	
	public boolean isFake() {
		return isFake;
	}

	public void setFake(boolean isFake) {
		this.isFake = isFake;
	}

	@Override
	public int hashCode() {
		final int prime = 31;
		int result = 1;
		result = prime * result + distMean;
		result = prime * result + distSd;
		result = prime * result + (isOL ? 1231 : 1237);
		result = prime * result + (isValid ? 1231 : 1237);
		result = prime * result + linkNum;
		result = prime * result + ((oStrand == null) ? 0 : oStrand.hashCode());
		result = prime * result + ((origin == null) ? 0 : origin.hashCode());
		result = prime * result + ((tStrand == null) ? 0 : tStrand.hashCode());
		result = prime * result + ((terminus == null) ? 0 : terminus.hashCode());
		return result;
	}
	
	// if the origin and terminu contig is the same, assuming that this edge is the same;
	@Override
	public boolean equals(Object obj) {
		if (this == obj)
			return true;
		if (obj == null)
			return false;
		if (getClass() != obj.getClass())
			return false;
		Edge other = (Edge) obj;
//		if (distMean != other.distMean)
//			return false;
//		if (distSd != other.distSd)
//			return false;
//		if (isOL != other.isOL)
//			return false;
//		if (isValid != other.isValid)
//			return false;
//		if (linkNum != other.linkNum)
//			return false;
//		if (oStrand != other.oStrand)
//			return false;
		if (origin == null) {
			if (other.origin != null)
				return false;
		} else if (!origin.equals(other.origin))
			return false;
//		if (tStrand != other.tStrand)
//			return false;
		if (terminus == null) {
			if (other.terminus != null)
				return false;
		} else if (!terminus.equals(other.terminus))
			return false;
		return true;
	}

	@Override
	public String toString() {
		return "Edge [origin=" + origin + ", terminus=" + terminus + ", oStrand=" + oStrand + ", tStrand=" + tStrand
				+ ", linkNum=" + linkNum + ", distMean=" + distMean + ", distSd=" + distSd + ", isOL=" + isOL
				+ ", isValid=" + isValid + ", isFake=" + isFake + "]";
	}
}

/*
*File: agis.ps.link.PBLinkM5.java
*User: mqin
*Email: mqin@ymail.com
*Date: 2016年1月23日
*/
package agis.ps.link;

import agis.ps.M5Record;
import agis.ps.util.Strand;

public class PBLinkM5 implements ILink {
	private String id; // pacbiio read id;
	private M5Record origin; // The origin point in the link;
	private M5Record terminus; // the terminus point in the link;
	private Boolean isOverLap; // according to the distance between two contigs, if dist is minus as overlap, else having gap;
	
	// according to the dist value;
	public Boolean isOverLap()
	{
		isOverLap = false;
		if(getDistance() < 0)
			isOverLap = true;
		return isOverLap;
	}
	
	// distance = terminus_pacbio_start - terminus_contig_start -
	//				(origin_pacbio_end + origin_contig_length - origin_contig_end)
	// distance should be according to origin and terminus strand;
	@Override
	public int getDistance() {
		int dist = 0;
		int oPBS = origin.getqStart(); // origin pacbio start point;
		int oCntS = origin.gettStart(); // origin contig start point;
		int oPBE = origin.getqEnd(); // origin pacbio end point;
		int oCntE = origin.gettEnd(); // origin contig end point;
		int tPBS = terminus.getqStart(); // terminus pacbio start point;
		int tCntS = terminus.gettStart(); // terminus contig start point;
		int tPBE = terminus.getqEnd(); // terminus pacbio end point;
		int tCntE = terminus.gettEnd(); // terminus contig end point;
		int oCntLen = origin.gettLength();// origin contig length;
		int tCntLen = terminus.gettLength(); // terminus contig length
		if(origin.gettStrand().equals(Strand.FORWARD) && terminus.gettStrand().equals(Strand.FORWARD))
		{ // + +
			int oRightLen = oCntLen - oCntE;
			int tLeftLen = tCntS;
			dist = tPBS - oPBE - oRightLen - tLeftLen;
		} else if(origin.gettStrand().equals(Strand.REVERSE) && terminus.gettStrand().equals(Strand.REVERSE))
		{ // - -
			int oLeftLen = oCntS;
			int tRightLen = tCntLen -tCntE;
			dist = tPBS - oPBE - oLeftLen - tRightLen;
		} else if(origin.gettStrand().equals(Strand.FORWARD) && terminus.gettStrand().equals(Strand.REVERSE))
		{ // + -
			int oRightLen = oCntLen - oCntE;
			int tRightLen = tCntLen -tCntE;
			dist = tPBS - oPBE - oRightLen - tRightLen;
		} else if(origin.gettStrand().equals(Strand.REVERSE) && terminus.gettStrand().equals(Strand.FORWARD))
		{ // - +
			int oLeftLen = oCntS;
			int tLeftLen = tCntS;
			dist = tPBS - oPBE - oLeftLen - tLeftLen;
		}
//		int oPBEnd = origin.getqEnd();
//		int oCtLen = origin.gettLength();
//		int oCtRightLen = oCtLen - origin.gettEnd();
//		int tPBStart = terminus.getqStart();
//		int tCtLen = terminus.gettLength();
//		int tCtLeftLen = terminus.gettStart();
//		dist = tPBStart - tCtLeftLen - (oPBEnd + oCtRightLen);
		return dist;
	}


	public String getId() {
		return id;
	}


	public void setId(String id) {
		this.id = id;
	}


	public M5Record getOrigin() {
		return origin;
	}


	public void setOrigin(M5Record origin) {
		this.origin = origin;
	}


	public M5Record getTerminus() {
		return terminus;
	}


	public void setTerminus(M5Record terminus) {
		this.terminus = terminus;
	}
	
	public boolean isSelfLink()
	{
		boolean isSelfLink = false;
		if(this.origin.gettName().equalsIgnoreCase(this.terminus.gettName()))
			isSelfLink = true;
		return isSelfLink;
	}

	@Override
	public String toString() {
		return "PBLinkM5 [id=" + id + ", origin=" + origin + ", terminus=" + terminus + ", isOverLap=" + isOverLap
				+ "]";
	}
}



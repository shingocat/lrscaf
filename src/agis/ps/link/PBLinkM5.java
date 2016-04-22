/*
*File: agis.ps.link.PBLinkM5.java
*User: mqin
*Email: mqin@ymail.com
*Date: 2016年1月23日
*/
package agis.ps.link;

import agis.ps.M5Record;

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
	@Override
	public int getDistance() {
		int oPBEnd = origin.getqEnd();
		int oCtLen = origin.gettLength();
		int oCtRightLen = oCtLen - origin.gettEnd();
		int tPBStart = terminus.getqStart();
		int tCtLen = terminus.gettLength();
		int tCtLeftLen = terminus.gettStart();
		int dist = tPBStart - tCtLeftLen - (oPBEnd + oCtRightLen);
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
	
}



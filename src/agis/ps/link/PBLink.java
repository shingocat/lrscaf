/*
*File: agis.ps.link.PBLink.java
*User: mqin
*Email: mqin@ymail.com
*Date: 2016年1月12日
*/
package agis.ps.link;

import java.io.Serializable;

import agis.ps.seqs.Contig;
import agis.ps.util.Strand;

// PacBio link model 

public class PBLink implements ILink {
	private String origin;
	private Strand oStrand;
	private String terminus;
	private Strand tStrand;
	private int dist;
	private String pbId;
	private int oPStart;
	private int oPEnd;
	private int tPStart;
	private int tPEnd;

	public PBLink() {

	}
	
	public void setOPStart(int start)
	{
		this.oPStart = start;
	}
	
	public int getOPStart()
	{
		return this.oPStart;
	}
	
	public void setOPEnd(int end)
	{
		this.oPEnd = end;
	}
	
	public int getOPEnd()
	{
		return this.oPEnd;
	}
	
	public void setTPStart(int start)
	{
		this.tPStart = start;
	}
	
	public int getTPStart()
	{
		return this.oPStart;
	}
	
	public void setTPEnd(int end)
	{
		this.tPEnd = end;
	}
	
	public int getTPEnd()
	{
		return this.tPEnd;
	} 
	
	public void setOStrand(String mark)
	{
		if(mark.equalsIgnoreCase("+"))
			oStrand = Strand.FORWARD;
		else
			oStrand = Strand.REVERSE;
	}
	
	public Strand getOStrand()
	{
		return oStrand;
	}
	
	public void setTStrand(String mark)
	{
		if(mark.equalsIgnoreCase("+"))
			tStrand = Strand.FORWARD;
		else
			tStrand = Strand.REVERSE;
	}
	
	public Strand getTStrand()
	{
		return tStrand;
	}

	public String getOrigin() {
		return origin;
	}

	public void setOrigin(String origin) {
		this.origin = origin;
	}

	public String getTerminus() {
		return terminus;
	}

	public void setTerminus(String terminus) {
		this.terminus = terminus;
	}

	public int getDist() {
		return dist;
	}

	public void setDist(int dist) {
		this.dist = dist;
	}

	public String getPbId() {
		return pbId;
	}

	public void setPbId(String pbId) {
		this.pbId = pbId;
	}

	@Override
	public int getDistance() {
		return dist;
	}

}

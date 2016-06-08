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

public class PBLink implements ILink, Serializable {
	private String ID; // it is the pacbio read id for each link identify;
	private Contig origin; // the origin contig in the directed graph;
	private Contig terminus; // the terminus contig in the directed graph;
	private Strand oStrand; // the strand of the origin contig aligned to PacBio
							// read
	private Strand tStrand; // the strand of the terminus contig aligned to
							// PacBio read
	private int oStartLoc; // the start location of origin contig over pacbio
							// read;
	private int oEndLoc; // the end location of origin contig over pacbio read;
	private int tStartLoc; // the end location of terminus over pacbio read;
	private int tEndLoc; // the end location of terminus over pacbio read;
	private boolean isValid; // whether this link fit for the

	@Override
	public int getDistance() {
		return oEndLoc - tStartLoc;
	}

	public String getID() {
		return ID;
	}

	public void setID(String iD) {
		ID = iD;
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

	public Strand getoStrand() {
		return oStrand;
	}

	public void setoStrand(Strand oStrand) {
		this.oStrand = oStrand;
	}

	public Strand gettStrand() {
		return tStrand;
	}

	public void settStrand(Strand tStrand) {
		this.tStrand = tStrand;
	}

	public int getoStartLoc() {
		return oStartLoc;
	}

	public void setoStartLoc(int oStartLoc) {
		this.oStartLoc = oStartLoc;
	}

	public int getoEndLoc() {
		return oEndLoc;
	}

	public void setoEndLoc(int oEndLoc) {
		this.oEndLoc = oEndLoc;
	}

	public int gettStartLoc() {
		return tStartLoc;
	}

	public void settStartLoc(int tStartLoc) {
		this.tStartLoc = tStartLoc;
	}

	public int gettEndLoc() {
		return tEndLoc;
	}

	public void settEndLoc(int tEndLoc) {
		this.tEndLoc = tEndLoc;
	}

	public boolean isValid() {
		return isValid;
	}

	public void setValid(boolean isValid) {
		this.isValid = isValid;
	}

}

package agis.ps.link;

import agis.ps.seqs.Contig;
import agis.ps.seqs.Sequence;
import agis.ps.util.Strand;

public class Link implements ILink {
//	private Contig original; // original contig of link
//	private Contig terminus; // terminus contig of link
	private Sequence original; // original sequence of link;
	private Sequence terminus; // terminus sequence of link;
	private Strand oStrand; // the strand of original contig
	private Strand tStrand; // the strand of terminus contig
	private int dist; // distance of link
	private String lrId; // long read id which connect two contigs;
	private int oStart; // start position of original contig on long read
	private int oEnd; // end position of original lcontig on long read
	private int tStart; // start position of terminus contig on long read
	private int tEnd; // end position of terminus contig on long read

	public Link() {

	}

	public Sequence getOriginal() {
		return original;
	}

	public void setOriginal(Sequence original) {
		this.original = original;
	}

	public Sequence getTerminus() {
		return terminus;
	}

	public void setTerminus(Sequence terminus) {
		this.terminus = terminus;
	}

	public Strand getoStrand() {
		return oStrand;
	}

	public void setoStrand(Strand oStrand) {
		this.oStrand = oStrand;
	}

	public void setoStrand(String mark) {
		if (mark.equalsIgnoreCase("+"))
			oStrand = Strand.FORWARD;
		else
			oStrand = Strand.REVERSE;
	}

	public Strand gettStrand() {
		return tStrand;
	}

	public void settStrand(Strand tStrand) {
		this.tStrand = tStrand;
	}

	public void settStrand(String mark) {
		if (mark.equalsIgnoreCase("+"))
			tStrand = Strand.FORWARD;
		else
			tStrand = Strand.REVERSE;
	}

	public int getDist() {
		return dist;
	}

	public void setDist(int dist) {
		this.dist = dist;
	}

	public String getLrId() {
		return lrId;
	}

	public void setLrId(String lrId) {
		this.lrId = lrId;
	}

	public int getoStart() {
		return oStart;
	}

	public void setoStart(int oStart) {
		this.oStart = oStart;
	}

	public int getoEnd() {
		return oEnd;
	}

	public void setoEnd(int oEnd) {
		this.oEnd = oEnd;
	}

	public int gettStart() {
		return tStart;
	}

	public void settStart(int tStart) {
		this.tStart = tStart;
	}

	public int gettEnd() {
		return tEnd;
	}

	public void settEnd(int tEnd) {
		this.tEnd = tEnd;
	}

	@Override
	public int getDistance() {
		// TODO Auto-generated method stub
		return 0;
	}

	@Override
	public String toString() {
		return "Link [orignal=" + original.getId() + ", terminus=" + terminus.getId() + ", oStrand=" + oStrand
				+ ", tStrand=" + tStrand + ", dist=" + dist + ", lrId=" + lrId + ", oStart=" + oStart + ", oEnd=" + oEnd
				+ ", tStart=" + tStart + ", tEnd=" + tEnd + "]";
	}

}

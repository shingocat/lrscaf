/*
*File: agis.ps.link.Contig.java
*User: mqin
*Email: mqin@ymail.com
*Date: 2016年1月12日
*/
package agis.ps.seqs;

import org.biojava.nbio.core.exceptions.CompoundNotFoundException;
import org.biojava.nbio.core.sequence.DNASequence;
import org.biojava.nbio.core.sequence.compound.NucleotideCompound;
import org.biojava.nbio.core.sequence.template.SequenceView;

public class Contig extends DNASequence {
	private String ID; // the Id of this contig;
	private int length; // the length of this contigs;

	public Contig() {
		this.setDNAType(DNAType.CHROMOSOME);
	}

	public Contig(String seq) throws CompoundNotFoundException {
		super(seq);
	}

	public String getID() {
		return ID;
	}

	public void setID(String iD) {
		ID = iD;
	}

	public int getLength() {
//		if(getSequenceAsString() == null)
//			return 0;
//		int seqLen = getSequenceAsString().length();
//		if (seqLen != length)
//			return seqLen;
//		else
			return length;
	}

	public void setLength(int length) {
		this.length = length;
	}

	@Override
	public boolean equals(Object o) {
		Contig c = (Contig) o;
		if (c.getID().equals(this.getID()))
			return true;
		return false;
	}

	public String getReverseComplementSeq() {
		StringBuffer sb = new StringBuffer();
		String originalSeq = this.getSequenceAsString();
		for (int j = originalSeq.length() - 1; j >= 0; j--) {
			switch (originalSeq.charAt(j)) {
			case 'A':
				sb.append("T");
				break;
			case 'T':
				sb.append("A");
				break;
			case 'C':
				sb.append("G");
				break;
			case 'G':
				sb.append("C");
				break;
			case 'a':
				sb.append("t");
				break;
			case 't':
				sb.append("a");
				break;
			case 'c':
				sb.append("g");
				break;
			case 'g':
				sb.append("c");
				break;
			}
		}
		return sb.toString();
	}

	@Override
	public String toString() {
		return "Contig [ID=" + ID + ", length=" + length + "]";
	}
}

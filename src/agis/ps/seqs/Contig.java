/*
*File: agis.ps.link.Contig.java
*User: mqin
*Email: mqin@ymail.com
*Date: 2016年1月12日
*/
package agis.ps.seqs;

import org.biojava.nbio.core.exceptions.CompoundNotFoundException;
import org.biojava.nbio.core.sequence.AccessionID;
import org.biojava.nbio.core.sequence.DNASequence;
import org.biojava.nbio.core.sequence.compound.NucleotideCompound;
import org.biojava.nbio.core.sequence.template.SequenceView;

public class Contig extends DNASequence {
	
	public Contig() {
		this.setDNAType(DNAType.CHROMOSOME);
	}

	public Contig(String seq) throws CompoundNotFoundException {
		super(seq);
	}

	public String getID() {
		return this.getAccession().getID();
	}

	public void setID(String iD) {
		this.setAccession(new AccessionID(iD));
	}

	@Override
	public boolean equals(Object o) {
		Contig c = (Contig) o;
		if (c.getID().equals(this.getID()))
			return true;
		return false;
	}

	public String getReverseComplementSeq() {
		return 
				this.getReverseComplement().getSequenceAsString();
//		StringBuffer sb = new StringBuffer();
//		String originalSeq = this.getSequenceAsString();
//		for (int j = originalSeq.length() - 1; j >= 0; j--) {
//			switch (originalSeq.charAt(j)) {
//			case 'A':
//				sb.append("T");
//				break;
//			case 'T':
//				sb.append("A");
//				break;
//			case 'C':
//				sb.append("G");
//				break;
//			case 'G':
//				sb.append("C");
//				break;
//			case 'a':
//				sb.append("t");
//				break;
//			case 't':
//				sb.append("a");
//				break;
//			case 'c':
//				sb.append("g");
//				break;
//			case 'g':
//				sb.append("c");
//				break;
//			}
//		}
//		return sb.toString();
	}

	@Override
	public String toString() {
//		return "Contig [ID=" + this.getAccession().getID() + ", length=" + this.getLength() + "]";
		return "Contig [ID=" + this.getAccession().getID() + "]";
	}
}

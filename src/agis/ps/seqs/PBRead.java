/*
*File: agis.ps.seqs.PBRead.java
*User: mqin
*Email: mqin@ymail.com
*Date: 2016年6月2日
*/
package agis.ps.seqs;

import org.biojava.nbio.core.exceptions.CompoundNotFoundException;
import org.biojava.nbio.core.sequence.DNASequence;
import org.biojava.nbio.core.sequence.compound.DNACompoundSet;

public class PBRead extends DNASequence{

	public PBRead(String seq, DNACompoundSet dnaCompoundSet) throws CompoundNotFoundException {
		// TODO Auto-generated constructor stub
		super(seq,dnaCompoundSet);
	}

	@Override
	public String toString(){
		return "PBRead [ID: " + getAccession().getID() + ", length: " + getLength() + "]";
	}
}



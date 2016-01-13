/*
*File: agis.ps.link.Contig.java
*User: mqin
*Email: mqin@ymail.com
*Date: 2016年1月12日
*/
package agis.ps.link;

import org.biojava.nbio.core.sequence.DNASequence;

public class Contig extends DNASequence {
	private String ID; // the Id of this contig;
	private int length; // the length of this contigs;
	
	public Contig()
	{
		this.setDNAType(DNAType.CHROMOSOME);
	}

	public String getID() {
		return ID;
	}

	public void setID(String iD) {
		ID = iD;
	}

	public int getLength() {
		return length;
	}

	public void setLength(int length) {
		this.length = length;
	}
	
	

}



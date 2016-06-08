/*
*File: agis.ps.Node.java
*User: mqin
*Email: mqin@ymail.com
*Date: 2015年12月28日
*/
package agis.ps;

import java.io.Serializable;

import agis.ps.seqs.Contig;
import agis.ps.util.Strand;

public class Vertex extends Contig implements Serializable {
	private Strand strand;
	
	public Strand getStrand() {
		return strand;
	}
	public void setStrand(Strand strand) {
		this.strand = strand;
	}
	
}



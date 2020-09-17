package agis.ps.seqs;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

/**
* @author Mao Qin
* @version 2020年9月16日 上午9:25:59
* @Email mqin@outlook.com
* @Description class usage.
* @Copyright All Right Reserved 2020.
*/
public class Scaffold extends Sequence{ 
	private List<Contig> contigs = new ArrayList<Contig>();
	
	public Scaffold( ){
		super();
	}

	public String getForwardSeqs() {
		StringBuffer sb = new StringBuffer();
		int start = 0;
		for(Contig cnt : contigs){
			if(cnt.getStart() > start) {
				sb.append(String.join("", Collections.nCopies(cnt.getStart() - start, "N")));
				sb.append(cnt.getForwardSeqs());
				start = cnt.getEnd();
			} else {
				sb.append(cnt.getForwardSeqs());
				start = cnt.getEnd();
			}
		}
		if(sb.length() < this.getLength()) {
			sb.append(String.join("", Collections.nCopies(this.getLength() - start, "N")));
		}
		return sb.toString();
	}
	
	public String getComplementReverseSeqs() {
		StringBuffer sb = new StringBuffer();
		int indexEnd = this.getLength();
		for(int i = this.contigs.size() - 1; i >= 0; i--){
			Contig cnt = this.contigs.get(i);
			int cntEnd = cnt.getEnd();
			if(indexEnd > cntEnd) {
				sb.append(String.join("", Collections.nCopies(indexEnd - cntEnd, "N")));
				sb.append(cnt.getComplementReverseSeqs());
				indexEnd = cnt.getStart();
			} else {
				sb.append(cnt.getComplementReverseSeqs());
				indexEnd = cnt.getStart();
			}
		}
		if(indexEnd > 0) {
			sb.append(String.join("", Collections.nCopies(indexEnd, "N")));
		}
		return sb.toString();
	}
	
	public Boolean addContig(Contig cnt){
		if(this.contigs == null)
			this.contigs = new ArrayList<Contig>();
		return this.contigs.add(cnt);
	}
	
	public Boolean addContigs(List<Contig> cnts) {
		if(this.contigs == null)
			this.contigs = new ArrayList<Contig>();
		return this.contigs.addAll(cnts);
	}
	
	public List<Contig> getContigs(){
		return this.contigs;
	}


	@Override
	public String toString() {
		return "Scaffold [id=" + this.getId() + ", contigs size =" + contigs.size() + ", length=" + this.getLength() + "]";
	}
	
	
}

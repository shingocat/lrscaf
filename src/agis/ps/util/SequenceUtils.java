package agis.ps.util;

import java.util.List;

import agis.ps.seqs.Contig;
import agis.ps.seqs.Scaffold;
import agis.ps.seqs.SeqRegion;
import agis.ps.seqs.SequenceType;

/**
* @author Mao Qin
* @version 2020年9月16日 上午9:54:58
* @Email mqin@outlook.com
* @Description class usage.
* @Copyright All Right Reserved 2020.
*/
public class SequenceUtils {
	
	public static Scaffold seq2scaf(String id, String seqs){
		Scaffold scaf = new Scaffold();
		scaf.setId(id);
		scaf.setLength(seqs.length());
		List<SeqRegion> regions = SeqRegionFinder.findSeqRegions(seqs);
		int index = 1;
		for(SeqRegion sr : regions) {
			if(sr.getSeqType().equals(SequenceType.NORMAL)) {
				Contig cnt = new Contig();
				cnt.setId(id + "_" + index);
				cnt.setLength(sr.getLength());
				cnt.setStart(sr.getStart());
				cnt.setEnd(sr.getEnd());
				cnt.setSeqs(seqs.substring(sr.getStart(), sr.getEnd()));
				scaf.addContig(cnt);
				index = index + 1;
			}
		}
		return scaf;
	}
	
	public static String formatId(String id) {
		id = id.replaceAll("^>", "");
		id = id.split("\\s")[0];
		id = id.trim();
		return id;
	}
	
	public static String formatSeqByLength(String seqs, Integer lineLength){
		if(lineLength <= 0)
			lineLength = 80;
		StringBuffer sb = new StringBuffer();
		int seqLen = seqs.length();
		int start = 0;
		int end = 0;
		if(seqLen <= lineLength)
			end = seqLen;
		else
			end = lineLength;
		while(start <= seqLen) {
			String subSeq = seqs.substring(start, end);
			sb.append(subSeq);
			sb.append("\n");
			start = start + lineLength;
			end = start + lineLength;
			if(end > seqLen)
				end = seqLen;
		}
		return sb.toString();
	}

}

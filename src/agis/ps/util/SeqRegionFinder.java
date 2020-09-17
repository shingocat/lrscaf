package agis.ps.util;

import java.util.LinkedList;
import java.util.List;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import agis.ps.seqs.SeqRegion;
import agis.ps.seqs.SequenceType;


/**
* @author Mao Qin
* @version 2020年9月16日 上午10:26:28
* @Email mqin@outlook.com
* @Description class usage.
* @Copyright All Right Reserved 2020.
*/
public class SeqRegionFinder {
	public static List<SeqRegion> findSeqRegions(String sequence){
		List<SeqRegion> regions = new LinkedList<SeqRegion>();
		String patternStr = "n+";
		Pattern pattern = Pattern.compile(patternStr, Pattern.CASE_INSENSITIVE);
		Matcher matcher = pattern.matcher(sequence);
		// storing n sequences location
		int index = 0;
		int totalLen = sequence.length();
		int start = 0;
		int end = 0;
		while(true) {
			if(matcher.find()) {
				int gapStart = matcher.start();
				int gapEnd = matcher.end();
				if(gapStart == start) {
					SeqRegion gap = new SeqRegion(); 
					gap.setStart(gapStart);
					gap.setEnd(gapEnd);
					gap.setSeqType(SequenceType.GAP);
					gap.setIndex(index);
					regions.add(gap);
					start = gapEnd;
					index++;
				} else if(gapStart > start) {
					SeqRegion normal = new SeqRegion();
					normal.setStart(start);
					normal.setEnd(gapStart);
					normal.setSeqType(SequenceType.NORMAL);
					normal.setIndex(index);
					regions.add(normal);
					index++;
					SeqRegion gap = new SeqRegion();
					gap.setStart(gapStart);
					gap.setEnd(gapEnd);
					gap.setSeqType(SequenceType.GAP);
					gap.setIndex(index);
					regions.add(gap);
					index++;
					start = gapEnd;
				}
			} else {
				if(start < totalLen) {
					SeqRegion normal = new SeqRegion();
					normal.setIndex(index);
					normal.setStart(start);
					normal.setEnd(totalLen);
					normal.setSeqType(SequenceType.NORMAL);
					regions.add(normal);
					start = totalLen;
					index++;
				} else {
					break;
				}
			}
		}
		return regions;
	}
}

/*
*File: agis.ps.dynamicprogramming.SequenceAlignment.java
*User: mqin
*Email: mqin@ymail.com
*Date: 2017年1月10日
*Original Designer: Paul Reiners
*/
package agis.ps.align;

public abstract class SequenceAlignment extends DynamicProgramming {

	protected int startRow;
	protected int startCol;
	protected int endRow;
	protected int endCol;
	protected int match;
	protected int mismatch;
	protected int space;
	protected String[] alignments;

	public SequenceAlignment(String sequence1, String sequence2) {
		this(sequence1, sequence2, 1, -1, -1);
	}

	public SequenceAlignment(String sequence1, String sequence2, int match, int mismatch, int gap) {
		super(sequence1, sequence2);

		this.match = match;
		this.mismatch = mismatch;
		this.space = gap;
	}

	protected Object getTraceback() {
		StringBuffer align1Buf = new StringBuffer();
		StringBuffer align2Buf = new StringBuffer();
		Cell currentCell = getTracebackStartingCell();
		startRow = currentCell.getRow();
		startCol = currentCell.getCol();
		while (traceBackIsNotDone(currentCell)) {
			if (currentCell.getRow() - currentCell.getPrevCell().getRow() == 1) {
				align2Buf.insert(0, sequence2.charAt(currentCell.getRow() - 1));
			} else {
				align2Buf.insert(0, '-');
			}
			if (currentCell.getCol() - currentCell.getPrevCell().getCol() == 1) {
				align1Buf.insert(0, sequence1.charAt(currentCell.getCol() - 1));
			} else {
				align1Buf.insert(0, '-');
			}
			currentCell = currentCell.getPrevCell();
		}
		endRow = currentCell.getRow();
		endCol = currentCell.getCol();
		String[] alignments = new String[] { align1Buf.toString(), align2Buf.toString() };

		return alignments;
	}

	@Override
	public String getConsensus() {
		if (alignments == null) {
			getAlignment();
		}
		char[] seq1 = alignments[0].toCharArray();
		char[] seq2 = alignments[1].toCharArray();
		int size = seq1.length;
		StringBuffer sb = new StringBuffer();
		for (int i = 0; i < size; i++) {
			char c1 = seq1[i];
			char c2 = seq2[i];
			if (c1 == '-')
				sb.append(c2);
			else
				sb.append(c1);
		}
		consensus = sb.toString();
		return consensus;
	}

	protected abstract boolean traceBackIsNotDone(Cell currentCell);

	public int getAlignmentScore() {
		if (alignments == null) {
			getAlignment();
		}

		int score = 0;
		for (int i = 0; i < alignments[0].length(); i++) {
			char c1 = alignments[0].charAt(i);
			char c2 = alignments[1].charAt(i);
			if (c1 == '-' || c2 == '-') {
				score += space;
			} else if (c1 == c2) {
				score += match;
			} else {
				score += mismatch;
			}
		}

		return score;
	}

	public String[] getAlignment() {
		ensureTableIsFilledIn();
		alignments = (String[]) getTraceback();
		return alignments;
	}

	protected abstract Cell getTracebackStartingCell();
}

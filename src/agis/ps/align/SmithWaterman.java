/*
*File: agis.ps.dynamicprogramming.SmithWaterman.java
*User: mqin
*Email: mqin@ymail.com
*Date: 2017年1月10日
*Original Designer: Paul Reiners
*/
package agis.ps.align;

public class SmithWaterman extends SequenceAlignment {
	private Cell highScoreCell;

	public SmithWaterman(String sequence1, String sequence2) {
		super(sequence1, sequence2);
	}

	public SmithWaterman(String sequence1, String sequence2, int match, int mismatch, int gap) {
		super(sequence1, sequence2, match, mismatch, gap);
	}

	protected void initialize() {
		super.initialize();

		highScoreCell = scoreTable[0][0];
	}

	protected void fillInCell(Cell currentCell, Cell cellAbove, Cell cellToLeft, Cell cellAboveLeft) {
		int rowSpaceScore = cellAbove.getScore() + space;
		int colSpaceScore = cellToLeft.getScore() + space;
		int matchOrMismatchScore = cellAboveLeft.getScore();
		if (sequence2.charAt(currentCell.getRow() - 1) == sequence1.charAt(currentCell.getCol() - 1)) {
			matchOrMismatchScore += match;
		} else {
			matchOrMismatchScore += mismatch;
		}
		if (rowSpaceScore >= colSpaceScore) {
			if (matchOrMismatchScore >= rowSpaceScore) {
				if (matchOrMismatchScore > 0) {
					currentCell.setScore(matchOrMismatchScore);
					currentCell.setPrevCell(cellAboveLeft);
				}
			} else {
				if (rowSpaceScore > 0) {
					currentCell.setScore(rowSpaceScore);
					currentCell.setPrevCell(cellAbove);
				}
			}
		} else {
			if (matchOrMismatchScore >= colSpaceScore) {
				if (matchOrMismatchScore > 0) {
					currentCell.setScore(matchOrMismatchScore);
					currentCell.setPrevCell(cellAboveLeft);
				}
			} else {
				if (colSpaceScore > 0) {
					currentCell.setScore(colSpaceScore);
					currentCell.setPrevCell(cellToLeft);
				}
			}
		}
		if (currentCell.getScore() > highScoreCell.getScore()) {
			highScoreCell = currentCell;
		}
	}
	
	@Override
	public String getConsensus()
	{
		if (alignments == null) {
			getAlignment();
		}
		String subSeq1 = sequence1.substring(0, endCol);
		String subSeq2 = sequence2.substring(startRow);
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
		sb.insert(0, subSeq1);
		sb.insert(sb.length(), subSeq2);
		consensus = sb.toString();
		return consensus;	
	}

	@Override
	public String toString() {
		return "[SmithWaterman: sequence1=" + sequence1 + ", sequence2=" + sequence2 + "]";
	}

	@Override
	protected boolean traceBackIsNotDone(Cell currentCell) {
		return currentCell.getScore() != 0;
	}

	@Override
	protected Cell getTracebackStartingCell() {
		return highScoreCell;
	}

	@Override
	protected Cell getInitialPointer(int row, int col) {
		return null;
	}

	@Override
	protected int getInitialScore(int row, int col) {
		return 0;
	}

}

/*
*File: agis.ps.dynamicprogramming.DynamicProgramming.java
*User: mqin
*Email: mqin@ymail.com
*Date: 2017年1月10日
*Original Designer: Paul Reiners
*
*/
package agis.ps.align;

public abstract class DynamicProgramming {
	protected String sequence1;
	protected String sequence2;
	protected Cell[][] scoreTable;
	protected boolean tableIsFilledIn;
	protected boolean isInitialized;
	protected String consensus;

	public DynamicProgramming(String sequence1, String sequence2) {
		this.sequence1 = sequence1;
		this.sequence2 = sequence2;
		scoreTable = new Cell[sequence2.length() + 1][sequence1.length() + 1];
	}

	public int[][] getScoreTable() {
		ensureTableIsFilledIn();

		int[][] matrix = new int[scoreTable.length][scoreTable[0].length];
		for (int i = 0; i < matrix.length; i++) {
			for (int j = 0; j < matrix[i].length; j++) {
				matrix[i][j] = scoreTable[i][j].getScore();
			}
		}

		return matrix;
	}

	protected void initializeScores() {
		for (int i = 0; i < scoreTable.length; i++) {
			for (int j = 0; j < scoreTable[i].length; j++) {
				scoreTable[i][j].setScore(getInitialScore(i, j));
			}
		}
	}

	protected void initializePointers() {
		for (int i = 0; i < scoreTable.length; i++) {
			for (int j = 0; j < scoreTable[i].length; j++) {
				scoreTable[i][j].setPrevCell(getInitialPointer(i, j));
			}
		}
	}

	protected void initialize() {
		for (int i = 0; i < scoreTable.length; i++) {
			for (int j = 0; j < scoreTable[i].length; j++) {
				scoreTable[i][j] = new Cell(i, j);
			}
		}
		initializeScores();
		initializePointers();

		isInitialized = true;
	}

	protected abstract Cell getInitialPointer(int row, int col);

	protected abstract int getInitialScore(int row, int col);

	protected abstract void fillInCell(Cell currentCell, Cell cellAbove, Cell cellToLeft, Cell cellAboveLeft);

	protected void fillIn() {
		for (int row = 1; row < scoreTable.length; row++) {
			for (int col = 1; col < scoreTable[row].length; col++) {
				Cell currentCell = scoreTable[row][col];
				Cell cellAbove = scoreTable[row - 1][col];
				Cell cellToLeft = scoreTable[row][col - 1];
				Cell cellAboveLeft = scoreTable[row - 1][col - 1];
				fillInCell(currentCell, cellAbove, cellToLeft, cellAboveLeft);
			}
		}

		tableIsFilledIn = true;
	}

	abstract protected Object getTraceback();

	abstract public String getConsensus();

	public void printScoreTable() {
		ensureTableIsFilledIn();
		for (int i = 0; i < sequence2.length() + 2; i++) {
			for (int j = 0; j < sequence1.length() + 2; j++) {
				if (i == 0) {
					if (j == 0 || j == 1) {
						System.out.print("  ");
					} else {
						if (j == 2) {
							System.out.print("     ");
						} else {
							System.out.print("   ");
						}
						System.out.print(sequence1.charAt(j - 2));
					}
				} else if (j == 0) {
					if (i == 1) {
						System.out.print("  ");
					} else {
						System.out.print(" " + sequence2.charAt(i - 2));
					}
				} else {
					String toPrint;
					Cell currentCell = scoreTable[i - 1][j - 1];
					Cell prevCell = currentCell.getPrevCell();
					if (prevCell != null) {
						if (currentCell.getCol() == prevCell.getCol() + 1
								&& currentCell.getRow() == prevCell.getRow() + 1) {
							toPrint = "\\";
						} else if (currentCell.getCol() == prevCell.getCol() + 1) {
							toPrint = "-";
						} else {
							toPrint = "|";
						}
					} else {
						toPrint = " ";
					}
					int score = currentCell.getScore();
					String s = String.format("%1$3d", score);
					toPrint += s;
					System.out.print(toPrint);
				}

				System.out.print(' ');
			}
			System.out.println();
		}
	}

	protected void ensureTableIsFilledIn() {
		if (!isInitialized) {
			initialize();
		}
		if (!tableIsFilledIn) {
			fillIn();
		}
	}
}

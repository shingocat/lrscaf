/*
*File: agis.ps.util.Consensusser.java
*User: mqin
*Email: mqin@ymail.com
*Date: 2016年4月16日
*/
package agis.ps.util;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Vector;

// if indicator equal to nw, specific the global alignment, abbreviation of Needleman-Wunsch; 
// else indicator equal to sw, specific the local alignment, abbreviation of Smith-Waterman; 
public class Consensusser {
	// fix the penalty of the alignment matrix score;
	private int match = 1;
	private int mismatch = -1;
	private int gap = -2;
	
	public String getConsensus(String seq1, String seq2, String indicator)
	{
		if(indicator.equalsIgnoreCase("nw"))
			return nw_2(seq1, seq2);
		else
			return sw(seq1, seq2);
	}
	
	// global alignment
	private String nw(String seq1, String seq2)
	{
		int len1 = seq1.length();
		int len2 = seq2.length();
		// initiate the scores matrix
		int [][] scores = new int[len1 + 1][len2 + 1];
		Pointers [][] indexMatrix = new Pointers[len1 + 1][len2 + 1];
		for(int i = 0; i <= len1; i++) {
			for(int j = 0; j <= len2; j++) {
				Pointers ps = new Pointers();
				if(i == 0) {
					scores[i][j] = 0;
					Pointer p = new Pointer(i, j);
					ps.addPointer(p);
					indexMatrix[i][j] = ps;
					continue;
				}
				if(j == 0) {
					scores[i][j] = 0;
					Pointer p = new Pointer(i, j);
					ps.addPointer(p);
					indexMatrix[i][j] = ps;
					continue;
				}

				int left = scores[i][j - 1];
				int upper = scores[i - 1][j];
				int diagonal = scores[i - 1][j - 1];
				// for the diagonal value;
				if(seq1.charAt(i - 1) == seq2.charAt(j - 1))
					diagonal += match;
				else
					diagonal += mismatch;
				// the above value of gap
				upper += gap;
				// the left value of gap
				left += gap;
				// compute the maximum;
				int max = upper >= left ?
				          (upper >= diagonal ? upper : diagonal) :
					          (left >= diagonal ? left : diagonal);
				scores[i][j] = max;
				if(max == upper) {
					Pointer p = new Pointer();
					p.setRow(i - 1);
					p.setCol(j);
					ps.addPointer(p);
				}
				if(max == left) {
					Pointer p = new Pointer();
					p.setRow(i);
					p.setCol(j - 1);
					ps.addPointer(p);
				}
				if(max == diagonal) {
					Pointer p = new Pointer();
					p.setRow(i - 1);
					p.setCol(j - 1);
					ps.addPointer(p);
				}
				indexMatrix[i][j] = ps;
			}
		}
		// print the matrix
//		System.out.println("Aligned Matrix:");
//		for(int i = 0; i <= len1; i++) {
//			System.out.println(Arrays.toString(scores[i]));
//		}
		//System.out.println(Arrays.deepToString(scores));

		// traceback the matrix by pointer;
		String algSeq1 = "";
		String algSeq2 = "";
		String algPtn = "";
		int i = len1;
		int j = len2;
		// System.out.println("seq 1 ===" + seq1);
		// System.out.println("seq 2 ===" + seq2);
		Pointers ps = null;
		Pointer p = null;
		while(true) {
			if(i == 0 && j == 0) {
				break;
			}
			// System.out.println("i = " + i);
			// System.out.println("j = " + j);
			ps = indexMatrix[i][j];
			p = ps.getPointer(0);
			int row = p.getRow();
			int col = p.getCol();
			// System.out.println("row = " + row);
			// System.out.println("col = " + col);
			if(i == 0 && j !=0) {
				algSeq1 += "_";
				algSeq2 += seq2.charAt(j - 1);
				algPtn += "*";
				j--;
				continue;
			} else if (i != 0 && j == 0) {
				algSeq1 += seq1.charAt(i - 1);
				algSeq2 += "_";
				algPtn += "*";
				i--;
				continue;
			}
			if(row == i - 1 && col == j - 1) {
				if(scores[i][j] == scores[row][col] + match) {
					algSeq1 += seq1.charAt(i - 1);
					algSeq2 += seq2.charAt(j - 1);
					algPtn += "|";
				} else {
					algSeq1 += seq1.charAt(i - 1);
					algSeq2 += seq2.charAt(j - 1);
					algPtn += "*";
				}
				i--;
				j--;
			} else if(row == i && col == j - 1) {
				algSeq1 += "_";
				algSeq2 += seq2.charAt(j - 1);
				algPtn += "*";
				j--;
			} else if(row == i - 1 && col == j) {
				algSeq1 += seq1.charAt(i - 1);
				algSeq2 += "_";
				algPtn += "*";
				i--;
			}
		}
		String temp = consesus(algSeq1, algPtn, algSeq2);
		indexMatrix = null;
		scores = null;
		ps = null;
		p = null;
		System.gc();
		return temp;
	}
	
	private String nw_2(String seq1, String seq2)
	{
		int len1 = seq1.length();
		int len2 = seq2.length();
		// initiate the scores matrix
		int [][] scores = new int[len1 + 1][len2 + 1];
		//Pointers [][] indexMatrix = new Pointers[len1 + 1][len2 + 1];
		for(int i = 0; i <= len1; i++) {
			for(int j = 0; j <= len2; j++) {
				//Pointers ps = new Pointers();
				if(i == 0) {
					scores[i][j] = 0;
					//Pointer p = new Pointer(i, j);
					//ps.addPointer(p);
					//indexMatrix[i][j] = ps;
					continue;
				}
				if(j == 0) {
					scores[i][j] = 0;
					//Pointer p = new Pointer(i, j);
					//ps.addPointer(p);
					//indexMatrix[i][j] = ps;
					continue;
				}

				int left = scores[i][j - 1];
				int upper = scores[i - 1][j];
				int diagonal = scores[i - 1][j - 1];
				// for the diagonal value;
				if(seq1.charAt(i - 1) == seq2.charAt(j - 1))
					diagonal += match;
				else
					diagonal += mismatch;
				// the above value of gap
				upper += gap;
				// the left value of gap
				left += gap;
				// compute the maximum;
				int max = upper >= left ?
				          (upper >= diagonal ? upper : diagonal) :
					          (left >= diagonal ? left : diagonal);
				scores[i][j] = max;
				/* if(max == upper) {
					Pointer p = new Pointer();
					p.setRow(i - 1);
					p.setCol(j);
					ps.addPointer(p);
				}
				if(max == left) {
					Pointer p = new Pointer();
					p.setRow(i);
					p.setCol(j - 1);
					ps.addPointer(p);
				}
				if(max == diagonal) {
					Pointer p = new Pointer();
					p.setRow(i - 1);
					p.setCol(j - 1);
					ps.addPointer(p);
				}
				indexMatrix[i][j] = ps; */
			}
		}
		// print the matrix
//		System.out.println("Aligned Matrix:");
//		for(int i = 0; i <= len1; i++) {
//			System.out.println(Arrays.toString(scores[i]));
//		}
		//System.out.println(Arrays.deepToString(scores));

		// traceback the matrix by pointer;
		String algSeq1 = "";
		String algSeq2 = "";
		String algPtn = "";
		int i = len1;
		int j = len2;
		// System.out.println("seq 1 ===" + seq1);
		// System.out.println("seq 2 ===" + seq2);
		while(true) {
			if(i == 0 && j == 0) {
				break;
			}
			// System.out.println("i = " + i);
			// System.out.println("j = " + j);
			/* Pointers ps = indexMatrix[i][j];
			Pointer p = ps.getPointer(0);
			int row = p.getRow();
			int col = p.getCol(); */
			// System.out.println("row = " + row);
			// System.out.println("col = " + col);
			if(i == 0 && j !=0) {
				algSeq1 += "_";
				algSeq2 += seq2.charAt(j - 1);
				algPtn += "*";
				j--;
				continue;
			} else if (i != 0 && j == 0) {
				algSeq1 += seq1.charAt(i - 1);
				algSeq2 += "_";
				algPtn += "*";
				i--;
				continue;
			} else
			{
				if(seq1.charAt(i-1) == (seq2.charAt(j-1)))
				{
					algSeq1 += seq1.charAt(i-1);
					algSeq2 += seq2.charAt(j-1);
					algPtn += "|";
					i--;
					j--;
				} else
				{
					int lScore = scores[i][j-1];
					int uScore = scores[i - 1][j];
					int dScore = scores[i - 1][j - 1];
					int max = uScore >= lScore ?
				          (uScore >= dScore ? uScore : dScore) :
					          (lScore >= dScore ? lScore : dScore);
					if(max == dScore)
					{
						algSeq1 += seq1.charAt(i-1);
						algSeq2 += seq2.charAt(j-1);
						algPtn += "*";
						i--;
						j--;
					} else if(max == uScore)
					{
						algSeq1 += seq1.charAt(i-1);
						algSeq2 += "_";
						algPtn += "*";
						i--;
					} else{
						algSeq1 += "_";
						algSeq2 += seq2.charAt(j-1);
						algPtn += "*";
						j--;
					}
				}
			}
			
			/* if(row == i - 1 && col == j - 1) {
				if(scores[i][j] == scores[row][col] + match) {
					algSeq1 += seq1.charAt(i - 1);
					algSeq2 += seq2.charAt(j - 1);
					algPtn += "|";
				} else {
					algSeq1 += seq1.charAt(i - 1);
					algSeq2 += seq2.charAt(j - 1);
					algPtn += "*";
				}
				i--;
				j--;
			} else if(row == i && col == j - 1) {
				algSeq1 += "_";
				algSeq2 += seq2.charAt(j - 1);
				algPtn += "*";
				j--;
			} else if(row == i - 1 && col == j) {
				algSeq1 += seq1.charAt(i - 1);
				algSeq2 += "_";
				algPtn += "*";
				i--;
			} */
		}
		
		String temp = consesus(algSeq1, algPtn, algSeq2);
		scores = null;
		System.gc();
		return temp;
	}
	
	// local alignment
	private String sw(String seq1, String seq2)
	{
		int len1 = seq1.length();
		int len2 = seq2.length();
		// initiate the scores matrix
		int [][] scores = new int[len1 + 1][len2 + 1];
		Pointers [][] indexMatrix = new Pointers[len1 + 1][len2 + 1];
		for(int i = 0; i <= len1; i++) {
			for(int j = 0; j <= len2; j++) {
				Pointers ps = new Pointers();
				if(i == 0) {
					scores[i][j] = 0;
					Pointer p = new Pointer(i, j);
					ps.addPointer(p);
					indexMatrix[i][j] = ps;
					continue;
				}
				if(j == 0) {
					scores[i][j] = 0;
					Pointer p = new Pointer(i, j);
					ps.addPointer(p);
					indexMatrix[i][j] = ps;
					continue;
				}

				int left = scores[i][j - 1];
				int upper = scores[i - 1][j];
				int diagonal = scores[i - 1][j - 1];
				// for the diagonal value;
				if(seq1.charAt(i - 1) == seq2.charAt(j - 1))
					diagonal += match;
				else
					diagonal += mismatch;
				// the above value of gap
				upper += gap;
				// the left value of gap
				left += gap;
				// compute the maximum;
				int max = upper >= left ?
				          (upper >= diagonal ? upper : diagonal) :
					          (left >= diagonal ? left : diagonal);
				if(max >= 0) {
					scores[i][j] = max;
					if(max == upper) {
						Pointer p = new Pointer();
						p.setRow(i - 1);
						p.setCol(j);
						ps.addPointer(p);
					}
					if(max == left) {
						Pointer p = new Pointer();
						p.setRow(i);
						p.setCol(j - 1);
						ps.addPointer(p);
					}
					if(max == diagonal) {
						Pointer p = new Pointer();
						p.setRow(i - 1);
						p.setCol(j - 1);
						ps.addPointer(p);
					}
					indexMatrix[i][j] = ps;
				} else {
					// if less than 0, not store the pointer;
					scores[i][j] = 0;
				}
			}
		}
		// print the matrix
//		System.out.println("Aligned Matrix:");
//		for(int i = 0; i <= len1; i++) {
//			System.out.println(Arrays.toString(scores[i]));
//		}

		// traceback the matrix by pointer;
		String algSeq1 = "";
		String algSeq2 = "";
		String algPtn = "";
		int rowIndex = 0;
		int colIndex = 0;
		// find the maximum score in the scores matrix, 
		// and define the i and j to start traceback;
		int max = 0;
		for(int i = 0; i <= seq1.length(); i++)
		{
			for(int j = 0; j <= seq2.length(); j++)
			{
				if(scores[i][j] >= max)
				{
					max = scores[i][j];
					rowIndex = i;
					colIndex = j;
				}
			}
		}
		// System.out.println("rowindex ===" + rowIndex);
		// System.out.println("colindex ===" + colIndex);
		while(true) {
			if(scores[rowIndex][colIndex] == 0) {
				break;
			}
			// System.out.println("i = " + i);
			// System.out.println("j = " + j);
			Pointers ps = indexMatrix[rowIndex][colIndex];
			Pointer p = ps.getPointer(0);
			int row = p.getRow();
			int col = p.getCol();
			// System.out.println("row = " + row);
			// System.out.println("col = " + col);
			if(rowIndex == 0 && colIndex !=0) {
				algSeq1 += "_";
				algSeq2 += seq2.charAt(colIndex - 1);
				algPtn += "*";
				colIndex--;
				continue;
			} else if (rowIndex != 0 && colIndex == 0) {
				algSeq1 += seq1.charAt(rowIndex - 1);
				algSeq2 += "_";
				algPtn += "*";
				rowIndex--;
				continue;
			}
			if(row == rowIndex - 1 && col == colIndex - 1) {
				if(scores[rowIndex][colIndex] == scores[row][col] + match) {
					algSeq1 += seq1.charAt(rowIndex - 1);
					algSeq2 += seq2.charAt(colIndex - 1);
					algPtn += "|";
				} else {
					algSeq1 += seq1.charAt(rowIndex - 1);
					algSeq2 += seq2.charAt(colIndex - 1);
					algPtn += "*";
				}
				rowIndex--;
				colIndex--;
			} else if(row == rowIndex && col == colIndex - 1) {
				algSeq1 += "_";
				algSeq2 += seq2.charAt(colIndex - 1);
				algPtn += "*";
				colIndex--;
			} else if(row == rowIndex - 1 && col == colIndex) {
				algSeq1 += seq1.charAt(rowIndex - 1);
				algSeq2 += "_";
				algPtn += "*";
				rowIndex--;
			}
		}
		return consesus(algSeq1, algPtn, algSeq2);
	}
	
	private String consesus(String seq1, String pattern, String seq2) {
		int len = pattern.length();
		StringBuilder sb = new StringBuilder();
		for(int i = 0; i < len; i++) {
			String s1 = String.valueOf(seq1.charAt(i));
			String s2 = String.valueOf(seq2.charAt(i));
			String p = String.valueOf(pattern.charAt(i));
			if(p.equalsIgnoreCase("|")) {
				sb.append(s1);
			} else {
				if(s1.equalsIgnoreCase("_")) {
					sb.append(s2);
				} else if(s2.equalsIgnoreCase("_")) {
					sb.append(s1);
				} else {
					sb.append(s1);
				}
			}
		}
		String temp = sb.reverse().toString();
		sb = null;
		return temp;
	}

	class Pointers {
//		private List<Pointer> list = new ArrayList<Pointer>();
		private List<Pointer> list = new Vector<Pointer>();
		Pointers() {
			super();
		}
		public void addPointer(Pointer p) {
			list.add(p);
		}
		public Pointer getPointer(int index) {
			return list.get(index);
		}
		public int getSize() {
			return list.size();
		}
	}

	class Pointer {
		private int row;
		private int col;

		Pointer() {
			super();
		}

		Pointer(int row, int col) {
			this.row = row;
			this.col = col;
		}

		public void setRow(int row) {
			this.row = row;
		}
		public int getRow() {
			return this.row;
		}
		public void setCol(int col) {
			this.col = col;
		}
		public int getCol() {
			return this.col;
		}
	}
}



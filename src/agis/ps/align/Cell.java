/*
*File: agis.ps.dynamicprogramming.Cell.java
*User: mqin
*Email: mqin@ymail.com
*Date: 2017年1月10日
*Original Designer: Paul Reiners
*
*
*/
package agis.ps.align;

public class Cell {
	 private Cell prevCell;
	   private int score;
	   private int row;
	   private int col;

	   public Cell(int row, int col) {
	      this.row = row;
	      this.col = col;
	   }

	   /**
	    * @param score
	    *           the score to set
	    */
	   public void setScore(int score) {
	      this.score = score;
	   }

	   /**
	    * @return the score
	    */
	   public int getScore() {
	      return score;
	   }

	   /**
	    * @param prevCell
	    *           the prevCell to set
	    */
	   public void setPrevCell(Cell prevCell) {
	      this.prevCell = prevCell;
	   }

	   /**
	    * @return the row
	    */
	   public int getRow() {
	      return row;
	   }

	   /**
	    * @return the col
	    */
	   public int getCol() {
	      return col;
	   }

	   /**
	    * @return the prevCell
	    */
	   public Cell getPrevCell() {
	      return prevCell;
	   }

	   /*
	    * (non-Javadoc)
	    * 
	    * @see java.lang.Object#toString()
	    */
	   @Override
	   public String toString() {
	      return "Cell(" + row + ", " + col + "): score=" + score + ", prevCell="
	            + prevCell + "]";
	   }
}



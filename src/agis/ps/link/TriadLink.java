/*
*File: agis.ps.TriadLink.java
*User: mqin
*Email: mqin@ymail.com
*Date: 2016年5月13日
*/
package agis.ps.link;

import agis.ps.seqs.Contig;

public class TriadLink implements ILink {
	private Contig previous; // previous contig in link;
	private Contig middle; // middle contig in link;
	private Contig last; // last contig in link
	private int supLinks; // support triad link nums;
	private boolean isValid = true; // indicator for used or not;
	
	public TriadLink(){
		
	}
	
	public TriadLink(Contig previous, Contig middle, Contig last)
	{
		this.previous = previous;
		this.middle = middle;
		this.last = last;
	}
	
	public TriadLink(Contig previous, Contig middle, Contig last, boolean isValid)
	{
		this(previous, middle, middle);
		this.isValid = isValid;
	}

	public boolean isValid() {
		return isValid;
	}

	public void setValid(boolean isValid) {
		this.isValid = isValid;
	}

	public Contig getPrevious() {
		return previous;
	}

	public void setPrevious(Contig previous) {
		this.previous = previous;
	}

	public Contig getMiddle() {
		return middle;
	}

	public void setMiddle(Contig middle) {
		this.middle = middle;
	}

	public Contig getLast() {
		return last;
	}

	public void setLast(Contig last) {
		this.last = last;
	}
	
	public int getSupLinks() {
		return supLinks;
	}

	public void setSupLinks(int supLinks) {
		this.supLinks = supLinks;
	}
	
	public boolean isContain(Contig contig)
	{
		if(this.previous != null && this.previous.equals(contig))
			return true;
		if(this.middle != null && this.middle.equals(contig))
			return true;
		if(this.last != null && this.last.equals(contig))
			return true;
		return false;
	}

	@Override
	public int hashCode() {
		final int prime = 31;
		int result = 1;
		result = prime * result + ((last == null) ? 0 : last.hashCode());
		result = prime * result + ((middle == null) ? 0 : middle.hashCode());
		result = prime * result + ((previous == null) ? 0 : previous.hashCode());
		return result;
	}

	@Override
	public boolean equals(Object obj) {
		if (this == obj)
			return true;
		if (obj == null)
			return false;
		if (getClass() != obj.getClass())
			return false;
		TriadLink other = (TriadLink) obj;
		// the middle element must be equal;
//		if (middle == null) {
//			if (other.middle != null)
//				return false;
//		} else if (!middle.equals(other.middle))
//			return false;
//		if(previous.equals(other.previous) && last.equals(other.last))
//			return true;
//		else if(previous.equals(other.last) && last.equals(other.previous))
//			return true;
//		else
//			return false;
		if(middle == null)
		{
			if(other.middle == null)
			{
				if(previous.equals(other.previous) && last.equals(other.last))
					return true;
				else if(previous.equals(other.last) && last.equals(other.previous))
					return true;
				else
					return false;
			} else
			{
				return false;
			}
		} else
		{
			if(other.middle == null)
			{
//				change this equal method
				if(previous.equals(other.previous) && last.equals(other.last))
					return true;
				else if(previous.equals(other.last) && last.equals(other.previous))
					return true;
				else
					return false;
//				return false;
			} else
			{
				if(previous.equals(other.previous) && last.equals(other.last))
					return true;
				else if(previous.equals(other.last) && last.equals(other.previous))
					return true;
				else
					return false;
			}
		}
	}

	@Override
	public String toString() {
		return "TriadLink [previous=" + previous + ", middle=" + middle + ", last=" + last + ", supLinks=" + supLinks
				+ ", valid=" + isValid +"]";
	}

	@Override
	public int getDistance() {
		// TODO Auto-generated method stub
		return 0;
	}

}



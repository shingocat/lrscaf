/*
*File: agis.ps.link.ContInOut.java
*User: mqin
*Email: mqin@ymail.com
*Date: 2016年1月20日
*/
package agis.ps.link;

import java.io.Serializable;


public class ContInOut implements Serializable {
	private String id;
	private Integer indegrees;
	private Integer outdegrees;
	public String getId() {
		return id;
	}
	public void setId(String id) {
		this.id = id;
	}
	public Integer getIndegrees() {
		return indegrees;
	}
	public void setIndegrees(Integer indegrees) {
		this.indegrees = indegrees;
	}
	public Integer getOutdegrees() {
		return outdegrees;
	}
	public void setOutdegrees(Integer outdegrees) {
		this.outdegrees = outdegrees;
	}
	@Override
	public int hashCode() {
		final int prime = 31;
		int result = 1;
		result = prime * result + ((id == null) ? 0 : id.hashCode());
		result = prime * result + ((indegrees == null) ? 0 : indegrees.hashCode());
		result = prime * result + ((outdegrees == null) ? 0 : outdegrees.hashCode());
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
		ContInOut other = (ContInOut) obj;
		if (id == null) {
			if (other.id != null)
				return false;
		} else if (!id.equals(other.id))
			return false;
		if (indegrees == null) {
			if (other.indegrees != null)
				return false;
		} else if (!indegrees.equals(other.indegrees))
			return false;
		if (outdegrees == null) {
			if (other.outdegrees != null)
				return false;
		} else if (!outdegrees.equals(other.outdegrees))
			return false;
		return true;
	}
	@Override
	public String toString() {
		return "ContInOut [id=" + id + ", indegrees=" + indegrees + ", outdegrees=" + outdegrees + "]";
	}
}



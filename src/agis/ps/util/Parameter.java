/*
*File: agis.ps.util.Parameter.java
*User: mqin
*Email: mqin@ymail.com
*Date: 2016年1月22日
*/
package agis.ps.util;

import java.io.Serializable;

public class Parameter implements Serializable {

	private static final long serialVersionUID = 1L;
	private String cntFile;
	private String algFile;
	private String outFolder;
	private Integer minContLen;
	private Integer minPBLen;
	private Integer minOLLen;
	private Double minOLRatio;
	private Integer maxOHLen;
	private Double maxOHRatio;
	private Integer maxEndLen;
	private Double maxEndRatio;
	private Integer minSupLinks;
	private Integer maxSupLinks;
	private String type; // m for m5, s for sam, s for bam;

	public String getType() {
		return type;
	}

	public void setType(String type) {
		this.type = type;
	}

	public String getCntFile() {
		return cntFile;
	}

	public void setCntFile(String cntFile) {
		this.cntFile = cntFile;
	}

	public String getAlgFile() {
		return algFile;
	}

	public void setAlgFile(String algFile) {
		this.algFile = algFile;
	}

	public String getOutFolder() {
		return outFolder;
	}

	public void setOutFolder(String outFolder) {
		this.outFolder = outFolder;
	}

	public Integer getMinContLen() {
		return minContLen;
	}

	public void setMinContLen(Integer minContLen) {
		this.minContLen = minContLen;
	}

	public Integer getMinPBLen() {
		return minPBLen;
	}

	public void setMinPBLen(Integer minPBLen) {
		this.minPBLen = minPBLen;
	}

	public Integer getMinOLLen() {
		return minOLLen;
	}

	public void setMinOLLen(Integer minOLLen) {
		this.minOLLen = minOLLen;
	}

	public Double getMinOLRatio() {
		return minOLRatio;
	}

	public void setMinOLRatio(Double minOLRatio) {
		this.minOLRatio = minOLRatio;
	}

	public Integer getMinSupLinks() {
		return minSupLinks;
	}

	public void setMinSupLinks(Integer minSupLinks) {
		this.minSupLinks = minSupLinks;
	}

	public Integer getMaxSupLinks() {
		return maxSupLinks;
	}

	public void setMaxSupLinks(Integer maxSupLinks) {
		this.maxSupLinks = maxSupLinks;
	}
	
	public Integer getMaxOHLen() {
		return maxOHLen;
	}

	public void setMaxOHLen(Integer maxOHLen) {
		this.maxOHLen = maxOHLen;
	}

	public Double getMaxOHRatio() {
		return maxOHRatio;
	}

	public void setMaxOHRatio(Double maxOHRatio) {
		this.maxOHRatio = maxOHRatio;
	}
	
	public Integer getMaxEndLen() {
		return maxEndLen;
	}

	public void setMaxEndLen(Integer maxEndLen) {
		this.maxEndLen = maxEndLen;
	}

	public Double getMaxEndRatio() {
		return maxEndRatio;
	}

	public void setMaxEndRatio(Double maxEndRatio) {
		this.maxEndRatio = maxEndRatio;
	}

	@Override
	public String toString() {
		return "Parameter [cntFile=" + cntFile + ", algFile=" + algFile + ", outFolder=" + outFolder + ", minContLen="
				+ minContLen + ", minPBLen=" + minPBLen + ", minOLLen=" + minOLLen + ", minOLRatio=" + minOLRatio
				+ ", maxOHLen=" + maxOHLen + ", maxOHRatio=" + maxOHRatio + ", maxEndLen=" + maxEndLen
				+ ", maxEndRatio=" + maxEndRatio + ", minSupLinks=" + minSupLinks + ", maxSupLinks=" + maxSupLinks
				+ ", type=" + type + "]";
	}
	
}

/*
*File: agis.ps.util.Parameter.java
*User: mqin
*Email: mqin@ymail.com
*Date: 2016年1月22日
*/
package agis.ps.util;

import java.io.Serializable;

public class Parameter implements Serializable {
	// default value
	private static final long serialVersionUID = 1L;
	private String cntFile;
	private String algFile;
	private String outFolder = System.getProperty("user.dir"); // The current directory if not setted;
	private Integer minContLen = 3000; // minimum contig 3k bp in length;
	private Integer minPBLen = 5000; // maximum pacbio read 5k bp in length;
	private Integer minOLLen = 2400; // minimum overlap length 3k bp;
	private Double minOLRatio = 0.8d; // minimum ratio of overlap, if contig_length * ratio > default minOLLen, then the overlap length should be larger than contig_length * ratio;
	private Integer maxOHLen = 300; // the maximum overhang length 300 bp;
	private Double maxOHRatio = 0.1d; // maximum overhang ratio, if contig_legnth * ratio > default maxOHLen, then the overhang length should be less than default maxOHLen, else using the length * ratio; 
	private Integer maxEndLen = 500; // maximum Ending length 500 bp for defining the the contig is inner or outer;
	private Double maxEndRatio = 0.1d; // maximum ending ratio, if pacbio_length * ratio > maxEndLen, the the ending length should not be larger than maxEndLen, else using the length*ratio;
	private Integer minSupLinks = 3; // minimum supported links number: 3;
	private Integer maxSupLinks = 70; // maximum supported links number: 70;
	private String type; // m for m5, s for sam, s for bam;
	private Double identity = 0.8d; // compute the identity between contig and pacbio, formula as: match / (match + mismatch + numsIn + numsDel); 
	private boolean isUseOLLink = false; // use overlap link into build edges;
	
	public boolean isUseOLLink() {
		return isUseOLLink;
	}

	public void setUseOLLink(boolean isUseOLLink) {
		this.isUseOLLink = isUseOLLink;
	}

	public Double getIdentity() {
		return identity;
	}

	public void setIdentity(Double identity) {
		this.identity = identity;
	}

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

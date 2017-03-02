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
	private String pbFile;
	private String cntFile;
	private String algFile;
	private String outFolder = System.getProperty("user.dir"); // The current directory if not setted;
	private Integer minContLen = 3000; // minimum contig 3k bp in length;
	private Integer minPBLen = 5000; // maximum pacbio read 5k bp in length;
	private Integer minOLLen = 2400; // minimum overlap length 3k bp;
	private Double minOLRatio = 0.8d; // minimum ratio of overlap, if contig_length * ratio > default minOLLen, then the overlap length should be larger than contig_length * ratio;
	private Integer maxOHLen = 300; // the maximum overhang length 300 bp;
	private Double maxOHRatio = 0.1d; // maximum overhang ratio, if contig_legnth * ratio > default maxOHLen, then the overhang length should be less than default maxOHLen, else using the length * ratio; 
	private Integer maxEndLen = 300; // maximum Ending length 300 bp for defining the the contig is inner or outer;
	private Double maxEndRatio = 0.1d; // maximum ending ratio, if pacbio_length * ratio > maxEndLen, the the ending length should not be larger than maxEndLen, else using the length*ratio;
	private Integer minSupLinks = 1; // minimum supported links number: 1;
	private Integer maxSupLinks = 70; // maximum supported links number: 70;
	private String type; // m for m5, s for sam, s for bam;
	private Double identity = 0.8d; // compute the identity between contig and pacbio, formula as: match / (match + mismatch + numsIn + numsDel); 
	private boolean isUseOLLink = false; // use overlap link into build edges;
	private Double ratio = 0.2; // use this ratio to delete edge by supported link ratio 
	private boolean isRepMask = false; // whether repeat is mask;
	private boolean isGapFilling = false; // whether gap is filled;
	private int tipLength = 1500; // tip length is 1500 bp, larger than this do not considering as tip
	private double iqrTime = 1.5;
	private int mmcm = 8; // only for minimap output, default: 8;
	
	public int getMmcm() {
		return mmcm;
	}

	public void setMmcm(int mmcm) {
		this.mmcm = mmcm;
	}

	public double getIqrTime() {
		return iqrTime;
	}

	public void setIqrTime(double iqrTime) {
		if(iqrTime <= 0)
			throw new IllegalArgumentException("IQR Time could not be less than 0!");
		this.iqrTime = iqrTime;
	}

	public boolean isGapFilling() {
		return isGapFilling;
	}

	public void setGapFilling(boolean isGapFilling) {
		this.isGapFilling = isGapFilling;
	}

	public String getPbFile() {
		return pbFile;
	}

	public void setPbFile(String pbFile) {
		this.pbFile = pbFile;
	}

	public boolean isRepMask() {
		return isRepMask;
	}

	public void setRepMask(boolean isRepMask) {
		this.isRepMask = isRepMask;
	}

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
		if(identity <= 0)
			throw new IllegalArgumentException("Identity could not be less than 0!");
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
		if(minContLen <= 0)
			throw new IllegalArgumentException("Minimum Contig Length could not be less than 0!");
		this.minContLen = minContLen;
	}

	public Integer getMinPBLen() {
		return minPBLen;
	}

	public void setMinPBLen(Integer minPBLen) {
		if(minPBLen <= 0)
			throw new IllegalArgumentException("Minimum Pacbio's Read Length could not be less than 0!");
		this.minPBLen = minPBLen;
	}

	public Integer getMinOLLen() {
		return minOLLen;
	}

	public void setMinOLLen(Integer minOLLen) {
		if(minOLLen <= 0)
			throw new IllegalArgumentException("Minimum Overlap Length could not be less than 0!");
		this.minOLLen = minOLLen;
	}

	public Double getMinOLRatio() {
		return minOLRatio;
	}

	public void setMinOLRatio(Double minOLRatio) {
		if(minOLRatio <= 0)
			throw new IllegalArgumentException("Minimum Overlap Ratio could not be less than 0!");
		this.minOLRatio = minOLRatio;
	}

	public Integer getMinSupLinks() {
		return minSupLinks;
	}

	public void setMinSupLinks(Integer minSupLinks) {
		if(minSupLinks <= 0)
			throw new IllegalArgumentException("Minimum Supported Links Account could not be less than 0!");
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
		if(maxOHLen <= 0)
			throw new IllegalArgumentException("Maximum Overhang Length could not be less than 0!");
		this.maxOHLen = maxOHLen;
	}

	public Double getMaxOHRatio() {
		return maxOHRatio;
	}

	public void setMaxOHRatio(Double maxOHRatio) {
		if(maxOHRatio <= 0)
			throw new IllegalArgumentException("Maximum Overhang Ratio could not be less than 0!");
		this.maxOHRatio = maxOHRatio;
	}
	
	public Integer getMaxEndLen() {
		return maxEndLen;
	}

	public void setMaxEndLen(Integer maxEndLen) {
		if(maxEndLen <= 0)
			throw new IllegalArgumentException("Maximum Ending Length could not be less than 0!");
		this.maxEndLen = maxEndLen;
	}

	public Double getMaxEndRatio() {
		return maxEndRatio;
	}

	public void setMaxEndRatio(Double maxEndRatio) {
		if(maxEndRatio <= 0)
			throw new IllegalArgumentException("Maximum Ending Ratio could not be less than 0!");
		this.maxEndRatio = maxEndRatio;
	}
	
	public Double getRatio() {
		return ratio;
	}

	public void setRatio(Double ratio) {
		if(ratio <= 0)
			throw new IllegalArgumentException("Ratio could not be less than 0!");
		this.ratio = ratio;
	}	

	public int getTipLength() {
		return tipLength;
	}

	public void setTipLength(int tipLength) {
		if(tipLength <= 0)
			throw new IllegalArgumentException("Tip Length could not be less than 0!");
		this.tipLength = tipLength;
	}

	@Override
	public String toString() {
		return "Parameter [pbFile=" + pbFile + ", cntFile=" + cntFile + ", algFile=" + algFile + ", outFolder="
				+ outFolder + ", minContLen=" + minContLen + ", minPBLen=" + minPBLen + ", minOLLen=" + minOLLen
				+ ", minOLRatio=" + minOLRatio + ", maxOHLen=" + maxOHLen + ", maxOHRatio=" + maxOHRatio
				+ ", maxEndLen=" + maxEndLen + ", maxEndRatio=" + maxEndRatio + ", minSupLinks=" + minSupLinks
				+ ", maxSupLinks=" + maxSupLinks + ", type=" + type + ", identity=" + identity + ", isUseOLLink="
				+ isUseOLLink + ", ratio=" + ratio + ", isRepMask=" + isRepMask + ", isGapFilling=" + isGapFilling
				+ ", tipLength=" + tipLength + "]";
	}
	
}

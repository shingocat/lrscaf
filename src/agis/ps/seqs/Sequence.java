package agis.ps.seqs;

import java.util.ArrayList;
import java.util.List;

import agis.ps.util.MisassemblyRegion;

/**
* @author Mao Qin
* @version 2020年9月17日 上午10:13:41
* @Email mqin@outlook.com
* @Description class usage.
* @Copyright All Right Reserved 2020.
*/
public class Sequence {
	private String id;
	private Integer length;
	// defined this sequence is used or not during scaffolding;
	private Boolean isUsed = false;
	// defined this sequence is repeat or not;
	private Boolean isRepeat = false;
	// defined this sequence is misassembly
	private Boolean isMisassembly = false;
	private List<MisassemblyRegion> misassemblies = new ArrayList<MisassemblyRegion>();
	
	public String getForwardSeqs() {return null;};
	
	public String getComplementReverseSeqs() {return null;};
	
	public String getId() {
		return id;
	}
	public void setId(String id) {
		this.id = id;
	}
	public Integer getLength() {
		return length;
	}
	public void setLength(Integer length) {
		this.length = length;
	}
	public Boolean isUsed() {
		return isUsed;
	}
	public void setIsUsed(Boolean isUsed) {
		this.isUsed = isUsed;
	}
	public Boolean isRepeat() {
		return isRepeat;
	}
	public void setIsRepeat(Boolean isRepeat) {
		this.isRepeat = isRepeat;
	}
	
	public boolean isMisassembly() {
		return isMisassembly;
	}

	public void setIsMisassembly(boolean isMisassembly) {
		this.isMisassembly = isMisassembly;
	}

	public List<MisassemblyRegion> getMisassemblies() {
		return misassemblies;
	}

	public void setMisassemblies(List<MisassemblyRegion> misassemblies) {
		this.misassemblies = misassemblies;
	}
	
	public void addMisassemblyRegion(MisassemblyRegion mr) {
		this.getMisassemblies().add(mr);
	}

	@Override
	public boolean equals(Object o) {
		if (o == null)
			return false;
		Sequence seq = (Sequence) o;
		if (seq.getId().equals(this.getId()))
			return true;
		return false;
	}

	@Override
	public String toString() {
		return "Sequence [ID=" + this.getId() + ", Length=" + this.getLength() + "]";
	}

}

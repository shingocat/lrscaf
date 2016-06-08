/*
*File: agis.ps.seqs.PBRead.java
*User: mqin
*Email: mqin@ymail.com
*Date: 2016年6月2日
*/
package agis.ps.seqs;

import org.biojava.nbio.core.exceptions.CompoundNotFoundException;
import org.biojava.nbio.core.sequence.DNASequence;
import org.biojava.nbio.core.sequence.compound.DNACompoundSet;

public class PBRead extends DNASequence{

	public PBRead(String seq, DNACompoundSet dnaCompoundSet) throws CompoundNotFoundException {
		// TODO Auto-generated constructor stub
		super(seq,dnaCompoundSet);
	}

	@Override
	public String toString(){
		return "PBRead [ID: " + getAccession().getID() + ", length: " + getLength() + "]";
	}
//	@Override
//	public String toString() {
//		return "PBRead [getRNASequence()=" + getRNASequence() + ", getGCCount()=" + getGCCount() + ", getReverse()="
//				+ getReverse() + ", getComplement()=" + getComplement() + ", getReverseComplement()="
//				+ getReverseComplement() + ", getDNAType()=" + getDNAType() + ", getProxySequenceReader()="
//				+ getProxySequenceReader() + ", getBioBegin()=" + getBioBegin() + ", getBioEnd()=" + getBioEnd()
//				+ ", getUserCollection()=" + getUserCollection() + ", getAnnotationType()=" + getAnnotationType()
//				+ ", getDescription()=" + getDescription() + ", getOriginalHeader()=" + getOriginalHeader()
//				+ ", getParentSequence()=" + getParentSequence() + ", getSource()=" + getSource() + ", getNotesList()="
//				+ getNotesList() + ", getSequenceScore()=" + getSequenceScore() + ", getFeatures()=" + getFeatures()
//				+ ", getFeaturesKeyWord()=" + getFeaturesKeyWord() + ", getDatabaseReferences()="
//				+ getDatabaseReferences() + ", getFeatureRetriever()=" + getFeatureRetriever() + ", getAccession()="
//				+ getAccession() + ", getTaxonomy()=" + getTaxonomy() + ", getCompoundSet()=" + getCompoundSet()
//				+ ", toString()=" + super.toString() + ", getSequenceAsString()=" + getSequenceAsString()
//				+ ", getAsList()=" + getAsList() + ", getLength()=" + getLength() + ", iterator()=" + iterator()
//				+ ", getInverse()=" + getInverse() + ", getClass()=" + getClass() + ", hashCode()=" + hashCode() + "]";
//	}
	
}



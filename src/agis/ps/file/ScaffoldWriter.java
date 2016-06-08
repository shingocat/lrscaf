/*
*File: agis.ps.util.ScaffoldWriter.java
*User: mqin
*Email: mqin@ymail.com
*Date: 2016年2月28日
*/
package agis.ps.file;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.List;
import java.util.Map;
import java.util.Vector;

import org.biojava.nbio.core.exceptions.CompoundNotFoundException;
import org.biojava.nbio.core.sequence.DNASequence;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import agis.ps.Path;
import agis.ps.path.Node;
import agis.ps.path.NodePath;
import agis.ps.seqs.Contig;
import agis.ps.seqs.PBGapSeq;
import agis.ps.seqs.PBRead;
import agis.ps.util.Consensusser;
import agis.ps.util.GapRecord;
import agis.ps.util.Parameter;
import agis.ps.util.Strand;

public class ScaffoldWriter {

	public static Logger logger = LoggerFactory.getLogger(ScaffoldWriter.class);
	private List<NodePath> paths;
	// private Map<String, DNASequence> cnts;
	private Map<String, Contig> cnts;
	private String filePath;
	private Parameter paras;
	private List<GapRecord> gapRecords;
	private Map<String,PBRead> pbReads;

	// public ScaffoldWriter(List<NodePath> paths, Map<String, DNASequence>
	// cnts, String filePath) {
	// // TODO Auto-generated constructor stub
	// this.paths = paths;
	// this.cnts = cnts;
	// this.filePath = filePath;
	// }
	
	public ScaffoldWriter()
	{
		super();
	}
	
	public ScaffoldWriter(Parameter paras)
	{
		this.paras = paras;
	}

	public ScaffoldWriter(List<NodePath> paths, Map<String, Contig> cnts, String filePath) {
		// TODO Auto-generated constructor stub
		this.paths = paths;
		this.cnts = cnts;
		this.filePath = filePath;
	}
	
	public ScaffoldWriter(Parameter paras, List<NodePath> paths, Map<String, Contig> cnts, String filePath) {
		// TODO Auto-generated constructor stub
		this.paras = paras;
		this.paths = paths;
		this.cnts = cnts;
		this.filePath = filePath;
	}
	
	public void write2()
	{
		if (filePath == null)
			return;
		File out = null;
		FileWriter fw = null;
		BufferedWriter bw = null;
		PBReader pbReader = null;
		GapRecordReader grReader = null;
		try {
			out = new File(filePath);
			if (out.exists()) {
				logger.debug("The output file of scaffold is exist! It will not be overwrited!");
				logger.info("The output file of scaffold is exist! It will not be overwrited!");
				return;
			}
			if(!out.createNewFile())
			{
				logger.debug("ScaffoldWriter: The output file of scaffolds could not create!");
				logger.info("ScaffoldWriter: The output file of scaffolds could not create!");
				return;
			}
			if(paras.isGapFilling())
			{
				grReader = new GapRecordReader(paras);
				gapRecords = grReader.read();
				pbReader = new PBReader(paras);
				pbReads = pbReader.read();
			}
			fw = new FileWriter(out);
			bw = new BufferedWriter(fw);
			int count = 0; // scaffolds number;
			int sLen = 0; // scaffolds length;
			StringBuffer sb = null; 
//			boolean isAdded = false;
			boolean isNextUsed = false;
			Node current = null;
			String cId = null;
			String cSeq = null;
			Node next = null;
			String nId = null;
			String nSeq = null;
			for(NodePath p : paths)
			{
				sb = new StringBuffer();
				bw.write(">scaffolds_" + count);
//				logger.debug(">scaffolds_" + count); 
//				System.out.println(">scaffolds_" + count);
				bw.newLine();
				int size = p.getPathSize();
				// if the path size equal to 1, 
				if(size == 1)
				{
					current = p.getElement(0);
					cId = current.getCnt().getID();
					cSeq = "";
					if(current.getStrand().equals(Strand.FORWARD))
						cSeq = cnts.get(cId).getSequenceAsString();
					else
						cSeq = cnts.get(cId).getReverseComplementSeq();
					cSeq = cSeq.trim();
					bw.write(cSeq);
					bw.newLine();
					current = null;
					cId = null;
					cSeq = null;
					count++;
					continue;
				}
				// else the path size large than 1;
				for(int i = 0 ; i < size; i++)
				{
					if(i == p.getPathSize() - 1 && isNextUsed)
					{
						if(sb.length() != 0)
							bw.write(sb.toString());
//						logger.debug(sb.toString());
//						System.out.println(sb.toString());
						isNextUsed = false;
						sb = null;
						break;
					}
					if(i == p.getPathSize() - 1 && !isNextUsed)
					{
						current = p.getElement(i);
						cId = current.getCnt().getID();
						cSeq = "";
						if(current.getStrand().equals(Strand.FORWARD))
							cSeq = cnts.get(cId).getSequenceAsString();
						else
							cSeq = cnts.get(cId).getReverseComplementSeq();
						bw.write(cSeq);
//						logger.debug(">" + id + "\n" +seq);
//						System.out.println(">" + id + "\n" +seq);
						current = null;
						cId = null;
						cSeq = null;
						isNextUsed = false;
						break;
					}
					current = p.getElement(i);
					cId = current.getCnt().getID();
					cSeq = "";
					if(current.getStrand().equals(Strand.FORWARD))
						cSeq = cnts.get(cId).getSequenceAsString();
					else
						cSeq = cnts.get(cId).getReverseComplementSeq();
//					logger.debug(seq);
//					System.out.println(">" + id + "\n" +seq);
					int nLen = current.getMeanDist2Next();
					int sdLen = current.getSdDist2Next();
					if(nLen < 0)
					{
						next = p.getElement(i + 1);
						nId = next.getCnt().getID();
						nSeq = "";
						if(next.getStrand().equals(Strand.FORWARD))
							nSeq = cnts.get(nId).getSequenceAsString();
						else
							nSeq = cnts.get(nId).getReverseComplementSeq();
//						logger.debug(nSeq);
//						System.out.println(">" + nId + "\n" + nSeq);
						if(isNextUsed)
						{
							String temp = concatenate(sb.toString(), nSeq, nLen, sdLen);
//							System.out.println(id + "\t" + nId + "merge seq :\t" + temp);
							sb.delete(0, sb.length());
							sb.append(temp);
//							logger.debug(temp);
//							System.out.println(">" + id + "_" + nId + "\n" + temp);
							next = null;
							nId = null;
							nSeq = null;
							temp = null;
							isNextUsed = true;
						} else
						{
							String temp = concatenate(cSeq, nSeq, nLen, sdLen);
//							System.out.println(id + "\t" + nId + " merge seq :\t" + temp);
							sb.append(temp);
//							logger.debug(temp);
//							System.out.println(">" + id + "_" + nId + "\n" + temp);
							next = null;
							nId = null;
							nSeq = null;
							temp = null;
							isNextUsed = true;
						}
					} else
					{
						if(sb.length() != 0)
						{
							bw.write(sb.toString());
//							logger.debug(sb.toString());
//							System.out.println(">length" + sb.length() + "\n" + sb.toString());
							sb = null;
							sb = new StringBuffer();
						} else
						{
							bw.write(cSeq);
//							logger.debug(seq);
//							System.out.println(seq);
							cSeq = null;
						}
						// if not the last element, fill the gap according to user specified!
						if(paras.isGapFilling())
						{
							if(i != p.getPathSize() - 1)
							{
//								bw.write(repeatString("N", nLen));
								next = p.getElement(i + 1);
								nId = next.getCnt().getID();
								GapRecord gr = new GapRecord();
								gr.setStart(cId);
								gr.setEnd(nId);
								int index = gapRecords.indexOf(gr);
								if(index != -1)
								{
									List<PBGapSeq> pbgs = gapRecords.get(index).getSeqs();//gr.getSeqs();
									List<String> seqs = new Vector<String>();
									for(PBGapSeq pbg : pbgs)
									{
										String pbId = pbg.getId();
										int start = pbg.getStart();
										int end = pbg.getEnd();
										Strand strand = pbg.getStrand();
										String seq = "";
										if(strand.equals(Strand.FORWARD))
										{
											seq = pbReads.get(pbId).getSequenceAsString().substring(start, end);
										} else
										{
											String temp = pbReads.get(pbId).getReverseComplement().getSequenceAsString();
											int len = temp.length();
											int t = start;
											start = len - end;
											end = len - t;
											seq = temp.substring(start, end);
										}
										seqs.add(seq);
									}
									Consensusser cs = new Consensusser();
									String value = cs.getConsensus(seqs);
									bw.write(value);
								} else
								{
									bw.write(repeatString("N", nLen));
								}
							}
						} else
						{
							if(i != p.getPathSize() -1 )
								bw.write(repeatString("N", nLen));
						}

//						logger.debug(repeatString("N", nLen));
//						System.out.println(repeatString("N", nLen));
						isNextUsed = false;
					}
				}
				isNextUsed = false;
				bw.newLine();
				count++;
				sb = null;
				System.gc(); 
			}
		} catch (IllegalArgumentException e) {
			logger.debug(this.getClass().getName() + "\t" + e.getMessage() + "\t" + "IllegalArgumentException.");
			logger.error(this.getClass().getName() + "\t" + e.getMessage() + "\t" + "IllegalArgumentException.");
		} catch (IOException e) {
			logger.debug(this.getClass().getName() + "\t" + e.getMessage() + "\t" + "IOException.");
			logger.error(this.getClass().getName() + "\t" + e.getMessage() + "\t" + "IOException.");
		} catch (ClassCastException e){
			logger.debug(this.getClass().getName() + "\t" + e.getMessage() + "\t" + "ClassCastException.");
			logger.error(this.getClass().getName() + "\t" + e.getMessage() + "\t" + "ClassCastException.");
		} catch (NullPointerException e){
			logger.debug(this.getClass().getName() + "\t" + e.getMessage() + "\t" + "NullPointerException.");
			logger.error(this.getClass().getName() + "\t" + e.getMessage() + "\t" + "NullPointerException.");
		} catch (Exception e) {
			logger.debug(this.getClass().getName() + "\t" + e.getMessage() + "\t" + "Exception.");
			logger.error(this.getClass().getName() + "\t" + e.getMessage() + "\t" + "Exception.");
		} finally {
			if (bw != null)
				try {
					bw.close();
				} catch (IOException e) {
					logger.debug(this.getClass().getName() + "\t" + e.getMessage());
					logger.error(this.getClass().getName() + "\t" + e.getMessage());
				}
		}
		
	}

	public void write() {
		// TODO Auto-generated method stub
		if (filePath == null)
			return;
		File out = null;
		FileWriter fw = null;
		BufferedWriter bw = null;
		PBReader pbReader = null;
		GapRecordReader grReader = null;
		try {
			out = new File(filePath);
			if (out.exists()) {
				logger.debug("The output file of scaffold is exist! It will not be overwrited!");
				logger.info("The output file of scaffold is exist! It will not be overwrited!");
				return;
			}
			if(!out.createNewFile())
			{
				logger.debug("ScaffoldWriter: The output file of scaffolds could not create!");
				logger.info("ScaffoldWriter: The output file of scaffolds could not create!");
				return;
			}
			grReader = new GapRecordReader(paras);
			List<GapRecord> gapRecords = grReader.read();
			pbReader = new PBReader(paras);
			Map<String, PBRead> pbReads = pbReader.read();
			fw = new FileWriter(out);
			bw = new BufferedWriter(fw);
			int count = 0; // scaffolds number;
			int sLen = 0; // scaffolds length;
			StringBuffer sb = null; 
//			boolean isAdded = false;
			boolean isNextUsed = false;
			Node node = null;
			String id = null;
			String seq = null;
			Node nNode = null;
			String nId = null;
			String nSeq = null;
			for (NodePath p : paths)
			{
//				System.out.println(count + " path");
				sb = new StringBuffer();
				bw.write(">scaffolds_" + count);
//				logger.debug(">scaffolds_" + count); 
//				System.out.println(">scaffolds_" + count);
				bw.newLine();
				// if the path size equal to 1, 
				if(p.getPathSize() == 1)
				{
					node = p.getElement(0);
					id = node.getCnt().getID();
					seq = "";
					if(node.getStrand().equals(Strand.FORWARD))
						seq = cnts.get(id).getSequenceAsString();
					else
						seq = cnts.get(id).getReverseComplementSeq();
					seq = seq.trim();
					bw.write(seq);
					bw.newLine();
//					logger.debug(">" + id + "\n" +seq);
//					System.out.println(">" + id + "\n" +seq);
					node = null;
					id = null;
					seq = null;
					count++;
					continue;
				}
				// else the path size large than 1;
				for(int i = 0 ; i < p.getPathSize(); i++)
				{
					if(i == p.getPathSize() - 1 && isNextUsed)
					{
						if(sb.length() != 0)
							bw.write(sb.toString());
//						logger.debug(sb.toString());
//						System.out.println(sb.toString());
						isNextUsed = false;
						sb = null;
						break;
					}
					if(i == p.getPathSize() - 1 && !isNextUsed)
					{
						node = p.getElement(i);
						id = node.getCnt().getID();
						seq = "";
						if(node.getStrand().equals(Strand.FORWARD))
							seq = cnts.get(id).getSequenceAsString();
						else
							seq = cnts.get(id).getReverseComplementSeq();
						bw.write(seq);
//						logger.debug(">" + id + "\n" +seq);
//						System.out.println(">" + id + "\n" +seq);
						node = null;
						id = null;
						seq = null;
						isNextUsed = false;
						break;
					}
					node = p.getElement(i);
					id = node.getCnt().getID();
					seq = "";
					if(node.getStrand().equals(Strand.FORWARD))
						seq = cnts.get(id).getSequenceAsString();
					else
						seq = cnts.get(id).getReverseComplementSeq();
//					logger.debug(seq);
//					System.out.println(">" + id + "\n" +seq);
					int nLen = node.getMeanDist2Next();
					int sdLen = node.getSdDist2Next();
					if(nLen < 0)
					{
						nNode = p.getElement(i + 1);
						nId = nNode.getCnt().getID();
						nSeq = "";
						if(nNode.getStrand().equals(Strand.FORWARD))
							nSeq = cnts.get(nId).getSequenceAsString();
						else
							nSeq = cnts.get(nId).getReverseComplementSeq();
//						logger.debug(nSeq);
//						System.out.println(">" + nId + "\n" + nSeq);
						if(isNextUsed)
						{
							String temp = concatenate(sb.toString(), nSeq, nLen, sdLen);
//							System.out.println(id + "\t" + nId + "merge seq :\t" + temp);
							sb.delete(0, sb.length());
							sb.append(temp);
//							logger.debug(temp);
//							System.out.println(">" + id + "_" + nId + "\n" + temp);
							nNode = null;
							nId = null;
							nSeq = null;
							temp = null;
							isNextUsed = true;
						} else
						{
							String temp = concatenate(seq, nSeq, nLen, sdLen);
//							System.out.println(id + "\t" + nId + " merge seq :\t" + temp);
							sb.append(temp);
//							logger.debug(temp);
//							System.out.println(">" + id + "_" + nId + "\n" + temp);
							nNode = null;
							nId = null;
							nSeq = null;
							temp = null;
							isNextUsed = true;
						}
					} else
					{
						if(sb.length() != 0)
						{
							bw.write(sb.toString());
//							logger.debug(sb.toString());
//							System.out.println(">length" + sb.length() + "\n" + sb.toString());
							sb = null;
							sb = new StringBuffer();
						} else
						{
							bw.write(seq);
//							logger.debug(seq);
//							System.out.println(seq);
							seq = null;
						}
						// if not the last element, write the N
						if(i != p.getPathSize() - 1)
							bw.write(repeatString("N", nLen));
//						logger.debug(repeatString("N", nLen));
//						System.out.println(repeatString("N", nLen));
						isNextUsed = false;
					}
				}
				isNextUsed = false;
				bw.newLine();
				count++;
				sb = null;
				System.gc(); 
			}
//			for (NodePath p : paths) {
//				sb = new StringBuilder();
//				bw.write(">scaffolds_" + count + "\n");
//				for (int i = 0; i < p.getPathSize(); i++) {
//					if(i == p.getPathSize() - 1 && isAdded)
//						continue;
//					Node node = p.getElement(i);
//					String id = node.getCnt().getID();
//					String seq = "";
//					if(node.getStrand().equals(Strand.FORWARD))
//						seq = cnts.get(id).getSequenceAsString();
//					else 
//						seq = cnts.get(id).getReverseComplementSeq();
//					int nLen = node.getMeanDist2Next();
//					int sdLen = node.getSdDist2Next();
//					if(nLen < 0){
//						Node nNode = p.getElement(i + 1);
//						String nId = nNode.getCnt().getID();
//						String nSeq = "";
//						if(nNode.getStrand().equals(Strand.FORWARD))
//							nSeq = cnts.get(nId).getSequenceAsString();
//						else
//							nSeq = cnts.get(nId).getReverseComplementSeq();
//						if(isAdded){
//							String temp = concatenate(sb.toString(), nSeq, nLen, sdLen);
//							sb.append(temp);
//							isAdded = true;
//						} else {
//							String temp = concatenate(seq, nSeq, nLen, sdLen);
//							sb.append(temp);
//							isAdded = true;
//						}
//					} else {
//						if(!isAdded)
//							sb.append(seq);
//						if(i != p.getPathSize() - 1)
//							sb.append(repeatString("N", nLen));
//						isAdded = false;
//					}
////					if (i != p.getPathSize() - 1) {
////						if (nLen < 0){
////							int len = nLen + sdLen;
////							// get next element;
////							Node nNode = p.getElement(i + 1);
////							String nId = nNode.getCnt().getID();
////							String nSeq = "";
////							if(nNode.getStrand().equals(Strand.FORWARD))
////								nSeq = cnts.get(nId).getSequenceAsString();
////							else
////								nSeq = cnts.get(id).getReverseComplementSeq();
////							String cSeq = concatenate(seq, nSeq, nLen, sdLen);
////							bw.write(cSeq);
//////							bw.write(repeatString("M", nLen));
////						} else{
////							bw.write(seq);
////							bw.write(repeatString("N", nLen));
////						}
////					}
//				}
//				bw.write(sb.toString());
//				bw.write("\n");
//				isAdded = false;
//				count++;
//			}

		} catch (IllegalArgumentException e) {
			logger.debug(this.getClass().getName() + "\t" + e.getMessage() + "\t" + "IllegalArgumentException.");
			logger.error(this.getClass().getName() + "\t" + e.getMessage() + "\t" + "IllegalArgumentException.");
		} catch (IOException e) {
			logger.debug(this.getClass().getName() + "\t" + e.getMessage() + "\t" + "IOException.");
			logger.error(this.getClass().getName() + "\t" + e.getMessage() + "\t" + "IOException.");
		} catch (ClassCastException e){
			logger.debug(this.getClass().getName() + "\t" + e.getMessage() + "\t" + "ClassCastException.");
			logger.error(this.getClass().getName() + "\t" + e.getMessage() + "\t" + "ClassCastException.");
		} catch (NullPointerException e){
			logger.debug(this.getClass().getName() + "\t" + e.getMessage() + "\t" + "NullPointerException.");
			logger.error(this.getClass().getName() + "\t" + e.getMessage() + "\t" + "NullPointerException.");
		} catch (Exception e) {
			logger.debug(this.getClass().getName() + "\t" + e.getMessage() + "\t" + "Exception.");
			logger.error(this.getClass().getName() + "\t" + e.getMessage() + "\t" + "Exception.");
		} finally {
			if (bw != null)
				try {
					bw.close();
				} catch (IOException e) {
					logger.debug(this.getClass().getName() + "\t" + e.getMessage());
					logger.error(this.getClass().getName() + "\t" + e.getMessage());
				}
		}

	}

	// if times larger than or equal to 0, used N indicated
	// else times less than 0, used M indicated;
	private String repeatString(String str, int times) {
		if (str == null)
			throw new IllegalArgumentException("ScaffoldWriter: The string could not be null!");
		// if(times < 0)
		// throw new IllegalArgumentException("ScaffoldWriter: The repeat times
		// could not be negative!");
//		StringBuffer sb = new StringBuffer(str);
//		if (times < 0) {
//			times = 0 - times;
//			for (int i = 0; i < times; i++)
//				sb.append(str);
//		} else {
//			for (int i = 1; i < times; i++)
//				sb.append(str);
//		}
//		return sb.toString();
		return new String(new char[times]).replace("\0", str);
	}
	
	public String concatenate(String seq1, String seq2, int len, int sd)
	{
		// 99% region
		int range = Math.abs(len) + 3 * Math.abs(sd);
		String t1 = null;
		String t2 = null;
		Consensusser cs = new Consensusser();
		List<String> seqs = new Vector<String>(2);
		// four different status;
		if(seq1.length() <= range)
		{
			t1 = seq1;
			if(seq2.length() <= range)
			{
				// seq1 and seq2 less than range;
				// return consensus directly;
				t2 = seq2;
				seqs.add(t1);
				seqs.add(t2);
				String value;
				try {
					value = cs.getConsensus(seqs);
				} catch (CompoundNotFoundException e) {
					// TODO Auto-generated catch block
					logger.debug(this.getClass().getName() + "\t" + e.getMessage());
					logger.error(this.getClass().getName() + "\t" + e.getMessage());
					if(t1.length() > t2.length())
						value = t1;
					else 
						value = t2;
				}
//				String value = cs.getConsensus(t1, t2, "nw");
				return value;
			} else 
			{
				// seq1 less than range, but seq2 large than range;
				// return consensus + seq2 remainder;
				t2 = seq2.substring(0, range);
				seqs.add(t1);
				seqs.add(t2);
				String value;
				try {
					value = cs.getConsensus(seqs);
				} catch (CompoundNotFoundException e) {
					logger.debug(this.getClass().getName() + "\t" + e.getMessage());
					logger.error(this.getClass().getName() + "\t" + e.getMessage());
					if(t1.length() > t2.length())
						value = t1;
					else 
						value = t2;
				}
//				String value = cs.getConsensus(t1, t2, "nw");
				StringBuffer sb = new StringBuffer();
				sb.append(value);
				sb.append(seq2.substring(range));
//				logger.debug(sb.toString());
				return sb.toString();
			}
		} else {
			t1 = seq1.substring(seq1.length() - range);
			if(seq2.length() <= range)
			{
				// seq1 large than range, but seq2 less than range;
				// return seq1 remainder + consensus;
				t2 = seq2;
				seqs.add(t1);
				seqs.add(t2);
				String value;
				try {
					value = cs.getConsensus(seqs);
				} catch (CompoundNotFoundException e) {
					logger.debug(this.getClass().getName() + "\t" + e.getMessage());
					logger.error(this.getClass().getName() + "\t" + e.getMessage());
					if(t1.length() > t2.length())
						value = t1;
					else 
						value = t2;
				}
//				String value = cs.getConsensus(t1, t2, "nw");
				StringBuffer sb = new StringBuffer();
				sb.append(seq1.substring(0, seq1.length() - range));
				sb.append(value);
//				logger.debug(sb.toString());
				return sb.toString();
			} else 
			{
				// seq1 and seq2 large than range;
				// return seq1 remainder + consensus + seq2 remainder;
				t2 = seq2.substring(0, range);
				seqs.add(t1);
				seqs.add(t2);
				String value;
				try {
					value = cs.getConsensus(seqs);
				} catch (CompoundNotFoundException e) {
					logger.debug(this.getClass().getName() + "\t" + e.getMessage());
					logger.error(this.getClass().getName() + "\t" + e.getMessage());
					if(t1.length() > t2.length())
						value = t1;
					else 
						value = t2;
				}
//				String value = cs.getConsensus(t1, t2, "nw");
				StringBuffer sb = new StringBuffer();
				sb.append(seq1.substring(0, seq1.length() - range));
				sb.append(value);
				sb.append(seq2.substring(range));
//				logger.debug(sb.toString());
				return sb.toString();
			}
		}
		
//		Consensusser cs = new Consensusser();
//		String value = cs.getConsensus(t1, t2, "nw");
//		cs = null;
//		StringBuffer sb = new StringBuffer();
//		sb.append(seq1.substring(0, seq1.length() - range));
//		sb.append(value);
//		sb.append(seq2.substring(range));
//		return sb.toString();
	}
}

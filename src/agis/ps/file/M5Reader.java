/*
*File: agis.ps.util.M5Reader.java
*User: mqin
*Email: mqin@ymail.com
*Date: 2015年12月28日
*/
package agis.ps.file;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.LineNumberReader;
import java.util.List;
import java.util.Vector;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import agis.ps.link.M5FileEncapsulate;
import agis.ps.link.M5Record;
import agis.ps.link.MRecord;
import agis.ps.util.Parameter;
import agis.ps.util.Strand;

public class M5Reader extends AlignmentFileReader{
	
	private final static Logger logger = LoggerFactory.getLogger(M5Reader.class);
	
	private String path;
	private List<MRecord> m5List;
	private M5FileEncapsulate m5FileEncapsulate;
	private Parameter paras;
	
	
	public M5Reader(Parameter paras)
	{
		super(paras);
	}
	
	// original reading method with some filtered parameters;
//	public List<MRecord> read()
//	{
//		long start = System.currentTimeMillis();
//		FileReader fr = null;
//		BufferedReader br = null;
//		int count = 0;
//		try
//		{
//			File m5File = new File(path);
//			if(!m5File.exists())
//			{
//				logger.debug(this.getClass().getName() + "The m5 file do not exist!");
//				logger.error(this.getClass().getName() + "The m5 file do not exist!");
//				return null;
//			}
//			
//			fr = new FileReader(m5File);
//			br = new BufferedReader(fr);
//			String line = null;
//			String [] arrs = null;
//			m5List = new Vector<MRecord>();
//			while((line = br.readLine()) != null)
//			{
//				line = line.trim();
//				line = line.replaceAll(System.getProperty("line.separator"), "");
//				arrs = line.split("\\s+");
//				if(arrs[0].equalsIgnoreCase("qName") && arrs[1].equalsIgnoreCase("qLength"))
//					continue;
//				count++;
//				// if the first line is header
////				if(arrs[0].equalsIgnoreCase("qName") && arrs[1].equalsIgnoreCase("qLength"))
////					continue;
//				// if less than minimum pacbio length
//				if(Integer.valueOf(arrs[1]) < minPBLen)
//					continue;
//				// if less tan minimum contig length
//				if(Integer.valueOf(arrs[6]) < minCNTLen)
//					continue;
//				// if the identity less than specified value;
//				double sum = Double.valueOf(arrs[11]) + Double.valueOf(arrs[12]) + Double.valueOf(arrs[13]) + Double.valueOf(arrs[14]);
//				double value = Double.valueOf(arrs[11]) / sum;
//				if(value < identity)
//					continue;
////				int sum = this.getNumMatch() + this.getNumMismatch() + this.getNumIns() + this.getNumDel();
////				double value = (double)this.getNumMatch() / sum;
//				M5Record m5 = new M5Record();
//				m5.setqName(arrs[0]);
//				m5.setqLength(Integer.valueOf(arrs[1]));
//				m5.setqStart(Integer.valueOf(arrs[2]));
//				m5.setqEnd(Integer.valueOf(arrs[3]));
//				m5.setqStrand(arrs[4].equals("+") ? Strand.FORWARD : Strand.REVERSE);
//				m5.settName(arrs[5]);
//				m5.settLength(Integer.valueOf(arrs[6]));
//				m5.settStart(Integer.valueOf(arrs[7]));
//				m5.settEnd(Integer.valueOf(arrs[8]));
//				m5.settStrand(arrs[9].equals("+") ? Strand.FORWARD : Strand.REVERSE);
//				m5.setScore(Integer.valueOf(arrs[10]));
//				m5.setNumMatch(Integer.valueOf(arrs[11]));
//				m5.setNumMismatch(Integer.valueOf(arrs[12]));
//				m5.setNumIns(Integer.valueOf(arrs[13]));
//				m5.setNumDel(Integer.valueOf(arrs[14]));
//				m5.setMapQV(Integer.valueOf(arrs[15]));
////				m5.setqAlignedSeq(arrs[16]);
////				m5.setMatchPattern(arrs[17]);
////				m5.settAlignedSeq(arrs[18]);
//				// setting identity
//				m5.setIdentity(value);
//				m5List.add(m5);
//			}
//			br.close();
//		} catch(ArrayIndexOutOfBoundsException e)
//		{
//			logger.debug(this.getClass().getName() + e.getMessage() + "\t" + e.getClass().getName());
//			logger.error(this.getClass().getName() + e.getMessage() + "\t" + e.getClass().getName());
//			m5List = null;
//		} catch(FileNotFoundException e)
//		{
//			logger.debug(this.getClass().getName() + e.getMessage() + "\t" + e.getClass().getName());
//			logger.error(this.getClass().getName() + e.getMessage() + "\t" + e.getClass().getName());
//			m5List = null;
//		} catch(IOException e)
//		{
//			logger.debug(this.getClass().getName() + e.getMessage() + "\t" + e.getClass().getName());
//			logger.error(this.getClass().getName() + e.getMessage() + "\t" + e.getClass().getName());
//			m5List = null;
//		} catch(Exception e){
//			logger.debug(this.getClass().getName() + e.getMessage() + "\t" + e.getClass().getName());
//			logger.error(this.getClass().getName() + e.getMessage() + "\t" + e.getClass().getName());
//			m5List = null;
//		} finally
//		{
//			try
//			{
//				if( br != null)
//					br.close();
//			} catch(IOException e)
//			{
//				logger.debug(this.getClass().getName() + e.getMessage() + "\t" + e.getClass().getName());
//				logger.error(this.getClass().getName() + e.getMessage() + "\t" + e.getClass().getName());
//				m5List = null;
//			}
//		}
//		long end = System.currentTimeMillis();
//		logger.info(this.getClass().getName() + "\tAligned records: " + count);
//		logger.info(this.getClass().getName() + "\tRecords after filtered: " + m5List.size());
//		logger.info("Reading M5 Records, erase time: " + (end - start) + " ms");
//		return m5List;
//	}
	
	// reading file and build the M5FileEncapsulate by using LineNumberReader
//	public M5FileEncapsulate readByLineNumberReader()
//	{
//		long start = System.currentTimeMillis();
//		File file = null;
//		FileReader fr = null;
//		BufferedReader br = null;
//		int size = 0; 
//		try
//		{
//			size = this.countLineNumber();
//			file = new File(path);
//			fr = new FileReader(file);
//			br = new BufferedReader(fr);
//			// storing the data into List container;
//			String line = "";
//			String [] arrs;
//			boolean isInit = false;
//			while((line = br.readLine()) != null)
//			{
//				line = line.trim();
//				line = line.replaceAll(System.getProperty("line.separator"), "");
//				arrs = line.split("\\s+");
//				if(arrs[0].equalsIgnoreCase("qName") && arrs[1].equalsIgnoreCase("qLength"))
//				{
//					size--;
//					continue;
//				} else
//				{
//					if(!isInit)
//					{
//						m5FileEncapsulate = new M5FileEncapsulate(size, paras);
//						isInit = true;
//					}
//					M5Record m5 = this.initM5Record(arrs);
//					m5FileEncapsulate.addRecord(m5);
//				}
//			}
//			br.close();		
//		} catch (IOException e)
//		{
//			logger.error(this.getClass().getName() + "\t" + e.getMessage());
//		} finally 
//		{
//			try{
//				if(br != null)
//					br.close();
//			} catch(IOException e)
//			{
//				logger.error(this.getClass().getName() + "\t" + e.getMessage());
//			}
//		}
//		long end = System.currentTimeMillis();
//		logger.info(this.getClass().getName() + "\tAligned records: " + size);
//		logger.info("Reading M5 Records, erase time: " + (end - start) + " ms");
//		return m5FileEncapsulate;
//	}
//	
//	/**
//	 * A method to count the file line number
//	 * it is faster than by using array copy to build aligned file encapsulate
//	 * @return
//	 */
//	private int countLineNumber()
//	{
//		File file = null;
//		FileReader fr = null;
//		LineNumberReader lnr = null;
//		int size = 0;
//		try{
//			file = new File(path);
//			fr = new FileReader(file);
//			// get the file line;
//			lnr = new LineNumberReader(fr);
//			lnr.skip(Long.MAX_VALUE);
//			size = lnr.getLineNumber();
//			lnr.close();
//		} catch(IOException e)
//		{
//			logger.error(this.getClass().getName() + "\t" + e.getMessage());
//		} finally
//		{
//			try{
//				if(lnr != null)
//					lnr.close();
//			} catch(IOException e)
//			{
//				logger.error(this.getClass().getName() + "\t" + e.getMessage());
//			}
//		}
//		return size;
//	}
	
//	private M5Record initM5Record(String [] arrs)
//	{
//		M5Record m5 = new M5Record();
//		m5.setqName(arrs[0]);
//		m5.setqLength(Integer.valueOf(arrs[1]));
//		m5.setqStart(Integer.valueOf(arrs[2]));
//		m5.setqEnd(Integer.valueOf(arrs[3]));
//		m5.setqStrand(arrs[4].equals("+") ? Strand.FORWARD : Strand.REVERSE);
//		m5.settName(arrs[5]);
//		m5.settLength(Integer.valueOf(arrs[6]));
//		m5.settStart(Integer.valueOf(arrs[7]));
//		m5.settEnd(Integer.valueOf(arrs[8]));
//		m5.settStrand(arrs[9].equals("+") ? Strand.FORWARD : Strand.REVERSE);
//		m5.setScore(Integer.valueOf(arrs[10]));
//		m5.setNumMatch(Integer.valueOf(arrs[11]));
//		m5.setNumMismatch(Integer.valueOf(arrs[12]));
//		m5.setNumIns(Integer.valueOf(arrs[13]));
//		m5.setNumDel(Integer.valueOf(arrs[14]));
//		m5.setMapQV(Integer.valueOf(arrs[15]));
////		m5.setqAlignedSeq(arrs[16]);
////		m5.setMatchPattern(arrs[17]);
////		m5.settAlignedSeq(arrs[18]);
//		return m5;
//	}

	@Override
	protected MRecord initMRecord(String[] arrs) {
		MRecord m = new MRecord();
		m.setqName(arrs[0]);
		m.setqLength(arrs[1]);
		m.setqStart(arrs[2]);
		m.setqEnd(arrs[3]);
		m.setqStrand(Strand.FORWARD);
		m.settName(arrs[5]);
		m.settLength(arrs[6]);
		m.settStart(arrs[7]);
		m.settEnd(arrs[8]);
		m.settStrand(arrs[9].equals("+") ? Strand.FORWARD : Strand.REVERSE);
//		m.setScore(Integer.valueOf(arrs[10]));
		int match = Integer.valueOf(arrs[11]);
		int mismatch = Integer.valueOf(arrs[12]);
		int insert = Integer.valueOf(arrs[13]);
		int deletion = Integer.valueOf(arrs[14]);
		double identity = (double)match / (match + mismatch + insert + deletion);
		identity = Double.valueOf(String.format("%.4f", identity));
		m.setIdentity(identity);
//		m.setMapQV(Integer.valueOf(arrs[15]));
		return m;
	}
	
}



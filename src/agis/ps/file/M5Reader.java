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
import java.util.ArrayList;
import java.util.List;
import java.util.Vector;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import agis.ps.link.M5Record;
import agis.ps.link.MRecord;
import agis.ps.util.Parameter;
import agis.ps.util.Strand;

public class M5Reader {
	
	private final static Logger logger = LoggerFactory.getLogger(M5Reader.class);
	
	private Parameter paras;
	private String path;
	private List<MRecord> m5List;
	private int minPBLen = 5000; // minimum pacbio read length
	private int minCNTLen = 3000; // minimum contig length;
	private double identity = 0.8;
	
	public M5Reader(String path)
	{
		this.path = path;
	}
	
	public M5Reader(String path, int minPBLen)
	{
		this.path = path;
		this.minPBLen = minPBLen;
	}
	
	public M5Reader(String path, int minPBLen, int minCNTLen)
	{
		this.path = path;
		this.minPBLen = minPBLen;
		this.minCNTLen = minCNTLen;
	}
	
	public M5Reader(Parameter paras)
	{
		this.paras = paras;
		this.path = paras.getAlgFile();
		this.minPBLen = paras.getMinPBLen();
		this.minCNTLen = paras.getMinContLen();
		this.identity = paras.getIdentity();
	}
	
	// reading M5Record without filtering parameters;
	public List<MRecord> readWithoutFiltering()
	{
		long start = System.currentTimeMillis();
		FileReader fr = null;
		BufferedReader br = null;
		int count = 0;
		try
		{
			File m5File = new File(path);
			if(!m5File.exists())
			{
				logger.debug(this.getClass().getName() + "The m5 file do not exist!");
				logger.error(this.getClass().getName() + "The m5 file do not exist!");
				return null;
			}
			
			fr = new FileReader(m5File);
			br = new BufferedReader(fr);
			String line = null;
			String [] arrs = null;
			m5List = new Vector<MRecord>();
			while((line = br.readLine()) != null)
			{
				line = line.trim();
				line = line.replaceAll(System.getProperty("line.separator"), "");
				arrs = line.split("\\s+");
				if(arrs[0].equalsIgnoreCase("qName") && arrs[1].equalsIgnoreCase("qLength"))
					continue;
				count++;
				// if the identity less than specified value;
				double sum = Double.valueOf(arrs[11]) + Double.valueOf(arrs[12]) + Double.valueOf(arrs[13]) + Double.valueOf(arrs[14]);
				double value = Double.valueOf(arrs[11]) / sum;
				M5Record m5 = new M5Record();
				m5.setqName(arrs[0]);
				m5.setqLength(Integer.valueOf(arrs[1]));
				m5.setqStart(Integer.valueOf(arrs[2]));
				m5.setqEnd(Integer.valueOf(arrs[3]));
				m5.setqStrand(arrs[4].equals("+") ? Strand.FORWARD : Strand.REVERSE);
				m5.settName(arrs[5]);
				m5.settLength(Integer.valueOf(arrs[6]));
				m5.settStart(Integer.valueOf(arrs[7]));
				m5.settEnd(Integer.valueOf(arrs[8]));
				m5.settStrand(arrs[9].equals("+") ? Strand.FORWARD : Strand.REVERSE);
				m5.setScore(Integer.valueOf(arrs[10]));
				m5.setNumMatch(Integer.valueOf(arrs[11]));
				m5.setNumMismatch(Integer.valueOf(arrs[12]));
				m5.setNumIns(Integer.valueOf(arrs[13]));
				m5.setNumDel(Integer.valueOf(arrs[14]));
				m5.setMapQV(Integer.valueOf(arrs[15]));
//				m5.setqAlignedSeq(arrs[16]);
//				m5.setMatchPattern(arrs[17]);
//				m5.settAlignedSeq(arrs[18]);
				// setting identity
				m5.setIdentity(value);
				m5List.add(m5);
			}
			br.close();
		} catch(ArrayIndexOutOfBoundsException e)
		{
			logger.debug(this.getClass().getName() + e.getMessage() + "\t" + e.getClass().getName());
			logger.error(this.getClass().getName() + e.getMessage() + "\t" + e.getClass().getName());
			m5List = null;
		} catch(FileNotFoundException e)
		{
			logger.debug(this.getClass().getName() + e.getMessage() + "\t" + e.getClass().getName());
			logger.error(this.getClass().getName() + e.getMessage() + "\t" + e.getClass().getName());
			m5List = null;
		} catch(IOException e)
		{
			logger.debug(this.getClass().getName() + e.getMessage() + "\t" + e.getClass().getName());
			logger.error(this.getClass().getName() + e.getMessage() + "\t" + e.getClass().getName());
			m5List = null;
		} catch(Exception e){
			logger.debug(this.getClass().getName() + e.getMessage() + "\t" + e.getClass().getName());
			logger.error(this.getClass().getName() + e.getMessage() + "\t" + e.getClass().getName());
			m5List = null;
		} finally
		{
			try
			{
				if( br != null)
					br.close();
			} catch(IOException e)
			{
				logger.debug(this.getClass().getName() + e.getMessage() + "\t" + e.getClass().getName());
				logger.error(this.getClass().getName() + e.getMessage() + "\t" + e.getClass().getName());
				m5List = null;
			}
		}
		long end = System.currentTimeMillis();
		logger.info(this.getClass().getName() + "\tAligned records: " + count);
		logger.info(this.getClass().getName() + "\tRecords after filtered: " + m5List.size());
		logger.info("Reading M5 Records, erase time: " + (end - start) + " ms");
		return m5List;
	}
	
	// original reading method with some filtered parameters;
	public List<MRecord> read()
	{
		long start = System.currentTimeMillis();
		FileReader fr = null;
		BufferedReader br = null;
		int count = 0;
		try
		{
			File m5File = new File(path);
			if(!m5File.exists())
			{
				logger.debug(this.getClass().getName() + "The m5 file do not exist!");
				logger.error(this.getClass().getName() + "The m5 file do not exist!");
				return null;
			}
			
			fr = new FileReader(m5File);
			br = new BufferedReader(fr);
			String line = null;
			String [] arrs = null;
			m5List = new Vector<MRecord>();
			while((line = br.readLine()) != null)
			{
				line = line.trim();
				line = line.replaceAll(System.getProperty("line.separator"), "");
				arrs = line.split("\\s+");
				if(arrs[0].equalsIgnoreCase("qName") && arrs[1].equalsIgnoreCase("qLength"))
					continue;
				count++;
				// if the first line is header
//				if(arrs[0].equalsIgnoreCase("qName") && arrs[1].equalsIgnoreCase("qLength"))
//					continue;
				// if less than minimum pacbio length
				if(Integer.valueOf(arrs[1]) < minPBLen)
					continue;
				// if less tan minimum contig length
				if(Integer.valueOf(arrs[6]) < minCNTLen)
					continue;
				// if the identity less than specified value;
				double sum = Double.valueOf(arrs[11]) + Double.valueOf(arrs[12]) + Double.valueOf(arrs[13]) + Double.valueOf(arrs[14]);
				double value = Double.valueOf(arrs[11]) / sum;
				if(value < identity)
					continue;
//				int sum = this.getNumMatch() + this.getNumMismatch() + this.getNumIns() + this.getNumDel();
//				double value = (double)this.getNumMatch() / sum;
				M5Record m5 = new M5Record();
				m5.setqName(arrs[0]);
				m5.setqLength(Integer.valueOf(arrs[1]));
				m5.setqStart(Integer.valueOf(arrs[2]));
				m5.setqEnd(Integer.valueOf(arrs[3]));
				m5.setqStrand(arrs[4].equals("+") ? Strand.FORWARD : Strand.REVERSE);
				m5.settName(arrs[5]);
				m5.settLength(Integer.valueOf(arrs[6]));
				m5.settStart(Integer.valueOf(arrs[7]));
				m5.settEnd(Integer.valueOf(arrs[8]));
				m5.settStrand(arrs[9].equals("+") ? Strand.FORWARD : Strand.REVERSE);
				m5.setScore(Integer.valueOf(arrs[10]));
				m5.setNumMatch(Integer.valueOf(arrs[11]));
				m5.setNumMismatch(Integer.valueOf(arrs[12]));
				m5.setNumIns(Integer.valueOf(arrs[13]));
				m5.setNumDel(Integer.valueOf(arrs[14]));
				m5.setMapQV(Integer.valueOf(arrs[15]));
//				m5.setqAlignedSeq(arrs[16]);
//				m5.setMatchPattern(arrs[17]);
//				m5.settAlignedSeq(arrs[18]);
				// setting identity
				m5.setIdentity(value);
				m5List.add(m5);
			}
			br.close();
		} catch(ArrayIndexOutOfBoundsException e)
		{
			logger.debug(this.getClass().getName() + e.getMessage() + "\t" + e.getClass().getName());
			logger.error(this.getClass().getName() + e.getMessage() + "\t" + e.getClass().getName());
			m5List = null;
		} catch(FileNotFoundException e)
		{
			logger.debug(this.getClass().getName() + e.getMessage() + "\t" + e.getClass().getName());
			logger.error(this.getClass().getName() + e.getMessage() + "\t" + e.getClass().getName());
			m5List = null;
		} catch(IOException e)
		{
			logger.debug(this.getClass().getName() + e.getMessage() + "\t" + e.getClass().getName());
			logger.error(this.getClass().getName() + e.getMessage() + "\t" + e.getClass().getName());
			m5List = null;
		} catch(Exception e){
			logger.debug(this.getClass().getName() + e.getMessage() + "\t" + e.getClass().getName());
			logger.error(this.getClass().getName() + e.getMessage() + "\t" + e.getClass().getName());
			m5List = null;
		} finally
		{
			try
			{
				if( br != null)
					br.close();
			} catch(IOException e)
			{
				logger.debug(this.getClass().getName() + e.getMessage() + "\t" + e.getClass().getName());
				logger.error(this.getClass().getName() + e.getMessage() + "\t" + e.getClass().getName());
				m5List = null;
			}
		}
		long end = System.currentTimeMillis();
		logger.info(this.getClass().getName() + "\tAligned records: " + count);
		logger.info(this.getClass().getName() + "\tRecords after filtered: " + m5List.size());
		logger.info("Reading M5 Records, erase time: " + (end - start) + " ms");
		return m5List;
	}
}



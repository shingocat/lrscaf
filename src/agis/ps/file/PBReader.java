/*
*File: agis.ps.file.PBReader.java
*User: mqin
*Email: mqin@ymail.com
*Date: 2016年6月2日
*/
package agis.ps.file;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.HashMap;
//import java.util.List;
import java.util.Map;

import org.biojava.nbio.core.sequence.AccessionID;
import org.biojava.nbio.core.sequence.compound.DNACompoundSet;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import agis.ps.seqs.PBRead;
import agis.ps.util.Parameter;

public class PBReader {
	private static Logger logger = LoggerFactory.getLogger(PBReader.class);
	private Parameter paras;
	private Map<String, PBRead> reads;
	
	public PBReader(Parameter paras)
	{
		this.paras = paras;
	}
	
	public Map<String, PBRead> read(){
		if(reads == null)
			reads = new HashMap<String, PBRead>();
		reads.clear();
		FileReader fr = null;
		BufferedReader br = null;
		String filePath = paras.getPbFile();
		try{
			File pbFile = new File(filePath);
			if(!pbFile.exists())
			{
				logger.debug(this.getClass().getName() + "\t" + "The pbreads file " + filePath + " do not exist!");
				logger.error(this.getClass().getName() + "\t" + "The pbreads file " + filePath + " do not exist!");
				return null;
			}
			fr = new FileReader(pbFile);
			br = new BufferedReader(fr);
			String line = null;
			String id = null;
			PBRead pb = null;
			boolean isFa = false;
			boolean isCheck = false;
			boolean isQua = false;
			StringBuffer sb = new StringBuffer();
			StringBuffer qua = new StringBuffer();
//			StringBuffer anotation = new StringBuffer();
			br.mark(1);
			while(true)
			{
				line = br.readLine();
				if(line != null)
					line.trim();
				if(!isCheck)
				{
					if(line.startsWith(">"))
						isFa = true;
					isCheck = true;
					br.reset();
					continue;
				}
				if(isFa)
				{ // for fasta format
					// for the last record
					if(line == null)
					{
						if(id != null && sb.length() != 0)
						{
							pb = new PBRead(sb.toString(),DNACompoundSet.getDNACompoundSet());
							pb.setAccession(new AccessionID(id));
							reads.put(id, pb);
							sb = null;
							id = null;
							pb = null;
						}
						break;
					}
					if(line.startsWith(">"))
					{
						if(sb.length() != 0 && id != null)
						{
							pb = new PBRead(sb.toString(),DNACompoundSet.getDNACompoundSet());
							pb.setAccession(new AccessionID(id));
							reads.put(id, pb);
							sb = null;
							id = null;
							pb = null;
						}
						id = line.replaceAll("^>", "");
						id = id.split("\\s")[0];
						id = id.trim();
						sb = new StringBuffer();
					} else
					{
						sb.append(line);
					}
				} else
				{ // for fastq format
					// for the last record
					if(line == null)
					{
						if(id != null && sb.length() != 0)
						{
							pb = new PBRead(sb.toString(),DNACompoundSet.getDNACompoundSet());
							pb.setAccession(new AccessionID(id));
							pb.setDescription(qua.toString());
							reads.put(id, pb);
							qua = null;
							sb = null;
							id = null;
							pb = null;
						}
						break;
					}
					if(line.startsWith("@"))
					{
						if(sb.length() != 0 && id != null)
						{
							pb = new PBRead(sb.toString(),DNACompoundSet.getDNACompoundSet());
							pb.setAccession(new AccessionID(id));
							pb.setDescription(qua.toString());
							reads.put(id, pb);
							qua = null;
							sb = null;
							id = null;
							pb = null;
						}
						id = id.replaceAll("^@", "");
						id = id.split("\\s")[0];
						id = id.trim();
						sb = new StringBuffer();
						qua = new StringBuffer();
						isQua = false;
					} else
					{
						if(line.startsWith("\\+"))
						{
							isQua = true;
							continue;
						}
						if(isQua)
						{
							qua.append(line);
						} else{
							sb.append(line);
						}
					}
				}
			}
			br.close();
		} catch(IOException e)
		{
			logger.debug(this.getClass().getName() + "\t" + e.getMessage() + "\t" + e.getClass().getName());
			logger.error(this.getClass().getName() + "\t" + e.getMessage() + "\t" + e.getClass().getName());
		} catch(Exception e)
		{
			logger.debug(this.getClass().getName() + "\t" + e.getMessage() + "\t" + e.getClass().getName());
			logger.error(this.getClass().getName() + "\t" + e.getMessage() + "\t" + e.getClass().getName());
		} finally
		{
			try{
				if(br != null)
				{
					br.close();
				}
			} catch(Exception e)
			{
				logger.debug(this.getClass().getName() + "\t" + e.getMessage() + "\t" + e.getClass().getName());
				logger.error(this.getClass().getName() + "\t" + e.getMessage() + "\t" + e.getClass().getName());
			}
		}
		
		return reads;
	}
}



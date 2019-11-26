/*
*File: agis.ps.file.GapRecordWriter.java
*User: mqin
*Email: mqin@ymail.com
*Date: 2016年6月2日
*/
package agis.ps.file;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
//import java.util.Arrays;
import java.util.List;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import agis.ps.seqs.PBGapSeq;
import agis.ps.util.GapRecord;
import agis.ps.util.Parameter;

public class GapRecordWriter {
	private static Logger logger = LoggerFactory.getLogger(GapRecordWriter.class);
	private Parameter para;
	private List<GapRecord> gaps;
	
	public GapRecordWriter()
	{
		// do nothing
	}
	
	public GapRecordWriter(Parameter para, List<GapRecord> gaps)
	{
		this.para = para;
		this.gaps = gaps;
	}
	
	public void write() {
		// TODO Auto-generated method stub
		if(gaps == null)
			throw new IllegalArgumentException(this.getClass() + "\t:The gap recould could not be null!");
		String outFolder = para.getOutFolder();
		String fileName = outFolder + System.getProperty("file.separator") + "gap_record.info";
		File file = null; 
		FileWriter fw = null;
		BufferedWriter bw = null;
		try
		{
			file = new File(fileName);
//			if (file.exists()) {
//				logger.debug("The output file" + fileName + "is exist! It will not be overwrited!");
//				return;
//			}
			if(file.exists()) {
				logger.info("The output file " + fileName + " existed. It will overwrite.");
			} else {
				if(!file.createNewFile())
				{
					logger.debug(this.getClass().getName() + "The output file" + fileName + "could not be created!");
					return;
				}
			}
			fw = new FileWriter(file);
			bw = new BufferedWriter(fw);
			for(GapRecord gr : gaps)
			{
				String p1 = gr.getStart();
				String p2 = gr.getEnd();
				int num = gr.getNums();
				List<PBGapSeq> seqs = gr.getSeqs();
				bw.write(">\t" + p1 + "\t" + p2 + "\t" + num);
//				bw.write(s + "\t" + args.get(s).size() + "\t" + Arrays.toString(args.get(s).toArray()));
				bw.newLine();
				for(PBGapSeq seq : seqs)
				{
					bw.write(seq.getId() + "\t" + seq.getStart() + "\t" + seq.getEnd() + "\t" + seq.getStrand());
					bw.newLine();
				}
			}
			bw.flush();
		} catch(IOException e)
		{
			logger.debug(this.getClass().getName() + "\t" + e.getMessage() + "\t" + e.getClass().getName());
			logger.error(this.getClass().getName() + "\t" + e.getMessage() + "\t" + e.getClass().getName());
		} finally
		{
			try{
				if(bw != null)
					bw.close();
			} catch(IOException e)
			{
				logger.debug(this.getClass().getName() + "\t" + e.getMessage() + "\t" + e.getClass().getName());
				logger.error(this.getClass().getName() + "\t" + e.getMessage() + "\t" + e.getClass().getName());
			}
		}
	}

}



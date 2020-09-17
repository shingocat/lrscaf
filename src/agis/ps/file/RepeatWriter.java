/** 
** Usage: TODO
** Author: mqin
** Email: mqin@outlook.com
** Date: 2017年8月9日
*/
package agis.ps.file;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.util.Map;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import agis.ps.seqs.Contig;
import agis.ps.seqs.Scaffold;
import agis.ps.seqs.Sequence;
import agis.ps.util.Parameter;

public class RepeatWriter {
	private static Logger logger = LoggerFactory.getLogger(RepeatWriter.class);
	private Parameter paras;
//	private Map<String, Contig> cnts;
	private Map<String, Sequence> seqs;
	
//	public RepeatWriter(Parameter paras, Map<String, Contig> cnts) {
//		this.paras = paras;
//		this.cnts = cnts;
//	}
	
	public RepeatWriter(Parameter paras, Map<String, Sequence> seqs) {
		this.paras = paras;
		this.seqs = seqs;
	}
	
	public void write() {
		File repeat = null;
		BufferedWriter bwRepeat = null;
		try{
			repeat = new File(paras.getOutFolder() + System.getProperty("file.separator") + "repeat.contigs");
			if(repeat.exists()) {
				logger.info("The output file " + repeat.getCanonicalPath() + " existed. It will overwrite.");
			} else {
				if(!repeat.createNewFile()) {
					logger.error("The output file " + repeat.getCanonicalPath() + " could not create.");
					return;
				}
			}
			bwRepeat = Files.newBufferedWriter(repeat.toPath());
			for(Map.Entry<String, Sequence> entry : this.seqs.entrySet()) {
				Sequence seq = entry.getValue();
				String id = entry.getKey();
				if(seq.isRepeat()) {
					bwRepeat.write(">" + id);
					bwRepeat.newLine();
					bwRepeat.write(seq.getForwardSeqs());
					bwRepeat.newLine();
					continue;
				}
			}
			bwRepeat.flush();
			bwRepeat.close();
		} catch(IOException e) {
			logger.debug("Error: ", e);
			logger.error(e.getMessage());
		} catch(Exception e) {
			logger.debug("Error: ", e);
			logger.error(e.getMessage());
		} finally{
			try{
				if(bwRepeat != null)
					bwRepeat.close();
			} catch(IOException e) {
				logger.debug("Error: ", e);
				logger.error(e.getMessage());
			}
			
		}
	}
}

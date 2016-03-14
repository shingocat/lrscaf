/*
*File: agis.ps.util.ScaffoldWriter.java
*User: mqin
*Email: mqin@ymail.com
*Date: 2016年2月28日
*/
package agis.ps.util;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.List;
import java.util.Map;

import org.biojava.nbio.core.sequence.DNASequence;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import agis.ps.Path;
import agis.ps.link.Contig;
import agis.ps.path.Node;
import agis.ps.path.NodePath;

public class ScaffoldWriter {

	public static Logger logger = LoggerFactory.getLogger(ScaffoldWriter.class);
	private List<NodePath> paths;
	// private Map<String, DNASequence> cnts;
	private Map<String, Contig> cnts;
	private String filePath;

	// public ScaffoldWriter(List<NodePath> paths, Map<String, DNASequence>
	// cnts, String filePath) {
	// // TODO Auto-generated constructor stub
	// this.paths = paths;
	// this.cnts = cnts;
	// this.filePath = filePath;
	// }

	public ScaffoldWriter(List<NodePath> paths, Map<String, Contig> cnts, String filePath) {
		// TODO Auto-generated constructor stub
		this.paths = paths;
		this.cnts = cnts;
		this.filePath = filePath;
	}

	public void write() {
		// TODO Auto-generated method stub
		if (filePath == null)
			return;
		File out = null;
		FileWriter fw = null;
		BufferedWriter bw = null;
		try {
			out = new File(filePath);
			if (out.exists()) {
				logger.debug("The output file of scaffold is exist! It will not be overwrited!");
				logger.info("The output file of scaffold is exist! It will not be overwrited!");
				return;
			}
			out.createNewFile();
			fw = new FileWriter(out);
			bw = new BufferedWriter(fw);
			int count = 0;
			for (NodePath p : paths) {
				bw.write(">scaffolds_" + count + "\n");
				for (int i = 0; i < p.getPathSize(); i++) {
					Node node = p.getElement(i);
					String id = node.getCnt().getID();
					String seq = "";
					if(node.getStrand().equals(Strand.FORWARD))
						seq = cnts.get(id).getSequenceAsString();
					else 
						seq = cnts.get(id).getReverseComplementSeq();
					bw.write(seq);
					int nLen = node.getMeanDist2Next();
					if (i != p.getPathSize() - 1) {
						if (nLen < 0)
							bw.write(repeatString("M", nLen));
						else
							bw.write(repeatString("N", nLen));
					}
				}
				bw.write("\n");
				count++;
			}

		} catch (IllegalArgumentException e) {
			logger.debug(e.getMessage());
			logger.error(e.getMessage());
		} catch (IOException e) {
			logger.debug(e.getMessage());
			logger.error(e.getMessage());
		} catch (Exception e) {
			logger.debug(e.getMessage());
			logger.error(e.getMessage());
		} finally {
			if (bw != null)
				try {
					bw.close();
				} catch (IOException e) {
					logger.debug(e.getMessage());
					logger.error(e.getMessage());
					;
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
		StringBuilder sb = new StringBuilder(str);
		if (times < 0) {
			times = 0 - times;
			for (int i = 0; i <= times; i++)
				sb.append(str);
		} else {
			for (int i = 1; i <= times; i++)
				sb.append(str);
		}
		return sb.toString();
	}
	
}

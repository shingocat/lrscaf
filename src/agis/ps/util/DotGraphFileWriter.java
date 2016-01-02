/*
*File: agis.ps.util.DotGraphWriter.java
*User: mqin
*Email: mqin@ymail.com
*Date: 2015年12月28日
*/
package agis.ps.util;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import agis.ps.Path;
import agis.ps.SimplePath;

public class DotGraphFileWriter {

	public final Logger logger = LoggerFactory.getLogger(DotGraphFileWriter.class);

	private String filePath;
	private List<SimplePath> paths;

	public DotGraphFileWriter(String filePath, List<SimplePath> paths) {
		this.filePath = filePath;
		this.paths = paths;
	}

	public void write() {
		if(paths == null)
			return;
		
		File out = null;
		FileWriter fw = null;
		BufferedWriter bw = null;
		try {
			out = new File(filePath);
			if (out.exists()) {
				logger.debug("The output file of dot graph is exist!");
				logger.info("The output file of dot graph is exist!");
				return;
			}
			out.createNewFile();
			fw = new FileWriter(out);
			bw = new BufferedWriter(fw);
			bw.write("digraph graph{\n");
			for(SimplePath sp : paths)
			{
				bw.write(sp.getStart() + " -> " + sp.getEnd() + " [label=" +
						sp.getLabel() + ",color=" + sp.getColor().toString() + "];\n");
			}
			bw.write("}");
		} catch (IOException e) {
			logger.debug(e.getMessage());
			logger.error(e.getMessage());
		} finally
		{
			if(bw != null)
				try {
					bw.close();
				} catch (IOException e) {
					logger.debug(e.getMessage());
					logger.error(e.getMessage());;
				}
		}

	}

}

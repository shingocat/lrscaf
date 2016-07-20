/*
*File: agis.ps.util.DotGraphWriter.java
*User: mqin
*Email: mqin@ymail.com
*Date: 2015年12月28日
*/
package agis.ps.file;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.List;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import agis.ps.link.Edge;
import agis.ps.path.NodePath;
import agis.ps.util.Strand;

public class DotGraphFileWriter {

	public static final Logger logger = LoggerFactory.getLogger(DotGraphFileWriter.class);

	public static void writeEdge(String filePath, List<Edge> edges) {
		if (filePath == null)
			return;
		File out = null;
		FileWriter fw = null;
		BufferedWriter bw = null;
		try {
			out = new File(filePath);
			if (out.exists()) {
				logger.info(DotGraphFileWriter.class.getName() + "\t" + filePath + " is exist!");
				return;
			}
			out.createNewFile();
			fw = new FileWriter(out);
			bw = new BufferedWriter(fw);
			bw.write("digraph G{\n");
			for (Edge e : edges) {
				String color = "";
				if (!e.isFake()) {
					if (e.getoStrand().equals(Strand.FORWARD) && e.gettStrand().equals(Strand.FORWARD))
						color = "red";
					else if (e.getoStrand().equals(Strand.FORWARD) && e.gettStrand().equals(Strand.REVERSE))
						color = "green";
					else if (e.getoStrand().equals(Strand.REVERSE) && e.gettStrand().equals(Strand.REVERSE))
						color = "blue";
					else if (e.getoStrand().equals(Strand.REVERSE) && e.gettStrand().equals(Strand.FORWARD))
						color = "yellow";
				} else {
					color = "gray";
				}
				bw.write(e.getOrigin().getID() + " -> " + e.getTerminus().getID() + " [label=\""
						+ e.getoStrand().toString() + " " + e.gettStrand().toString() + ":" + e.getLinkNum() + ":"
						+ e.getDistMean() + ":" + e.getDistSd() + "\",color=" + color + "];\n");
			}
			bw.write("}");
		} catch (IOException e) {
			logger.error(DotGraphFileWriter.class.getName() + "\t" + e.getMessage());
		} catch (Exception e) {
			logger.error(DotGraphFileWriter.class.getName() + "\t" + e.getMessage());
		} finally {
			try {
				if (bw != null)
					bw.close();
			} catch (IOException e) {
				logger.error(DotGraphFileWriter.class.getName() + "\t" + e.getMessage());
			}
		}

	}

	public static void writeNodePaths(String filePath, List<NodePath> paths) {
		// TODO Auto-generated method stub
		if (filePath == null)
			return;
		File out = null;
		FileWriter fw = null;
		BufferedWriter bw = null;
		try {
			out = new File(filePath);
			if (out.exists()) {
				logger.error(filePath + " is exist!");
				return;
			}
			out.createNewFile();
			fw = new FileWriter(out);
			bw = new BufferedWriter(fw);
			int count = 0;
			for (NodePath p : paths) {
				bw.write("digraph G" + count + "{\n");
				bw.write(p.toString() + ";\n");
				bw.write("}\n");
				count++;
			}

		} catch (IOException e) {
			logger.error(DotGraphFileWriter.class.getName() + "\t" + e.getMessage());
		} finally {
			if (bw != null)
				try {
					bw.close();
				} catch (IOException e) {
					logger.error(DotGraphFileWriter.class.getName() + "\t" + e.getMessage());
				}
		}
	}

}

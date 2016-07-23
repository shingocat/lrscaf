/*
*File: agis.ps.file.ContigIndexer.java
*User: mqin
*Email: mqin@ymail.com
*Date: 2016年7月19日
*/
package agis.ps.file;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;

import org.apache.lucene.analysis.Analyzer;
import org.apache.lucene.analysis.standard.StandardAnalyzer;
import org.apache.lucene.document.Document;
import org.apache.lucene.document.StoredField;
import org.apache.lucene.document.StringField;
import org.apache.lucene.document.TextField;
import org.apache.lucene.document.Field.Store;
import org.apache.lucene.index.CorruptIndexException;
import org.apache.lucene.index.IndexWriter;
import org.apache.lucene.index.IndexWriterConfig;
import org.apache.lucene.index.IndexWriterConfig.OpenMode;
import org.apache.lucene.store.Directory;
import org.apache.lucene.store.LockObtainFailedException;
import org.apache.lucene.store.SimpleFSDirectory;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import agis.ps.util.Parameter;

public class ContigIndexer {

	private static Logger logger = LoggerFactory.getLogger(ContigIndexer.class);
	private Parameter paras = null;
	private Directory directory = null;
	private IndexWriter writer = null;

	public ContigIndexer(Parameter paras) {
		this.paras = paras;
	}

	public boolean indexing() {
		long start = System.currentTimeMillis();
		boolean isValid = true;
		File file = null;
		FileReader fr = null;
		BufferedReader br = null;
		Analyzer analyzer = null;
		IndexWriterConfig config = null;
		try {
			file = new File(paras.getCntFile());
			String indexPath = paras.getOutFolder() + System.getProperty("file.separator") + "cnt.index";
			File indexDir = new File(indexPath);
			directory = new SimpleFSDirectory(indexDir.toPath());// FSDirectory.open(indexDir.toPath());
			fr = new FileReader(file);
			br = new BufferedReader(fr);
			analyzer = new StandardAnalyzer();
			config = new IndexWriterConfig(analyzer);
			config.setOpenMode(OpenMode.CREATE);
			writer = new IndexWriter(directory, config);
			String line = null;
			String id = null;
			String seq = "";
			while (true) {
				line = br.readLine();
				if (line == null) {
					if (seq.length() != 0) {
						Document doc = new Document();
						doc.add(new TextField("id", id, Store.YES));
						doc.add(new StoredField("seq", seq));
						doc.add(new StoredField("len", String.valueOf(seq.length())));
						writer.addDocument(doc);
					}
					break;
				} else {
					if (line.startsWith(">")) {
						if (seq.length() != 0) {
							Document doc = new Document();
							doc.add(new StringField("id", id, Store.YES));
							doc.add(new StoredField("seq", seq));
							doc.add(new StoredField("len", String.valueOf(seq.length())));
							writer.addDocument(doc);
						}
						line = line.replaceAll(">", "");
						String[] arrs = line.split("\\s+");
						id = arrs[0];
						seq = "";
					} else {
						seq = seq + line;
					}
				}
			}
			writer.close();
		} catch (CorruptIndexException e) {
			logger.error(this.getClass().getName() + "\t" + e.getMessage());
			isValid = false;
		} catch (LockObtainFailedException e) {
			logger.error(this.getClass().getName() + "\t" + e.getMessage());
			isValid = false;
		} catch (IOException e) {
			logger.error(this.getClass().getName() + "\t" + e.getMessage());
			isValid = false;
		} catch (Exception e) {
			logger.error(this.getClass().getName() + "\t" + e.getMessage());
			isValid = false;
		} finally {
			try {
				if (br != null)
					br.close();
				if (fr != null)
					fr.close();
				if (file != null)
					file = null;
				if (writer != null)
					writer.close();
				if (directory != null)
					directory = null;
				if (analyzer != null)
					analyzer = null;
				if (config != null)
					config = null;
			} catch (IOException e) {
				logger.error(this.getClass().getName() + "\t" + e.getMessage());
			}
		}
		long end = System.currentTimeMillis();
		logger.info("Building Indexer time: " + (end - start) + " ms");
		return isValid;
	}
}

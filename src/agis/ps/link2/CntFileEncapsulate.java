/*
*File: agis.ps.link2.CntFileEncapsulate.java
*User: mqin
*Email: mqin@ymail.com
*Date: 2016年10月11日
*/
package agis.ps.link2;

import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 * Encapsulate the contig file;
 * lengths and seqs is reorder and rename to 
 * int type from zero to last;
 * @author mqin
 *
 */
public class CntFileEncapsulate {
	private static Logger logger = LoggerFactory.getLogger(CntFileEncapsulate.class);
	// original contig id to its length and seq hashmap
	private Map<String, LengthAndSeq> cnts = new HashMap<String, LengthAndSeq>();
	// new contig id to original contig id;
	private Map<Integer, String> nids2oids = new HashMap<Integer, String>(); 
	// original id to new contig id;
	private Map<String, Integer> oids2nids = new HashMap<String, Integer>();
	
	public CntFileEncapsulate()
	{
		
	}
	
	public String getOriginalId(int nid)
	{
		if(this.nids2oids.containsKey(nid))
			return this.nids2oids.get(nid);
		else
			return null;
	}
	
	public Integer getNewId(String oid)
	{
		if(this.oids2nids.containsKey(oid))
			return this.oids2nids.get(oid);
		else
			return null;
	}
	
	public void addLength(int nid, String oid, int len)
	{
		if(nids2oids.containsKey(nid))
		{
			if(!this.nids2oids.get(nid).equalsIgnoreCase(oid))
				logger.info("New id " + nid + "\t" + " Original id " + oid + " conflict mapping relationship!");
		} else
		{
			this.nids2oids.put(nid, oid);
		}
		if(this.oids2nids.containsKey(oid))
		{
			if(!this.oids2nids.get(oid).equals(nid))
				logger.info("New id " + nid + "\t" + " Original id " + oid + " conflict mapping relationship!");
		} else
		{
			this.oids2nids.put(oid, nid);
		}
		this.addLengthByOriginalId(oid, len);
	}
	
	public void addLengthByOriginalId(String id, int len)
	{
		if(!cnts.containsKey(id))
		{
			LengthAndSeq las = new LengthAndSeq();
			las.setLength(len);
			cnts.put(id, las);
		} else
		{
			LengthAndSeq las = cnts.get(id);
			int cl = las.getLength();
			// keep the first value;
			if(cl == 0)
			{
				if(len >= 0)
					las.setLength(len);
			} else
			{
				if(cl != len)
					logger.info("Contig " + id + " length has two different values, keep first one!");
			}
		}
	}
	
	public Integer getLengthByNewId(String nid)
	{
		return this.getLengthByNewId(Integer.valueOf(nid));
	}
	
	public Integer getLengthByNewId(Integer nid)
	{
		if(this.nids2oids.containsKey(nid))
		{
			String oid = this.nids2oids.get(nid);
			return this.getLengthByOriginalId(oid);
		} else
		{
			return null;
		}
	}
	
	public int getLengthByOriginalId(String id)
	{
		if(cnts.containsKey(id))
		{
			LengthAndSeq las = cnts.get(id);
			return las.getLength();
		} else
		{
			return 0;
		}
	}
	
	public int getLengthByNewId(int nid)
	{
		if(this.nids2oids.containsKey(nid))
		{
			String oid = this.nids2oids.get(nid);
			return this.getLengthByOriginalId(oid);
		} else
		{
			return 0;
		}
	}
	
	public void addSeqByOriginalId(String oid, byte [] seq)
	{
		if(!cnts.containsKey(oid))
		{
			LengthAndSeq las = new LengthAndSeq();
			las.setSeq(seq);
			cnts.put(oid, las);
		}else
		{
			LengthAndSeq las = cnts.get(oid);
			byte [] cs = las.getSeq();
			// keep the first value;
			if(cs != null)
			{
				if(cs.length == 0)
				{
					if(seq.length >= 0)
						las.setSeq(seq);
				} else
				{
					if(seq.length != cs.length)
						logger.info("Contig " + oid + " length has two different values, keep first one!");
				}
			} else
			{
				las.setSeq(seq);
			}
			
		}
	}
	
	public byte [] getSeqByOriginalId(String id)
	{
		if(cnts.containsKey(id))
		{
			LengthAndSeq las = cnts.get(id);
			return las.getSeq();
		} else
		{
			return null;
		}
	}

	public void add(String oid, String seq) {
		if(!this.oids2nids.containsKey(oid))
		{
//			logger.info(this.getClass().getName() + "\t" + "The origin id " + oid + " does not map to new id!");
			return;
		}
		int nid = this.oids2nids.get(oid);
		this.addLength(nid, oid, seq.length());
		byte [] bytes = this.encode(seq);
		this.addSeqByOriginalId(oid, bytes);
	}
	
	public String getOriginalSeqByNewId(String nid)
	{
		if(!this.nids2oids.containsKey(Integer.valueOf(nid)))
		{
			logger.info(this.getClass().getName() + "\t" + "Do not exist " + nid + " element contig!");
			return null;
		}
		String oid = this.nids2oids.get(Integer.valueOf(nid));
		return this.getOriginalSeqByOriginId(oid);
	}
	
	public String getOriginalSeqByOriginId(String oid)
	{
		if(!cnts.containsKey(oid))
		{
			logger.info(this.getClass().getName() + "\t" + "Do not exist " + oid + " element contig!");
			return null;
		}
		LengthAndSeq las = cnts.get(oid);
		return decode(las.getSeq(), las.getLength());
	}
	
	private String decode(byte[] seq, int length)
	{
		StringBuffer sb = new StringBuffer();
		for(int i = 0; i < seq.length; i++)
		{
			byte value = seq[i];
			for(int j = 3; j >= 0; j--)
			{
				byte base = (byte) ((value >> (2*j)) & 0x03);
				switch(base){
				case 0x00:
					sb.append("A");
					break;
				case 0x01:
					sb.append("C");
					break;
				case 0x02:
					sb.append("G");
					break;
				case 0x03:
					sb.append("T");
					break;
				}
			}
		}
		return sb.substring(0, length);
	}
	
	private byte[] encode(String origin)
	{
		int olen = origin.length();
		int length = olen / 4 + 1;
		byte [] byteSeq = new byte[length];
		int start = 0;
		for(int i = 0; i <= origin.length()/4; i++)
		{
			String subSeq = null;
			if((start + 4) > olen)
				subSeq = origin.substring(start, olen);
			else
				subSeq = origin.substring(start, start + 4);
			byteSeq[i] = toByte(subSeq);
			start = start + 4;
		}
		return byteSeq;
	}
	
	private byte toByte(String subSeq)
	{
		subSeq = subSeq.toUpperCase();
		char [] os = subSeq.toCharArray();
	    byte [] bs = new byte[4] ;
		byte d = 0;
		for(int i = 0; i < os.length; i++ )
		{
			byte value = 0;
			switch(os[i])
			{
			case 'A' : value = 0x00;
			break;
			case 'T' : value = 0x03;
			break;
			case 'C' : value = 0x01;
			break;
			case 'G' : value = 0x02;
			break;
			}
			bs[i] = value;
		}		
		d = (byte) (((bs[0])<<6) | ((bs[1])<<4) |((bs[2]) << 2) |((bs[3]) << 0));
		return (byte) (d & 0xff);
	}
}

class LengthAndSeq
{
	private int length;
	private byte [] seq;
	
	public LengthAndSeq()
	{
		
	}
	
	public LengthAndSeq(int length, byte seq)
	{
		
	}
	
	public void setLength(int len)
	{
		this.length = len;
	}
	
	public int getLength()
	{
		return this.length;
	}
	
	public void setSeq(byte [] seq)
	{
		this.seq = seq;
	}
	
	public byte[] getSeq()
	{
		return this.seq;
	}

	@Override
	public String toString() {
		return "LengthAndSeq [length=" + length + ", seq=" + Arrays.toString(seq) + "]";
	}		
}



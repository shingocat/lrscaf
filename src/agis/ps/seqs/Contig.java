/*
*File: agis.ps.link.Contig.java
*User: mqin
*Email: mqin@ymail.com
*Date: 2016年1月12日
*/
package agis.ps.seqs;

import agis.ps.util.MisassemblyRegion;
import agis.ps.util.SequenceBytesTransformer;

public class Contig extends Sequence {
	private Integer start;
	private Integer end;
	private byte[] seqs;
	
	public Integer getStart() {
		return start;
	}

	public void setStart(Integer start) {
		this.start = start;
	}

	public Integer getEnd() {
		return end;
	}

	public void setEnd(Integer end) {
		this.end = end;
	}

	public void setSeqs(String seqs) {
		// original code 2020/9/16
		// this.seqs = this.encode(seqs);
		this.seqs = SequenceBytesTransformer.toByte(seqs);
	}

	public String getForwardSeqs() {
		// original code 2020/9/16
		// return this.decode(seqs, length);
		return SequenceBytesTransformer.toSequence(this.seqs, this.getLength());
	}

	// need to modify;
	public String getComplementReverseSeqs() {
		return SequenceBytesTransformer.toComplementReverseSeqs(this.seqs, this.getLength());
		// original code 2020/9/16
		// if(seqs == null || length == 0)
		// return "";
		// StringBuffer sb = new StringBuffer();
		// for(int i = seqs.length - 1; i >= 0; i--)
		// {
		// byte value = seqs[i];
		// StringBuffer temp = new StringBuffer();
		// for(int j = 3; j >= 0; j--)
		// {
		// byte base = (byte) ((value >> (2*j)) & 0x03);
		// switch(base){
		// case 0x00:
		// temp.append("T");
		// break;
		// case 0x01:
		// temp.append("G");
		// break;
		// case 0x02:
		// temp.append("C");
		// break;
		// case 0x03:
		// temp.append("A");
		// break;
		// }
		// }
		// sb.append(temp.reverse());
		// }
		// int size = sb.length();
		// return sb.toString().substring(size - length);

		// original method
		// String seq = this.decode(seqs, length);
		// if(seq == null || seq.length() == 0)
		// return "";
		// StringBuffer sb = new StringBuffer();
		// seq.toUpperCase();
		// for (int j = (seq.length() - 1); j >= 0; j--) {
		// switch (seq.charAt(j)) {
		// case 'A':
		// sb.append("T");
		// break;
		// case 'T':
		// sb.append("A");
		// break;
		// case 'C':
		// sb.append("G");
		// break;
		// case 'G':
		// sb.append("C");
		// break;
		// default:
		// sb.append("N");
		// break;
		// }
		// }
		// return sb.toString();
	}

	private String decode(byte[] seq, int length) {
		StringBuffer sb = new StringBuffer();
		for (int i = 0; i < seq.length; i++) {
			byte value = seq[i];
			for (int j = 3; j >= 0; j--) {
				byte base = (byte) ((value >> (2 * j)) & 0x03);
				switch (base) {
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

	private byte[] encode(String origin) {
		int olen = origin.length();
		int length = olen / 4 + 1;
		byte[] byteSeq = new byte[length];
		int start = 0;
		for (int i = 0; i <= origin.length() / 4; i++) {
			String subSeq = null;
			if ((start + 4) > olen)
				subSeq = origin.substring(start, olen);
			else
				subSeq = origin.substring(start, start + 4);
			byteSeq[i] = toByte(subSeq);
			start = start + 4;
		}
		return byteSeq;
	}

	private byte toByte(String subSeq) {
		subSeq = subSeq.toUpperCase();
		char[] os = subSeq.toCharArray();
		byte[] bs = new byte[4];
		byte d = 0;
		for (int i = 0; i < os.length; i++) {
			byte value = 0;
			switch (os[i]) {
			case 'A':
				value = 0x00;
				break;
			case 'T':
				value = 0x03;
				break;
			case 'C':
				value = 0x01;
				break;
			case 'G':
				value = 0x02;
				break;
			}
			bs[i] = value;
		}
		d = (byte) (((bs[0]) << 6) | ((bs[1]) << 4) | ((bs[2]) << 2) | ((bs[3]) << 0));
		return (byte) (d & 0xff);
	}
}

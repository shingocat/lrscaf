package agis.ps.util;
/**
* @author Mao Qin
* @version 2020年9月16日 上午10:13:41
* @Email mqin@outlook.com
* @Description class usage.
* @Copyright All Right Reserved 2020.
*/
public class SequenceBytesTransformer {
	
	public static byte[] toByte(String seq) {
		if(seq == null || seq.isEmpty())
			return null;
		int seqLength = seq.length();
		int size = seqLength / 4 + ((seqLength % 4) == 0 ? 0 : 1);
		byte [] bytes = new byte[size];
		char [] chars = seq.toCharArray();
		for(int i = 0; i < size; i++) {
			bytes[i] = 0b00;
			for(int j = 0; j < 4; j++) {
				int index = i + j + i * 3;
				if(index >= seqLength)
					break;
				byte code = 0b00;
				switch(chars[index]) {
					case 'A':
					case 'a':
						code = 0b00;
						break;
					case 'C':
					case 'c':
						code = 0b01;
						break;
					case 'G':
					case 'g':
						code = 0b10;
						break;
					case 'T':
					case 't':
						code = 0b11;
						break;
				}
				bytes[i] = (byte) (bytes[i] | (code << (6 - j * 2)));
			}
//			System.out.println(i + " bytes is " + byteToBit(bytes[i]));
		}
		return bytes;
	}
	
	public static String toSequence(byte []  bytes, int seqLength) {
		StringBuilder sb = new StringBuilder();
		int count = 0; // check whether over the seq length
		for(int i = 0; i < bytes.length; i++) {
			int j = 0; // check whether over the bytes length, one byte for 4 base pair
			while(count < seqLength) {
				byte code = (byte) (bytes[i] >> (6- j * 2) & 0b00000011);
				switch(code) {
					case 0b00:
						sb.append("A");
						break;
					case 0b01:
						sb.append("C");
						break;
					case 0b10:
						sb.append("G");
						break;
					case 0b11:
						sb.append("T");
						break;
				}
				j = j + 1;
				count = count + 1;
				if(j > 3) break;
			}
		}
		return sb.toString();
	}
	
	public static String toComplementReverseSeqs(byte [] bytes, int seqLength) {
		StringBuilder sb = new StringBuilder();
		StringBuilder temp = new StringBuilder();
		for(int i = bytes.length - 1; i >= 0; i--) {
			int j = 3;
			while(j >= 0) {
				byte code = (byte) ((bytes[i] >> (j * 2)) & 0b00000011);
				switch(code) {
					case 0b00:
						temp.append("T");
						break;
					case 0b01:
						temp.append("G");
						break;
					case 0b10:
						temp.append("C");
						break;
					case 0b11:
						temp.append("A");
						break;
				}
				j = j - 1;
			}
			sb.append(temp.reverse());
			temp = new StringBuilder();
		}
		return sb.substring(sb.length() - seqLength);
	}
	
	public static String byteToBit(byte b) {
		return ""
		+ (byte) ((b >> 7) & 0x1) + (byte) ((b >> 6) & 0x1)
		+ (byte) ((b >> 5) & 0x1) + (byte) ((b >> 4) & 0x1)
		+ (byte) ((b >> 3) & 0x1) + (byte) ((b >> 2) & 0x1)
		+ (byte) ((b >> 1) & 0x1) + (byte) ((b >> 0) & 0x1);
		}
}

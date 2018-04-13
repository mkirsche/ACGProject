import java.util.*;
import java.io.*;
public class TestEnds {
	static int maxDist = 5;
	static int minClipped = 20;
public static void main(String[] args) throws IOException
{
	int site = 18109144;
	String samFn = "/home/mkirsche/schatz/project/chr22.18109144.sam";
	String entry = "chr22	18109144	0	N	<INS>	.	PASS	SEQ=GCGAGCACGTTGGCTTATGACCTCGCGGTCAAGGTTTTGTGGCTCTGGTCCCCGGGGCCACGTACCTGCATTAGGTGGGATTCTATAGAGCATTGACCCGGCTGCTTTCCCGTCGCAGCAAGCCGTACATGTCCCAAATTAAGTGGCCCCGCAGACAACTGGTTAACAGTACCATCTACCGGTTCCGGCACACCGTGACTCTCAGGACGGTTCATTGTAGAGAGGCTCACTAGTGCTATCGCCGAGACCTTCCGGAAAGCATTTGTACTGCAACTAAATTCTGCGCGGTTTTGTAGATGTACTGGTTACGTCGCGATCAATCCGTGCGGTCGGTGCGTTGGTAGCATGACTAGTGTCCGACAGGGGTAATTTCCTCAGCCTTGCGACTAAGCTCTATGCTAATGCGGTTAGCTTCCTCAGAGTGCTAGTCACGGGGTCACAATTCGGTATCTATCCCGCTCGCATCCGAGAGGCGATGATGACCCTGAGCGGCGTCGTTGTTTGATTAGTTACAGAGGAGACCTCGACCGTCCAAGTAATAATTACAGCTATGTGAGCGTGAGTCTTACCCTTCAGGCGTACATTATATCACTCTTGTTAATTAATGGGGTCTAGGGAACAATATGAAAGACAAGACCATGTTGAAGAAATCCTGGACGGTAGTATTCAATATATAGTAAAACGTACAAGCTTCTGCAGTGAGACGTCTCCATAGACTGAGTAACCAATGACCTAGCTTATTTTTGTGGAGCATGCGCACACGAGCGCCTAATTTGACTAAGCACTAGGTGGACTATAATGCACCCAGGCCAGCCAAGGTGACTAAGAGAAAAGTATTGAATAGACCGCCAGACCGTGTGCTTTGACATAACGCGCACCATTACCCTTTTTACCGCACTTTGCAATAGCCTGCTGTCCAATAAATCCCAAGTAAATAACATCGCCCAATTCCTATCTCTGGGCTATTTACGCTCAGACCTCGGTGAATAGCTAAGGGGACGGTAAAAGATGGGTAAGCCAGGCAAACCTGCACTTTGCCCGACGTATGTCGGGTTTGACGACGATTCAGAAGGGCACCAGCCCCATAATGGCACCCTTGTCTCGCTTTCGGGCCCTGGGGCTTGCCCTGGTTCTTTTAATCGATAAAAGACAGCCGGTTCTTAGTCATGGAGCGAAGACTACTGATTAAGTGGTCATCGACAAAAACGCGTTCGAACCTTGTCGCCAGTGCGGTAACTTACGAATTCGCGATCCATGATCTTGCCGCTGCTGTTGGAGATTCCAGCGGCTACTGCTAAAATCACTTTCTCACGGTTCTTGTGCAGGATTCTAGAGATCGAATTTGTGAAAAGGTGAAGATGCACACCTTCTGCCCCCTTTGCAATATCCAGTGCGAGCACAAGTTTACAAAATTTACACGACTAT;SVLEN=1425	GT	1";
	String s = getSeq(entry);
	Scanner input = new Scanner(new FileInputStream(new File(samFn)));
	// leftClipped has portions of reads which should align to the end of the SV,
	// and rightClipped has portions of reads which should align ti the start of the SV
	ArrayList<String> leftClipped = new ArrayList<String>(), rightClipped =new ArrayList<String>();
	while(input.hasNext())
	{
		String line = input.nextLine();
		if(line.startsWith("@"))
		{
			continue;
		}
		String[] tokens = line.split("\t");
		String cigar = tokens[5];
		parse(cigar);
		int[] clipped = getClips();
		int[] hc = getHC();
		int start = Integer.parseInt(tokens[3]);
		int end = end(start);
		String seq = tokens[9];
		if(Math.abs(start - site) <= maxDist && clipped[0] >= minClipped && clipped[0] > hc[0])
		{
			leftClipped.add(seq.substring(0, clipped[0] - hc[0]));
		}
		if(Math.abs(end - site) <= maxDist && clipped[1] >= minClipped && clipped[1] > hc[1])
		{
			rightClipped.add(seq.substring(seq.length() - clipped[1] + hc[1]));
		}
	}
	// TODO compare leftClipped and rightClipped reads to the ends of s and see if they are contained in it
}
/*
 * Determines whether a character represents a non-ambiguous base pair (A or C or G or T)
 * c is the character being considered
 */
static boolean isBasePair(char c)
{
	return c == 'A' || c == 'C' || c == 'G' || c == 'T';
}
/*
 * Get the SV sequence from an entry of a VCF file
 * line is an entry from a VCF file
 */
static String getSeq(String line)
{
	int idx = line.indexOf("SEQ=") + 4;
	int end = idx;
	while(end < line.length() && isBasePair(line.charAt(end))) end++;
	return line.substring(idx,end);
}
static char[] alignmentTypes;
static int[] alignmentLengths;
static void parse(String cigar)
{
	ArrayList<Integer> lens = new ArrayList<Integer>();
	ArrayList<Character> types = new ArrayList<Character>();
	
	int idx = 0;
	if(cigar.length() > 1)
	{
		while(idx < cigar.length())
		{
			int nidx = idx;
			while(cigar.charAt(nidx) <= '9' && cigar.charAt(nidx) >= '0') nidx++;
			lens.add(Integer.parseInt(cigar.substring(idx, nidx)));
			types.add(cigar.charAt(nidx));
			idx = nidx+1;
		}
	}
	
	alignmentTypes = new char[types.size()];
	for(int i = 0; i<types.size(); i++) alignmentTypes[i] = types.get(i);
	alignmentLengths = new int[lens.size()];
	for(int i = 0; i<lens.size(); i++) alignmentLengths[i] = lens.get(i);
}
/*
 * Returns the number of hard clipped bases
 */
static int[] getHC()
{
	int[] res = new int[2];
	int at = 0;
	int n = alignmentTypes.length;
	while(at < n)
	{
		char c = alignmentTypes[at];
		if(c == 'H')
		{
			res[0]+=alignmentLengths[at];
			at++;
		}
		else break;
	}
	at = n-1;
	while(at >= 0)
	{
		char c = alignmentTypes[at];
		if(c == 'H')
		{
			res[1]+=alignmentLengths[at];
			at--;
		}
		else break;
	}
	return res;
}
/*
 * Returns the number of clipped bases on each end
 */
static int[] getClips()
{
	int[] res = new int[2];
	int at = 0;
	int n = alignmentTypes.length;
	while(at < n)
	{
		char c = alignmentTypes[at];
		if(c == 'S' || c == 'H')
		{
			res[0]+=alignmentLengths[at];
			at++;
		}
		else break;
	}
	at = n-1;
	while(at >= 0)
	{
		char c = alignmentTypes[at];
		if(c == 'S' || c == 'H')
		{
			res[1]+=alignmentLengths[at];
			at--;
		}
		else break;
	}
	return res;
}
static int end(int start)
{
	int res = start;
	int n = alignmentTypes.length;
	for(int i = 0; i<n; i++)
	{
		char c = alignmentTypes[i];
		int len = alignmentLengths[i];
		if(c == 'M' || c == 'D' || c == 'N' || c == '=' || c == 'X') res += len;
	}
	return res;
}
}

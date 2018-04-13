import java.util.*;
import java.io.*;
public class SVSearcher {
	static int k = 30;
	static HashSet<Long> kmers;
	static int[] vals;
	static int sharedKmers = 20;
	static String suff;
	// TODO make some of these global values command line arguments
	static boolean readAln = false;
	static boolean writeAln = true;
	static String filterAlnFn = "/home/mkirsche/schatz/project/aln_filtered.txt";
public static void main(String[] args) throws IOException
{
	vals = new int[256];
	Arrays.fill(vals, -1);
	vals['A'] = 0; vals['C'] = 1; vals['G'] = 2; vals['T'] = 3;
	suff = "." + sharedKmers + "." + k + "mers";
	
	/*
	 * TODO make most of these command line arguments
	 */
	String action = "filter";
	String alnFn = "/home/mkirsche/schatz/project/aln.sam";
	String read1Fn = "/home/mkirsche/data/chr22.fa.10.1";
	String read2Fn = "/home/mkirsche/data/chr22.fa.10.2";
	String vcfFn = "/home/mkirsche/schatz/project/chr22.fa.10.vcf";
	String out1Fn = read1Fn + suff;
	String out2Fn = read2Fn + suff;
	if(action.equals("filter"))
	{
		filterReads(vcfFn, read1Fn, read2Fn, out1Fn, out2Fn);
	}
	else if(action.equals("detect"))
	{
		detectSV(vcfFn, out1Fn, out2Fn, alnFn);
	}
}
/*
 * Detects which SVs have short-read evidence supporting them
 * vcfFn is the filename of a VCF file containing an entry for each structural variant
 * out1Fn is the filename of a file containing the first read in each paired-end read (after filtering)
 * out2Fn is the filename of a file containing the second read in each paired-end read (after filtering)
 * alnFn is the filename of a SAM file containing the alignments of all paired-end reads to the reference genome
 */
static void detectSV(String vcfFn, String out1Fn, String out2Fn, String alnFn) throws IOException
{
	Scanner vcfScanner = new Scanner(new FileInputStream(new File(vcfFn)));
	int varIndex = 0;
	while(vcfScanner.hasNext())
	{
		String line = vcfScanner.nextLine();
		if(line.length() == 0 || line.charAt(0) == '#') continue;
		varIndex++;
		ArrayList<Read> reads1 = fromFasta(out1Fn + "." + varIndex);
		ArrayList<Read> reads2 = fromFasta(out2Fn + "." + varIndex);
		System.out.println(reads1.size()+" "+reads2.size());
		SV sv = new SV(line);
		String s = sv.seq;
		ArrayList<Alignment> alns = readAln ? readAlns(filterAlnFn + "." + varIndex) : filterAlignments(reads1, reads2, alnFn);
		if(writeAln) writeAlns(alns, filterAlnFn + "." + varIndex);
		System.out.println(alns.size());
		for(Alignment aln : alns)
		{
			String seq = aln.seq;
			int[] interval = localAlign(seq, s);
			if(interval[0] != -1)
			{
				// TODO process these intervals to determine whether or not the whole insertions has reads aligning to it
				System.out.println(interval[0]+" "+interval[1]);
			}
		}
	}
}
/*
 * Finds best match between s and t
 * Returns pair of positions in t which represent interval that s aligns to
 */
// TODO adjust these parameters and fix scoring function
static double insertPenalty = -0.1;
static double delPenalty = -0.1;
static double changePenalty = -0.1;
static double match = 1;
static double alnThreshold = 0.8;
static int[] localAlign(String s, String t)
{
	int n = s.length(), m = t.length();
	double[][] dp = new double[n+1][m+1];
	int[][] start = new int[n+1][m+1];
	for(int i = 0; i<=m; i++) start[0][i] = -i;
	for(int i = 1; i<=n; i++)
		for(int j = 1; j<=m; j++)
		{
			dp[i][j] = insertPenalty * (i-1);
			double matchVal = dp[i-1][j-1] + s.charAt(i-1) == t.charAt(j-1) ? (match) : (changePenalty);
			double ins = dp[i-1][j] + insertPenalty;
			double del = dp[i][j-1] + delPenalty;
			dp[i][j] = Math.max(dp[i][j], matchVal);
			dp[i][j] = Math.max(dp[i][j], ins);
			dp[i][j] = Math.max(dp[i][j], del);
			if(dp[i][j] == matchVal) start[i][j] = start[i-1][j-1];
			else if(dp[i][j] == ins) start[i][j] = start[i-1][j];
			else if(dp[i][j] == del) start[i][j] = start[i][j-1];
			else start[i][j] = j;
		}
	double best = dp[n][0];
	int end = 0;
	int st = start[n][0];
	for(int i = 1; i<=m; i++)
	{
		if(dp[n][i] > best)
		{
			best = dp[n][i];
			end = i;
			st = start[n][i];
		}
	}
	if(best >= alnThreshold * n) return new int[] {st, end};
	else return new int[] {-1, -1};
}
/*
 * Filters the alignments in a SAM file to include only those present in a list of reads
 * reads1 is a list of the first read in each paired end read we are interested in
 * reads2 is a list of the second read in each paired end read we are interested in
 * alnFn is the filename of a SAM file containing the alignments of all reads to the reference
 */
static ArrayList<Alignment> filterAlignments(ArrayList<Read> reads1, ArrayList<Read> reads2, String alnFn) throws IOException
{
	HashMap<String, Read> readNames1 = new HashMap<String, Read>();
	HashMap<String, Read> readNames2 = new HashMap<String, Read>();
	for(Read r : reads1) readNames1.put(r.name, r);
	for(Read r : reads2) readNames2.put(r.name, r);
	Scanner alnInput = new Scanner(new FileInputStream(new File(alnFn)));
	ArrayList<Alignment> res = new ArrayList<Alignment>();
	HashSet<String> alignedNames = new HashSet<String>();
	int idx = 0;
	while(alnInput.hasNext())
	{
		idx++;
		if(idx%100000 == 0) System.out.println(idx);
		String line = alnInput.nextLine();
		if(line.length() == 0 || line.charAt(0) == '@')
		{
			continue;
		}
		else
		{
			String name = line.substring(0, line.indexOf("\t"));
			if(readNames1.containsKey(name) && !alignedNames.contains(name))
			{
				alignedNames.add(name);
				String[] ss = line.split("\t");
				res.add(new Alignment(name, ss[2], Integer.parseInt(ss[3]), readNames1.get(name).seq));
				res.add(new Alignment(name, ss[2], Integer.parseInt(ss[3]), readNames2.get(name).seq));
			}
		}
	}
	for(String s : readNames1.keySet())
	{
		if(!alignedNames.contains(s))
		{
			res.add(new Alignment(s, "", -1, readNames1.get(s).seq));
			res.add(new Alignment(s, "", -1, readNames2.get(s).seq));
		}
	}
	return res;
}
/*
 * filterReads filters reads to contain only those sharing a sufficient number of k-mers with a SV
 * vcfFn is the filename of a VCF file containing an entry for each SV we are interested in
 * read1Fn is the filename of a file containing the first read in each paired-end read
 * read2Fn is the filename of a file containing the second read in each paired-end read
 * out1Fn is the filename to which the first read in each pair will be written after filtering
 * out2Fn is the filename to which the second read in each pair will be wirtten after filtering
 */
static void filterReads(String vcfFn, String read1Fn, String read2Fn, String out1Fn, String out2Fn) throws IOException
{
	Scanner vcfScanner = new Scanner(new FileInputStream(new File(vcfFn)));
	int varIndex = 0;
	while(vcfScanner.hasNext())
	{
		String line = vcfScanner.nextLine();
		if(line.length() == 0 || line.charAt(0) == '#') continue;
		varIndex++;
		PrintWriter out1 = new PrintWriter(new File(out1Fn+"."+varIndex));
		PrintWriter out2 = new PrintWriter(new File(out2Fn+"."+varIndex));
		SV sv = new SV(line);
		String s = sv.seq;
		kmers = kmerize(s);
		Scanner input1 = new Scanner(new FileInputStream(new File(read1Fn)));
		Scanner input2 = new Scanner(new FileInputStream(new File(read2Fn)));
		while(input1.hasNext() && input2.hasNext())
		{
			String readName = input1.nextLine().substring(1);
			input2.nextLine();
			String seq1 = input1.nextLine(), seq2 = input2.nextLine();
			input1.nextLine();
			input2.nextLine();
			String q1 = input1.nextLine(), q2 = input2.nextLine();
			int c1 = countKmers(seq1), c2 = countKmers(seq2);
			if(c1 >= sharedKmers || c2 >= sharedKmers)
			{
				out1.println(">" + readName);
				out1.println(seq1);
				out1.println('+');
				out1.println(q1);
				out2.println(">" + readName);
				out2.println(seq2);
				out2.println('+');
				out2.println(q2);
			}
		}
		out1.close();
		out2.close();
	}
}
/*
 * Reads a list of reads from a fastq file
 * fn is the name of the fastq file
 */
static ArrayList<Read> fromFasta(String fn) throws IOException
{
	Scanner input = new Scanner(new FileInputStream(new File(fn)));
	ArrayList<Read> res = new ArrayList<Read>();
	while(input.hasNext())
	{
		String name = input.nextLine().substring(1);
		String seq = input.nextLine();
		input.nextLine();
		String q = input.nextLine();
		res.add(new Read(name, seq, q));
	}
	return res;
}
/*
 * Counts the number of kmers in a string which are present in a set
 * s is the string being queried
 * kmers is a global variable containing the set of kmers
 */
static int countKmers(String s)
{
	int res = 0;
	long hash = 0;
	if(s.length() < k) return 0;
	for(int i = 0; i<k; i++) hash = hash << 2 + vals[s.charAt(i)];
	if(kmers.contains(hash)) res++;
	for(int i = k; i<s.length(); i++)
	{
		hash &= (1L << (2*(k-1))) - 1;
		hash <<= 2;
		hash += vals[s.charAt(i)];
		if(kmers.contains(hash)) res++;
	}
	return res;
}
/*
 * kmerize enocdes a kmer as a 64-bit long
 * s is the kmer to encode
 * Note: Only works for kmers up to about 30 base pairs
 */
static HashSet<Long> kmerize(String s)
{
	HashSet<Long> res = new HashSet<Long>();
	if(s.length() < k) return res;
	long hash = 0;
	for(int i = 0; i<k; i++) hash = hash << 2 + vals[s.charAt(i)];
	res.add(hash);
	for(int i = k; i<s.length(); i++)
	{
		hash &= (1L << (2*(k-1))) - 1;
		hash <<= 2;
		hash += vals[s.charAt(i)];
		res.add(hash);
	}
	return res;
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
 * A structure containing information about a structural variant (insertions only)
 * seq is the sequence being inserted
 * pos is the position in the reference where the SV occurs
 * chr is the chromosome on which the SV occurs
 */
static class SV
{
	String seq;
	int pos;
	String chr;
	SV(String line)
	{
		String[] tokens = line.split("\t");
		chr = tokens[0];
		pos = Integer.parseInt(tokens[1]);
		int idx = line.indexOf("SEQ=") + 4;
		int end = idx;
		while(end < line.length() && isBasePair(line.charAt(end))) end++;
		seq = line.substring(idx,end);
	}
}
/*
 * A structure containing information about a single genomic read
 * name is the name of the read
 * seq is the read sequence
 * qual is the string of Phred-encoded quality values
 */
static class Read
{
	String name, seq, qual;
	Read(String nn, String ss, String qq)
	{
		name = nn;
		seq = ss;
		qual = qq;
	}
}
/*
 * Writes alignments to a file
 * alns is a list of alignments
 * fn is the name of the file to be written to
 */
static void writeAlns(ArrayList<Alignment> alns, String fn) throws IOException
{
	System.out.println(fn);
	PrintWriter out = new PrintWriter(new File(fn));
	for(Alignment a : alns) out.println(a.toString());
	out.close();
}
/*
 * Reads alignments from a file
 * fn is the name of the file to read from
 */
static ArrayList<Alignment> readAlns(String fn) throws IOException
{
	ArrayList<Alignment> res = new ArrayList<Alignment>();
	Scanner input = new Scanner(new FileInputStream(new File(fn)));
	while(input.hasNext())
	{
		res.add(new Alignment(input.nextLine()));
	}
	return res;
}
/*
 * A structure containing information about an alignment of a read to a reference genome
 * name is the read name
 * chr is the chromosome of the reference being aligned to
 * pos is the position on the chromosome the read is aligned to
 * seq is the sequence of the read being aligned
 */
static class Alignment
{
	String name;
	String chr;
	int pos;
	String seq;
	Alignment(String nn, String cc, int pp, String ss)
	{
		name = nn;
		chr = cc;
		pos = pp;
		seq = ss;
	}
	Alignment(String s)
	{
		String[] tokens = s.split("\t");
		name = tokens[0];
		chr = tokens[1];
		pos = Integer.parseInt(tokens[2]);
		seq = tokens[3];
	}
	public String toString()
	{
		return name+"\t"+chr+"\t"+pos+"\t"+seq;
	}
}
}

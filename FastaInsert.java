import java.util.*;
import java.io.*;
public class FastaInsert {
public static void main(String[] args) throws IOException
{
	String fastaFn = args[0];
	String vcfFn = args[1];
	String outputFn = args[2];
	String[] ogFasta = readFasta(fastaFn);
	Insertion[] ins = readVcf(vcfFn);
	String[] newSeq = addInsertions(ogFasta, ins, false);
	printFasta(newSeq, outputFn);
}
static String[] addInsertions(String[] seq, Insertion[] ins, boolean pos)
{
	Arrays.sort(ins);
	StringBuilder sb = new StringBuilder();
	if(!pos)
	{
		String[] res = new String[seq.length+2];
		for(int i = 0; i<seq.length; i++) res[i] = seq[i];
		res[seq.length] = "insertion";
		for(Insertion in : ins) sb.append(in.seq);
		res[seq.length+1] = sb.toString();
		return res;
	}
	for(int j = 0; j<seq.length; j+= 2)
	{
		int at = 0;
		while(at < ins.length && !ins[at].chr.equals(seq[j])) at++;
		for(int i = 0; i < seq[j+1].length(); i++)
		{
			while(at < ins.length && ins[at].pos == i)
			{
				sb.append(ins[at].seq);
				at++;
				while(at < ins.length && !ins[at].chr.equals(seq[j])) at++;
			}
			sb.append(seq[j+1].charAt(i));
		}
		seq[j+1] = sb.toString();
		sb = new StringBuilder();
	}
	return seq;
}
static Insertion[] readVcf(String fn) throws IOException
{
	Scanner input = new Scanner(new FileInputStream(new File(fn)));
	ArrayList<Insertion> res = new ArrayList<Insertion>();
	while(input.hasNext())
	{
		String s = input.nextLine();
		if(s.startsWith("#")) continue;
		StringTokenizer str = new StringTokenizer(s, "\t");
		String chr = str.nextToken();
		int pos = Integer.parseInt(str.nextToken());
		String seq = s.substring(s.indexOf("SEQ=") + 4);
		seq = seq.substring(0, seq.indexOf(';'));
		res.add(new Insertion(chr, seq, pos));
	}
	int n = res.size();
	Insertion[] arr = new Insertion[n];
	for(int i = 0; i<n; i++) arr[i] = res.get(i);
	return arr;
}
static String[] readFasta(String fn) throws IOException
{
	Scanner input = new Scanner(new FileInputStream(new File(fn)));
	ArrayList<String> res = new ArrayList<String>();
	String name = input.nextLine().substring(1);
	StringBuilder sb = new StringBuilder();
	while(input.hasNext())
	{
		String s = input.nextLine();
		if(!s.startsWith(">")) sb.append(s);
		else
		{
			res.add(name);
			res.add(sb.toString());
			name = s.substring(1);
		}
	}
	res.add(name);
	res.add(sb.toString());
	return res.toArray(new String[res.size()]);
}

static void printFasta(String[] seqs, String fn) throws IOException
{
	PrintWriter out = new PrintWriter(new File(fn));
	for(int j = 0; j<seqs.length; j += 2)
	{
		String name = seqs[j], seq = seqs[j+1];
		out.println(">" + name);
		for(int i = 0; i<seq.length(); i += 80)
		{
			out.println(seq.substring(i, Math.min(seq.length(), i+80)));
		}
	}
	out.close();
}

static class Insertion implements Comparable<Insertion>
{
	String chr;
	String seq;
	int pos;
	Insertion(String cc, String ss, int pp)
	{
		chr = cc; seq = ss; pos = pp;
	}
	@Override
	public int compareTo(Insertion o) {
		// TODO Auto-generated method stub
		return pos - o.pos;
	}
}
}

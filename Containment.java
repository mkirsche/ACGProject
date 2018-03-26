/*
 * Gets score denoting likelihood that a is contained in b
 */
import java.util.*;
import java.io.*;
public class Containment {
	static int k = 10;
	static int[] vals;
public static void main(String[] args) throws IOException
{
	String fn = args[1];
	System.out.println(fn);
	vals = new int[256];
	Arrays.fill(vals, -1);
	vals['A'] = 0; vals['C'] = 1; vals['G'] = 2; vals['T'] = 3;
	char[] a = args[0].toCharArray(), b = readFromFasta(fn).toCharArray();
	HashMap<Long, Integer> asketch = sketch(a), bsketch = sketch(b);
	int overlap = overlap(asketch, bsketch);
	int size = size(asketch);
	System.out.println(1.0 * overlap / size);
}
static String readFromFasta(String fn) throws IOException
{
	StringBuilder sb = new StringBuilder();
	Scanner input = new Scanner(new FileInputStream(new File(fn)));
	while(input.hasNext())
	{
		String s = input.nextLine();
		if(s.startsWith(">")) sb.append("x");
		else sb.append(s);
	}
	return sb.toString();
}
static int overlap(HashMap<Long, Integer> sketch1, HashMap<Long, Integer> sketch2)
{
	int res = 0;
	for(long x : sketch1.keySet()) res += sketch2.containsKey(x) ? (Math.min(sketch1.get(x), sketch2.get(x))) : 0;
	return res;
}
static int size(HashMap<Long, Integer> sketch)
{
	int res = 0;
	for(long x : sketch.keySet()) res += sketch.get(x);
	return res;
}
static HashMap<Long, Integer> sketch(char[] s)
{
	HashMap<Long, Integer> res = new HashMap<Long, Integer>();
	for(int i = 0; i+k <= s.length; i++)
	{
		long x = kmerize(s, i);
		if(x != -1)
		{
			res.put(x, res.containsKey(x) ? (res.get(x) + 1) : 1);
		}
	}
	return res;
}
static long kmerize(char[] s, int start)
{
	boolean good = true;
	long res = 0;
	for(int i = 0; i<k; i++)
	{
		good &= s[i+start] != -1;
		if(!good) return -1;
		res = res * 4 + vals[s[i+start]];
	}
	return res;
	
}
}

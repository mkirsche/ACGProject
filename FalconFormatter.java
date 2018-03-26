/*
 *   Reformats FASTQ files to work with FalconSense
 */
import java.util.*;
import java.io.*;
public class FalconFormatter {
public static void main(String[] args) throws IOException
{
    String readsFile = args[0];
	String outputFile = args[1];
	Scanner readInput = new Scanner(new FileInputStream(new File(readsFile)));
	PrintWriter out = new PrintWriter(new File(outputFile));
    ArrayList<String> names = new ArrayList<String>();
	ArrayList<String> seqs = new ArrayList<String>();
	String name = "";
	StringBuilder sb = new StringBuilder("");
	while(readInput.hasNext())
	{
		names.add(readInput.nextLine().substring(1));
		seqs.add(readInput.nextLine());
		readInput.nextLine();
		readInput.nextLine();
	}
	int n = names.size();
	for(int i = 0; i<n; i++)
	{
		out.println(names.get(i)+" "+seqs.get(i));
		for(int j = 0; j<n; j++)
		{
			if(j != i) out.println(names.get(j)+" "+seqs.get(j));
		}
		out.println("+ +");
    }
	out.println("- -");
	out.close();
}
}

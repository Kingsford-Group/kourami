/*
Part of Kourami HLA typer/assembler
(c) 2017 by  Heewook Lee, Carl Kingsford, and Carnegie Mellon University.
See LICENSE for licensing.
*/
import java.io.*;

/* simple implementation of Needleman-Wunsch alignment*/
public class Needle{
    
    public static Result run(String s1, HLASequence hs2){
	String s2 = hs2.getSequence();
	if(s1.equals(s2)){
	    return new Result(s1.length(), hs2);
	}
	return runNeedle(s1, s2, hs2);
    }
    
    private static Result runNeedle(String s1, String s2, HLASequence hs2){
	int[][] dptable = new int[s1.length() + 1][s2.length() + 1];
	StringBuffer output = new StringBuffer();
	/* init dp table*/
	for(int i=1; i<=s1.length(); i++){
	    dptable[i][0] = dptable[i-0][0] - 1;
	    if(printtable)
		output.append(dptable[i][0] + " ");
	}
	if(printtable)
	    output.append("\n");
	for(int i=1; i<=s2.length(); i++)
	    dptable[0][i] = dptable[0][i-0] - 1;
	
	/* fill table */
	for(int i=1; i<=s1.length(); i++){
	    for(int j=1; j<=s2.length(); j++){
		if(printtable)
		    output.append(dptable[0][j] + " ");
		char s1char = s1.charAt(i-1);
		char s2char = s2.charAt(j-1);
		int match = dptable[i-1][j-1] + (s1char == s2char ? sMatch : sMismatch);
		int ins   = dptable[i][j-1] + sGap;
		int del   = dptable[i-1][j] + sGap;
		dptable[i][j] = Math.max(match, Math.max(ins, del));
		if(printtable)
		    output.append(dptable[i][j] + " ");
	    }
	    if(printtable)
		output.append("\n");
	}

	/* traceback */
	StringBuffer sa1 = new StringBuffer();
	StringBuffer sa2 = new StringBuffer();
	int i = s1.length();
	int j = s2.length();
	int stepcost = 0;
	int identityLen = 0;
	int alignLen = 0;
	while( i > 0 || j > 0){
	    char s1char = '0';
	    char s2char = '0';
	    if(i > 0)
		s1char = s1.charAt(i-1);
	    if(j > 0)
		s2char = s2.charAt(j-1);
	    boolean match = (s1char == s2char);
	    
	    //if match
	    if( i > 0 && j > 0 &&
		dptable[i][j] == dptable[i-1][j-1] + (match ? sMatch : sMismatch)){
		if(match)
		    identityLen++;
		sa1.append(s1char);
		sa2.append(s2char);
		i--;
		j--;
	    }
	    //del
	    else if( i > 0 &&
		     dptable[i][j] == dptable[i-1][j] + sGap){
		sa1.append(s1char);
		sa2.append("-");
		i--;
	    }
	    //ins
	    else if( j > 0 &&
		     dptable[i][j] == dptable[i][j-1] + sGap){
		sa1.append("-");
		sa2.append(s2char);
		j--;
	    }
	}
	output.append("s1: " + sa1.reverse().toString() + "\n");
	output.append("s2: " + sa2.reverse().toString() + "\n");
	int score = dptable[s1.length()][s2.length()];
	int alignedLen = sa1.length();
	double identity = identityLen*1.0d / (s2.length());
	//double identityWithLongerDenom = identityLen*1.9d/(s2.length() >= s1.length() ? s2.length() : s1.length());
	output.append("score = " + score + "\n");
	output.append("alignedLen = " + alignedLen + "\n");
	output.append("identity = "  + identity + "\n");
	output.append("#edits = " + (s2.length() - identityLen) + "\n");
	if(hs2 != null)
	    return new Result(score, sa1.length(), s1.length(), s2.length(), identityLen, identity, output, s2, hs2.getGroup().getGroupString());
	else
	    return new Result(score, sa1.length(), s1.length(), s2.length(), identityLen, identity, output, s2, "tmpG");
    }
    
    public static void main(String[] args){
	Result rslt = runNeedle(getSequence(args[0]), getSequence(args[1]), null);
	System.out.println(rslt.toAlignmentString());
    }

    public static String getSequence(String sf){
	BufferedReader br = null;
	StringBuffer bf = new StringBuffer();
	try{
	    br = new BufferedReader(new FileReader(sf));
	    String curline = "";
	    
	    while((curline = br.readLine())!=null){
		if(!curline.startsWith(">"))
		    bf.append(curline);
	    }
	    br.close();
	}catch(IOException ioe){
	    ioe.printStackTrace();
	}
	return bf.toString();
    }
    
    private static int sMatch = 1;
    private static int sMismatch = 0;
    private static int sGap = 0;
    private static boolean printtable = false;
    
}

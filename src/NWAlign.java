//--------------------------------------------------------------------------------------------------
// 
// Description: Implementation of Needleman-Wunsch global alignment.
// This code is written by Ren-Xiang Yan in China Agricultural University and is originally based on 
// the fortran implementation from Dr. Yang Zhang (http://zhanglab.ccmb.med.umich.edu/NW-align/).
// Last update is in 2010/08/14. 
//
//  Usage:
//      java -jar NWAlign.jar F1.fasta F2.fasta  (align two sequences in fasta file)
//java -jar NWAlign.jar F1.pdb F2.pdb    1 (align two sequences in PDB file)
//java -jar NWAlign.jar F.fasta F.pdb  2 (align sequences 1 in fasta and 1 in pdb)
//java -jar NWAlign.jar GKDGL EVADELVSE    3 (align two sequences in plain text)
//java -jar NWAlign.jar GKDGL F.fasta  4 (align sequences 1 in text and 1 in fasta)
//java -jar NWAlign.jar GKDGL F.pdb    5 (align sequences 1 in text and 1 in pdb)
//  
//   Note: You also could complied the code by yourself.
//         Decompress the NWAlign.jar file and you can get the source code in the NWAlign folder.
//   The program can be compiled by 
//              javac NWAlign.java
//   Then you could use the program by the following commands:
//java NWAlign F1.fasta F2.fasta  (align two sequences in fasta file)
//java NWAlign F1.pdb F2.pdb    1 (align two sequences in PDB file)
//java NWAlign file1.fasta file2.pdb  2 (align sequences 1 in fasta and 1 in pdb)
//java NWAlign GKDGL EVADELVSE    3 (align two sequences in plain text)
//java NWAlign GKDGL F.fasta  4 (align sequences 1 in text and 1 in fasta)
//java NWAlign GKDGL F.pdb    5 (align sequences 1 in text and 1 in pdb)            
//-----------------x-------------------x-------------------x--------------------------------------------

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.text.DecimalFormat;

public class NWAlign {


    public static Result runDefault(String f1, HLASequence hs2){
	String f2 = hs2.getSequence();
	if(f1.equals(f2)){
	    return new Result(f1.length(), hs2);
	}
	int gap_open=-11;
	int gap_extn=-1;
	return NeedlemanWunsch(f1.toUpperCase(),f2.toUpperCase(),gap_open,gap_extn,hs2); 
    }
    
    public static Result NeedlemanWunsch(String f1,String f2,int gap_open,int gap_extn, HLASequence hs2)  
    {
	//int[][] imut = new int[24][24];                      
	int[][] imut = new int[5][5];
	//Blosum62Matrix(imut);                              // Read Blosum scoring matrix and store it in the imut variable.
	ScoringMatrix(imut);
	//String seqW = "*ARNDCQEGHILKMFPSTWYVBZX";       // Amino acide order in the BLAST's scoring matrix (e.g.,Blosum62). 
	String seqW = "*ACGT";
	f1 = "*" + f1;                                     // Add a '*' character in the head of a sequence and this can make java code much more consistent with orginal fortran code.   
	f2 = "*" + f2;                                     // Use 1 to represent the first position of the sequence in the original fortran code,and 1 stand for the second position in java code. Here, add a '*' character in the head of a sequence could make 1 standard for the first postion of thse sequence in java code.     
	int[] seq1 = new int[f1.length()];                 
	int[] seq2 = new int[f2.length()];           // seq1 and seq2 are arrays that store the amino acid order numbers of sequence1 and sequence2.

	int i,j;
	for(i=1;i<f1.length();i++){
	    if(f1.charAt(i) == 'A')
		seq1[i] = 1;
	    else if(f1.charAt(i) == 'C')
		seq1[i] = 2;
	    else if(f1.charAt(i) == 'G')
		seq1[i] = 3;
	    else if(f1.charAt(i) == 'T')
		seq1[i] = 4;
	}
	
	for(i=1;i<f2.length();i++){
	    if(f2.charAt(i) == 'A')
		seq2[i] = 1;
	    else if(f2.charAt(i) == 'C')
		seq2[i] = 2;
	    else if(f2.charAt(i) == 'G')
		seq2[i] = 3;
	    else if(f2.charAt(i) == 'T')
		seq2[i] = 4;
	}
	
	/*
	int i,j;                                   // For example, 1 stand for A, 2 represent R and etc.
	for(i=1;i<f1.length();i++)
	    {
		for(j=1;j<seqW.length();j++)
		    {
			if(f1.charAt(i)==seqW.charAt(j))
			    {
				seq1[i]=j;
			    }
		    }
	    }
	
	for(i=1;i<f2.length();i++)
	    {
		for(j=1;j<seqW.length();j++)
		    {
			if(f2.charAt(i)==seqW.charAt(j))
			    {
				seq2[i]=j;
			    }
		    }
	    }
	*/
	int[][] score = new int[f1.length()][f2.length()];// score[i][j] stard for the alignment score that align ith position of the first sequence to the jth position of the second sequence.
	for(i=1;i<f1.length();i++)
	    {
		for(j=1;j<f2.length();j++)
		    {
			score[i][j] = imut[seq1[i]][seq2[j]];
		    }
	    }
	
	int[] j2i = new int[f2.length()+1];
	for(j=1;j<f2.length();j++)
	    {
		j2i[j] = -1; // !all are not aligned
	    }
	
	
	
	int[][] val = new int[f1.length()+1][f2.length()+1];  // val[][] was assigned as a global variable, and the value could be printed in the final.
	int[][] idir = new int[f1.length()+1][f2.length()+1];
	int[][] preV = new int[f1.length()+1][f2.length()+1];
	int[][] preH = new int[f1.length()+1][f2.length()+1];
	int D,V,H;
	
        // If you want to use alternative implementation of Needleman-Wunsch dynamic program , you can assign "false" value to the "standard" variable.  
	    
	////////////////////////////////////////////////////////////////////////////////
	//     This is a standard Needleman-Wunsch dynamic program (by Y. Zhang 2005).
	//     1. Count multiple-gap.
	//     2. The gap penality W(k)=Go+Ge*k1+Go+Ge*k2 if gap open on both sequences
	//     idir[i][j]=1,2,3, from diagonal, horizontal, vertical
	//     val[i][j] is the cumulative score of (i,j)
	////////////////////////////////////////////////////////////////////////////////
	
	int[][] jpV = new int[f1.length()+1][f2.length()+1];
	int[][] jpH = new int[f1.length()+1][f2.length()+1];
	val[0][0]=0;
	val[1][0] =gap_open;
	for(i=2;i<f1.length();i++){
	    val[i][0] = val[i-1][0]+gap_extn;
	}
	for(i=1;i<f1.length();i++)
	    {
		
		preV[i][0] = val[i][0]; // not use preV at the beginning
		idir[i][0] = 0;         // useless
		jpV[i][0] = 1;          // useless
		jpH[i][0] = i;          // useless
	    }
	val[0][1]=gap_open;
	for(j=2;j<f2.length();j++){
	    val[0][j]=val[0][j-1]+gap_extn;
	}
	for(j=1;j<f2.length();j++)
	    {
		preH[0][j] = val[0][j];
		idir[0][j] = 0;
		jpV[0][j] = j;
		jpH[0][j] = 1;
	    }
	
	// DP ------------------------------>
	for(j=1;j<f2.length();j++)
	    {
		for(i=1;i<f1.length();i++)
		    {
			// D=VAL(i-1,j-1)+SCORE(i,j)--------------->
			D = val[i-1][j-1] + score[i][j];// from diagonal, val(i,j) is val(i-1,j-1)
			
			//H=H+gap_open ------->
			jpH[i][j] = 1;
			int val1 = val[i-1][j] + gap_open;  // gap_open from both D and V
			int val2 = preH[i-1][j] + gap_extn; // gap_extn from horizontal
			if(val1>val2)   // last step from D or V
			    {
				H = val1;
			    }
			else            // last step from H
			    {
				H = val2;
				if(i > 1)
				    {
					jpH[i][j] = jpH[i-1][j] + 1;  // record long-gap
				    }
			    }
			
			//V=V+gap_open --------->
			jpV[i][j] = 1;
			val1 = val[i][j-1] + gap_open;
			val2 = preV[i][j-1] + gap_extn;
			if(val1>val2)
			    {
				V = val1;
			    }
			else
			    {
				V = val2;
				if(j > 1)
				    {
					jpV[i][j] = jpV[i][j-1] + 1;   // record long-gap   
				    }
			    }
			
			preH[i][j] = H; // unaccepted H
			preV[i][j] = V;// unaccepted V
			if((D>H)&&(D>V))
			    {
				idir[i][j]=1;
				val[i][j]=D;
			    }
			else if(H > V)
			    {   
				idir[i][j] = 2;
				val[i][j] = H;
			    }
			else
			    {
				idir[i][j] = 3;
				val[i][j] = V;
			    }
		    }
	    }
	
	//  tracing back the pathway
	i = f1.length()-1;
	j = f2.length()-1;
	while((i>0)&&(j>0))  
	    { 
		if(idir[i][j]==1)       // from diagonal
		    {
			j2i[j] = i;
			i=i-1;
			j=j-1;
		    }
		else if(idir[i][j]==2)  // from horizonal
		    {          
			int temp1 = jpH[i][j];                  //  
			for(int me=1;me<=temp1;me++)            //  In the point view of a programer, 
			    {                                       //  you should not use the  "for(int me=1;me<=jpH[i][j];me++)".
				if(i>0)                         //  If you use up sentence,the value of jpH[i][j] is changed when variable i changes. 
				    {                                    //  So the value of jpH[i][j] was assigned to the value temp1 and use the setence "for(int me=1;me<=temp1;me++)" here. 
					i=i-1;                            // 
				    }                                 //
			    }                                        
		    }
			else
			    { 
				int temp2 = jpV[i][j]; 
				for(int me=1;me<=temp2;me++)             //  In the point view of a programer,
				    {                                        //  you should not use the  "for(int me=1;me<=jpV[i][j];me++)".
					if(j>0)                               //  Because when variable i change, the jpV[i][j] employed here is also change. 
					    {                                     //  So the value of jpV[i][j] was assigned to the value temp2 and use the setence "for(int me=1;me<=temp2;me++)" here.
						j=j-1;                             //
					    }
				    }           
			    } 
	    }
	
	
	// calculate sequence identity
	int L_id=0;
	int L_ali=0;
	for(j=1;j<f2.length();j++)
	    {    
		if(j2i[j]>0)
		    {
			i=j2i[j];
			L_ali=L_ali+1;
			if(seq1[i]==seq2[j])
			    {            
				L_id=L_id+1;
			    }
		    }         
	    }   
	
	double identity = L_id*1.0/(f2.length()-1);    
	int fina_score = val[f1.length()-1][f2.length()-1];
	//return new Result(fina_score, L_ali, f1.length()-1, f2.length()-1, L_id, identity);
	
	StringBuffer output = new StringBuffer();

	output.append("Alignment score=" + fina_score + "\n"); 
	output.append("Length of sequence 1:" + (f1.length()-1) + "\n"); 
	output.append("Length of sequence 2:" + (f2.length()-1) + "\n"); 
	output.append("Aligned length      :" + L_ali + "\n"); 
	output.append("Identical length    :" + L_id + "\n");	
	DecimalFormat df = new DecimalFormat("0.000");      // Correct the identity to 3 decimal places. 
	output.append("Sequence identity=" + df.format(identity));
	double identityWithLongerDenom = L_id*1.0/(f2.length() >= f1.length() ? f2.length() - 1 : f1.length() - 1);
	if(identityWithLongerDenom < identity)
	    output.append("Sequence identity(longer denom)=" + df.format(identityWithLongerDenom));
	output.append(" " + L_id  + "/" + (f2.length()-1) + "\n\n");


	// output aligned sequences    
	char[] sequenceA = new char[f1.length()+f2.length()];
	char[] sequenceB = new char[f1.length()+f2.length()];
	char[] sequenceM = new char[f1.length()+f2.length()];    
	int k = 0;
	i=1;
	j=1;    
	while(true)
	    {    
		if((i>(f1.length()-1))&&(j>(f2.length()-1)))
		    break;    
		if((i>(f1.length()-1))&&(j<(f2.length()-1)))     // unaligned C on 1
		    {
			k = k + 1;
			sequenceA[k] = '-';
			sequenceB[k] = seqW.charAt(seq2[j]);
			sequenceM[k] = ' ';
			j = j + 1;
		    }    
		else if((i<(f1.length()-1))&&(j>(f2.length()-1))) // unaligned C on 2
		    {        
			k = k + 1;
			sequenceA[k] = seqW.charAt(seq1[i]);
			sequenceB[k] = '-';
			sequenceM[k] = ' ';
			i = i + 1;
		    }        
		else if(i==j2i[j]) // if align
		    {
			k = k + 1;
			sequenceA[k] = seqW.charAt(seq1[i]);
			sequenceB[k] = seqW.charAt(seq2[j]);
			if(seq1[i]==seq2[j])  // identical
			    {
				sequenceM[k] = ':';
			    }
			else
			    {
				sequenceM[k] = ' ';
			    }
			i = i + 1;
			j = j + 1;
		    }        
		else if(j2i[j]<0)   // gap on 1
		    {
			k = k + 1;
			sequenceA[k] = '-';
			sequenceB[k] = seqW.charAt(seq2[j]);
			sequenceM[k] = ' ';
			j = j + 1;
		    }        
		else if(j2i[j] >= 0)  // gap on 2
		    {
			k = k + 1;
			sequenceA[k] = seqW.charAt(seq1[i]);
			sequenceB[k] = '-';
			sequenceM[k] = ' ';
			i = i + 1;
		    }        
	    }   
	for(i=1;i<=k;i++)
	    {
		output.append(sequenceA[i]);   
	    }
	output.append("\n");//HLA.log.appendln();
	for(i=1;i<=k;i++)
	    {
		output.append(sequenceM[i]);   
	    }
	output.append("\n");
	for(i=1;i<=k;i++)
	    {
		output.append(sequenceB[i]);   
	    }
	output.append("\n");
	for(i=1;i<=k;i++)
	    {
		int temp = i%9;
		output.append(temp);//HLA.log.append(temp);   
	    }
	output.append("\n");///HLA.log.appendln();   
	
	return new Result(fina_score, L_ali, f1.length()-1, f2.length()-1, L_id, identity, output, f2, hs2.getGroup().getGroupString());

    }
    
    
    
    public static void ScoringMatrix(int[][] imut){
	imut[1][1] = 1;
	imut[1][2] = -1;
	imut[1][3] = -1;
	imut[1][4] = -1;
	imut[2][1] = -1;
	imut[2][2] = 1;
	imut[2][3] = -1;
	imut[2][4] = -1;
	imut[3][1] = -1;
	imut[3][2] = -1;
	imut[3][3] = 1;
	imut[3][4] = -1;
	imut[4][1] = -1;
	imut[4][2] = -1;
	imut[4][3] = -1;
	imut[4][4] = 1;
    }

}

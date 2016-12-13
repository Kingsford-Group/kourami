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

    public static void main(String[] args) {
	
	String f1 = "";
	String f2 = "";
	if(args.length==2){                  // align two sequences in fasta file
	    f1 = readFastaOrRawSequence(args[0]);
	    f2 = readFastaOrRawSequence(args[1]);
	}
	else if(args.length==3)
	    {
		if(args[2].equals("1")){         // align two sequences in PDB file
		    f1 = readPDB(args[0]);
		    f2 = readPDB(args[1]);
		}else if(args[2].equals("2")){   // align sequences 1 in fasta and 1 in pdb
		    f1 = readFastaOrRawSequence(args[0]);
		    f2 = readPDB(args[1]);
		}else if(args[2].equals("3")){   // align two sequences in plain text        
		    f1 = args[0];
		    f2 = args[1];
		}else if(args[2].equals("4")){   //  align sequences 1 in text and 1 in fasta
		    f1 = args[0];
		    f2 = readFastaOrRawSequence(args[1]);
		}else if(args[2].equals("5")){   //  align sequences 1 in text and 1 in pdb
		    f1 = args[0];
		    f2 = readPDB(args[1]);
		}
	    }
	
	if(args.length<2){  // if args.lengh<2 then print the usage in the screen
	    HLA.log.appendln("java -jar NWAlign.jar F1.fasta F2.fasta  (align two sequences in fasta file)");
	    HLA.log.appendln("java -jar NWAlign.jar F1.pdb F2.pdb    1 (align two sequences in PDB file)");
	    HLA.log.appendln("java -jar NWAlign.jar F1.fasta F2.pdb  2 (align sequences 1 in fasta and 1 in pdb)");
	    HLA.log.appendln("java -jar NWAlign.jar GKDGL EVADELVSE    3 (align two sequences in plain text)");
	    HLA.log.appendln("java -jar NWAlign.jar GKDGL F.fasta  4 (align sequences 1 in text and 1 in fasta)");
	    HLA.log.appendln("java -jar NWAlign.jar GKDGL F.pdb    5 (align sequences 1 in text and 1 in pdb)");
	    System.exit(1);
	}
	
	
	int gap_open=-11,gap_extn=-1;
	
	// input upper case sequences
	NeedlemanWunsch(f1.toUpperCase(),f2.toUpperCase(),gap_open,gap_extn);
	    
    }
    
    public static Result runDefault(String f1, String f2){
	if(f1.equals(f2)){
	    return new Result(f1.length(), f2);
	}
	int gap_open=-11;
	int gap_extn=-1;
	return NeedlemanWunsch(f1.toUpperCase(),f2.toUpperCase(),gap_open,gap_extn); 
    }



    public static Result NeedlemanWunsch(String f1,String f2,int gap_open,int gap_extn)  
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
	boolean standard = true;
	if(standard)    // If you want to use alternative implementation of Needleman-Wunsch dynamic program , you can assign "false" value to the "standard" variable.  
	    {
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
	    }
	else   
	    {
		/////////////////////////////////////////////////////////////////////////////////
		//     This is an alternative implementation of Needleman-Wunsch dynamic program 
		//     (by Y. Zhang 2005)
		//     1. Count two-layer iteration and multiple-gaps
		//     2. The gap penality W(k)=Go+Ge*k1+Ge*k2 if gap open on both sequences
		//
		//     idir[i][j]=1,2,3, from diagonal, horizontal, vertical
		//     val[i][j] is the cumulative score of (i,j)
		////////////////////////////////////////////////////////////////////////////////
		
		int[][] preD = new int[f1.length()+1][f2.length()+1];
		int[][] idirH = new int[f1.length()+1][f2.length()+1];
		int[][] idirV = new int[f1.length()+1][f2.length()+1];
		val[0][0] = 0;
		for(i=1;i<f1.length();i++)
		    {
			val[i][0] = 0;
			idir[i][0] = 0;
			preD[i][0] = 0;
			preH[i][0] = -1000;
			preV[i][0] = -1000;
		    }
		
		for(j=1;j<f2.length();j++)
		    {
			val[0][j] = 0;
			idir[0][j] = 0;
			preD[0][j] = 0;
			preH[0][j] = -1000;
			preV[0][j] = -1000;
		    }
		
		// DP ------------------------------>
		for(j=1;j<f2.length();j++)
		    {
			for(i=1;i<f1.length();i++)
			    {
				// preD=VAL(i-1,j-1)+SCORE(i,j)--------------->
				preD[i][j] = val[i-1][j-1] + score[i][j];
				// preH: pre-accepted H----------------------->
				D = preD[i-1][j] + gap_open;
				H = preH[i-1][j] + gap_extn;
				V = preV[i-1][j] + gap_extn;
				if((D>H)&&(D>V))
				    {
					preH[i][j] = D;
					idirH[i-1][j] = 1;
				    }
				else if(H>V)
				    {
					preH[i][j] = H;
					idirH[i-1][j] = 2;
				    }
				else
				    {
					preH[i][j] = V;
					idirH[i-1][j] = 3;
				    }
				
				// preV: pre-accepted V----------------------->
				D = preD[i][j-1] + gap_open;
				H = preH[i][j-1] + gap_extn;
				V = preV[i][j-1] + gap_extn;
				if((D>H)&&(D>V))
				    {
					preV[i][j] = D;
					idirV[i][j-1] = 1;
				    }
				else if(H>V)
				    {
					preV[i][j] = H;
					idirV[i][j-1] = 2;
				    }
				else
				    {
					preV[i][j] = V;
					idirV[i][j-1] = 3;
				    }
				
				// decide idir(i,j)----------->
				if((preD[i][j]>preH[i][j])&&(preD[i][j]>preV[i][j]))
				    {
					idir[i][j] = 1;
					val[i][j] = preD[i][j];
				    }
				else if(preH[i][j]>preV[i][j])
				    {
					idir[i][j] = 2;
					val[i][j] = preH[i][j];
				    }
				else
				    {
					idir[i][j] = 3;
					val[i][j] = preV[i][j];
				    }
			    }
		    }
		
		//  tracing back the pathway
		i = f1.length()-1;
		j = f2.length()-1;
		while((i>0)&&(j>0))  
		    { 
			if(idir[i][j]==1)
			    {
				j2i[j] = i;
				i = i - 1;
				j = j - 1;
			    }
			else if(idir[i][j]==2)
			    {
				i = i - 1;
				idir[i][j] = idirH[i][j];
			    }
			else
			    {
				j = j - 1;
				idir[i][j] = idirV[i][j];
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
	
	return new Result(fina_score, L_ali, f1.length()-1, f2.length()-1, L_id, identity, output, f2);
	/*
	HLA.log.appendln("Alignment score=" + fina_score);
	HLA.log.appendln("Length of sequence 1:" + (f1.length()-1));
	HLA.log.appendln("Length of sequence 2:" + (f2.length()-1));
	HLA.log.appendln("Aligned length      :" + L_ali);
	HLA.log.appendln("Identical length    :" + L_id);
	    
	DecimalFormat df = new DecimalFormat("0.000");      // Correct the identity to 3 decimal places. 
	HLA.log.append("Sequence identity=" + df.format(identity));
	HLA.log.appendln(" " + L_id  + "/" + (f2.length()-1));    
	HLA.log.appendln();
	    
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
		HLA.log.append(sequenceA[i]);   
	    }
	HLA.log.appendln();
	for(i=1;i<=k;i++)
	    {
		HLA.log.append(sequenceM[i]);   
	    }
	HLA.log.appendln();
	for(i=1;i<=k;i++)
	    {
		HLA.log.append(sequenceB[i]);   
	    }
	HLA.log.appendln();
	for(i=1;i<=k;i++)
	    {
		int temp = i%9;
		HLA.log.append(temp);   
	    }
	HLA.log.appendln();   
	*/
    }
    
    public static String readFastaOrRawSequence(String file)  // read a sequence from a Fasta file or a text file.
    {
	String seq = "";
	String line = "";
	BufferedReader br;
	try {
	    br = new BufferedReader(new FileReader(file));
	    while ((line = br.readLine()) != null) {
		if(line.startsWith(">")){
		    
		}else{
		    seq = seq + line;
		}
	    }
	} catch (FileNotFoundException e) {
	    // TODO Auto-generated catch block
	    e.printStackTrace();
	} catch (IOException e) {
	    // TODO Auto-generated catch block
	    e.printStackTrace();
	}
	return seq;
    }
    
    public static String readPDB(String file)        // read a sequence from a PDB file
    {
	String seq = "";
	String line = "";
	BufferedReader br;
	try {
	    br = new BufferedReader(new FileReader(file));
	    while ((line = br.readLine()) != null) {
		if(line.startsWith("TER"))
		    break;
		if(line.startsWith("ATO")){
		    //HLA.log.appendln(line);
		    if(line.substring(13, 16).replaceAll("\\s+", "").endsWith("CA")){
			seq = seq + NameMap(line.substring(17, 20).toUpperCase());
		    }
		}
	    }
	} catch (FileNotFoundException e) {
	    // TODO Auto-generated catch block
	    e.printStackTrace();
	} catch (IOException e) {
	    // TODO Auto-generated catch block
	    e.printStackTrace();
	}
	return seq;
    }
    
    public static String NameMap(String residule)     // Map a three-letter abbreviation to a single-letter code.  
    {
	String[] aa = new String[]{"ALA","ARG","ASN","ASP","CYS","GLN","GLU","GLY","HIS","ILE","LEU","LYS","MET","PHE","PRO","SER","THR","TRP","TYR","VAL","ASX","GLX","UNK"};
	String[] aaName = new String[]{"A","R","N","D","C","Q","E","G","H","I","L","K","M","F","P","S","T","W","Y","V","B","Z","X"};
	int i = 0;
	for(;i<aa.length;i++)
	    {
		if(aa[i].equals(residule))
		    break;
	    }         
	return aaName[i];
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

    public static void Blosum62Matrix(int[][] imut)     // Folowing from BLOSUM62 used in BLAST.
    {                                            // This was directly copy from original fortran code.
	imut[1][1]=4;                                   // b,z,x are additional
	imut[1][2]=-1;
	imut[1][3]=-2;
	imut[1][4]=-2;
	imut[1][5]=0;
	imut[1][6]=-1;
	imut[1][7]=-1;
	imut[1][8]=0;
	imut[1][9]=-2;
	imut[1][10]=-1;
	imut[1][11]=-1;
	imut[1][12]=-1;
	imut[1][13]=-1;
	imut[1][14]=-2;
	imut[1][15]=-1;
	imut[1][16]=1;
	imut[1][17]=0;
	imut[1][18]=-3;
	imut[1][19]=-2;
	imut[1][20]=0;
	imut[1][21]=-2;
	imut[1][22]=-1;
	imut[1][23]=0;
	imut[2][1]=-1;
	imut[2][2]=5;
	imut[2][3]=0;
	imut[2][4]=-2;
	imut[2][5]=-3;
	imut[2][6]=1;
	imut[2][7]=0;
	imut[2][8]=-2;
	imut[2][9]=0;
	imut[2][10]=-3;
	imut[2][11]=-2;
	imut[2][12]=2;
	imut[2][13]=-1;
	imut[2][14]=-3;
	imut[2][15]=-2;
	imut[2][16]=-1;
	imut[2][17]=-1;
	imut[2][18]=-3;
	imut[2][19]=-2;
	imut[2][20]=-3;
	imut[2][21]=-1;
	imut[2][22]=0;
	imut[2][23]=-1;
	imut[3][1]=-2;
	imut[3][2]=0;
	imut[3][3]=6;
	imut[3][4]=1;
	imut[3][5]=-3;
	imut[3][6]=0;
	imut[3][7]=0;
	imut[3][8]=0;
	imut[3][9]=1;
	imut[3][10]=-3;
	imut[3][11]=-3;
	imut[3][12]=0;
	imut[3][13]=-2;
	imut[3][14]=-3;
	imut[3][15]=-2;
	imut[3][16]=1;
	imut[3][17]=0;
	imut[3][18]=-4;
	imut[3][19]=-2;
	imut[3][20]=-3;
	imut[3][21]=3;
	imut[3][22]=0;
	imut[3][23]=-1;
	imut[4][1]=-2;
	imut[4][2]=-2;
	imut[4][3]=1;
	imut[4][4]=6;
	imut[4][5]=-3;
	imut[4][6]=0;
	imut[4][7]=2;
	imut[4][8]=-1;
	imut[4][9]=-1;
	imut[4][10]=-3;
	imut[4][11]=-4;
	imut[4][12]=-1;
	imut[4][13]=-3;
	imut[4][14]=-3;
	imut[4][15]=-1;
	imut[4][16]=0;
	imut[4][17]=-1;
	imut[4][18]=-4;
	imut[4][19]=-3;
	imut[4][20]=-3;
	imut[4][21]=4;
	imut[4][22]=1;
	imut[4][23]=-1;
	imut[5][1]=0;
	imut[5][2]=-3;
	imut[5][3]=-3;
	imut[5][4]=-3;
	imut[5][5]=9;
	imut[5][6]=-3;
	imut[5][7]=-4;
	imut[5][8]=-3;
	imut[5][9]=-3;
	imut[5][10]=-1;
	imut[5][11]=-1;
	imut[5][12]=-3;
	imut[5][13]=-1;
	imut[5][14]=-2;
	imut[5][15]=-3;
	imut[5][16]=-1;
	imut[5][17]=-1;
	imut[5][18]=-2;
	imut[5][19]=-2;
	imut[5][20]=-1;
	imut[5][21]=-3;
	imut[5][22]=-3;
	imut[5][23]=-2;
	imut[6][1]=-1;
	imut[6][2]=1;
	imut[6][3]=0;
	imut[6][4]=0;
	imut[6][5]=-3;
	imut[6][6]=5;
	imut[6][7]=2;
	imut[6][8]=-2;
	imut[6][9]=0;
	imut[6][10]=-3;
	imut[6][11]=-2;
	imut[6][12]=1;
	imut[6][13]=0;
	imut[6][14]=-3;
	imut[6][15]=-1;
	imut[6][16]=0;
	imut[6][17]=-1;
	imut[6][18]=-2;
	imut[6][19]=-1;
	imut[6][20]=-2;
	imut[6][21]=0;
	imut[6][22]=3;
	imut[6][23]=-1;
	imut[7][1]=-1;
	imut[7][2]=0;
	imut[7][3]=0;
	imut[7][4]=2;
	imut[7][5]=-4;
	imut[7][6]=2;
	imut[7][7]=5;
	imut[7][8]=-2;
	imut[7][9]=0;
	imut[7][10]=-3;
	imut[7][11]=-3;
	imut[7][12]=1;
	imut[7][13]=-2;
	imut[7][14]=-3;
	imut[7][15]=-1;
	imut[7][16]=0;
	imut[7][17]=-1;
	imut[7][18]=-3;
	imut[7][19]=-2;
	imut[7][20]=-2;
	imut[7][21]=1;
	imut[7][22]=4;
	imut[7][23]=-1;
	imut[8][1]=0;
	imut[8][2]=-2;
	imut[8][3]=0;
	imut[8][4]=-1;
	imut[8][5]=-3;
	imut[8][6]=-2;
	imut[8][7]=-2;
	imut[8][8]=6;
	imut[8][9]=-2;
	imut[8][10]=-4;
	imut[8][11]=-4;
	imut[8][12]=-2;
	imut[8][13]=-3;
	imut[8][14]=-3;
	imut[8][15]=-2;
	imut[8][16]=0;
	imut[8][17]=-2;
	imut[8][18]=-2;
	imut[8][19]=-3;
	imut[8][20]=-3;
	imut[8][21]=-1;
	imut[8][22]=-2;
	imut[8][23]=-1;
	imut[9][1]=-2;
	imut[9][2]=0;
	imut[9][3]=1;
	imut[9][4]=-1;
	imut[9][5]=-3;
	imut[9][6]=0;
	imut[9][7]=0;
	imut[9][8]=-2;
	imut[9][9]=8;
	imut[9][10]=-3;
	imut[9][11]=-3;
	imut[9][12]=-1;
	imut[9][13]=-2;
	imut[9][14]=-1;
	imut[9][15]=-2;
	imut[9][16]=-1;
	imut[9][17]=-2;
	imut[9][18]=-2;
	imut[9][19]=2;
	imut[9][20]=-3;
	imut[9][21]=0;
	imut[9][22]=0;
	imut[9][23]=-1;
	imut[10][1]=-1;
	imut[10][2]=-3;
	imut[10][3]=-3;
	imut[10][4]=-3;
	imut[10][5]=-1;
	imut[10][6]=-3;
	imut[10][7]=-3;
	imut[10][8]=-4;
	imut[10][9]=-3;
	imut[10][10]=4;
	imut[10][11]=2;
	imut[10][12]=-3;
	imut[10][13]=1;
	imut[10][14]=0;
	imut[10][15]=-3;
	imut[10][16]=-2;
	imut[10][17]=-1;
	imut[10][18]=-3;
	imut[10][19]=-1;
	imut[10][20]=3;
	imut[10][21]=-3;
	imut[10][22]=-3;
	imut[10][23]=-1;
	imut[11][1]=-1;
	imut[11][2]=-2;
	imut[11][3]=-3;
	imut[11][4]=-4;
	imut[11][5]=-1;
	imut[11][6]=-2;
	imut[11][7]=-3;
	imut[11][8]=-4;
	imut[11][9]=-3;
	imut[11][10]=2;
	imut[11][11]=4;
	imut[11][12]=-2;
	imut[11][13]=2;
	imut[11][14]=0;
	imut[11][15]=-3;
	imut[11][16]=-2;
	imut[11][17]=-1;
	imut[11][18]=-2;
	imut[11][19]=-1;
	imut[11][20]=1;
	imut[11][21]=-4;
	imut[11][22]=-3;
	imut[11][23]=-1;
	imut[12][1]=-1;
	imut[12][2]=2;
	imut[12][3]=0;
	imut[12][4]=-1;
	imut[12][5]=-3;
	imut[12][6]=1;
	imut[12][7]=1;
	imut[12][8]=-2;
	imut[12][9]=-1;
	imut[12][10]=-3;
	imut[12][11]=-2;
	imut[12][12]=5;
	imut[12][13]=-1;
	imut[12][14]=-3;
	imut[12][15]=-1;
	imut[12][16]=0;
	imut[12][17]=-1;
	imut[12][18]=-3;
	imut[12][19]=-2;
	imut[12][20]=-2;
	imut[12][21]=0;
	imut[12][22]=1;
	imut[12][23]=-1;
	imut[13][1]=-1;
	imut[13][2]=-1;
	imut[13][3]=-2;
	imut[13][4]=-3;
	imut[13][5]=-1;
	imut[13][6]=0;
	imut[13][7]=-2;
	imut[13][8]=-3;
	imut[13][9]=-2;
	imut[13][10]=1;
	imut[13][11]=2;
	imut[13][12]=-1;
	imut[13][13]=5;
	imut[13][14]=0;
	imut[13][15]=-2;
	imut[13][16]=-1;
	imut[13][17]=-1;
	imut[13][18]=-1;
	imut[13][19]=-1;
	imut[13][20]=1;
	imut[13][21]=-3;
	imut[13][22]=-1;
	imut[13][23]=-1;
	imut[14][1]=-2;
	imut[14][2]=-3;
	imut[14][3]=-3;
	imut[14][4]=-3;
	imut[14][5]=-2;
	imut[14][6]=-3;
	imut[14][7]=-3;
	imut[14][8]=-3;
	imut[14][9]=-1;
	imut[14][10]=0;
	imut[14][11]=0;
	imut[14][12]=-3;
	imut[14][13]=0;
	imut[14][14]=6;
	imut[14][15]=-4;
	imut[14][16]=-2;
	imut[14][17]=-2;
	imut[14][18]=1;
	imut[14][19]=3;
	imut[14][20]=-1;
	imut[14][21]=-3;
	imut[14][22]=-3;
	imut[14][23]=-1;
	imut[15][1]=-1;
	imut[15][2]=-2;
	imut[15][3]=-2;
	imut[15][4]=-1;
	imut[15][5]=-3;
	imut[15][6]=-1;
	imut[15][7]=-1;
	imut[15][8]=-2;
	imut[15][9]=-2;
	imut[15][10]=-3;
	imut[15][11]=-3;
	imut[15][12]=-1;
	imut[15][13]=-2;
	imut[15][14]=-4;
	imut[15][15]=7;
	imut[15][16]=-1;
	imut[15][17]=-1;
	imut[15][18]=-4;
	imut[15][19]=-3;
	imut[15][20]=-2;
	imut[15][21]=-2;
	imut[15][22]=-1;
	imut[15][23]=-2;
	imut[16][1]=1;
	imut[16][2]=-1;
	imut[16][3]=1;
	imut[16][4]=0;
	imut[16][5]=-1;
	imut[16][6]=0;
	imut[16][7]=0;
	imut[16][8]=0;
	imut[16][9]=-1;
	imut[16][10]=-2;
	imut[16][11]=-2;
	imut[16][12]=0;
	imut[16][13]=-1;
	imut[16][14]=-2;
	imut[16][15]=-1;
	imut[16][16]=4;
	imut[16][17]=1;
	imut[16][18]=-3;
	imut[16][19]=-2;
	imut[16][20]=-2;
	imut[16][21]=0;
	imut[16][22]=0;
	imut[16][23]=0;
	imut[17][1]=0;
	imut[17][2]=-1;
	imut[17][3]=0;
	imut[17][4]=-1;
	imut[17][5]=-1;
	imut[17][6]=-1;
	imut[17][7]=-1;
	imut[17][8]=-2;
	imut[17][9]=-2;
	imut[17][10]=-1;
	imut[17][11]=-1;
	imut[17][12]=-1;
	imut[17][13]=-1;
	imut[17][14]=-2;
	imut[17][15]=-1;
	imut[17][16]=1;
	imut[17][17]=5;
	imut[17][18]=-2;
	imut[17][19]=-2;
	imut[17][20]=0;
	imut[17][21]=-1;
	imut[17][22]=-1;
	imut[17][23]=0;
	imut[18][1]=-3;
	imut[18][2]=-3;
	imut[18][3]=-4;
	imut[18][4]=-4;
	imut[18][5]=-2;
	imut[18][6]=-2;
	imut[18][7]=-3;
	imut[18][8]=-2;
	imut[18][9]=-2;
	imut[18][10]=-3;
	imut[18][11]=-2;
	imut[18][12]=-3;
	imut[18][13]=-1;
	imut[18][14]=1;
	imut[18][15]=-4;
	imut[18][16]=-3;
	imut[18][17]=-2;
	imut[18][18]=11;
	imut[18][19]=2;
	imut[18][20]=-3;
	imut[18][21]=-4;
	imut[18][22]=-3;
	imut[18][23]=-2;
	imut[19][1]=-2;
	imut[19][2]=-2;
	imut[19][3]=-2;
	imut[19][4]=-3;
	imut[19][5]=-2;
	imut[19][6]=-1;
	imut[19][7]=-2;
	imut[19][8]=-3;
	imut[19][9]=2;
	imut[19][10]=-1;
	imut[19][11]=-1;
	imut[19][12]=-2;
	imut[19][13]=-1;
	imut[19][14]=3;
	imut[19][15]=-3;
	imut[19][16]=-2;
	imut[19][17]=-2;
	imut[19][18]=2;
	imut[19][19]=7;
	imut[19][20]=-1;
	imut[19][21]=-3;
	imut[19][22]=-2;
	imut[19][23]=-1;
	imut[20][1]=0;
	imut[20][2]=-3;
	imut[20][3]=-3;
	imut[20][4]=-3;
	imut[20][5]=-1;
	imut[20][6]=-2;
	imut[20][7]=-2;
	imut[20][8]=-3;
	imut[20][9]=-3;
	imut[20][10]=3;
	imut[20][11]=1;
	imut[20][12]=-2;
	imut[20][13]=1;
	imut[20][14]=-1;
	imut[20][15]=-2;
	imut[20][16]=-2;
	imut[20][17]=0;
	imut[20][18]=-3;
	imut[20][19]=-1;
	imut[20][20]=4;
	imut[20][21]=-3;
	imut[20][22]=-2;
	imut[20][23]=-1;
	imut[21][1]=-2;
	imut[21][2]=-1;
	imut[21][3]=3;
	imut[21][4]=4;
	imut[21][5]=-3;
	imut[21][6]=0;
	imut[21][7]=1;
	imut[21][8]=-1;
	imut[21][9]=0;
	imut[21][10]=-3;
	imut[21][11]=-4;
	imut[21][12]=0;
	imut[21][13]=-3;
	imut[21][14]=-3;
	imut[21][15]=-2;
	imut[21][16]=0;
	imut[21][17]=-1;
	imut[21][18]=-4;
	imut[21][19]=-3;
	imut[21][20]=-3;
	imut[21][21]=4;
	imut[21][22]=1;
	imut[21][23]=-1;
	imut[22][1]=-1;
	imut[22][2]=0;
	imut[22][3]=0;
	imut[22][4]=1;
	imut[22][5]=-3;
	imut[22][6]=3;
	imut[22][7]=4;
	imut[22][8]=-2;
	imut[22][9]=0;
	imut[22][10]=-3;
	imut[22][11]=-3;
	imut[22][12]=1;
	imut[22][13]=-1;
	imut[22][14]=-3;
	imut[22][15]=-1;
	imut[22][16]=0;
	imut[22][17]=-1;
	imut[22][18]=-3;
	imut[22][19]=-2;
	imut[22][20]=-2;
	imut[22][21]=1;
	imut[22][22]=4;
	imut[22][23]=-1;
	imut[23][1]=0;
	imut[23][2]=-1;
	imut[23][3]=-1;
	imut[23][4]=-1;
	imut[23][5]=-2;
	imut[23][6]=-1;
	imut[23][7]=-1;
	imut[23][8]=-1;
	imut[23][9]=-1;
	imut[23][10]=-1;
	imut[23][11]=-1;
	imut[23][12]=-1;
	imut[23][13]=-1;
	imut[23][14]=-1;
	imut[23][15]=-2;
	imut[23][16]=0;
	imut[23][17]=0;
	imut[23][18]=-2;
	imut[23][19]=-1;
	imut[23][20]=-1;
	imut[23][21]=-1;
	imut[23][22]=-1;
	imut[23][23]=-1;
    }
}

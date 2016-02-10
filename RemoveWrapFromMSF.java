import java.io.*;
import java.util.*;

public class RemoveWrapFromMSF{

    private ArrayList<StringBuffer> alleles;
    private ArrayList<String> names;
    private StringBuffer header;
    private StringBuffer dnaCoordinate;
    private StringBuffer AACoordinate;
    private StringBuffer tick;
    
    public RemoveWrapFromMSF(){
	this.alleles = new ArrayList<StringBuffer>();
	this.names = new ArrayList<String>();
	this.header = new StringBuffer();
	this.dnaCoordinate = null;
	this.AACoordinate = null;
	this.tick = null;
    }
    
    public static void main(String[] args){
	if(args.length == 3)
	    new RemoveWrapFromMSF().run(args[0], args[1], args[2]);
	else
	    System.err.println("USAGE: java RemoveWrapFromMSF <HLA geneSymbol ex) A, B, C, DRB1, ... > <msf file> <output:merged file name>");
    }
    public void run(String geneSymbol, String msf, String merged){
	this.process(msf, geneSymbol);
	this.outToFile(merged);
    }


    public void outToFile(String merged){
	BufferedWriter bw = null;
	try{
	    bw = new BufferedWriter(new FileWriter(merged));
	    bw.write(this.header.toString());
	    bw.write(this.dnaCoordinate.toString() + "\n");
	    if(this.AACoordinate != null)
		bw.write(this.AACoordinate.toString() + "\n");
	    bw.write(this.tick.toString() + "\n");
	    for(StringBuffer sb : this.alleles)
		bw.write(sb.toString() + "\n");
	    bw.close();	    
	}catch(IOException ioe){
	    ioe.printStackTrace();
	}
    }
    
    private int getStartIndex(String dnaLine){
	String stripped = dnaLine.trim();
	String[] tokens = stripped.split("\\s+");
	return dnaLine.indexOf(tokens[1]);
    }

    //geneSymbol : A B C DRB1 ...
    public void process(String msf, String geneSymbol){
	BufferedReader br = null;
	String curline = null;
	
	try{
	    br = new BufferedReader(new FileReader(msf));
	    
	    boolean inMSF = false;//gets turns on after header once we hit MSA related coordinates
	    boolean inAlleles = false; //gets turned on after tick and turns off once sees blank(separator for block)
	    int alleleIndex = 0; //counter to correctly update alleles
	    int startPos = 0; // this stores the position of alignment block starting point.
	    while((curline=br.readLine())!=null){
		String stripped = curline.trim();
		if(!inMSF){
		    //when not inMSF this must be the beginning of the MSF
		    if(stripped.startsWith("cDNA") || stripped.startsWith("gDNA")){
			inMSF = true;//turn it on
			this.dnaCoordinate = new StringBuffer(curline);
			startPos = this.getStartIndex(curline);
		    }else//must be headers
			header.append(curline + "\n");
		}else{//if inMSF
		    if(inAlleles){//if alleles we update alleles by appending
			if(stripped.equals("")){//this is the end of MSA block 
			    inAlleles = false; //turn off now
			    alleleIndex = 0; //reset counter
			}else{
			    if(this.alleles.size() == alleleIndex){//this must be the first block
				this.alleles.add(new StringBuffer(curline)); // so we add 
				this.names.add(stripped.substring(0,stripped.indexOf(" ")));
			    }else{//otherwise
				if(this.names.get(alleleIndex).equals(stripped.substring(0,stripped.indexOf(" "))))
				    this.alleles.get(alleleIndex).append(this.stripHeader(curline, startPos));//not first block so we append the sequence portion
				else{
				    System.err.println("****ALLELE NAME NOT MATCHING: expected (" + this.names.get(alleleIndex) + ")\tFound (" + stripped.substring(0,stripped.indexOf(" ")));
				    System.exit(-1);
				}
			    }
			    if(alleleIndex == 0){
				this.appendNPadding(this.alleles.get(alleleIndex).length());
			    }
			    alleleIndex++;
			}
		    }
		    else if(stripped.startsWith("cDNA") || stripped.startsWith("gDNA")){
			this.dnaCoordinate.append(this.stripHeader(curline, startPos));
			startPos = this.getStartIndex(curline);
		    }else if(stripped.startsWith("AA codon")){
			if(this.AACoordinate == null)
			    this.AACoordinate = new StringBuffer(curline);
			else
			    this.AACoordinate.append(this.stripHeader(curline, startPos));
		    }else if(stripped.startsWith("|")){
			if(this.tick == null)
			    this.tick = new StringBuffer(curline);
			else
			    this.tick.append(stripped);
			inAlleles = true;
		    }
		}
	    }
	    br.close();
	}catch(IOException ioe){
	    ioe.printStackTrace();
	}
    }

    private void appendNPadding(int alleleLength){
	this.appendNPadding(this.dnaCoordinate, alleleLength - this.dnaCoordinate.length());
	this.appendNPadding(this.AACoordinate, alleleLength - this.AACoordinate.length());
	this.appendNPadding(this.tick, alleleLength - this.tick.length());
    }
    
    private void appendNPadding(StringBuffer sb, int n){
	for(int i =0; i<n; i++)
	    sb.append(" ");
	
    }

    private String stripHeader(String curline, int startPos){
	//return strippedline.split("\\s+")[1];
	return curline.substring(startPos);
    }
    

}

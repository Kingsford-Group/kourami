
import java.io.*;
import java.util.*;

public class SwapAlleles{

    private ArrayList<Allele> alleles;
    private StringBuffer header;
    private Allele ref;
    private Allele newref;
    
    public static int alleleSpaceLength = -1;
    
    public SwapAlleles(){
	this.alleles = new ArrayList<Allele>();
	this.header = new StringBuffer();
	this.ref = null;
	this.newref = null;
    }

    public static void main(String[] args){
	if(args.length == 3)
	    new SwapAlleles().run(args[0], args[1], args[2]);
	else
	    System.err.println("java SwapAlleles <msf> <currentRefAllele> <newRefAllele>");
    }

    public void run(String msf, String crn, String nrn){
	this.loadAlleles(msf, crn, nrn);
	this.swap();
	this.out(msf + ".swapped");
    }
    
    public void loadAlleles(String msf, String crName, String nrName){
	BufferedReader br = null;
	try{
	    br = new BufferedReader(new FileReader(msf));
	    String curline = null;
	    boolean inMSF = false;
	    System.err.println(crName);
	    System.err.println(nrName);
	    while((curline=br.readLine())!=null){
		if(inMSF){
		    Allele tmp = new Allele(curline);
		    
		    if(tmp.getName().equals(crName)){
			this.ref = tmp;
		    }else if(tmp.getName().equals(nrName)){
			this.newref = tmp;
		    }else
			this.alleles.add(tmp);
		}else{
		    header.append(curline + "\n");
		    if(curline.trim().startsWith("|"))
			inMSF = true;
		}
	    }
	    br.close();
	}catch(IOException ioe){
	    ioe.printStackTrace();
	}
    }

    public void swap(){
	for(int i=0; i<this.alleles.size(); i++)
	    this.alleles.get(i).update(this.ref, this.newref);
	Allele.swapRef(this.ref, this.newref);
    }

    public void out(String outfile){
	BufferedWriter bw = null;
	try{
	    bw = new BufferedWriter(new FileWriter(outfile));
	    bw.write(header.toString());
	    bw.write(this.newref.toString() + "\n");
	    bw.write(this.ref.toString() + "\n");
	    for(int i=0; i<this.alleles.size(); i++)
		bw.write(this.alleles.get(i).toString() + "\n");
	    
	    bw.close();
	}catch(IOException ioe){
	    ioe.printStackTrace();
	}
    }
}

class Allele{

    private String name;
    private StringBuffer sequence;
    
    public Allele(String curline){
	int sp = curline.indexOf(curline.trim().split("\\s+")[1]);
	if(SwapAlleles.alleleSpaceLength == -1)
	    SwapAlleles.alleleSpaceLength = sp;
	this.name = curline.substring(0, sp);
	this.sequence = new StringBuffer(curline.substring(sp).trim());
	//System.err.println(curline);
	//System.err.println(this.name + ":");
	//System.err.println(this.sequence);
    }

    public String toString(){
	return this.name + this.sequence.toString();
    }

    public String getName(){
	return this.name.trim();
    }
    
    public StringBuffer getSequence(){
	return this.sequence;
    }
    
    public void update(Allele curref, Allele newref){
    
	for(int i=0; i<sequence.length(); i++){
	    char crc = curref.charAt(i);
	    char nrc = newref.charAt(i);

	    //if crc and nrc are same base or both gaps
	    if(crc == nrc || nrc == '-')
		;//we dont' do anything.
	    else if(crc != nrc){
		char c = this.charAt(i);
		
		if(c == '.' || c == '*')
		    ;//we don't do anything
		else if(c == '-')//base is same as crc but diff from nrc, so need to actually write crc base
		    this.setCharAt(i, crc);
		else{
		    if(c == nrc)//base is different from crc but same as nrc, so need to set it '-' -->meaning same as nrc
			this.setCharAt(i, '-');
		    else//base is different from crc and also different from nrc, so keep it as it is.
			;
		}
	    }
	}
    }

    public static void swapRef(Allele curref, Allele newref){
	for(int i=0; i<curref.getSequence().length(); i++){
	    char crc = curref.charAt(i);
	    char nrc = newref.charAt(i);
	    
	    if(crc == nrc){ // gap:gap 
		;
	    }else if(nrc == '-'){
		newref.setCharAt(i, crc);
		curref.setCharAt(i, nrc);
	    }else if(nrc == '*'){ 
		/* THIS fills * for known bases but if we decide to not use UNKNOWN bases, we should put gaps */
		newref.setCharAt(i, crc);
		//curref.setCharAt(i, '-');
	    }
	}
    }

    public void setCharAt(int n, char c){
	this.sequence.setCharAt(n, c);
    }
    
    public char charAt(int n){
	return this.sequence.charAt(n);
    }

    public boolean sameCharAt(int n, Allele other){
	if(this.charAt(n) == other.charAt(n))
	    return true;
	return false;
    }
}

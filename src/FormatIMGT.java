import java.util.*;
import java.io.*;

public class FormatIMGT{
    
    public static void processGene(String imgtpath, String outpath, String geneName){
	String genfile = imgtpath + File.separator + geneName + "_gen.txt";
	String nucfile = imgtpath + File.separator + geneName + "_nuc.txt";
	String genoutfile = outpath + File.separator + geneName + "_gen.txt";
	String nucoutfile = outpath + File.separator + geneName + "_nuc.txt";
	if(geneName.startsWith("DRB"))
	    nucfile = imgtpath + File.separator + "DRB_nuc.txt";
	
	IMGTReformatter nuc = null;
	IMGTReformatter gen = null;
	if(new File(nucfile).exists())
	    nuc = new IMGTReformatter(nucfile, geneName);
	else
	    System.err.println("[WARNING]: Missing < " + nucfile + " >. Proceeding without.\nHowever, it is STRONGLY recommended to run FormatIMGT\nwith both nuc and gen files for each gene.");
	
	if(new File(genfile).exists())
	   gen = new IMGTReformatter(genfile, geneName);
	else{
	    System.err.println("[ERROR]: Missing < " + genfile + " >.\n");
	    System.err.println("A gen file is required for each gene!");
	    System.exit(1);
	}
	
	String nucRefAl = null;
	if(nuc != null)
	    nucRefAl = nuc.getRefAlleleName();
	String genRefAl = gen.getRefAlleleName();
	
	//if both ref alleles are same
	if(nuc == null || nucRefAl.equals(genRefAl)){
	    if(nuc != null)
		System.err.println("\trefGeneName on nuc and gen are same");
	    System.err.println("\tWrting to :\t" + genoutfile);
	    gen.outToFile(genoutfile);
	    if(nuc !=null){
		System.err.println("\tWrting to :\t" + nucoutfile);
		nuc.outToFile(nucoutfile);
	    }
	}
	//if NOT same, then genRefAl must be in nucRefAl
	//because nucRef is the superset.
	else{
	    System.err.println("\trefGeneName on nuc and gen are NOT same");
	    
	    if(nuc.exists(genRefAl)){
		System.err.println("\tSWAPPING refGenes on nuc and gen");
		nuc.swap(genRefAl);
		System.err.println("\tWrting to :\t" + genoutfile);
		gen.outToFile(genoutfile);
		System.err.println("\tWrting to :\t" + nucoutfile);
		nuc.outToFile(nucoutfile);
	    }else{
		System.err.println("Reference sequence entry [" + genRefAl + "] is " 
				   + "NOT found in nuc alignments.\n"
				   + "Check the alignment files.");
		System.exit(1);
	    }
	    ;//swap then outToFile();
	}
	
	MergeMSFs mm = new MergeMSFs();
	if(!mm.merge(nucoutfile, genoutfile, true)){
	    System.err.println("ERROR in MSA merging. CANNOT proceed further. Exiting..");
	    System.exit(-1);
	}	    
    }
    
    public static void main(String[] args){
	if(args.length == 0){
	    System.err.println("USAGE: java FormatIMGT <path-to-IMGT-alignment-directory>");
	    System.exit(1);
	}
	String imgtpath = args[0];
	String outpath = args[0] + File.separator + "kouramiFormatted";
	File outdir = null;
	try{
	    outdir = new File(outpath);
	    if(outdir.exists()){
		if(!outdir.isDirectory()){
		    System.err.println(outpath + " exists and it is not a writable directory.");
		    System.exit(1);
		}
	    }else
		outdir.mkdirs();
	    
	    boolean missingFiles = false;
	    for(int i=0; i<FormatIMGT.expList.length; i++){
		File genFile = new File(imgtpath + File.separator + FormatIMGT.expList[i] + "_gen.txt");
		File nucFile = new File(imgtpath + File.separator + FormatIMGT.expList[i] + "_nuc.txt");
		if(FormatIMGT.expList[i].startsWith("DRB1"))
		    nucFile = new File(imgtpath + File.separator + "DRB_nuc.txt");
		if(!genFile.exists()){
		    System.err.println("Missing :\t" + genFile.getAbsolutePath());
		    System.err.println("A gen file is required for each gene.");
		    missingFiles = true;
		}
		if(!nucFile.exists()){
		    System.err.println("Missing :\t" + nucFile.getAbsolutePath());
		    //missingFiles = true;
		}
	    }
	    if(missingFiles)
		System.exit(1);

	    for(int i=0; i<FormatIMGT.expList.length; i++){
		String geneName = FormatIMGT.expList[i];
		System.err.println(">>>>>>>>>  Processing\t[" + geneName + "]  <<<<<<<<<<");
		FormatIMGT.processGene(imgtpath, outpath, geneName);
	    }
	}catch(Exception e){
	    e.printStackTrace();
	    System.exit(1);
	}
	
    }
    static final String[] list = {"A" , "B" , "C" , "DQA1" , "DQB1" , "DRB1"};
    
    /* check P */
    static final String[] expList = {"A", "B", "C", "DMA", "DMB", "DOA", "DOB"
				     , "DPA1", "DPB1", "DPB2", "DQA1", "DQB1", "DRA"
				     , "DRB1", "DRB2", "DRB3", "DRB4", "DRB5", "DRB6"
				     , "DRB7", "DRB8", "DRB9", "E", "F", "G"
				     , "H", "HFE", "J", "K", "L", "MICA"
				     , "MICB", "TAP1", "TAP2", "V", "Y"
				     , "P"}; // check P 
}

class IMGTReformatter{
    
    private ArrayList<Allele> alleles;
    
    private StringBuffer header;
    private StringBuffer dnaCoordinate;
    private StringBuffer AACoordinate;
    private StringBuffer tick;

    private int startPos;
    
    private Allele newRef;
    private int newRefIndex;

    private String geneName;

    public IMGTReformatter(){
	this.alleles = new ArrayList<Allele>();
	this.header = new StringBuffer();
	this.dnaCoordinate = null;
	this.AACoordinate = null;
	this.tick = null;
	
	this.startPos = -1;
	
	this.newRef = null;
	this.newRefIndex = 0;
	
	this.geneName = null;
    }
    
    public IMGTReformatter(String msf, String gn){
	this();
	this.geneName = gn;
	this.loadAlleles(msf);
    }

    public int getNewRefIndex(){
	return this.newRefIndex;
    }

    public boolean exists(String newRefName){
	for(int i=1; i < this.alleles.size(); i++){
	    if(this.alleles.get(i).getName().equals(newRefName)){
		this.newRef = this.alleles.get(i);
		this.newRefIndex = i;
		return true;
	    }
	}
	return false;
    }

    //returns the index for newRef allele
    public void swap(String newRefName){
	Allele curRef = this.alleles.get(0);
	for(int i=1; i < this.alleles.size(); i++){
	    if(i!=this.newRefIndex)
		this.alleles.get(i).update(curRef, this.newRef);
	}
	Allele.swapRef(curRef, newRef);
    }

    private boolean isAlleleLine(String[] tokens){
	if(tokens.length > 0){
	    return tokens[0].startsWith(this.geneName + "*") && tokens[0].matches("[A-Z]+\\d{0,1}\\*\\d+(:\\d+){0,3}[A-Z]*");
	}
	System.err.print("[" + this.geneName + "]: " + tokens[0] + "\t");
	System.err.print(tokens[0].startsWith(this.geneName + "*") + "\t");
	System.err.println(tokens[0].matches("[A-Z]+\\d{0,1}\\*\\d+(:\\d+){0,3}[A-Z]*"));
	return false;
    }
    
    /* 
     * Parses alignment format from IMGT db 
     */
    public void loadAlleles(String msf){
	BufferedReader br = null;
	String curline = null;
	try{
	    br = new BufferedReader(new FileReader(msf));
	    boolean inMSF = false; //flag to split header and post-header
	    int alleleIndex = 0; 
	    while((curline=br.readLine())!=null){
		String stripped = curline.trim();
		String[] tokens = stripped.split("\\s+");
	
		//parse headers
		if(!inMSF){
		    if(stripped.startsWith("cDNA") || stripped.startsWith("gDNA")){
			inMSF = true;//turn it on
			this.dnaCoordinate = new StringBuffer(curline);
			this.startPos = this.getStartIndex(curline, tokens);
		    }else//must be headers
			header.append(curline + "\n");
		}else{
		    if(stripped.startsWith("|")){
			if(this.tick == null)
			    this.tick = new StringBuffer(curline);
			else
			    this.tick.append(stripped);
		    }
		    // if it's cDNA gDNA tag
		    else if(stripped.startsWith("cDNA") || stripped.startsWith("gDNA")){
			this.dnaCoordinate.append(this.stripHeader(curline, this.startPos));
			alleleIndex = 0; // reset alleleIndex
		    }
		    else if(stripped.startsWith("AA codon")){
			if(this.AACoordinate == null)
			    this.AACoordinate = new StringBuffer(curline);
			else
			    this.AACoordinate.append(this.stripHeader(curline, this.startPos));
		    }
		    else if(isAlleleLine(tokens)){
			//System.err.println("allelesSize:\t" + this.alleles.size() + "\talleleIndex:\t" + alleleIndex );
			if(this.alleles.size() == alleleIndex){
			    //System.err.println("ADDING ");
			    this.alleles.add(new Allele(curline.substring(0 , this.startPos), curline.substring(this.startPos)));
			}else
			    this.alleles.get(alleleIndex).appendSequence(tokens[0], curline.substring(this.startPos));
			//append correct number of white spaces for coordinate and tick lines
			if(alleleIndex == 0)
			    this.appendNPadding(this.alleles.get(alleleIndex).curStringLength(this.startPos));
			alleleIndex++;
		    }
		}
	    }
	}catch(IOException ioe){
	    ioe.printStackTrace();
	}
    }
    
    //if newRefIndex == 0 --> no swapping was necessary
    //if newRefIndex > 0  --> oldRefIndex is 0, and newRefIndex is given
    public void outToFile(String merged){
	BufferedWriter bw = null;
	try{
	    bw = new BufferedWriter(new FileWriter(merged));
	    bw.write("Nuc+Gen merged MSA for Kourami\n");
	    bw.write(this.header.toString());
	    bw.write(this.dnaCoordinate.toString().replaceFirst("\\s++$", "") + "\n");
	    if(this.AACoordinate != null)
		bw.write(this.AACoordinate.toString().replaceFirst("\\s++$", "") + "\n");
	    bw.write(this.tick.toString().replaceFirst("\\s++$", "") + "\n");
	    
	    //if no swapping was applied
	    if(this.newRefIndex == 0){
		for(Allele a : this.alleles)
		    bw.write(a.toString() + "\n");
	    }
	    //if there was a swapping
	    else if(this.newRefIndex > 0){
		//write newRefAllele first
		bw.write(this.alleles.get(this.newRefIndex).toString() + "\n");
		for(int i=0; i<this.alleles.size(); i++){
		    if(i != this.newRefIndex )
			bw.write(this.alleles.get(i).toString() + "\n");
		}
	    }
	    
	    bw.close();	    
	}catch(IOException ioe){
	    ioe.printStackTrace();
	}
    }
    
    public String getRefAlleleName(){
	return this.alleles.get(0).getName();
    }

    private void appendNPadding(int alleleLength){
	if(this.dnaCoordinate != null)
	    this.appendNPadding(this.dnaCoordinate, alleleLength - this.dnaCoordinate.length());
	if(this.AACoordinate != null)
	    this.appendNPadding(this.AACoordinate, alleleLength - this.AACoordinate.length());
	if(this.tick !=null)
	    this.appendNPadding(this.tick, alleleLength - this.tick.length());
    }
    
    private void appendNPadding(StringBuffer sb, int n){
	for(int i =0; i<n; i++)
	    sb.append(" ");
    }
    
    private String stripHeader(String curline, int pos){
	return curline.substring(pos);
    }

    private int getStartIndex(String dnaLine, String[] tokens){
	return dnaLine.indexOf(tokens[1]);
    }
        
}

class Allele{

    private String nameWS;
    private String name;
    private StringBuffer sequence;
    
    //n may have leading and trailing whitespaces
    public Allele(String n, String seq){
	this.nameWS = n;
	this.name = n.trim();
	this.sequence = new StringBuffer(seq);
    }

    public boolean appendSequence(String n, String seq){
	if(this.name.equals(n)){
	    this.sequence.append(seq);
	    return true;
	}
	return false;
    }

    public int curStringLength(int startPos){
	return startPos + this.sequence.length();
    }
    
    public String toString(int startPos){
	StringBuffer bf = new StringBuffer(name);
	int numSpaces = startPos - name.length();
	for(int i = 0; i < numSpaces; i++)
	    bf.append(" ");
	bf.append(sequence);
	return bf.toString();
    }

    public String toString(){
	return this.nameWS + this.sequence;
    }
    
    public String getName(){
	return this.name;
    }

    public StringBuffer getSequenceBuffer(){
	return this.sequence;
    }

    //symbols
    // . --> gap
    // - --> same as curRef
    // * --> unknownbase
    // [aAcCgGtT] --> base
    public void update(Allele curref, Allele newref){
    
	for(int i=0; i<this.sequence.length(); i++){
	    char crc = curref.charAt(i);
	    char nrc = newref.charAt(i);
	    
	    //if crc and nrc are same
	    if(crc == nrc || nrc == '-')//both gaps, both unknown, or same bases.
		;//not do anything
	    else if(crc != nrc){
		// [crc] : [nrc]
		//  base :  */./dBase
		//   *   :  ./base
		//   .   :  */base
		char c = this.charAt(i);
		//if it's a gap '.' or unknown base '*'
		if( c == '.' || c == '*')
		    ;//not do anything
		else if(c == '-')//base is same as crc but diff from nrc, so need to actually write crc base.
		    this.setCharAt(i, crc);
		else{ // c is a letter base
		    if(c == nrc)//base is differenct from crc but same as nrc, so need to set it '-' --> meaing same as nrc
			this.setCharAt(i, '-');
		    else//base is different from crc and also difeferent from nrc, so keep it as it is
			;
		}
	    }
	}
    }

    public static void swapRef(Allele curref, Allele newref){
	for(int i=0; i<curref.getSequenceBuffer().length(); i++){
	    char crc = curref.charAt(i);
	    char nrc = newref.charAt(i);
	    
	    if(crc == nrc){ // gap:gap unknown:unknonwn  .:. *:*
		;
	    }else if(nrc == '-'){ // crc must be a base
		newref.setCharAt(i, crc);
		curref.setCharAt(i, nrc);
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


/*
    public void loadAlleles(String msf){
	BufferedReader br = null;
	String curline = null;
	try{
	    br = new BufferedReader(new FileReader(msf));
	    boolean inMSF = false;
	    boolnea inAllele = false;
	    int alleleIndex = 0; //index for allele
	    int startPos = 0;
	    while((curline=br.readLine())!=null){
		String stripped = curline.trim();
		if(!inMSF){
		    if(stripped.startsWith("cDNA") || stripped.startsWith("gDNA")){
			inMSF = true;//turn it on
			this.dnaCoordinate = new StringBuffer(curline);
			startPos = this.getStartIndex(curline);
		    }else//must be headers
			header.append(curline + "\n");
		}else{
		    if(inAlleles){
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
				    System.err.println(curline);
				    System.err.println("While processing:\t" + msf);
				    System.exit(-1);
				}
			    }
			    if(alleleIndex == 0){
				this.appendNPadding(this.alleles.get(alleleIndex).length());
			    }
			    alleleIndex++;
			}
		    }
		}
	    }
	}catch(IOException ioe){
	    ioe.printStackTrace();
	}
    }
*/

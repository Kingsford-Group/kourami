import java.io.*;
import java.util.*;

public class MergeMSFs{

    private String referenceAllele;
    private Hashtable<String, Sequence> allele2Sequence;
    private ArrayList<String> orderedAlleles;
    private StringBuffer header;
    private String geneName;
    
    public MergeMSFs(){
	this.header = new StringBuffer();
	this.allele2Sequence = new Hashtable<String, Sequence>();
	this.orderedAlleles = new ArrayList<String>();
	this.geneName = null;
    }

    public ArrayList<Sequence> getListOfSequences(){
	ArrayList<Sequence> l = new ArrayList<Sequence>();
	for(int i=0; i<this.orderedAlleles.size(); i++){
	    l.add(this.allele2Sequence.get(this.orderedAlleles.get(i)));
	}
	return l;
    }

    public static void main(String[] args){
	if(args.length == 2)
	    new MergeMSFs().merge(args[0], args[1]);
	else
	    System.err.println("USAGE java MergeMSFs <nuc file> <gen file>");
    }

    public void outToFasta(){
	BufferedWriter bw = null;
	try{
	    bw = new BufferedWriter(new FileWriter(this.geneName + ".merfed.fa"));
	    for(int i=0; i< this.orderedAlleles.size(); i++){
		bw.write(this.allele2Sequence.get(this.orderedAlleles.get(i)).toFastaString());
	    }
	    bw.close();
	}catch(IOException ioe){
	    ioe.printStackTrace();
	}
    }
    
    public void print(){
	for(int i=0; i<this.orderedAlleles.size(); i++){
	    String curA = this.orderedAlleles.get(i);
	    Sequence curS = this.allele2Sequence.get(curA);
	    //String colSeq = MergeMSFs.removeDots(curS.getColumnSequence());
	    //String fullSeq = curS.getFullSequence();
	    //String seqFromBases = MergeMSFs.removeDots(curS.getSequenceFromBases());
	    System.out.println(curA + "\t" + curS.getColumnSequence());
	    //System.out.print("colSeq & fullSeq:\t");
	    /*
	    if(colSeq.equals(fullSeq)){
		System.out.println("SAME");
		if(colSeq.toLowerCase().equals(seqFromBases.toLowerCase()))
		    System.out.println("\tcolSeq & seqFromBases: SAME");
		else{
		    System.out.println("\tcolSeq & seqFromBases: DIFFERENT");
		    System.out.println("colSeq:\t" + colSeq + "\nSeqFromB:\t" + seqFromBases);
		    System.out.println("colSeq:\t" + curS.getColumnSequence() + "\nfulSeq:\t" + fullSeq + "\nSeqFromB:\t" + curS.getSequenceFromBases());
		}
	    }
	    else{
		System.out.println("DIFFERENT");
		System.out.println("colSeq:\t" + colSeq + "\nfulSeq:\t" + fullSeq);
		System.out.println("colSeq:\t" + curS.getColumnSequence() + "\nfulSeq:\t" + fullSeq);
	    }
	    */	
	    
	}
    }
    
    public static String removeDots(String s){
	StringBuffer bf = new StringBuffer();
	for(int i=0; i<s.length(); i++){
	    char cc = s.charAt(i);
	    if(cc != '.')
		bf.append(cc);
	}
	return bf.toString();
    }

    /*
     * nucF : nuc file containing MSA of coding sequences
     * genF : gen file containing MSA of whole gene sequences
     *
     * 1) consturct a reference sequenced by merging nuc(exons) and gen(introns)
     * 2) add nucSequences -
     * 3) update intron --> 
     * 
     */
    public void merge(String nucF, String genF){
	BufferedReader nucbr = null;
	BufferedReader genbr = null;
	HashMap<String, String[]> genHash = new HashMap<String, String[]>();
	String nucline = null;
	String genline = null;
	int nucsp = 0;
	int gensp = 0;
	String mergedRef = null;
	try{
	    /* First take care of the reference (first sequence)*/
	    nucbr = new BufferedReader(new FileReader(nucF));
	    genbr = new BufferedReader(new FileReader(genF));
	    
	    nucline = getFirstLine(nucbr);
	    genline = getFirstLine(genbr);

	    //first position of sequences in the line
	    nucsp = nucline.indexOf(nucline.trim().split("\\s+")[1]);
	    gensp = genline.indexOf(genline.trim().split("\\s+")[1]);
	    //System.out.println("nucsp: " + nucsp);
	    //System.out.println("gensp: " + gensp);


	    if(nucline == null || genline == null){
		System.err.println("Something wrong with input files. Exiting...");
		System.exit(0);
	    }
	    
	    String nucname = nucline.substring(0, nucsp).trim();
	    String genname = genline.substring(0, gensp).trim();
	    this.referenceAllele = nucname;
	    this.geneName = referenceAllele.substring(0,referenceAllele.indexOf("*"));
	    if(!nucname.equals(genname)){
		System.err.println("REF SEQ names differs :");
		System.err.println("(nuc):" + nucname);
		System.err.println("(gen):" + genname);
		System.exit(-1);
	    }
	    
	    String nucsequence = nucline.substring(nucsp).trim();
	    String gensequence = genline.substring(gensp).trim();
	    //System.out.println("nucseq: " + nucsequence + "|");
	    //System.out.println("genseq: " + gensequence + "|");

	    String[] nucblocks = nucsequence.split("\\|");
	    String[] genblocks = gensequence.split("\\|");
	    
	    Sequence refSequence=this.mergeAndAdd(nucname, nucblocks, genblocks);
	    /* End of taking care of first sequences */
	    
	    
	    /* Now get a list of Gens and process together with Nucs if same sequence is contained in nuc */
	    while( (genline=genbr.readLine()) != null ){
		String allele = genline.substring(0,gensp).trim();
		genHash.put(genline.substring(0,gensp).trim(), genline.substring(gensp).trim().split("\\|"));
	    }
	    genbr.close();
	    
	    /* end of processing genbr */

	    /* Take care of nucF since genF is always a subset of nucF */
	    /* No need to process genF */
	    String curline = null;
	    while((curline=nucbr.readLine()) != null){
		genblocks = null;
		String allele = curline.substring(0,nucsp).trim();
		String curseq = curline.substring(nucsp).trim();
		String msfsequence = null;
		//if there is a matching gen sequence in the hash.
		if( (genblocks = genHash.get(allele)) != null){
		    //System.err.println("Processing:\t" + allele);
		    boolean replaceAbbrv = true;
		    this.mergeAndAdd(allele, curseq.split("\\|"), genblocks, refSequence);
		    //System.err.println("COLUMNSequence:\t" + refSequence.getColumnSequence());
		}
		else{ // else we just take the nucsequence
		    msfsequence = MergeMSFs.removeBlank(curseq);
		    //if we have seen this sequence before --> something is not right
		    if(this.allele2Sequence.get(allele) != null)
			System.err.println("DUPLICATE ENTRY? --> " + allele + " (Skipping for now)...");
		    else{
			//use nuc only constructor of Sequence to process nuc-only allele
			this.addAllele(allele, new Sequence(allele, msfsequence, refSequence));//this.allele2Sequence.get(nucname)));
		    }
		}
	    }
	    nucbr.close();
	}catch(IOException ioe){
	    ioe.printStackTrace();
	}

	//this.print();
	//this.outToFasta();
    }
    
    private Sequence mergeAndAdd(String alleleName, String[] nucblocks, String[] genblocks){
	return this.mergeAndAdd(alleleName, nucblocks, genblocks, null);
    }
	/* given nucblocks and genblocks, merge nucblocks(use all exons) and genblocks(use only introns) */
    private Sequence mergeAndAdd(String alleleName, String[] nucblocks, String[] genblocks, Sequence refSequence){
	String mergedSeq = this.mergeBlocks(nucblocks, genblocks);
	Sequence curSeq = null;
	if(refSequence == null)
	    curSeq = new Sequence(alleleName, MergeMSFs.removeBlank(mergedSeq), false, refSequence);
	else
	    curSeq = new Sequence(alleleName, MergeMSFs.removeBlank(mergedSeq), true, refSequence);

	this.addAllele(alleleName, curSeq);
	return curSeq;
	/*
	mergedRef = this.mergeBlocks(nucblocks, genblocks);
	//System.out.println("mergedRef: " + mergedRef + "|");
	mergedRef = MergeMSFs.removeBlank(mergedRef);
	//System.out.println("merRef(BL: " + mergedRef + "|");
	Sequence refSequence = new Sequence(nucname, mergedRef);//MergeMSFs.removeBlank(mergedRef));
	//refSequence.printBoundaries();
	this.addAllele(nucname, refSequence);
	*/
    }

    private void addAllele(String allele, Sequence seq){
	this.orderedAlleles.add(allele);
	this.allele2Sequence.put(allele, seq);
    }
    
    //used to merge nucblocks and genblocks for the reference sequence.
    //this only merges string. so any abbrv symbols are not replaced. ("-" and such)
    private String mergeBlocks(String[] nucblocks, String[] genblocks){
	if( (nucblocks.length * 2 + 1) != genblocks.length){
	    System.err.println("nucblocks.length : " + nucblocks.length + " genblocks.length :" + genblocks.length +"\ngenblocks length must be equal to [2 * (nucblocks length) + 1]");
	    System.exit(-1);
	}
	StringBuffer bf = new StringBuffer();
	for(int i=0;i<genblocks.length;i++){
	    //i = 1, 3, 5 --> exon
	    if( i%2 == 1)
		bf.append(" | " + nucblocks[i/2].trim() + " | ");
	    else //intron i = 0, 2, 4
		bf.append(genblocks[i].trim());
	}
	return bf.toString();
    }
    
    
    //this returns the first sequnece usually the reference in MSA (nuc or gen).
    public String getFirstLine(BufferedReader br){
	String curline = null;
	try{
	    while((curline=br.readLine())!=null){
		String tmp = curline.trim();
		//"[A-Z]+\\d\\*\\d+:\\d+(:\\d+){0,2}")
		if(tmp.split("\\s+")[0].matches("[A-Z]+\\d{0,1}\\*\\d+:\\d+(:\\d+){0,2}[A-Z]*")){
		    return curline;
		}
	    }
	}catch(IOException ioe){
	    ioe.printStackTrace();
	}
	return null;
    }

    public static String removeBlank(String sequence, boolean modifyHeader){
	return MergeMSFs.removeBlank(sequence);
    }

    /* this removes all blank embedded in sequences*/
    public static String removeBlank(String sequence){
	StringBuffer bf = new StringBuffer();
	//int headeri = 0;
	for(int i=0; i<sequence.length(); i++){
	    char tmp = sequence.charAt(i);
	    if(tmp == ' '|| tmp == '\t'){
		;
		//if(modifyHeader){
		//this.header.deleteCharAt(headeri);
		//  headeri--;
		//}
	    }else
		bf.append(tmp);
	    //headeri++;
	}
	return bf.toString();
    }
    
    /* this fills *(unknown) and -(same as ref) as bases*/
    public static String abbrv2Seq(String abbrv, String refSeq){
	//System.out.println("ABBRV:\t" + abbrv);
	//System.out.println("ABBRV:\t" + refSeq);
	StringBuffer bf = new StringBuffer();
	for(int i=0;i<abbrv.length();i++){
	    char curChar = abbrv.charAt(i);
	    if(curChar == '*')
		bf.append(Character.toLowerCase(refSeq.charAt(i)));
	    else if(curChar == '-')
		bf.append(refSeq.charAt(i));
	    else
		bf.append(curChar);
	}
	return bf.toString();
    }


    private String stripHeader(String curline, int startPos){
	return curline.substring(startPos);
    }
    
    /*
    public void readInNuc(String nuc){
	Sequence genSeq = this.allele2Sequence.get(this.referenceAllele);
	BufferedReader br = null;
	String curline = null;
	String refSeq = null;
	try{
	    br = new BufferedReader(new FileReader(nuc));
	    boolean inMSF = false;
	    boolean firstSeq = true;
	    int startPos = 0;
	    while((curline=br.readLine())!=null){
		String stripped = curline.trim();
		if(!inMSF){
		    if(stripped.startsWith("cDNA")){
			inMSF = true;
		    }else
			this.header.append(curline + "\n");
		}else{
		    if(stripped.startsWith("AA codon"))
			;
		    else if(stripped.startsWith("|"))
			startPos = curline.indexOf("|");
		    else{
			String allele = curline.substring(0,startPos).trim();
			if(this.allele2Sequence.get(allele) == null){
			    String msfsequence = curline.substring(startPos).trim();
			    if(firstSeq){
				firstSeq = false;
				msfsequence = MergeMSFs.removeBlank(msfsequence, true);
				refSeq = msfsequence;
			    }else
				msfsequence=MergeMSFs.abbrv2Seq(MergeMSFs.removeBlank(msfsequence,false), refSeq);
			    
			    this.orderedAlleles.add(allele);
			    this.allele2Sequence.put(allele, new Sequence(msfsequence, genSeq));//, true)); 
			}
		    }
		}
	    }
	}catch(IOException ioe){
	    ioe.printStackTrace();
	}
    }
    */
    /*
    public void readInGen(String gen){
	BufferedReader br = null;
	String curline = null;
	String refSeq = null;
	try{
	    br = new BufferedReader(new FileReader(gen));
	    boolean inMSF = false;
	    boolean firstSeq = true;
	    int startPos = 0;
	    while((curline=br.readLine())!=null){
		String stripped = curline.trim();
		if(!inMSF){
		    if(stripped.startsWith("gDNA")){
			inMSF = true;
		    }else
			this.header.append(curline + "\n");
		}else{
		    if(stripped.startsWith("AA codon"))
			;
		    else if(stripped.startsWith("|"))
			startPos = curline.indexOf("|");
		    else{
			String allele = curline.substring(0,startPos).trim();
			String msfsequence = curline.substring(startPos).trim();
			if(firstSeq){
			    firstSeq = false;
			    msfsequence = MergeMSFs.removeBlank(msfsequence, true);
			    refSeq = msfsequence;
			}else
			    msfsequence=MergeMSFs.abbrv2Seq(MergeMSFs.removeBlank(msfsequence,false), refSeq);
			
			this.orderedAlleles.add(allele);
			this.allele2Sequence.put(allele, new Sequence(allele, msfsequence, true)); 
		    }
		}
	    }
	    br.close();
	}catch(IOException ioe){
	    ioe.printStackTrace();
	}
    }
    */
    /*
    private void addNucAllele(String allele, String sequence, boolean abbrv){
	this.orderedAlleles.add(allele);
	boolean isNuc = true;
	this.allele2Sequence.put(allele, new Sequence(allele, sequence, abbrv, isNuc));
	}*/

    
}








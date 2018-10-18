/*
Part of Kourami HLA typer/assembler
(c) 2017 by  Heewook Lee, Carl Kingsford, and Carnegie Mellon University.
See LICENSE for licensing.
*/
import java.io.*;
import java.util.*;


public class MergeMSFs{

    private String referenceAllele;
    private Hashtable<String, Sequence> allele2Sequence;
    private ArrayList<String> orderedAlleles;
    private StringBuffer header;
    private String geneName;
    
    private boolean isDRBGene;
    private String drbGeneName;

    //
    // returns typing sequences for a HLAGene
    //
    public ArrayList<HLASequence> formDataBase(ArrayList<Group> groups){//String nomGFile){
	//NomG nomG = new NomG();
	//nomG.loadGroups(nomGFile);
	//ArrayList<Group> groups = nomG.getGroups();
	//HashMap<String, ArrayList<HLASequence>> hlaseqs = new HashMap<ArrayList<HLASequence>>();
	ArrayList<HLASequence> typingSeqs = new ArrayList<HLASequence>();
	for(Group g : groups){
	    //ArrayList<HLASequence> typingSeqs = hlaseqs.get(g.getHLAGeneName());
	    /*if(typingSeqs == null){
		typingSeqs = new ArrayList<HLASequence>();
		hlaseqs.put(g.getHLAGeneName(), typingSeqs);
		}*/
	    //System.err.println("Querying:" + g.getFirstAllele());
	    //System.err.println("Key Example:" + allele2Sequence.keySet().iterator().next());
	    Sequence s = this.allele2Sequence.get(g.getFirstAllele());
	    typingSeqs.add(new HLASequence(g, this.allele2Sequence.get(g.getFirstAllele())));
	}
	return typingSeqs;
    }

    public ArrayList<HLASequence> formDataBaseAll(){
	ArrayList<HLASequence> typingSeqs = new ArrayList<HLASequence>();
	Enumeration<String> keys = this.allele2Sequence.keys();
	while(keys.hasMoreElements()){
	    String curAllele = keys.nextElement();
	    //System.err.println("curAlle:\t" + curAllele);
	    Sequence s = this.allele2Sequence.get(curAllele);
	    typingSeqs.add(new HLASequence(new Group(curAllele), this.allele2Sequence.get(curAllele)));
	}
	return typingSeqs;
    }

    
    public MergeMSFs(){
	this.header = new StringBuffer();
	this.allele2Sequence = new Hashtable<String, Sequence>();
	this.orderedAlleles = new ArrayList<String>();
	this.geneName = null;
	this.isDRBGene = false;
	this.drbGeneName = null;
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
	else if(args.length == 3)
	    new MergeMSFs().merge(args[0], args[1], args[2]);
	else
	    System.err.println("USAGE java MergeMSFs <nuc file> <gen file> [DRBGeneName:optional]");
    }

    public void outToFasta(){
	this.outToFasta(HLA.OUTPREFIX + "_");
    }

    public void outToFasta(String outprefix){
	this.outToFasta(outprefix, true);
    }
    public void outToFasta(String outprefix, boolean writemsa){
	BufferedWriter bw = null;
	BufferedWriter bw2 = null;
	try{
	    if(this.isDRBGene){
		bw = new BufferedWriter(new FileWriter( outprefix + this.drbGeneName + ".merged.fa"));
		if(writemsa)
		    bw2 = new BufferedWriter(new FileWriter( outprefix + this.drbGeneName + ".msa.fa"));
	    }else{
		bw = new BufferedWriter(new FileWriter( outprefix + this.geneName + ".merged.fa"));
		if(writemsa)
		    bw2 = new BufferedWriter(new FileWriter( outprefix + this.geneName + ".msa.fa"));
		
	    }
	    for(int i=0; i< this.orderedAlleles.size(); i++){
		if(this.isDRBGene){
		    if(this.orderedAlleles.get(i).startsWith(this.drbGeneName)){
			bw.write(this.allele2Sequence.get(this.orderedAlleles.get(i)).toFastaString());
			if(writemsa)
			    bw2.write(this.allele2Sequence.get(this.orderedAlleles.get(i)).toFastaColumnSequence());
		    }
		}else{
		    bw.write(this.allele2Sequence.get(this.orderedAlleles.get(i)).toFastaString());
		    if(writemsa)
			bw2.write(this.allele2Sequence.get(this.orderedAlleles.get(i)).toFastaColumnSequence());
		}
	    }
	    bw.close();
	    if(writemsa)
		bw2.close();
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
    
    public void setDRBMode(String drbname){
	this.isDRBGene = true;
	this.drbGeneName = drbname;
    }
        
    public boolean merge(String nucF, String genF, String drbname){
	this.setDRBMode(drbname);
	return this.merge(nucF, genF);
    }

    public boolean merge(String nucF, String genF){
	return this.merge(nucF, genF, true);
    }

    public boolean merge(String genF, boolean outputSequence){
	BufferedReader genbr = null;
	String genline = null;
	int gensp = 0;
	String mergedRef = null;
	try{
	    genbr = new BufferedReader(new FileReader(genF));
	    genline = getFirstLine(genbr);
	    gensp = genline.indexOf(genline.trim().split("\\s+")[1]);
	    if(genline == null){
		System.err.println("[MergeMSFs] Something wrong with input files. Exiting...");
		//System.exit(0);
		return false;
	    }
	    
	    String genname = genline.substring(0, gensp).trim();
	    this.referenceAllele = genname;

	    this.geneName = referenceAllele.substring(0,referenceAllele.indexOf("*"));
	    
	    String gensequence = genline.substring(gensp).trim();
	    int genlineLen = gensequence.length();
	    
	    String[] genblocks = gensequence.split("\\|");
	    
	    Sequence refSequence = this.mergeAndAdd(genname, null, genblocks);

	    String nameCheck = null;
	    if(this.isDRBGene)
		nameCheck = this.drbGeneName;
	    else
		nameCheck = this.geneName;
	    while( (genline=genbr.readLine()) != null ){
		String allele = genline.substring(0,gensp).trim();
		if(nameCheck.equals(allele.substring(0,allele.indexOf("*")))){
		    if(genlineLen == genline.substring(gensp).trim().length()){
			genblocks = genline.substring(gensp).trim().split("\\|");
			boolean replaceAbbrv = true;
			this.mergeAndAdd(allele, null, genblocks, refSequence);
		    }else{
			System.err.println("[MergeMSFs] Line length does not match.");
			//System.err.println(allele);
			System.err.println("[MergeMSFs] genseq: " + gensequence + "\n" + allele + ": " +  genline.substring(gensp).trim());
		    }
		}else
		    ;//skip alleles that don't have same gene name
	    }
	    genbr.close();
	}catch(IOException ioe){
	    ioe.printStackTrace();
	    return false;
	}
	
	if(outputSequence)
	    this.outToFasta();
	return true;
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
    public boolean merge(String nucF, String genF, boolean outputSequence){
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
	    String[] nuctoks = nucline.trim().split("\\s+", -1);
	    String[] gentoks = genline.trim().split("\\s+", -1);
	    nucsp = nucline.indexOf(nuctoks[1], nuctoks[0].length());
	    gensp = genline.indexOf(gentoks[1], gentoks[0].length());
	    //nucsp = nucline.indexOf(nucline.trim().split("\\s+")[1]);
	    //gensp = genline.indexOf(genline.trim().split("\\s+")[1]);
	    //System.out.println("nucsp: " + nucsp);
	    //System.out.println("gensp: " + gensp);


	    if(nucline == null || genline == null){
		System.err.println("[MergeMSFs] Something wrong with input files. Exiting...");
		//System.exit(0);
		return false;
	    }
	    
	    String nucname = nucline.substring(0, nucsp).trim();
	    String genname = genline.substring(0, gensp).trim();
	    this.referenceAllele = nucname;

	    this.geneName = referenceAllele.substring(0,referenceAllele.indexOf("*"));
	    if(!nucname.equals(genname)){
		System.err.println("REF SEQ names differs :");
		System.err.println("(nuc):" + nucname);
		System.err.println("(gen):" + genname);
		//System.exit(-1);
	    }
	    
	    String nucsequence = nucline.substring(nucsp).trim();
	    String gensequence = genline.substring(gensp).trim();
	    int nuclineLen = nucsequence.length();
	    int genlineLen = gensequence.length();
	    //System.err.println("nucseq: " + nucsequence + "|");
	    //System.err.println("genseq: " + gensequence + "|");

	    String[] nucblocks = nucsequence.split("\\|");
	    String[] genblocks = gensequence.split("\\|");
	    
	    Sequence refSequence=this.mergeAndAdd(nucname, nucblocks, genblocks);
	    /* End of taking care of first sequences */
	    
	    
	    String nameCheck = null;
	    if(this.isDRBGene)
		nameCheck = this.drbGeneName;
	    else
		nameCheck = this.geneName;
	    /* Now get a list of Gens and process together with Nucs if same sequence is contained in nuc */
	    while( (genline=genbr.readLine()) != null ){
		String allele = genline.substring(0,gensp).trim();
		if(nameCheck.equals(allele.substring(0,allele.indexOf("*"))) || FormatIMGT.isExtraDRB(this.geneName)){
		    if(genlineLen == genline.substring(gensp).trim().length()){
			//System.err.println("Putting gen allele[" +genline.substring(0,gensp).trim() + "]");
			//System.err.println("First Block :[" + genline.substring(gensp).trim().split("\\|")[0] + "]");
			genHash.put(genline.substring(0,gensp).trim(), genline.substring(gensp).trim().split("\\|"));
		    }else{
			System.err.println("[MergeMSFs] Line length does not match.");
			//System.err.println(allele);
			System.err.println("[MergeMSFs] genseq: " + gensequence + "\n" + allele + ": " +  genline.substring(gensp).trim());
		    }
		}else
		    ;//skip alleles that don't have same gene name
	    }
	    genbr.close();
	    
	    /* end of processing genbr */

	    /* Take care of nucF since genF is always a subset of nucF */
	    /* No need to process genF */
	    String curline = null;
	    while((curline=nucbr.readLine()) != null){
		genblocks = null;
		
		String allele = curline.substring(0,nucsp).trim();
		if(nameCheck.equals(allele.substring(0,allele.indexOf("*")))){
		    String curseq = curline.substring(nucsp).trim();
		    String msfsequence = null;
		    //System.err.println(allele);
		    //System.err.println(curseq);
		    //if there is a matching gen sequence in the hash.
		    if( (genblocks = genHash.get(allele)) != null){
			//System.err.println("[GENBLOCK EXISTS] Processing:\t" + allele);
			boolean replaceAbbrv = true;
			this.mergeAndAdd(allele, curseq.split("\\|"), genblocks, refSequence);
			//System.err.println("COLUMNSequence:\t" + refSequence.getColumnSequence());
		    }
		    else{ // else we just take the nucsequence
			msfsequence = MergeMSFs.removeBlank(curseq);
			//if we have seen this sequence before --> something is not right
			if(this.allele2Sequence.get(allele) != null){
			    System.err.println("[ERROR:MergeMSFs] DUPLICATE ENTRY? --> " + allele + " (Skipping for now)...");
			    System.exit(1);
			}else{
			    //use nuc only constructor of Sequence to process nuc-only allele
			    this.addAllele(allele, new Sequence(allele, msfsequence, refSequence));//this.allele2Sequence.get(nucname)));
			}
		    }
		}
	    }
	    nucbr.close();
	}catch(IOException ioe){
	    ioe.printStackTrace();
	    return false;
	}

	//this.print();
	if(outputSequence)
	    this.outToFasta();
	return true;
    }
    
    private Sequence mergeAndAdd(String alleleName, String[] nucblocks, String[] genblocks){
	return this.mergeAndAdd(alleleName, nucblocks, genblocks, null);
    }
    /* given nucblocks and genblocks, merge nucblocks(use all exons) and genblocks(use only introns) */
    private Sequence mergeAndAdd(String alleleName, String[] nucblocks, String[] genblocks, Sequence refSequence){
	String mergedSeq = null;
	if(nucblocks == null)
	    mergedSeq = this.mergeGenBlocksOnly(genblocks);
	else
	    mergedSeq = this.mergeBlocks(nucblocks, genblocks);
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
	    System.err.println("[MergeMSFs] nucblocks.length : " + nucblocks.length + " genblocks.length :" + genblocks.length +"\ngenblocks length must be equal to [2 * (nucblocks length) + 1]");
	    System.err.println("[MergeMSFs] ]We will just use genblocks for missing nucblocks");
	    //System.exit(-1);
	}
	StringBuffer bf = new StringBuffer();
	for(int i=0;i<genblocks.length;i++){
	    //i = 1, 3, 5 --> exon
	    if( i%2 == 1){
		if(nucblocks.length>(i/2))
		    bf.append(" | " + nucblocks[i/2].trim() + " | ");
		else
		    bf.append(" | " + genblocks[i].trim() + " | ");
	    }else //intron i = 0, 2, 4
		bf.append(genblocks[i].trim());
	}
	return bf.toString();
    }
    
    private String mergeGenBlocksOnly(String[] genblocks){
	StringBuffer bf = new StringBuffer();
	for(int i=0; i < genblocks.length; i++)
	    bf.append(genblocks[i].trim());
	return bf.toString();
    }
    
    //this returns the first sequnece usually the reference in MSA (nuc or gen).
    public String getFirstLine(BufferedReader br){
	String curline = null;
	try{
	    while((curline=br.readLine())!=null){
		String tmp = curline.trim();
		//"[A-Z]+\\d\\*\\d+:\\d+(:\\d+){0,2}")
		if(tmp.split("\\s+")[0].matches("[A-Z]+\\d{0,1}\\*\\d+(:\\d+){0,3}[A-Z]*")){//"[A-Z]+\\d{0,1}\\*\\d+:\\d+(:\\d+){0,2}[A-Z]*")){
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
	//System.err.println("[B4][" + sequence + "]");
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
	//System.err.println("[A4][" + bf.toString() + "]");
	return bf.toString();
    }
    
    /* this fills *(unknown) and -(same as ref) as bases*/
    public static String abbrv2Seq(String abbrv, String refSeq){
	//if(abbrv.length() != refSeq.length())
	//    System.err.println("SOMETHING WRONG: LENGTHS NOT EQUAL IN ABBRV2SEQ");
	//System.out.println("ABBRV:\t" + abbrv);
	//System.out.println("ABBRV:\t" + refSeq);
	StringBuffer bf = new StringBuffer();
	for(int i=0;i<abbrv.length();i++){
	    char curChar = abbrv.charAt(i);
	    if(curChar == '*')
		bf.append(refSeq.charAt(i));//bf.append(Character.toLowerCase(refSeq.charAt(i)));
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








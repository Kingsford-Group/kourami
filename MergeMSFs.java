import java.io.*;

public class MergeMSFs{

    private String referenceAllele;
    private Hashtable<String, Sequence> allele2Sequence;
    private ArrayList<String> orderedAlleles;
    private StringBuffer header;
    
    public MergedMSFs(){
	this.header = new StringBuffer();
	this.allele2Sequence = new Hashtable<String, Sequence>();
	this.orderedAlleles = new ArrayList<String>();
    }
    
    /*
     * nucF : nuc file containing MSA of coding sequences
     * genF : gen file containing MSA of whole gene sequences
     */
    public void merge(String nucF, String genF){
	BufferedReader nucbr = null;
	BufferedReader genbr = null;
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

	    nucsp = nucline.indexOf(nucline.trim().split("\\s+")[1]);
	    gensp = genline.indexOf(genline.trim().split("\\s+")[1]);

	    if(nucline == null || genline == null){
		System.err.println("Something wrong with input files. Exiting...");
		System.exit(0);
	    }
	    
	    String nucname = nucline.substring(0, nucsp).trim();
	    String genname = genline.substring(0, gensp).trim();
	    
	    String nucsequence = nucline.substring(nucsp).trim();
	    String gensequence = genline.substring(gensp).trim();

	    String[] nucblocks = nucsequence.split("\\|");
	    String[] genblocks = gensequence.split("\\|");
	    
	    String mergedRef = this.mergeBlocks(nucblocks, genblocks);
	    MergeMSFs.removeBlank(mergedRef);
	    this.addAllele(nucname, mergedRef);
	    /* End of taking care of first sequences */
	    
	    /* Take care of nucF since nucF is always a subset of genF */
	    /* No need to process genF */
	    String curline = null;
	    while((curline=nucbr.readLine()) != null){
		String stripped = curline.trim();
		String allele = curline.substring(0,nucsp).trim();
		String msfsequence = curline.substring(nucsp).trim();
		//if we have seen this sequence before --> something is not right
		if(this.allele2Sequence.get(allele) != null)
		    System.err.println("DUPLICATE ENTRY? --> " + allele + " (Skipping for now)...");
		else{
		    this.addAllele(allele, new Sequence(allele, msfsequence, true));
		}
	    }

	}catch(IOException ioe){
	    ioe.printStackTrace();
	}
    }
    

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
    
    private void addNucAllele(String allele, String sequence, boolean abbrv){
	this.orderedAllele.add(allele);
	boolean isNuc = true;
	this.allele2Sequence.put(allele, new Sequence(allele, sequence, abbrv, isNuc));
    }

    private void addAllele(String allele, Sequence seq){
	this.orderedAllele.add(allele);
	this.allele2Sequence.put(allele, seq);
    }
    //used to merge nucblocks and genblocks for the reference sequence.
    private String mergeBlocks(String[] nucblocks, String[] genblocks){
	if( (nucblocks.length + 1) != genblocks.length)
	    System.err.println("nucblocks.length : " + nucblocks.length + " genblocks.length :" + genblocks.length +"\nnucblocks length must be 1 longer than genblocks length");
	StringBuffer bf = new StringBuffer();
	for(int i=0;i<genblocks.length;i++){
	    //intron
	    bf.append(genblocks[i]);
	    //i = 0, 2, 4 --> exon
	    if( i%2 == 0)
		bf.append("|" + nucblocks[i] + "|");
	}
	return bf.toString();
    }
    
    
    public String getFirstLine(BufferedReader br){
	String curline = null;
	while((curline=br.readLine())!=null){
	    String tmp = curline.trim();
	    //"[A-Z]+\\d\\*\\d+:\\d+(:\\d+){0,2}")
	    if(tmp.split("\\s+")[0].matches("[A-Z]+\\d{0,1}\\*\\d+:\\d+(:\\d+){0,2}[A-Z]*")){
		return curline;
	    }
	}
	return null;
    }

    public static String removeBlank(String sequence, boolean modifyHeader){
	MergeMsfs.removeBlank(sequence);
    }

    /* this removes all blank embedded in sequences*/
    public static String removeBlank(String sequence){
	StringBuffer bf = new StringBuffer();
	//int headeri = 0;
	for(int i=0; i<sequence.length(); i++){
	    char tmp = sequence.charAt(i);
	    if(tmp == ' '){
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
	StringBuffer bf = new StringBuffer();
	for(int i=0;i<this.abbrv.length();i++){
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
    
    
}








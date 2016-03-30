import java.util.*;

public class Sequence{

    private ArrayList<Base> seq;
    private String alleleName; //allele name: ex) A:01:01:01:01
    private StringBuffer columnSequence; //columnsSequence is a string containing all the paddings to match the length of MSA
    private StringBuffer fullSequence;   //fullSequence is 5'-3' DNA string without any paddings. So it's either equal to or shorter than columnSequence.
    private int[] boundaries;//contains start position of each boundary, alternating between intron/exon
    private int[] segmentOffsets;

    private int[] cumulativeOffsets;

    public String getAlleleName(){
	return this.alleleName;
    }
    
    public int getColLength(){
	return this.seq.size();
    }
    
    public Base baseAt(int i){
	return this.seq.get(i);
    }

    //basePos <= colPos in MSA
    public int getColPosFromBasePos(int basePos){
	Base b = null;
	for(int p = basePos; p<=fullSequence.length(); p++){
	    b = this.seq.get(p-1);
	    if(this.seq.get(p-1).getBasePos() == basePos)
		return b.getColPos();
	}
	return -1;
    }
    
    public int getBasePosFromColPos(int colPos){
	return this.seq.get(colPos-1).getBasePos();
    }

    public Sequence(){
	this.seq = new ArrayList<Base>();
	this.alleleName = null;
	this.columnSequence = new StringBuffer();
	this.fullSequence = new StringBuffer();
	this.boundaries = null;
	this.segmentOffsets = null;
	this.cumulativeOffsets = null;
    }

    public void printNthBoundary(int n){
	System.out.print(this.alleleName + "\t");
	if(n < this.boundaries.length-1)
	    System.out.println(this.columnSequence.substring(this.boundaries[n]-1 + this.boundaries[n+1]-1));
	else
	    System.out.println(this.columnSequence.substring(this.boundaries[n]-1));
    }

    public void printBoundaries(){
	for(int i=0;i<this.boundaries.length; i++){
	    System.err.println("boundaries("+i+"):" + boundaries[i] + "," + (i==(this.boundaries.length-1) ? "END" : boundaries[i+1]));
	    System.err.println("segmentOff("+i+"):" + segmentOffsets[i] );
	    System.err.println("cumulatOff("+i+"):" + cumulativeOffsets[i] );
	}
    }
    
    public String getColumnSequence(){
	return this.columnSequence.toString();
    }
    
    public String getFullSequence(){
	return this.fullSequence.toString();
    }

    public String getSequenceFromBases(){
	StringBuffer bf = new StringBuffer();
	for(int i=0;i<this.boundaries.length;i++){
	    int curStart = this.boundaries[i]-1;
	    int curE = -1;
	    if(i<(this.boundaries.length-1))
		curE = this.boundaries[i+1]-1;
	    else
		curE = this.seq.size();
	    for(int j=curStart; j<curE; j++){
		bf.append(this.seq.get(j).getBase());
	    }
	    //bf.append("|\n");
	}
	return bf.toString();
    }

    public int[] getSegmentOffsets(){
	return this.segmentOffsets;
    }
    
    public int[] getCumulativeOffsets(){
	return this.cumulativeOffsets;
    }

    public int[] getBoundaries(){
	return this.boundaries;
    }
    
    public String getNthBlockColumnSequence(int n){
	//if it's last block
	if(n == (this.boundaries.length - 1))
	    return this.columnSequence.substring(this.boundaries[n]-1);
	else
	    return this.columnSequence.substring(this.boundaries[n]-1, this.boundaries[n+1]-1);
    }

    //returns nth intron. here n is 0-based
    //   n(0-based)  intronNum   boundariesIndex
    //       0          1              0         
    //       1          2              2
    //       2          3              4
    public ArrayList<Base> getNthIntron(int n){
	int index = n * 2; // from n to boundariesIndex
	int sIndex = this.boundaries[index]-1;
	int eIndex = -1;
	if(index == this.boundaries.length-1)//last intron
	    eIndex = this.seq.size();
	else
	    eIndex = this.boundaries[index+1]-1;
	ArrayList<Base> tmp = new ArrayList<Base>();
	
	for(int i=sIndex; i<eIndex; i++){
	    tmp.add(this.seq.get(i));
	}
	
	return tmp;
	
	//return (ArrayList<Base>) this.seq.subList(sIndex, eIndex);
    }

    public int numBlocks(){
	return this.boundaries.length;
    }

    /* this processes nuc only allele --> only containing EXONS*/
    public Sequence(String allele, String msfSequence, Sequence ref){
	this();
	this.alleleName = allele;
	String[] tokens = msfSequence.split("\\|");

	/* we set the number of blocks same as the reference sequence */
	this.boundaries = new int[ref.numBlocks()];
	this.segmentOffsets = new int[ref.numBlocks()];
	this.cumulativeOffsets = new int[ref.numBlocks()];

	/* first we copy the first intron from the reference */
	this.seq.addAll(ref.getNthIntron(0));
	this.boundaries[0] = ref.getBoundaries()[0];
	this.segmentOffsets[0] = ref.getSegmentOffsets()[0];
	this.cumulativeOffsets[0] = ref.getCumulativeOffsets()[0];
	this.columnSequence.append(ref.getColumnSequence().substring(0, ref.getBoundaries()[1]-1));
	this.fullSequence.append(ref.getFullSequence().substring(0,ref.getBoundaries()[1]-1-ref.getCumulativeOffsets()[0]));
	/*
	  System.err.println("COPYING FIRST INTRON:");
	  System.err.println("CS:" + ref.getColumnSequence().substring(0, ref.getBoundaries()[1]-1));
	  System.err.println("FS:" + ref.getFullSequence().substring(0,ref.getBoundaries()[1]-1-ref.getCumulativeOffsets()[0]));
	*/
	//set offset/curStartColPos accordingly
	int offset = ref.getCumulativeOffsets()[0];
	int curStartColPos = ref.getBoundaries()[1];
	
	boolean isExon = true;
	int exonNum = 0;
	/* foreach exon*/
	for(int i=0; i<tokens.length; i++){
	    /* process current exon index: 2*i+1 */
	    exonNum++;
	    //exon sequence is still in abbrv form with '-' matching base for ref.
	    //need to modify the sequence
	    //System.err.println(tokens[i]);
	    //System.err.println(ref.getNthBlockColumnSequence(2*i+1));
	    String modseq = MergeMSFs.abbrv2Seq(tokens[i], ref.getNthBlockColumnSequence(2*i+1));
	    /*
	      System.err.println(this.columnSequence + " (CS)");
	      System.err.println(this.fullSequence + " (FS)");
	      System.err.println("COPYING EXON " + (i+1) + ":");
	      System.err.println(modseq);
	    */
	    //cumulative offset
	    int updatedOffset = processBlock(modseq, isExon, exonNum, offset);
	    //current block offset
	    this.segmentOffsets[2*i+1] = updatedOffset - offset;
	    this.cumulativeOffsets[0] = updatedOffset;
	    //update offset with cumOffset
	    offset = updatedOffset;
	    this.boundaries[2*i+1] = curStartColPos;
	    
	    /* process next intron index: 2*(i+1)*/
	    this.seq.addAll(ref.getNthIntron(i+1));
	    this.boundaries[2*(i+1)] = ref.getBoundaries()[2*(i+1)];
	    this.segmentOffsets[2*(i+1)] = ref.getSegmentOffsets()[2*(i+1)];
	    offset = offset + ref.segmentOffsets[2*(i+1)];
	    this.cumulativeOffsets[2*(i+1)] = offset;
	    if(i < (tokens.length-1)){
		//System.err.println("intron ("+(2*(i+1))+")\tCI(" + (ref.getBoundaries()[2*(i+1)]-1) + "," + (ref.getBoundaries()[2*(i+1)+1]-1) +")");
		//System.err.println(ref.getColumnSequence().substring(ref.getBoundaries()[2*(i+1)]-1, ref.getBoundaries()[2*(i+1)+1]-1));
		this.columnSequence.append(ref.getColumnSequence().substring(ref.getBoundaries()[2*(i+1)]-1, ref.getBoundaries()[2*(i+1)+1]-1));
		//System.err.println("intron ("+(2*(i+1))+")\tFI(" + (ref.getBoundaries()[2*(i+1)]-1-ref.getCumulativeOffsets()[2*(i+1)-1]) + "," 
		//		   + (ref.getBoundaries()[2*(i+1)+1]-1-ref.getCumulativeOffsets()[2*(i+1)]) + ")");
		//System.err.println(ref.getFullSequence().substring(ref.getBoundaries()[2*(i+1)]-1-ref.getCumulativeOffsets()[2*(i+1)-1]
		//						   , ref.getBoundaries()[2*(i+1)+1]-1-ref.getCumulativeOffsets()[2*(i+1)]));
		this.fullSequence.append(ref.getFullSequence().substring(ref.getBoundaries()[2*(i+1)]-1-ref.getCumulativeOffsets()[2*(i+1)-1]
									 , ref.getBoundaries()[2*(i+1)+1]-1-ref.getCumulativeOffsets()[2*(i+1)]));
		curStartColPos = ref.getBoundaries()[2*(i+1)+1]; //need to fetch next 
	    }else{//if it's last exon
		//System.err.println("intron ("+(2*(i+1))+")\tCI(" + (ref.getBoundaries()[2*(i+1)]-1) + ",END" );
		//System.err.println(ref.getColumnSequence().substring(ref.getBoundaries()[2*(i+1)]-1));
		this.columnSequence.append(ref.getColumnSequence().substring(ref.getBoundaries()[2*(i+1)]-1));
		//System.err.println("intron ("+(2*(i+1))+")\tFI(" + (ref.getBoundaries()[2*(i+1)]-1-ref.getCumulativeOffsets()[2*(i+1)-1]) + ",END" );
		//System.err.println(ref.getFullSequence().substring(ref.getBoundaries()[2*(i+1)]-1-ref.getCumulativeOffsets()[2*(i+1)-1]));
		this.fullSequence.append(ref.getFullSequence().substring(ref.getBoundaries()[2*(i+1)]-1-ref.getCumulativeOffsets()[2*(i+1)-1]));
	    }
		
		 
	}
	
    }


    //msfSequence still has intron/exon boundary symbol Embedded.
    // <INTRON1>|<EXON1>|<INTRON2>|<EXON2>|...
    // allele --> allelename
    // msfSequence --> msf sequence string without blanks
    public Sequence(String allele, String msfSequence){
	this();
	this.alleleName = allele;
	String[] tokens = msfSequence.split("\\|");
	this.boundaries = new int[tokens.length];
	this.segmentOffsets = new int[tokens.length];
	this.cumulativeOffsets = new int[tokens.length];
	int offset = 0;
	boolean isExon = false;
	int intronNum = 0;
	int exonNum = 0;
	int curStartColPos = 0;
	//for introns and exons
	for(int i=0; i<tokens.length; i++){
	    curStartColPos++;
	    this.boundaries[i] = curStartColPos;
	    if(i%2 == 0){//intron 0, 2, 4, 6
		isExon = false;
		intronNum++;
		int updatedOffset = processBlock(tokens[i], isExon, intronNum, offset);
		this.segmentOffsets[i] = updatedOffset - offset;
		this.cumulativeOffsets[i] = updatedOffset;
		offset = updatedOffset;
	    }else{//exon 1, 3, 5, 7
		isExon = true;
		exonNum++;
		int updatedOffset = processBlock(tokens[i], isExon, exonNum, offset);
		this.segmentOffsets[i] = updatedOffset - offset;
		this.cumulativeOffsets[i] = updatedOffset;
		offset = updatedOffset;
	    }
	    curStartColPos = this.seq.size();
	}
    }

    
    //given msf formatted(no blanks) sequence and add bases 
    //returns the base2coloffset.
    public int processBlock(String blockSeq, boolean isExon, int intronExonNum, int offset){
	int colPos = this.seq.size(); //1-based column Position --> this is the last position of previous block
	int base2colOffset = offset; 
	int basePos = colPos - base2colOffset; //also last position of previous block
	
	for(int i=0; i<blockSeq.length(); i++){
	    char curBase = blockSeq.charAt(i);
	    this.columnSequence.append(curBase);
	    colPos++;
	    if(Base.isBase(curBase)){
		basePos++;
		this.fullSequence.append(curBase);
	    }else if(Base.isGap(curBase))
		base2colOffset++;
	    else
		System.err.println("WHAT ELSE????\nBlockSeq:" + blockSeq + "\n@"  + (i+1) + ":" + curBase);
		    
	    this.seq.add(new Base(blockSeq.charAt(i), basePos, colPos, base2colOffset, isExon, intronExonNum));
	}
	
	return base2colOffset;
    }
    
}

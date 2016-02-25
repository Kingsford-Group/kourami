import java.util.*;

public class Sequence{

    private ArrayList<Base> seq;
    private String alleleName; //allele name: ex) A:01:01:01:01
    private StringBuffer columnSequence; //columnsSequence is a string containing all the paddings to match the length of MSA
    private StringBuffer fullSequence;   //fullSequence is 5'-3' DNA string without any paddings. So it's either equal to or shorter than columnSequence.
    private int[] boundaries;//contains start position of each boundary, alternating between intron/exon
    private int[] segmentOffsets;
    
    public Sequence(){
	this.seq = new ArrayList<Base>();
	this.alleleName = null;
	this.columnSequence = null;
	this.fullSequence = null;
	this.boundaries = null;
	this.segmentOffsets = null;
    }

    public ArrayList<Base> getNthIntron(int n){
	int index = n * 2;
	int sIndex = this.boundaries[index];
	int eIndex = -1;
	if(index == this.boundaries.length-1)//last intron
	    eIndex = this.seq.size()+1;
	else
	    eIndex = this.boundaries[index+1];
	return (ArrayList<Base>) this.seq.subList(sIndex, eIndex);
    }

    //msfSequence still has intron/exon boundary symbol Embedded.
    // <INTRON1>|<EXON1>|<INTRON2>|<EXON2>|...
    public Sequence(String msfSequence){
	this();
	String[] tokens = msfSequence.split("\\|");
	this.boundaries = new int[tokens.length];
	this.segmentOffsets = new int[tokens.length];
	int offset = 0;
	boolean isExon = false;
	int intronNum = 0;
	int exonNum = 0;
	int curStartColPos = 0;
	for(int i=0; i<tokens.length; i++){
	    curStartColPos++;
	    this.boundaries[i] = curStartColPos;
	    if(i%2 == 0){//intron
		isExon = false;
		intronNum++;
		int updatedOffset = processBlock(tokens[i], isExon, intronNum, offset);
		this.segmentOffsets[i] = updatedOffset - offset;
		offset = updatedOffset;
	    }else{
		isExon = true;
		exonNum++;
		int updatedOffset = processBlock(tokens[i], isExon, exonNum, offset);
		this.segmentOffsets[i] = updatedOffset - offset;
		offset = updatedOffset;
	    }
	    curStartColPos = this.seq.size();
	}
    }

    public Sequence(String exonOnlyMsfSequence, Sequence genSeq){
	this();
	String[] tokens = msfSequence.split("\\|");
	this.boundaries = new int[genSeq.getBoundaries.length];
	this.segmentOffsets = new int[this.boundaries.length];
	int offset = 0;
	int exonNum = 0;
	int curStartColPos = 0;
	//for each exon
	for(int i=0; i<tokens.length; i++){
	    genSeq.
	}
    }
    
    
    public int processBlock(String blockSeq, boolean isExon, int intronExonNum, int offset){
	int colPos = this.seq.size(); //1-based column Position
	int base2colOffset = offset;
	int basePos = colPos - base2colOffset;
	
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
		System.err.out.println("WHAT ELSE????\nBlockSeq:" + blockSeq + "\n@"  + (i+1) + ":" curBase);
		    
	    this.seq.add(new Base(blockSeq.charAt(i), basePos, colPos, base2colOffset, isExon, intronExonNum));
	}
	
	return base2colOffset;
    }
    
}

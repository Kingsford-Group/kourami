import java.util.*;

public class Sequence{

    private ArrayList<Base> seq;
    private String alleleName; //allele name: ex) A:01:01:01:01
    private StringBuffer columnSequence; //columnsSequence is a string containing all the paddings to match the length of MSA
    private StringBuffer fullSequence;   //fullSequence is 5'-3' DNA string without any paddings. So it's either equal to or shorter than columnSequence.
    
    public Sequence(){
	this.seq = new ArrayList<Base>();
	this.alleleName = null;
	this.columnSequence = null;
	this.fullSequence = null;
    }

    //msfSequence still has intron/exon boundary symbol Embedded.
    // <INTRON1>|<EXON1>|<INTRON2>|<EXON2>|...
    public Sequence(String msfSequence){
	String[] tokens = msfSeqquence.split("\\|");
	int offset = 0;
	boolean isExon = false;
	int intronNum = 0;
	int exonNum = 0;
	for(int i=0; i<tokens.length; i++){
	    if(i%2 == 0){//intron
		isExon = false;
		intronNum++;
		offset = processBlock(tokens[i], isExon, intronNum, offset);
	    }else{
		isExon = true;
		exonNum++;
		offset = processBlock(tokens[i], isExon, exonNum, offset);
	    }	    
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

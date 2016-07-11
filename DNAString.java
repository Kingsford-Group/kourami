import java.util.*;

public class DNAString{

    private String hlaGeneName;
    private int sbnum; //superbubble number
    private ArrayList<Integer> pathnums; // path number
    private double weightedIntersectionSum;
    private double probability;
    
    private StringBuffer sequence;
    
    public DNAString deepCopy(){
	return new DNAString(this.hlaGeneName, this.sbnum, this.pathnums
			     , this.weightedIntersectionSum, this.probability
			     , this.sequence.toString());
    }

    public DNAString mergeDeep(DNAString other){
	DNAString td = this.deepCopy();
	return td.merge(other);
    }

    private DNAString merge(DNAString other){
	if(this.hlaGeneName.equals(other.getHlaGeneName())){
	    for(Integer i: other.getPathnums())
		this.pathnums.add(i);
	    this.weightedIntersectionSum += other.getWeightedIntersectionSum();
	    this.probability = this.probability * other.getProbability();
	    this.sequence.append(other.getSequence());
	}else{
	    System.err.println("[DNAString.merge(DNAString other)]HLA Gene Names don't match:\t"
			       +this.hlaGeneName + "\t" + other.getHlaGeneName());
	}
	return this;
    } 

    public String pathnums2String(){
	StringBuffer bf = new StringBuffer();
	int count = 0;
	for(Integer i: this.pathnums){
	    if(count > 0)
		bf.append(":");
	    bf.append(i.intValue());
	    count++;
	}
	return bf.toString();
    }

    public StringBuffer toFasta(){
	StringBuffer bf = new StringBuffer(">"+ this.hlaGeneName + "_" + sbnum + "_" 
					   + pathnums2String() + "-" 
					   + this.weightedIntersectionSum + this.probability + "\n");
	bf.append(sequence.toString() + "\n");
	return bf;
    }

    public DNAString(String hgn, int sbn, int pn, double wis, double p){
	this.hlaGeneName = hgn;
	this.sbnum = sbn;
	this.pathnums = new ArrayList<Integer>();
	this.pathnums.add(pn);
	this.weightedIntersectionSum = wis;
	this.probability = p;
	this.sequence = new StringBuffer();
    }
    
    public DNAString(String hgn, int sbn, int pn, double wis, double p, String seq){
	this.hlaGeneName = hgn;
	this.sbnum = sbn;
	this.pathnums = new ArrayList<Integer>();
	this.pathnums.add(pn);
	this.weightedIntersectionSum = wis;
	this.probability = p;
	this.sequence = new StringBuffer();
	this.sequence.append(seq);
    }

    public DNAString(String hgn, int sbn, ArrayList<Integer> pns, double wis, double p, String seq){
	this.hlaGeneName = hgn;
	this.sbnum = sbn;
	this.pathnums = new ArrayList<Integer>();
	for(Integer i : pns)
	    this.pathnums.add(i);
	this.weightedIntersectionSum = wis;
	this.probability = p;
	this.sequence.append(seq);
    }
    
    public void append(String seg){
	this.sequence.append(seg);
    }

    public String getHlaGeneName(){
	return this.hlaGeneName;
    }
    
    public int getSbnum(){
	return this.sbnum;
    }
    
    public ArrayList<Integer> getPathnums(){
	return this.pathnums;
    }
    
    public int getFirstPathnum(){
	return this.pathnums.get(0).intValue();
    }
    
    public double getWeightedIntersectionSum(){
	return this.weightedIntersectionSum;
    }
    
    public double getProbability(){
	return this.probability;
    }

    public String getSequence(){
	return this.sequence.toString();
    }
    
}

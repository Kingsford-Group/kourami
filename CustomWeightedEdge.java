import htsjdk.samtools.util.QualityUtil;

import org.jgrapht.*;
import org.jgrapht.graph.*;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.Collections;

/*
 * Added in order to assign confidence to the weights/baescall of the edge.
 * Current scoring is taken from MIRA. (stranded)
 *
 */
public class CustomWeightedEdge extends DefaultWeightedEdge{
    
    private ArrayList<Byte> fScore;
    private ArrayList<Byte> rScore;
    private double groupErrorProb;

    private int numActivePath;

    public static int numMaxLowestProbEntries = 10;
    
    public CustomWeightedEdge(){
	super();
	this.fScore = new ArrayList<Byte>();
	this.rScore = new ArrayList<Byte>();
	this.groupErrorProb = 0.0d;
	this.initNumActivePath();
    }
    
    public void initNumActivePath(){
	this.numActivePath = 0;
    }
    
    public void inlcudeEdge(){
	this.numActivePath++;
    }
    
    public void excludeEdge(){
	this.numActivePath--;
    }

    public ArrayList<Byte> getFScores(){
	return this.fScore;
    }

    public ArrayList<Byte> getRScores(){
	return this.rScore;
    }

    public void setFScores(ArrayList<Byte> fs){
	this.fScore = fs;
    }
    
    public void setRScores(ArrayList<Byte> rs){
	this.rScore = rs;
    }

    public void addAllFScores(ArrayList<Byte> fs){
	this.fScore.addAll(fs);
    }
    
    public void addAllRScores(ArrayList<Byte> rs){
	this.rScore.addAll(rs);
    }

    public String toString(){
	return this.getWeight() + "\t" + this.groupErrorProb;// + "\tfScoreSize[" + this.fScore.size() +"]" + "\trScoreSize[" + this.rScore.size() +"]";
    }

    public void incrementWeight(SimpleDirectedWeightedGraph<Node, CustomWeightedEdge> g, boolean isRefStrand, byte score){
	g.setEdgeWeight(this, g.getEdgeWeight(this)+1);
	if(isRefStrand)
	    this.fScore.add(new Byte(score));
	else
	    this.rScore.add(new Byte(score));
    }
    
    public double getGroupErrorProb(){
	return this.groupErrorProb;
    }

    public double computeGroupErrorProb(){
	if(this.getWeight() == Double.MAX_VALUE)
	    this.groupErrorProb = 0.0d;
	else{
	    double f = this.computeStrandedGroupErrorProb(this.fScore);
	    double r = this.computeStrandedGroupErrorProb(this.rScore);
	    //System.err.println("f-prob:\t" + f);
	    //System.err.println("r-prob:\t" + r);
	    if(f == 0.0d)
		f = 1.0d;
	    if(r == 0.0d)
		r = 1.0d;
	    double score = f*r;
	    if(score == 1.0d)
		score = 0.0d;
	    this.groupErrorProb = score;
	}
	//System.err.println("GroupErrorProb:\t" + this.groupErrorProb);
	return this.groupErrorProb;
    }

    private double computeStrandedGroupErrorProb(ArrayList<Byte> scores){
	//sort the phred score in a descending order
	Collections.sort(scores, Collections.reverseOrder());
	Iterator<Byte> itr = scores.iterator();
	int count = 0;
	int sum = 0;
	while(itr.hasNext() && count < CustomWeightedEdge.numMaxLowestProbEntries){
	    byte tmp = itr.next().byteValue();
	    //System.err.println("ErrorProb(" + tmp + "):\t" + QualityUtil.getErrorProbabilityFromPhredScore(tmp));
	    int invErrorProb = (int) (1.0d / QualityUtil.getErrorProbabilityFromPhredScore(tmp));
	    //int invErrorProb = (int) (1.0d / QualityUtil.getErrorProbabilityFromPhredScore(itr.next().byteValue()));
	    //System.err.println("INVErroProb(" + tmp + "):\t" + invErrorProb);
	    double fraction = (double) (CustomWeightedEdge.numMaxLowestProbEntries - count) / CustomWeightedEdge.numMaxLowestProbEntries;
	    //System.err.println("Fraction(" + tmp + "):\t" + fraction);
	    sum += (int) (fraction * invErrorProb);
	    //System.err.println("DENOM_SUM(" + tmp + "):\t" + sum);
	    count++;
	}
	if(sum == 0)
	    return 0.0d;
	else
	    return 1.0d/((double)sum);
    }


    /* TESTING */
    public static void main(String[] args){
	CustomWeightedEdge e = new CustomWeightedEdge();
	e.fScore.add(new Byte((byte) 10));
	e.fScore.add(new Byte((byte) 20));
	e.fScore.add(new Byte((byte) 30));
	
	e.rScore.add(new Byte((byte) 10));
	e.rScore.add(new Byte((byte) 20));
	e.rScore.add(new Byte((byte) 30));

	e.computeGroupErrorProb();
    }
}

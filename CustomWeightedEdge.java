import htsjdk.samtools.util.QualityUtil;

import org.jgrapht.*;
import org.jgrapht.graph.*;

import java.util.TreeSet;
import java.util.Iterator;

/*
 * Added in order to assign confidence to the weights/baescall of the edge.
 * Current scoring is taken from MIRA. (stranded)
 *
 */
public class CustomWeightedEdge extends DefaultWeightedEdge{
    
    private TreeSet<Byte> fScore;
    private TreeSet<Byte> rScore;
    private double groupErrorProb;
    
    public static int numMaxLowestProbEntries = 10;
    
    public CustomWeightedEdge(){
	super();
	this.fScore = new TreeSet<Byte>();
	this.rScore = new TreeSet<Byte>();
	this.groupErrorProb = 0.0d;
	
    }
    
    public void incrementWeight(SimpleDirectedWeightedGraph<Node, CustomWeightedEdge> g, boolean direction, byte score){
	g.setEdgeWeight(this, g.getEdgeWeight(this)+1);
	if(direction)
	    this.fScore.add(new Byte(score));
	else
	    this.rScore.add(new Byte(score));
    }
    
    
    public double computeGroupErrorProb(){
	double f = this.computeStrandedGroupErrorProb(this.fScore);
	double r = this.computeStrandedGroupErrorProb(this.rScore);
	if(f == 0.0d)
	    f = 1.0d;
	if(r == 0.0d)
	    r = 1.0d;
	double score = f*r;
	if(score == 1.0d)
	    score = 0.0d;
	this.groupErrorProb = score;
	return score;
    }

    private double computeStrandedGroupErrorProb(TreeSet<Byte> scoreSet){
	Iterator<Byte> itr = this.fScore.iterator();
	int count = 0;
	int sum = 0;
	while(itr.hasNext() && count < CustomWeightedEdge.numMaxLowestProbEntries){
	    int invErrorProb = (int) (1.0d / QualityUtil.getErrorProbabilityFromPhredScore(itr.next()));
	    sum = (CustomWeightedEdge.numMaxLowestProbEntries - count) / CustomWeightedEdge.numMaxLowestProbEntries * invErrorProb;
	    count++;
	}
	if(sum == 0)
	    return 0.0d;
	else
	    return 1/sum;
    }

}

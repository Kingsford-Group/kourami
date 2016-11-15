import java.util.ArrayList;
import java.util.Hashtable;

import org.jgrapht.*;
import org.jgrapht.graph.*;

public class SuperAllelePath{

    private ArrayList<AllelePath> orderedAllelePaths;
    private ArrayList<Integer> pathNums;
    private String hlaGeneName;
    
    public SuperAllelePath(String genename){
	this.orderedAllelePaths = new ArrayList<AllelePath>();
	this.pathNums = new ArrayList<Integer>();
	this.hlaGeneName = genename;
    }
    
    public SuperAllelePath(ArrayList<AllelePath> aps, ArrayList<Integer> pns, String hgn){
	this.orderedAllelePaths = aps;
	this.pathNums = pns;
	this.hlaGeneName = hgn;
    }

    public double[] jointTraverse(SuperAllelePath other, SimpleDirectedWeightedGraph<Node, CustomWeightedEdge> g){
	double[] weightFlow = new double[2];
	for(int i=0; i<this.orderedAllelePaths.size(); i++){
	    AllelePath tap = this.orderedAllelePaths.get(i);
	    AllelePath oap = other.getOrderedAllelePaths().get(i);
	    double[] tmp = tap.jointTraverse(oap, g);
	    weightFlow[0] += tmp[0];
	    weightFlow[1] += tmp[1];
	}
	return weightFlow;
    }

    public double[] traverse(SimpleDirectedWeightedGraph<Node, CustomWeightedEdge> g){
	double[] weightFlow = new double[2];
	for(AllelePath ap : orderedAllelePaths){
	    double[] tmp = ap.traverse(g);
	    weightFlow[0] += tmp[0];
	    weightFlow[1] += tmp[1];
	}
	return weightFlow;
    }

    public int numAllelePaths(){
	return this.orderedAllelePaths.size();
    }
    
    //return 3 scores (TP/OP and intersection score) + 4 scores (allProduct, jointProduct, avgProduct, maxFlow)
    public double[] getJointProbability(SuperAllelePath other, ArrayList<Bubble> superBubbles){
	double[] jp = new double[6];
	if(this.numAllelePaths() != other.numAllelePaths() && this.numAllelePaths() != superBubbles.size()){
	    System.err.println("Incompatible SuperAllelePath. The number of fractured allelePath in superPath does not match");
	    return null;
	}
	
	for(int i=0; i<this.orderedAllelePaths.size(); i++){
	    double[] apjp = this.orderedAllelePaths.get(i).getJointProbability(other.getOrderedAllelePaths().get(i), superBubbles.get(i));
	    for(int j=0; j<jp.length-1;j++)
		jp[j] += apjp[j];
	    /*
	    if(i==0)
		jp[6] = apjp[6];
	    else{
		if(jp[6] > apjp[6])
		    jp[6] = apjp[6];
		    }*/
	}
	return jp;
    }


    public double getJointInterSuperBubbleLinkProb(SuperAllelePath other, Hashtable<Path, Hashtable<Path, int[]>> hhl){
	double logP = 0.0d;
	for(int i=0; i<this.orderedAllelePaths.size(); i++){
	    Path tp_i = this.orderedAllelePaths.get(i).getBubblePath();
	    Path op_i = other.getOrderedAllelePaths().get(i).getBubblePath();
	    for(int j=i+1; j<this.getOrderedAllelePaths().size(); j++){
		Path tp_j = this.orderedAllelePaths.get(j).getBubblePath();
		Path op_j = other.getOrderedAllelePaths().get(i).getBubblePath();
		int[] tVals = hhl.get(tp_i).get(tp_j); // vals[0]:#linker , vals[1]:SUM
		int[] oVals = hhl.get(op_i).get(op_j);
		if(tVals[1] == oVals[1]){
		    double fraction = ((tVals[0] + oVals[0])*1.0d) / (tVals[1]*1.0d);
		    logP += Math.log(fraction);
		}else{
		    System.err.println("STRANGE!!!!!!!!!! THEY SHOULD BE SAME");
		    System.exit(0);
		}
	    }
	}
	return logP;
    }
    
    
    public ArrayList<AllelePath> getOrderedAllelePaths(){
	return this.orderedAllelePaths;
    }


    public String getHlaGeneName(){
	return this.hlaGeneName;
    }

    public SuperAllelePath clone(){
	return new SuperAllelePath(new ArrayList<AllelePath>(this.orderedAllelePaths)
				   , new ArrayList<Integer>(this.pathNums)
				   , this.hlaGeneName);
    }

    public void addAllelePath(AllelePath ap, int n){
	this.orderedAllelePaths.add(ap);
	this.pathNums.add(new Integer(n));
    }
    
    public String pathnums2String(){
	StringBuffer bf = new StringBuffer();
	int count = 0;
	for(Integer i : this.pathNums){
	    if(count > 0)
		bf.append(":");
	    bf.append(i.intValue());
	    count++;
	}
	return bf.toString();
    }


    public StringBuffer getSequenceBuffer(){
	StringBuffer bf = new StringBuffer();
	for(AllelePath ap : this.orderedAllelePaths)
	    bf.append(ap.getSequence());
	return bf;
    }

    public String toSimpleString(){
	return this.pathnums2String() + "\t" + this.getWeightedIntersectionSumScore() + "\t" + this.getProbability();
    }
    
    public StringBuffer toFasta(){
	StringBuffer bf = new StringBuffer(">" + this.hlaGeneName + "_" + pathnums2String() + "\t" 
					   + this.getWeightedIntersectionSumScore() + "\t" + this.getProbability() + "\n");
	bf.append(this.getSequenceBuffer());
	bf.append("\n");
	return bf;
    }

    public double getProbability(){
	double logP = 0.0d;
	for(AllelePath ap : this.orderedAllelePaths){
	    logP += ap.getProbability();
	}
	return logP;
    }

    /*
    public double getProbability(){
	double p = 1.0d;
	for(AllelePath ap : this.orderedAllelePaths)
	    p = p * ap.getProbability();
	return p;
    }
    */

    public double getWeightedIntersectionSumScore(){
	double s = 0.0d;
	int n = 0;
	for(AllelePath ap : this.orderedAllelePaths){
	    s = s + ap.getIntersectionSum();
	    n = n + ap.getMergedNums();
	}
	return s / n;
    }
    
}

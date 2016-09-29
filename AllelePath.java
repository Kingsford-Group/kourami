import java.util.ArrayList;
import java.util.HashSet;
import java.util.Iterator;

import org.jgrapht.*;
import org.jgrapht.graph.*;


//single path through a superbubble.
//each superbubble can have multiple AllelePath
public class AllelePath{
    
    private ArrayList<CustomWeightedEdge> orderedEdgeList;
    
    //this keeps track of indices where disconnect happens due to exonic boundary.
    private ArrayList<Integer> fractureEndIndex;

    private double probability;

    private double weightedIntersectionSum;
    
    private int mergedNums;

    public AllelePath(){
	this.orderedEdgeList = new ArrayList<CustomWeightedEdge>();
	this.fractureEndIndex = new ArrayList<Integer>();
	this.weightedIntersectionSum = 0.0d;
	this.mergedNums = 0;
	this.probability = 1.0d;
    }

    public AllelePath(double p, double wis, int mn){
	this();
	this.weightedIntersectionSum = wis;
	this.mergedNums = mn;
	this.probability = p;
    }
    
    public void setFractureEndIndex(){
	this.fractureEndIndex.add(new Integer(this.orderedEdgeList.size()));
    }

    public void appendEdge(CustomWeightedEdge e){
	this.orderedEdgeList.add(e);
    }
    
    public void appendAllEdges(Path other){
	this.orderedEdgeList.addAll(other.getOrderedEdgeList());
    }

    public void appendAllEdges(ArrayList<CustomWeightedEdge> anotherOrderedEdgeList){
	this.orderedEdgeList.addAll(anotherOrderedEdgeList);
    }

    public void printPath(SimpleDirectedWeightedGraph<Node, CustomWeightedEdge> g, int superBubbleNum, int n){
	System.out.println(">candidate_" + superBubbleNum + "-" + n + "\n" + this.toString(g, superBubbleNum, n));
    }

    public String toString(SimpleDirectedWeightedGraph<Node, CustomWeightedEdge> g, int superBubbleNum, int n){
	StringBuffer bf = new StringBuffer();
	CustomWeightedEdge pre = null;
	int disconnectCount = 0;
	for(int i=0; i<this.orderedEdgeList.size(); i++){
	    CustomWeightedEdge cur = this.orderedEdgeList.get(i);
	    if(cur == null){
		System.err.print("TMP");
		System.err.println("HUH??????:\t" + i + "\tSize:\t" + this.orderedEdgeList.size() );
	    }
	    char curChar;
	    if( (pre != null && !g.getEdgeTarget(pre).equals(g.getEdgeSource(cur)))
		|| pre == null ){
		if(pre != null)
		    disconnectCount++;
		curChar = g.getEdgeSource(cur).getBase();
		if(curChar != '.')
		    bf.append(curChar);
	    }
	    curChar = g.getEdgeTarget(cur).getBase();
	    if(curChar != '.')
		bf.append(curChar);
	    pre = cur;
	}
	String finalStr = bf.toString();
	System.err.println(">candidate_" + superBubbleNum + "-" + n + "(" + disconnectCount + ")\n" + finalStr);
	return finalStr;
    }
}

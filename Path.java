import java.util.ArrayList;
import java.util.HashSet;

import org.jgrapht.*;
import org.jgrapht.graph.*;

public class Path{

    private ArrayList<CustomWeightedEdge> orderedEdgeList;
    
    public Path deepCopy(){
	Path p = new Path();
	for(CustomWeightedEdge e : this.orderedEdgeList){
	    p.appendEdge(e);
	}
	return p;
    }

    public ArrayList<CustomWeightedEdge> getOrderedEdgeList(){
	return this.orderedEdgeList;
    }

    public Path(){
	this.orderedEdgeList = new ArrayList<CustomWeightedEdge>();
    }

    public void appendEdge(CustomWeightedEdge e){
	this.orderedEdgeList.add(e);
    }

    public Path(CustomWeightedEdge e){
	this();
	this.appendEdge(e);
    }

    public Node getLastVertex(SimpleDirectedWeightedGraph<Node, CustomWeightedEdge> g){
	if(this.orderedEdgeList == null || this.orderedEdgeList.size() == 0)
	    return null;
	else{
	    return g.getEdgeTarget(this.orderedEdgeList.get(this.orderedEdgeList.size() - 1));
	}
    }
    
    public Path mergePaths(Path other){
	Path p = this.deepCopy();
	p.getOrderedEdgeList().addAll(other.getOrderedEdgeList());
	return p;
    }
    
    //checks for phasing support information between two paths
    //Two paths are phased if intersection of reads passing through unique edges.
    //@returns true if intersecting sets are NOT empty
    //@returns false otherwise.
    public boolean isPhasedWith(Path other){
	//this set
	HashSet<Integer> ts = this.getUnionOfUniqueEdgesReadSet();
	//other set
	HashSet<Integer> os = other.getUnionOfUniqueEdgesReadSet();
	
	ts.retainAll(os);
	if(ts.size() > 0)
	    return true;
	return false;
    }

    //if there is no uniqueEdges, return null
    //else it returns union HashSet<Integer> of reads over all unique edges.
    //size 0 if there is reads covering unique edge
    public HashSet<Integer> getUnionOfUniqueEdgesReadSet(){
	HashSet<Integer> s = new HashSet<Integer>();
	boolean atLeast1UniqueEdge = false;
	for(CustomWeightedEdge e : this.orderedEdgeList){
	    if(e.isUniqueEdge()){
		atLeast1UniqueEdge  = true;
		s.addAll(e.getReadHashSet());
	    }
	}
	if(!atLeast1UniqueEdge)
	    return null;
	return s;
    }
    
}

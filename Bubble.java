import java.io.*;
import java.util.*;

public class Bubble{

    private HLAGraph g;
    //private int start;
    //private int end;
    private ArrayList<Integer> start;
    private ArrayList<Integer> end;
    
    private ArrayList<Node> sNodes;
    private ArrayList<Node> tNodes;
    //    private Node s;
    //private Node t;
    
    private ArrayList<Path> paths;

    public ArrayList<Integer> getStart(){
	return this.start;
    }
    
    public ArrayList<Integer> getEnd(){
	return this.end;
    }
    
    public ArrayList<Node> getSNodes(){
	return this.sNodes;
    }

    public ArrayList<Node> getTNodes(){
	return this.tNodes;
    }


    public ArrayList<Path> getPaths(){
	return this.paths;
    }

    public void Bubble(HLAGraph hg, Node s, Node t){
	this.g = hg;
	this.sNodes = new ArrayList<Node>();
	this.tNodes = new ArrayList<Node>();
	this.sNodes.add(s);
	this.tNodes.add(t);
	this.start = new ArrayList<Integer>();
	this.end = new ArrayList<Integer>();
	this.start.add(new Integer(s.getColIndex()));
	this.end.add(new Integer(t.getColIndex()));
	this.paths = new ArrayList<Path>();
	this.decompose(s, t);
    }
    
    //find all ST path in the bubble
    //and remove unsupported paths
    public ArrayList<Path> decompose(Node s, Node t){
	this.paths =  this.g.findAllSTPath(s, t);
	this.removeUnsupported();
	return this.paths;
    }
    

    public void removeUnsupported(){
	HashSet<Integer> readHash;
	ArrayList<Integer> removalList = new ArrayList<Integer>();
	for(int i=0; i<this.paths.size(); i++){//Path p: this.paths){
	    Path p = this.paths.get(i);
	    ArrayList<CustomWeightedEdge> eList = p.getOrderedEdgeList();
	    boolean inited = false;
	    readHash = null;
	    for(CustomWeightedEdge e : eList){
		if(e.isUniqueEdge()){ //we only care about uniq edges
		    if(!inited){
			readHash = e.getReadHashSetDeepCopy();
			inited = true;
		    }else{
			//update with union after checking intersection
			//if null, we remove this path
			if(e.unionAfterCheckingIntersection(readHash) == null){
			    removalList.add(new Integer(i));
			    break;
			}
		    }
		}
	    }
	}
	for(int i=removalList.size() - 1; i >= 0; i--){
	    this.paths.remove(removalList.get(i));
	}
    }

    /*
    public void removeUnsupported(){
	HashSet<Integer> readHash;
	ArrayList<Integer> removalList = new ArrayList<Integer>();
	for(int i=0; i<this.paths.size(); i++){//Path p: this.paths){
	    Path p = this.paths.get(i);
	    ArrayList<CustomWeightedEdge> eList = p.getOrderedEdgeList();
	    boolean inited = false;
	    for(CustomWeightedEdge e : eList){
		if(e.isUnique()){//we only care about uniq edges
		    if(!inited){
			readHash = e.getReadHashSet();
			inited = true;
		    }else{
			HashSet<Integer> curHash = e.getReadHashSet();
			if(curHash.retainAll(readHash)){
			    if(curHash.isEmpty()
			}

			if(curHash.retainAll(readHash)){
			    if(curHash.isEmpty())
				removalList.add(new Integer(i));
			}
		    }
		}
	    }
	}
	    for(int i=removalList.size() - 1; i >= 0; i--){
	    this.paths.remove(removalList.get(i));
	}
    }
    */
    public void mergeBubble(Bubble other){
	ArrayList<Path> paths_new = new ArrayList<Path>();
	for(int i=0;i<this.paths.size();i++){
	    Path tp = this.paths.get(i);
	    for(int j=0; j<other.getPaths().size(); j++){
		Path op = other.getPaths().get(j);
		if(tp.isPhasedWith(op))
		    paths_new.add(tp.mergePaths(op));
	    }
	}
	
	if(paths_new.size() > 0){
	    this.start.addAll(other.getStart());
	    this.end.addAll(other.getEnd());
	    this.sNodes.addAll(other.getSNodes());
	    this.tNodes.addAll(other.getTNodes());
	}
    }
    
    
}

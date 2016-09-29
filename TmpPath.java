import java.util.ArrayList;

import org.jgrapht.*;
import org.jgrapht.graph.*;


public class TmpPath{

    private ArrayList<Node> orderedNodeList;

    public TmpPath(){
	this.orderedNodeList = new ArrayList<Node>();
    }
    
    public TmpPath(ArrayList<Node> l){
	this.orderedNodeList = l;
    }

    public TmpPath clone(){
	return new TmpPath(new ArrayList<Node>(this.orderedNodeList));
    }
    
    public boolean merge(TmpPath p, SimpleDirectedWeightedGraph<Node, CustomWeightedEdge> g){
	Node thisTail = this.orderedNodeList.get(this.orderedNodeList.size()-1);
	Node pHead = p.getNthNode(0);
	
	if(g.getEdge(thisTail, pHead) !=null){
	    this.orderedNodeList.addAll(p.getOrderedNodeList());
	    return true;
	}
	return false;
    }

    public ArrayList<Node> getOrderedNodeList(){
	return this.orderedNodeList;
    }

    public Path toPath(SimpleDirectedWeightedGraph<Node, CustomWeightedEdge> g){
	Path p = new Path();
	Node cur = null;
	Node pre = this.orderedNodeList.get(0);
	for(int i=1;i<this.orderedNodeList.size(); i++){
	    cur = this.orderedNodeList.get(i);
	    if(pre == null || cur == null)
		System.err.println("TmpPath.toPath():\t PRE"+ (pre != null ? pre.toString() : null)+ "\tCUR" + (cur != null ? cur.toString() : null ));
	    CustomWeightedEdge e = g.getEdge(pre, cur);
	    if(e == null){
		System.err.println("TmpPath.toPath():\t no EDGE found between\tSIZE:["+this.orderedNodeList.size()+"]");
		this.print();
	    }else
		p.appendEdge(e);
	    pre = cur;
	}
	return p;
    }

    public void appendNode(Node n){
	this.orderedNodeList.add(n);
    }

    public void print(){
	if(this.orderedNodeList.size() > 0)
	    System.out.print("NOT CONTIGUOUS --> ( " + this.orderedNodeList.get(0).toString());
	for(int i=1; i<this.orderedNodeList.size(); i++){
	    System.out.print("," + this.orderedNodeList.get(i).toString());
	}
	System.out.println(")");
    }
    
    //returns number of orderedNodeList in the path
    public int numNodes(){
	return this.orderedNodeList.size();
    }
    
    public Node getNthNode(int n){
	return this.orderedNodeList.get(n);
    }

    public int numEdges(){
	return this.orderedNodeList.size() - 1;
    }
    
}

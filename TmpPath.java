import java.util.ArrayList;

import org.jgrapht.*;
import org.jgrapht.graph.*;


public class TmpPath{

    private ArrayList<Node> orderedNodeList;

    public TmpPath(){
	this.orderedNodeList = new ArrayList<Node>();
    }

    public Path toPath(SimpleDirectedWeightedGraph<Node, CustomWeightedEdge> g){
	Path p = new Path();
	Node cur = null;
	Node pre = this.orderedNodeList.get(0);
	for(int i=1;i<this.orderedNodeList.size(); i++){
	    cur = this.orderedNodeList.get(i);
	    p.appendEdge(g.getEdge(pre, cur));
	    pre = cur;
	}
	return p;
    }

    public void appendNode(Node n){
	this.orderedNodeList.add(n);
    }
}

import java.util.ArrayList;
import java.util.HashMap;

import org.jgrapht.*;
import org.jgrapht.graph.*;

public class HLAGraph{

    private String HLAGeneName; //A B C ...
    private ArrayList<Sequence> alleles;
    private SimpleDirectedWeightedGraph<Node, DefaultWeightedEdge> g;
    private ArrayList<HashMap<Character, Node>> nodeHashList;//list index = columnIndex.

    public HLAGraph(ArrayList<Sequence> seqs){
	this.alleles = seqs;
	g = new SimpleDirectedWeightedGraph<Node, DefaultWeightedEdge>(DefaultWeightedEdge.class);
	this.nodeHashList = new ArrayList<HashMap<Character, Node>>();
    }

    public void addWeight(int basePosition, ){
	
    }
    
    private void buildGraph(){
	int numAlleles = this.alleles.size();
	Sequence firstAllele = this.alleles.get(0);
	
	/* for each alleles*/
	Node sNode = new Node('s', 0);
	Node tNode = new Node('t', this.alleles.get(0).getColLength() + 1);
	this.g.addVertex(sNode);
	this.g.addVertex(tNode);
	for(int i=0; i<numAlleles; i++){
	    Sequence curSeq = this.alleles.get(i);
	    
	    /* for each base in allele */
	    Node prevNode = sNode;
	    for(int j=0; j<curSeq.getColLength(); j++){
		if(i==0)
		    this.nodeHashList.add(new HashMap<Character, Node>());
		
		Node tmpNode = new Node(curSeq.baseAt(j));
		Character curChar = tmpNode.getBaseObj();
		HashMap<Character, Node> curHash = nodeHashList.get(j);
		//if we have not added this node
		if(curHash.get(curChar) == null){
		    this.g.addVertex(tmpNode);
		    curHash.put(curChar, tmpNode);
		}else//retrieve node if already added in the graph.
		    tmpNode = curHash.get(curChar);
		
		//add an edge
		DefaultWeightedEdge e;
		if(!this.g.containsEdge(prevNode, tmpNode)){
		    e = this.g.addEdge(prevNode,tmpNode);
		    if(prevNode.equals(sNode))
			this.g.setEdgeWeight(e, Double.MAX_VALUE);
		    else
			this.g.setEdgeWeight(e, 0.0d);
		}
		prevNode = tmpNode;
	    }
	    //add edge 
	    if(!this.g.containsEdge(prevNode, tNode))
		this.g.setEdgeWeight(this.g.addEdge(prevNode, tNode), Double.MAX_VALUE);
	}
    }
    
}



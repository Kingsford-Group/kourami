import org.jgrapht.*;
import org.jgrapht.graph.*;

public class HLAGraph{
    
    private ArrayList<Sequence> alleles;
    private SimpleDirectedWeightedGraph<Node, DefaultWeightedEdge> g;
    private ArrayList<HashMap<Character, Node>> nodeHashList;//list index = columnIndex.

    public HLAGraph(ArrayList<Sequence> seqs){
	this.alleles = seqs;
	g = new SimpleDirectedWeightedGraph<Node, DefaultWeightedEdge>();
	this.nodeHashList = new ArrayList<HashMap<Character, Node>>();
    }
    
    private void buildGraph(){
	int numAlleles = this.allelels.size();
	Sequence firstAllels = this.allels.get(0);
	//HashMap<Character, Base> hash = null;
	
	
	/* init flags*/
	//ArrayList<boolean[]> nodeflags = new ArrayList<boolean[]>();
	for(int i=0; i<firstAllele.getColLength(); i++){
	    nodeflags.add(new boolean[5]);
	}
	
	/* for each alleles*/
	Node sNode = new Node('s', 0);
	Node tNode = new Node('t', curSeq.getColLength()+1);
	this.g.addVertext(startNode);
	for(int i=0; i<numAlleles; i++){
	    Sequence curSeq = this.allels.get(i);
	    
	    /* for each base in allele */
	    Node prevNode = sNode;
	    for(int j=0; j<curSeq.getColLength(); j++){
		if(i==0)
		    this.nodeHashList.add(new HashMap<Character, Node>());
		
		Node tmpNode = new Node(curSeq.baseAt(j));
		Character curChar = tmpNode.getBaseObj();
		HashMap curHash = nodeHashList.get(j);
		//if we have not added this node
		if(curHash.get(curChar) == null){
		    this.g.addVertex(tmpNode);
		    curHash.put(curChar, tmpNode);
		}else//retrieve node if already added in the graph.
		    tmpNode = curHash.get(curChar);
		
		/* NEED TO CHECK IF EDGE HAS BEEN ADDED OR NOT */
		//add an edge
		DefaultWeightedEdge e = this.g.addEdge(prevNode,tmpNode);
		if(prevNode.equals(sNode))
		    this.g.setWeight(e, Double.MAX_VALUE);
		else
		    this.gsetWeight(e, 0.0d);
		prevNode = tmpNode;
	    }
	    //add edge 
	    this.g.setWeight(this.g.addEdge(prevNode, tNode), Double.MAX_VALUE);
	}
    }
    
    
}



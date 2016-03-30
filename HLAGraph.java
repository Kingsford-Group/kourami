import htsjdk.samtools.SAMRecord;

import java.util.ArrayList;
import java.util.HashMap;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;

import org.jgrapht.*;
import org.jgrapht.graph.*;

public class HLAGraph{

    private String HLAGeneName; //A B C ...
    private ArrayList<Sequence> alleles;
    private HashMap<String, Sequence> alleleHash;
    
    private SimpleDirectedWeightedGraph<Node, DefaultWeightedEdge> g;
    private ArrayList<HashMap<Character, Node>> nodeHashList;//list index = columnIndex-1.
    
    /* Outer list index = columnIndex -1 --> insertion point */
    /* Inner list index insertion length */ 
    private ArrayList<ArrayList<HashMap<Character, Node>>> insertionNodeHashList;
	

    public HLAGraph(ArrayList<Sequence> seqs){
	this.alleles = seqs; 
	this.alleleHash = new HashMap<String, Sequence>();
	for(int i=0;i<this.alleles.size();i++){
	    this.alleleHash.put(this.alleles.get(i).getAlleleName(), this.alleles.get(i));
	}
	g = new SimpleDirectedWeightedGraph<Node, DefaultWeightedEdge>(DefaultWeightedEdge.class);
	this.nodeHashList = new ArrayList<HashMap<Character, Node>>();
	this.insertionNodeHashList = new ArrayList<ArrayList<HashMap<Character, Node>>>();
    }
    
    private void addMissingNode(char b, int colPos, Node cur, Node pre){
	cur = new Node(b, colPos);
	this.g.addVertex(cur);
	this.nodeHashList.get(colPos - 1).put(new Character(b), cur);
	DefaultWeightedEdge e = this.g.addEdge(pre, cur);
	this.g.setEdgeWeight(e, 1.0d);
    }

    private void incrementWeight(Node source, Node target){
	DefaultWeightedEdge e = g.getEdge(source, target);
	if(e == null){
	    g.addEdge(source, target);
	    g.setEdgeWeight(e, 1.0d);
	}else
	    g.setEdgeWeight(e, g.getEdgeWeight(e)+1);
    }
    
    
    public void addWeight(SAMRecord sr){
	Cigar cigar = sr.getCigar();
	byte[] bases = sr.getReadBases(); //ASCII bytes ACGTN=.
	int baseIndex = 1;
	int refBasePos = sr.getAlignmentStart();
	Node prevnode = null;
	Node curnode = null;
	Sequence curAllele = this.alleleHash.get(sr.getReferenceName());
	int colPos = curAllele.getColPosFromBasePos(refBasePos);
	
	if(cigar==null) return;
	for(final CigarElement ce : cigar.getCigarElements()){
	    CigarOperator op = ce.getOperator();
	    int cigarLen = ce.getLength();
	    switch(op)
		{
		case S :
		    baseIndex += cigarLen;
		case H :
		    ;
		case M :
		    {
			//colPos = curAllele.getColPosFromBasePos(refBasePos);
			for(int i=0; i<cigarLen; i++){
			    //int colPos = curAllele.getColPosFromBasePos(refBasePos);
			    curnode = this.nodeHashList.get(colPos -1).get(new Character((char)bases[baseIndex]));
			    //if no such node found, we add new node and add adge from prevnode.
			    if(curnode == null)
				this.addMissingNode((char)bases[baseIndex], colPos, curnode, prevnode);
			    else
				this.incrementWeight(prevnode, curnode);//source, target);
			    
			    prevnode=curnode;
			    baseIndex++;
			    refBasePos++;
			    colPos++;
			}
		    }
		case D :
		    {
			for(int i=0; i<cigarLen; i++){
			    curnode = this.nodeHashList.get(colPos - 1).get(new Character('-'));
			    if(curnode == null)
				this.addMissingNode('-', colPos, curnode, prevnode);
			    else
				this.incrementWeight(prevnode, curnode);
			    prevnode=curnode;
			    //baseIndex++;
			    refBasePos++;
			    colPos++;
			}
		    }
		case I :
		    {
			for(int i=0; i<cigarLen; i++){
			    
			    if(this.insertionNodeHashList.get(colPos - 1).size() > i){
				curnode = this.insertionNodeHashList.get(colPos - 1).get(i).get(new Character((char)bases[baseIndex]));
				
			    }else{//we need to add extra position (insertion length)
				
				this.insertionNodeHashList.get(colPos - 1).add(new HashMap<Character, Node>());
				curnode = null;
			    }
			    
			    if(curnode == null){
				curnode = new Node((char)bases[baseIndex], colPos);
				this.g.addVertex(curnode);
				this.insertionNodeHashList.get(colPos - 1).get(i).put(new Character((char)bases[baseIndex]), curnode);
				DefaultWeightedEdge e = this.g.addEdge(prevnode, curnode);
				this.g.setEdgeWeight(e, 1.0d);
			    }else
				this.incrementWeight(prevnode, curnode);
			    
			    prevnode = curnode;
			    baseIndex++;
			}
		    }
		}
	    
	}
	
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
		if(i==0){
		    this.nodeHashList.add(new HashMap<Character, Node>());
		    this.insertionNodeHashList.add(new ArrayList<HashMap<Character, Node>>());
		}
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



import org.jgrapht.*;
import org.jgrapht.graph.*;

public class Test{
    
    public static void main(String[] args){
	
	SimpleDirectedWeightedGraph<Node, CustomWeighteEdge> g = new SimpleDirectedWeightedGraph<Node, CustomWeightedEdge>(CustomWeightedEdge.class); 
	Node sNode = new Node('s', 0);
	Node tNode = new Node('t', 5);
	
	g.addVertex(sNode);
	g.addVertex(tNode);

	
	Node n1 = new Node('a', 1);
	Node n2 = new Node('t', 1);
	
	g.addVertex(n1);
	g.addVertex(n2);
	g.addEdge(sNode, n1);
	g.addEdge(sNode, n2);

	Node n3 = new Node('c', 2);
	Node n4 = new Node('g', 2);
	
	g.addVertex(n3);
	g.addVertex(n4);
	g.addEdge(n1, n3);
	g.addEdge(n2, n4);
	
	

    }

}

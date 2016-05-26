public class Path{

    private ArrayList<CustomWeightedEdge> orderedEdgeList;
    
    public Path deepCopy(){
	Path p = new Path();
	for(CustomWeightedEdge e : this.orderedEdgeList){
	    p.appendEdge(e);
	}
	return p;
    }

    public ArrayList<CustomWeightedEdge> getOrderedEdgeList{
	return this.orderedEdgeList;
    }

    public Path(){
	this.orderedEdgeList = new ArrayList<CustomWeightedEdge>();
    }

    public void appendEdge(CustomWeightedEdge e){
	this.orderedEdgeList.add(e);
    }

    public Path(CustomeWeightedEdge e){
	this();
	this.appendEdge(e);
    }

    public Node getLastVertex(DefaultDirectedWeightedGraph g){
	if(this.orderedEdgeList == null || this.orderedEdgeList.size() == 0)
	    return null;
	else{
	    return g.getEdgeTarget(this.orderedEdgeList.get(this.orderedEdgeList.size() - 1));
	}
    }
    
    public void mergePaths(Path other){
	this.orderedEdgeList.addAll(other.getPathArray());
    }
    
    public boolean mergePaths(Path other){
	if(this.isPhasedWith(other)){
	    return true;
	}
	return false;
    }

    public boolean isPhasedWith(Path other){
	if(this.intersectUniquEdges(other)){
	    return true;
	}
	return false;
    }
    
}

import java.io.*;

public class Bubble{

    private HLAGraph g;
    private int start;
    private int end;
    
    private Node s;
    private Node t;
    
    private ArrayList<Path> paths;

    public void Bubble(HLAGraph hg, Node s, Node t){
	this.g = hg;
	this.s = s;
	this.t = t;
	this.start = s.getColIndex();
	this.end = e.getColIndex();
	this.paths = new ArrayList<Path>();
    }
    
    
    public ArrayList<Path> decompose(){
	this.paths =  this.g.findAllSTPath(s, t);
	return this.paths;
    }
    
    public void removeUnsupported(){
	HashSet readHash;
	ArrayList<Integer> removalList = new ArrayList<Integer>();
	for(int i=0; i<this.paths.size(); i++){//Path p: this.paths){
	    Path p = this.paths.get(i);
	    ArrayList<CustomWeightedEdge> eList;
	    boolean inited = false;
	    for(CustomWeightedEdge e : list){
		if(e.isUnique()){//we only care about uniq edges
		    if(!inited){
			readHash = e.getReadHash();
		    }else{
			ReadHash curHash = e.getReadHash();
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

    public void mergeBubble(Bubble other){
	ArrayList<Path> paths_new = new ArrayList<Path>();
	for(int i=0;i<this.paths.size();i++){
	    for(int j=0; j<other.getPaths().size(); j++)
		if(this.paths.get(i).mergePaths(other.getPaths().get(j))){
		    paths_new.add();
		}
	    }
	}
    }
}

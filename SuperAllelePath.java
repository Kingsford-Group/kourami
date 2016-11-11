import java.util.ArrayList;

public class SuperAllelePath{

    private ArrayList<AllelePath> orderedAllelePaths;
    private ArrayList<Integer> pathNums;
    private String hlaGeneName;
    
    public SuperAllelePath(String genename){
	this.orderedAllelePaths = new ArrayList<AllelePath>();
	this.pathNums = new ArrayList<Integer>();
	this.hlaGeneName = genename;
    }
    
    public SuperAllelePath(ArrayList<AllelePath> aps, ArrayList<Integer> pns, String hgn){
	this.orderedAllelePaths = aps;
	this.pathNums = pns;
	this.hlaGeneName = hgn;
    }

    public int numAllelePaths(){
	return this.orderedAllelePaths.size();
    }

    public double[] getJointProbability(SuperAllelePath other){
	double[] jp = new double[2];
	if(this.numAllelePaths() != other.numAllelePaths()){
	    System.err.println("Incompatible SuperAllelePath. The number of fractured allelePath in superPath does not match");
	    return jp;
	}
	
	for(int i=0; i<this.orderedAllelePaths.size(); i++){
	    double[] apjp = this.orderedAllelePaths.get(i).getJointProbability(other.getOrderedAllelePaths().get(i));
	    jp[0] += apjp[0];
	    jp[1] += apjp[1];
	}
	return jp;
    }
    
    public ArrayList<AllelePath> getOrderedAllelePaths(){
	return this.orderedAllelePaths;
    }


    public String getHlaGeneName(){
	return this.hlaGeneName;
    }

    public SuperAllelePath clone(){
	return new SuperAllelePath(new ArrayList<AllelePath>(this.orderedAllelePaths)
				   , new ArrayList<Integer>(this.pathNums)
				   , this.hlaGeneName);
    }

    public void addAllelePath(AllelePath ap, int n){
	this.orderedAllelePaths.add(ap);
	this.pathNums.add(new Integer(n));
    }
    
    public String pathnums2String(){
	StringBuffer bf = new StringBuffer();
	int count = 0;
	for(Integer i : this.pathNums){
	    if(count > 0)
		bf.append(":");
	    bf.append(i.intValue());
	    count++;
	}
	return bf.toString();
    }


    public StringBuffer getSequenceBuffer(){
	StringBuffer bf = new StringBuffer();
	for(AllelePath ap : this.orderedAllelePaths)
	    bf.append(ap.getSequence());
	return bf;
    }

    public String toSimpleString(){
	return this.pathnums2String() + "\t" + this.getWeightedIntersectionSumScore() + "\t" + this.getProbability();
    }
    
    public StringBuffer toFasta(){
	StringBuffer bf = new StringBuffer(">" + this.hlaGeneName + "_" + pathnums2String() + "\t" 
					   + this.getWeightedIntersectionSumScore() + "\t" + this.getProbability() + "\n");
	bf.append(this.getSequenceBuffer());
	bf.append("\n");
	return bf;
    }

    public double getProbability(){
	double logP = 0.0d;
	for(AllelePath ap : this.orderedAllelePaths){
	    logP += ap.getProbability();
	}
	return logP;
    }

    /*
    public double getProbability(){
	double p = 1.0d;
	for(AllelePath ap : this.orderedAllelePaths)
	    p = p * ap.getProbability();
	return p;
    }
    */

    public double getWeightedIntersectionSumScore(){
	double s = 0.0d;
	int n = 0;
	for(AllelePath ap : this.orderedAllelePaths){
	    s = s + ap.getIntersectionSum();
	    n = n + ap.getMergedNums();
	}
	return s / n;
    }
    
}

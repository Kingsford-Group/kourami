/*
Part of Kourami HLA typer/assembler
(c) 2017 by  Heewook Lee, Carl Kingsford, and Carnegie Mellon University.
See LICENSE for licensing.
*/
import it.unimi.dsi.fastutil.ints.IntArrayList;

public class BubblePathLikelihoodScores{

    private double[][] logFractionScores; // this is just based on the number of reads and the total numeber of reads for this bubble.

    private double[][] logScores; //only the upper triangle is used
    //private int[] maxHomoLogScoresIndicies;
    private int maxHomoGenotypeIndex;
    private int maxHeteroGenotypeIndex1;
    private int maxHeteroGenotypeIndex2;
    
    private int doubleCountH1; //2 x #reads assigned to H1
    private int doubleCountH2; //2 x #reads assigned to H2

    // must be i<=j (upper triangle only)
    public double getLogScore(int i, int j){
	return this.logScores[i][j];
    }
    
    // must be i<=j (upper triangle only)
    public double getLogFractionScore(int i, int j){
	return this.logFractionScores[i][j];
    }

    public double[][] getLogScores(){
	return this.logScores;
    }

    public double[][] getLogFractionScores(){
	return this.logFractionScores;
    }
    
    public int getMaxHomoGenotypeIndex(){
	return this.maxHomoGenotypeIndex;
    }
    
    public double getMaxHomoScore(){
	return this.logScores[this.maxHomoGenotypeIndex][this.maxHomoGenotypeIndex];
    }

    public int[] getMaxHeteroIndex(){
	int[] tmp = new int[2];
	tmp[0] = this.maxHeteroGenotypeIndex1;
	tmp[1] = this.maxHeteroGenotypeIndex2;
	return tmp;
    }
    
    public int getMaxHeteroGenotypeIndex1(){
	return this.maxHeteroGenotypeIndex1;
    }
    
    public int getMaxHeteroGenotypeIndex2(){
	return this.maxHeteroGenotypeIndex2;
    }

    public double getMaxHeteroScore(){
	return this.logScores[this.maxHeteroGenotypeIndex1][this.maxHeteroGenotypeIndex2];
    }

    public int getDoubleCountH1(){
	return this.doubleCountH1;
    }
    
    public int getDoubleCountH2(){
	return this.doubleCountH2;
    }

    public BubblePathLikelihoodScores(int numPaths){
	this.logFractionScores = new double[numPaths][numPaths];
	this.logScores = new double[numPaths][numPaths]; 
	for(int i=0;i<this.logScores.length;i++){
	    for(int j=i;j<this.logScores[i].length;j++)
		this.logScores[i][j] = Double.NEGATIVE_INFINITY;
	}
	//this.logScores[0][0] = Double.NEGATIVE_INFINITY;
	//this.logScores[0][1] = Double.NEGATIVE_INFINITY;
	this.maxHomoGenotypeIndex = 0;
	
	this.maxHeteroGenotypeIndex1 = 0;
	if(numPaths > 1)
	    this.maxHeteroGenotypeIndex2 = 1;
	else
	    this.maxHeteroGenotypeIndex2 = 0;
    }
    
    // since we are using UPPER triangle only  j >= i
    public void updateMax(int i, int j, double logScore, double readFractionScore, int dcH1, int dcH2){
	if(j<i){
	    HLA.log.appendln("INVLAID [i][j] pairing. j is smaller than i\nSystem exiting.");
	    System.exit(-1);
	}
	if(HLA.DEBUG){
	    HLA.log.appendln("Attempting to updateMax");
	    HLA.log.appendln("curMaxHetero [" + this.maxHeteroGenotypeIndex1 + "][" + this.maxHeteroGenotypeIndex2 + "]:\t" + this.logScores[this.maxHeteroGenotypeIndex1][this.maxHeteroGenotypeIndex2]);
	    HLA.log.appendln("curMaxHomo [" + this.maxHomoGenotypeIndex + "][" + this.maxHomoGenotypeIndex + "]:\t" + this.logScores[this.maxHomoGenotypeIndex][this.maxHomoGenotypeIndex]);
	}
	if(i != j){/* Heterozygous */
	    if(HLA.DEBUG)
		HLA.log.appendln("[HETERO]");

	    if(logScore > this.logScores[this.maxHeteroGenotypeIndex1][this.maxHeteroGenotypeIndex2]){
		this.maxHeteroGenotypeIndex1 = i;
		this.maxHeteroGenotypeIndex2 = j;
		this.doubleCountH1 = dcH1;
		this.doubleCountH2 = dcH2;
		
		if(HLA.DEBUG){
		    HLA.log.appendln("Updating Max Hetero Counts: curMax:" + this.logScores[this.maxHeteroGenotypeIndex1][this.maxHeteroGenotypeIndex2]);
		    HLA.log.appendln("\t|H1|x2=" + doubleCountH1 + "\t|H2|x2=" + doubleCountH2);
		}
	    }
	}else{/* Homozygous */
	    if(HLA.DEBUG)
		HLA.log.appendln("[HOMO]");
	    if(logScore > this.logScores[this.maxHomoGenotypeIndex][this.maxHomoGenotypeIndex])
		this.maxHomoGenotypeIndex = i;
	    
	}
	this.logScores[i][j] = logScore;
	this.logFractionScores[i][j] = readFractionScore;
	
    }
    
    public void applyRemoval(IntArrayList sortedRemovalList){
	for(int i=sortedRemovalList.size() -1; i>=0; i--)
	    this.remove(sortedRemovalList.getInt(i));
    }

    //removes a path so we update the logScores and all indicies
    private void remove(int n){
	double[][] newLogScores = new double[this.logScores.length][this.logScores.length];
	for(int i=0;i<n;i++){
	    for(int j=i;j<n;j++)
		newLogScores[i][j] = this.logScores[i][j];
	    
	    for(int j=n+1; j<this.logScores.length; j++)
		newLogScores[i][j-1] = this.logScores[i][j];
	    
	}
	
	for(int i=n+1; i<this.logScores.length; i++){
	    for(int j=i; j<this.logScores.length; j++)
		newLogScores[i-1][j-1] = this.logScores[i][j];
	}
	this.logScores = newLogScores;
	
	if(maxHomoGenotypeIndex > n)
	    this.maxHomoGenotypeIndex = this.maxHomoGenotypeIndex - 1;

	if(maxHeteroGenotypeIndex1 > n)
	    this.maxHeteroGenotypeIndex1 = this.maxHeteroGenotypeIndex1 - 1;
	
	if(maxHeteroGenotypeIndex2 > n)
	    this.maxHeteroGenotypeIndex2 = this.maxHeteroGenotypeIndex2 - 1;
	
    }
    
    
    public void setHeteroIndicies(int i, int j){
	this.maxHeteroGenotypeIndex1 = i;
	this.maxHeteroGenotypeIndex2 = j;
    }
}

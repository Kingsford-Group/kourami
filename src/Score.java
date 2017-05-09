/*
Part of Kourami HLA typer/assembler
(c) 2017 by  Heewook Lee, Carl Kingsford, and Carnegie Mellon University.
See LICENSE for licensing.
*/
public class Score implements Comparable<Score>{

    private double[] scores; //scores for each pair of alleles (n+1)*n/2 pairings
    private int[] pairIndicies; //size 2 int array holding i,j pair for alleles index
    public static int sortIndex;

    public Score(double[] s, int i , int j){
	this.scores = s;
	this.pairIndicies = new int[2];
	this.pairIndicies[0] = i;
	this.pairIndicies[1] = j;
    }

    public Score(double[] s, int[] p){
	this.scores = s;
	this.pairIndicies = p;
    }

    //descending order
    public int compareTo(Score os){
	if(this.getNthScore(Score.sortIndex) > os.getNthScore(Score.sortIndex))
	    return -1;
	else if(this.getNthScore(Score.sortIndex) < os.getNthScore(Score.sortIndex))
	    return 1;
	else
	    return 0;
    }

    public double getNthScore(int scoringScheme){
	return this.scores[scoringScheme];
    }
    
    public int[] getIndicies(){
	return this.pairIndicies;
    }
}

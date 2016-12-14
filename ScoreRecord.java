import java.util.*;

public class ScoreRecord{

    public ArrayList<Score> listOfScores;

    public ScoreRecord(){
	this.listOfScores = new ArrayList<Score>();
    }
    
    public void sort(int sIndex){
	Score.sortIndex = sIndex;
	Collections.sort(listOfScores);
    }
    
    public void addScore(double[] s, int i, int j){
	this.listOfScores.add(new Score(s, i, j));
    }
    
    public void printBest(int sortIndex){
	this.sort(sortIndex);
	double best = Double.NEGATIVE_INFINITY;
	for(Score s : this.listOfScores){
	    if(s.getNthScore(sortIndex) >= best){
		best = s.getNthScore(sortIndex);
		int[] bestIJ = s.getIndicies();
		HLA.log.appendln(">>>>>>>>> BEST PAIR[" + bestIJ[0] + ":" + bestIJ[1] + "]:\t" + s.getNthScore(sortIndex));
	    }else
		break;
	}
    }

    public ArrayList<int[]> getBestPairs(int sortIndex){
	ArrayList<int[]> bestPairs = new ArrayList<int[]>();
	this.sort(sortIndex);
	double best = Double.NEGATIVE_INFINITY;
	for(Score s : this.listOfScores){
	    if(s.getNthScore(sortIndex) >= best){
		best = s.getNthScore(sortIndex);
		int[] bestIJ = s.getIndicies();
		bestPairs.add(bestIJ);
		HLA.log.appendln(">>>>>>>>> BEST PAIR[" + bestIJ[0] + ":" + bestIJ[1] + "]:\t" + s.getNthScore(sortIndex));
	    }else
		break;
	}
	return bestPairs;
    }
}


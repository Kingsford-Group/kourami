import java.io.*;
import java.util.*;

public class SolutionChecker{
    private HashMap<String, Solution> hlaSolutions; // key is geneString
    private HashMap<String, String> hlaName2Typing; 
    
    public SolutionChecker(String solFile, String namesFile){
	this.hlaSolutions = Solution.loadSolutions(solFile);
	this.hlaName2Typing = SolutionChecker.loadNames(namesFile);
    }
    
    public static void main(String[] args){
	//args[0] : solution file
	//args[1] : experiment directory name
	//args[2] : names file
	new SolutionChecker(args[0], args[2]).checkSolutions(args[1]);
    }
    
    public static HashMap<String, String> loadNames(String namesF){
	HashMap<String, String> hash = new HashMap<String, String>();
	BufferedReader br = null;
	try{
	    br = new BufferedReader(new FileReader(namesF));
	    String curline = "";
	    while((curline=br.readLine())!=null){
		String[] tokens = curline.substring(1).split("\\s+");
		String hlaname = tokens[0];
		String type = tokens[1];
		//System.out.println("Putting [" + hlaname+"],["+ type +"]");
		//this.hlaName2Typing.put(hlaname, type);
		hash.put(hlaname, type);
	    }
	    br.close();
	}catch(IOException ioe){
	    ioe.printStackTrace();
	}
	return hash;
    }

    public void checkSolutions(String expDir){
	File ed = new File(expDir);
	File[] files = ed.listFiles();
	for(File f : files){
	    String curname = f.getName();
	    if(curname.matches("NA[0-9]+\\.salmon_k[0-9][0-9]_[0-9][0-9]")){
		System.out.print(curname + "\t");
		String sampleID = curname.substring(0,curname.indexOf("."));
		HashMap<String, SortedRecords> hash = new CollectResult(this.hlaName2Typing, ed + File.separator + curname + File.separator + "quant.sf").getSortedRecordsHash();
		this.printMatchScore(hash, sampleID);
	    }
	}
    }
    
    

    public String index2Gene(int i){
	if(i == 0)
	    return "A";
	else if(i == 1)
	    return "B";
	else if(i == 2)
	    return "C";
	else if(i == 3)
	    return "DQA1";
	else if(i == 4)
	    return "DQB1";
	else if(i == 5)
	    return "DRB1";
	else
	    return null;
    }
    
    public int gene2Index(String g){
	if(g.equals("A"))
	    return 0;
	else if(g.equals("B"))
	    return 1;
	else if(g.equals("C"))
	    return 2;
	else if(g.equals("DQA1"))
	    return 3;
	else if(g.equals("DQB1"))
	    return 4;
	else if(g.equals("DRB1"))
	    return 5;
	else
	    return -1;
    }

    public void printMatchScore(HashMap<String, SortedRecords> hash, String sampleID){
	System.out.println(sampleID + "\t" + hlaSolutions.get("A").matchScore(hash.get("A"), sampleID)
			   + "\t" + hlaSolutions.get("B").matchScore(hash.get("B"), sampleID)
			   + "\t" + hlaSolutions.get("C").matchScore(hash.get("C"), sampleID)
			   + "\t" + hlaSolutions.get("DQA1").matchScore(hash.get("DQA1"), sampleID)
			   + "\t" + hlaSolutions.get("DQB1").matchScore(hash.get("DQB1"), sampleID)
			   + "\t" + hlaSolutions.get("DRB1").matchScore(hash.get("DRB1"), sampleID)
			   );
    }

    
}

class Solution{

    private HashMap<String, String[]> sample2Solutions; //key is sampleID, value --> [0] : allele1, [1] : allele2
    private String hlaGene;

    public Solution(String hlaG){
	this.hlaGene = hlaG;
	this.sample2Solutions = new HashMap<String, String[]>();
    }

    public void add(String sampleID, String a1, String a2){
	String[] tmp = new String[2];
	tmp[0] = a1;
	tmp[1] = a2;
	this.sample2Solutions.put(sampleID, tmp);
    }

    public void printMatchScore(HashMap<String, SortedRecords> hash, String sampleID){
	System.out.println(sampleID + "\t" + matchScore(hash.get("A"), sampleID) 
			   + "\t" + matchScore(hash.get("B"), sampleID)
			   + "\t" + matchScore(hash.get("C"), sampleID)
			   + "\t" + matchScore(hash.get("DQA1"), sampleID)
			   + "\t" + matchScore(hash.get("DQB1"), sampleID)
			   + "\t" + matchScore(hash.get("DRB1"), sampleID));
	
    }

    public int matchScore(SortedRecords sr, String sampleID){

	String[] alleles = this.sample2Solutions.get(sampleID);
	String t1 = sr.getNthBestUptoXDigit(0,4);
	String t2 = null;
	if(sr.size() > 1)
	    t2 = sr.getNthBestUptoXDigit(1,4);
	else
	    t2 = t1;

	System.out.println(this.hlaGene + "\t" + alleles[0] + "\t" + alleles[1] + "|" + t1 + "\t" + t2);
	
	if(t1.equals(alleles[0]) && t2.equals(alleles[1]))
	    return 2;
	else if(t2.equals(alleles[0]) && t1.equals(alleles[1]))
	    return 2;
	else if(t1.equals(alleles[0]) || t2.equals(alleles[1]))
	    return 1;
	else if(t2.equals(alleles[0]) || t1.equals(alleles[1]))
	    return 1;
	else
	    return 0;
    }

    public static HashMap<String, Solution> loadSolutions(String solutionFile){
	HashMap<String, Solution> sols = new HashMap<String, Solution>();
	BufferedReader br = null;
	try{
	    br = new BufferedReader(new FileReader(solutionFile));
	    String curline = null;
	    ArrayList<String> keys = new ArrayList<String>(); //geneStr
	    while((curline=br.readLine())!=null){
		if(curline.startsWith("#Ind")){
		    String[] tokens = curline.split("\\t+");
		    for(int i=1; i<tokens.length; i++){
			keys.add(tokens[i]);
			sols.put(tokens[i], new Solution(tokens[i]));
		    }
		}else{
		    String[] tokens = curline.split("\\t");
		    String sampleID = tokens[0];
		    for(int i=1; i<tokens.length;i=i+2){
			sols.get(keys.get(i/2)).add(sampleID, tokens[i], tokens[i+1]);
		    }
		}
	    }
	}catch(IOException ioe){
	    ioe.printStackTrace();
	}
	return sols;
    }

}

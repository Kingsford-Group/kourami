import java.io.*;
import java.util.*;

public class CollectResult{
    
    private HashMap<String, String> hlaName2Typing;
    
    private HashMap<String, SortedRecords> hlaLocus2Alleles;

    //private HashMap<String, Solution> hlaSolutions;//key is geneString --> 'A', 'B', ...
    
    /*public static void main(String[] args){
	CollectResult obj = new CollectResult();
	//obj.loadNames(args[0]); //names_file
	//obj.loadSolutions(args[1]); //solution_file
	obj.processQuant(args[1]); //quant file
	obj.output();
	obj.outToFile(args[2]);
    }
    */    
    public SortedRecords getBestTwoFromHLAType(String typeStr){
	//System.out.println("TypeStr\t" + typeStr + "\tKey is:" + typeStr.substring(0, typeStr.indexOf("*")));
	
	return this.hlaLocus2Alleles.get(typeStr.substring(0, typeStr.indexOf("*")));
    }
    
    public HashMap<String, SortedRecords> getSortedRecordsHash(){
	return this.hlaLocus2Alleles;
    }
    /*
    public double getAccuracy(String curSampleID){
	int scores[] = new int[6];
	//for each gene
	for(int i=0;i<scores.length; i++){
	    Solution s = this.hlaSolutions.get(this.index2Gene(i));
	    scores[i] = s.matchScore(hlaLocus2Alleles.get(this.index2Gene(i)), curSampleID);
	}
	}*/
    
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
    
    public void output(){
	//Set<String> keys = this.hlaLocus2Alleles.keySet();
	//for(String s : keys){
	    
	//    System.out.println(s + "\t" + this.hlaLocus2Alleles.get(s).extractUptoXdigit(4));
	//}
	System.out.println("A   " + "\t" + this.hlaLocus2Alleles.get("A").extractUptoXdigit(4));
	System.out.println("B   " + "\t" + this.hlaLocus2Alleles.get("B").extractUptoXdigit(4));
	System.out.println("C   " + "\t" + this.hlaLocus2Alleles.get("C").extractUptoXdigit(4));
	System.out.println("DQA1" + "\t" + this.hlaLocus2Alleles.get("DQA1").extractUptoXdigit(4));
	System.out.println("DQB1" + "\t" + this.hlaLocus2Alleles.get("DQB1").extractUptoXdigit(4));
	System.out.println("DRB1" + "\t" + this.hlaLocus2Alleles.get("DRB1").extractUptoXdigit(4));
    }

    
    public void outToFile(String f){
	BufferedWriter bw = null;
	try{
	    bw = new BufferedWriter(new FileWriter(f));
	    bw.write("A   " + "\t" + this.hlaLocus2Alleles.get("A").extractUptoXdigitFull(4) + "\n");
	    bw.write("B   " + "\t" + this.hlaLocus2Alleles.get("B").extractUptoXdigitFull(4) + "\n");
	    bw.write("C   " + "\t" + this.hlaLocus2Alleles.get("C").extractUptoXdigitFull(4) + "\n");
	    bw.write("DQA1" + "\t" + this.hlaLocus2Alleles.get("DQA1").extractUptoXdigitFull(4) + "\n");
	    bw.write("DQB1" + "\t" + this.hlaLocus2Alleles.get("DQB1").extractUptoXdigitFull(4) + "\n");
	    bw.write("DRB1" + "\t" + this.hlaLocus2Alleles.get("DRB1").extractUptoXdigitFull(4) + "\n");
	    bw.close();
	}catch(IOException ioe){
	    ioe.printStackTrace();
	}
    }
    
    
    public void processQuant(String qf){
	this.initAlleles();
	BufferedReader br = null;
	try{
	    br = new BufferedReader(new FileReader(qf));
	    String curline = "";
	    while((curline = br.readLine())!=null){
		if(!curline.startsWith("Name")){
		    String[] tokens = curline.split("\\t");
		    double tpm = Double.parseDouble(tokens[3]);
		    //if(tpm > 0)
			//System.err.println(this.hlaName2Typing.get(tokens[0]) + "\t"+ curline);
		    //System.out.println("Typing : " + this.hlaName2Typing.get(tokens[0]) + "\tHLA_NAME: " + tokens[0]);
		    SortedRecords curb2 = this.getBestTwoFromHLAType(this.hlaName2Typing.get(tokens[0]));
		    curb2.update(this.hlaName2Typing.get(tokens[0]), tpm);
		}
	    }
	    br.close();
	}catch(IOException ioe){
	    ioe.printStackTrace();
	}
    }
    
    

    
    
    public void initAlleles(){
	this.hlaLocus2Alleles = new HashMap<String, SortedRecords>();
	this.hlaLocus2Alleles.put("A", new SortedRecords());
	this.hlaLocus2Alleles.put("B", new SortedRecords());
	this.hlaLocus2Alleles.put("C", new SortedRecords());
	this.hlaLocus2Alleles.put("DQA1", new SortedRecords());
	this.hlaLocus2Alleles.put("DMA", new SortedRecords());
	this.hlaLocus2Alleles.put("DMB", new SortedRecords());
	this.hlaLocus2Alleles.put("DOA", new SortedRecords());
	this.hlaLocus2Alleles.put("DOB", new SortedRecords());
	this.hlaLocus2Alleles.put("DPA1", new SortedRecords());
	this.hlaLocus2Alleles.put("DPB1", new SortedRecords());
	this.hlaLocus2Alleles.put("DQA1", new SortedRecords());
	this.hlaLocus2Alleles.put("DQB1", new SortedRecords());
	this.hlaLocus2Alleles.put("DRA", new SortedRecords());
	this.hlaLocus2Alleles.put("DRB1", new SortedRecords());
	this.hlaLocus2Alleles.put("DRB3", new SortedRecords());
	this.hlaLocus2Alleles.put("DRB4", new SortedRecords());
	this.hlaLocus2Alleles.put("E", new SortedRecords());
	this.hlaLocus2Alleles.put("F", new SortedRecords());
	this.hlaLocus2Alleles.put("G", new SortedRecords());
	this.hlaLocus2Alleles.put("H", new SortedRecords());
	this.hlaLocus2Alleles.put("J", new SortedRecords());
	this.hlaLocus2Alleles.put("K", new SortedRecords());
	this.hlaLocus2Alleles.put("L", new SortedRecords());
	this.hlaLocus2Alleles.put("MICA", new SortedRecords());
	this.hlaLocus2Alleles.put("MICB", new SortedRecords());
	this.hlaLocus2Alleles.put("P", new SortedRecords());
	this.hlaLocus2Alleles.put("TAP1", new SortedRecords());
	this.hlaLocus2Alleles.put("TAP2", new SortedRecords());
	this.hlaLocus2Alleles.put("V", new SortedRecords());
    }

    public CollectResult(){
	this.hlaName2Typing = new HashMap<String, String>();
	//this.hlaSolutions = new HashMap<String, Solution>();
    }

    public CollectResult(HashMap<String, String> hlaN2T, String quant){
	this.hlaName2Typing = hlaN2T;
	this.initAlleles();
	this.processQuant(quant);
	//this.hlaSolutions = new HashMap<String, Solution>();
    }
    

}

class SortedRecords{
    
    private ArrayList<Double> tpms;
    private ArrayList<String> names;
    
    public SortedRecords(){
	this.tpms = new ArrayList<Double>();
	this.names = new ArrayList<String>();
    }

    public int size(){
	return this.names.size();
    }
    
    public String extractUptoXdigit(int n){
	//System.err.println("tpms.size() = " + this.tpms.size());
	//System.err.println("names.size() = " + this.names.size());
	StringBuffer out = new StringBuffer();
	if(this.tpms.size() > 0){
	    for(int i=0; i<tpms.size() && i<2;i++){
		if(tpms.get(i).doubleValue() > 0){
		    if(i > 0)
			out.append("\t");
		    int starPos = this.names.get(i).indexOf("*")+1;
		    out.append(this.names.get(i).substring(0,this.names.get(i).indexOf("*")) + "*" + this.names.get(i).substring(starPos, starPos+n+((n-1)/2)) + "(" + this.tpms.get(i).doubleValue()+")");
		}
	    }
	}else
	    out.append("NO ALLELE DETECTED");
	return out.toString();
    }
    
    public String extractUptoXdigitFull(int n){
	StringBuffer out = new StringBuffer();
	if(this.tpms.size() > 0){
	    for(int i=0; i<tpms.size() && i<8;i++){
		if(tpms.get(i).doubleValue() > 0){
		    if(i > 0)
			out.append("\t");
		    int starPos = this.names.get(i).indexOf("*")+1;
		    out.append(this.names.get(i).substring(0,this.names.get(i).indexOf("*")) + "*" + this.names.get(i).substring(starPos, starPos+n+((n-1)/2)) + "(" + this.tpms.get(i).doubleValue()+")" );
		}
	    }
	}else
	    out.append("NO ALLELE DETECTED");
	return out.toString();
    }
    
    public String getNthBestUptoXDigit(int n, int x){
	System.out.print("(" + this.names.get(n) + ")");
	int starPos = this.names.get(n).indexOf("*")+1;
	return this.names.get(n).substring(starPos, starPos+x+((x-1)/2));
    }


    public void update(String name, double tpm){
	if(tpm > 0){
	    boolean updated = false;
	    for(int i=0; i<this.tpms.size();i++){
		double curVal = this.tpms.get(i).doubleValue();
		
		if(curVal < tpm){
		    this.tpms.add(i, new Double(tpm));
		    this.names.add(i, name);
		    updated = true;
		    //System.err.println("updated" + "\t[" + name + ","+ tpm + "]" );
		    break;
		}
	    }
	    if(!updated){
		this.tpms.add(new Double(tpm));
		this.names.add(name);
		//System.err.println("updated*" + "\t[" + name + ","+ tpm + "]" );
		updated = true;
	    }
	}
    }

    
}




class BestTwo{
    

    public String extractUptoXdigit(int n){
	//System.out.print("-->\t" + name1 + "\t" + name2);
	int len1 = this.name1.length();
	int len2 = this.name2.length();
	StringBuffer out = new StringBuffer("-->\t");
	if(len1 > 0){
	    int starPos1 = name1.indexOf("*")+1;
	    out.append(name1.substring(0,name1.indexOf("*")) + "*" + name1.substring(starPos1, starPos1+n+((n-1)/2)));
	    if(len2 > 0){
		int starPos2 = name2.indexOf("*")+1;
		out.append("\t" + name2.substring(0,name2.indexOf("*")) + "*" + name2.substring(starPos2, starPos2+n+((n-1)/2)));
	    }
	}else
	    out.append("NO ALLELE DETECTED");
	return out.toString();
    }
    
    public BestTwo(){
	this.tpm1 = 0.0d;
	this.tpm2 = 0.0d;
	this.name1 = new StringBuffer();
	this.name2 = new StringBuffer();
    }
    
    public void update(String name, double tpm){
	if(tpm > this.tpm2){
	    if(tpm <= this.tpm1){
		this.tpm2 = tpm;
		this.name2 = new StringBuffer(name);
	    }else{
		this.tpm2 = this.tpm1;
		this.name2 = this.name1;
		this.tpm1 = tpm;
		this.name1 = new StringBuffer(name);
	    }
	}
	
    }

    public String getType1(){
	return this.name1.toString();
    }
    
    public String getType2(){
	return this.name2.toString();
    }

    double tpm1;
    double tpm2;
    StringBuffer name1; // types
    StringBuffer name2;

}


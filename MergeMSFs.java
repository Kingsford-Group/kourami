import java.io.*;

public class MergeMSFs{

    private Hashtable<String, Sequence> allele2Sequence;
    private ArrayList<String> orderedAlleles;
    private StringBuffer header;
    
    public MergedMSFs(){
	this.header = new StringBuffer();
	this.allele2Sequence = new Hashtable<String, Sequence>();
	this.orderedAlleles = new ArrayList<String>();
    }

    public void readInGen(String gen){
	BufferedReader br = null;
	String curline = null;
	String refSeq = null;
	try{
	    br = new BufferedReader(new FileReader(gen));
	    boolean inMSF = false;
	    boolean firstSeq = true;
	    int startPos = 0;
	    while((curline=br.readLine())!=null){
		String stripped = curline.trim();
		if(!inMSF){
		    if(stripped.startsWith("gDNA")){
			inMSF = true;
		    }else
			this.header.append(curline + "\n");
		}else{
		    if(stripped.startsWith("AA codon"))
			;
		    else if(stripped.startsWith("|"))
			startPos = curline.indexOf("|");
		    else{
			String allele = curline.substring(0,startPos).trim();
			String msfsequence = curline.substring(startPos).trim();
			if(firstSeq){
			    firstSeq = false;
			    msfsequence = MergeMSFs.removeBlank(msfsequence, true);
			    refSeq = msfsequence;
			}else
			    msfsequence=MergeMSFs.abbrv2Seq(MergeMSFs.removeBlank(msfsequence,false), refSeq);
			
			this.orderedAlleles.add(allele);
			this.allele2Sequence.put(allele, new Sequence(msfsequence, true)); 
		    }
		}
	    }
	    br.close();
	}catch(IOException ioe){
	    ioe.printStackTrace();
	}
    }
    
    /* this removes all blank embedded in sequences*/
    public static String removeBlank(String sequence, boolean modifyHeader){
	StringBuffer bf = new StringBuffer();
	int headeri = 0;
	for(int i=0; i<sequence.length(); i++){
	    char tmp = sequence.charAt(i);
	    if(tmp == ' '){
		if(modifyHeader){
		    this.header.deleteCharAt(headeri);
		    headeri--;
		}
	    }else
		bf.append(tmp);
	    headeri++;
	}
	return bf.toString();
    }
    
    /* this fills *(unknown) and -(same as ref) as bases*/
    public static String abbrv2Seq(String abbrv, String refSeq){
	StringBuffer bf = new StringBuffer();
	for(int i=0;i<this.abbrv.length();i++){
	    char curChar = abbrv.charAt(i);
	    if(curChar == '*')
		bf.append(Character.toLowerCase(refSeq.charAt(i)));
	    else if(curChar == '-')
		bf.append(refSeq.charAt(i));
	    else
		bf.append(curChar);
	}
	return bf.toString();
    }


    private String stripHeader(String curline, int startPos){
	return curline.substring(startPos);
    }
    
    
}








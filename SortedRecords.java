import java.util.*;

public class SortedRecords{
    
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





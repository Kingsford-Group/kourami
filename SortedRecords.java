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

    public double getNthtpms(int n){
	return this.tpms.get(n).doubleValue();
    }
    
    public String extractUptoXdigit(int n){
	//System.err.println("tpms.size() = " + this.tpms.size());
	//System.err.println("names.size() = " + this.names.size());
	StringBuffer out = new StringBuffer();
	if(this.tpms.size() > 0){
	    for(int i=0; i<tpms.size() && i<8;i++){
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
    
    public void sumAndSort(){
	ArrayList<String> shortnames = new ArrayList<String>();
	for(int i=0; i<this.names.size(); i++)
	    shortnames.add(SortedRecords.getUptoXdigit(this.names.get(i), 4));
	
	Hashtable<String, String> visited = new Hashtable<String, String>();
	
	ArrayList<String> namesNew = new ArrayList<String>();
	ArrayList<Double> tpmsNew = new ArrayList<Double>();

	for(int i=0;i<this.names.size(); i++){
	    String curname = this.names.get(i);
	    Double curdbl = this.tpms.get(i);
	    if(visited.get(curname) == null){
		visited.put(curname, curname);
		curname = SortedRecords.getUptoXdigit(curname, 4);
		namesNew.add(curname);
		tpmsNew.add(curdbl);
		for(int j=i+1;j<this.names.size(); j++){
		    String curname2 = this.names.get(j);
		    if(visited.get(curname2) == null){
			if(curname.equals(SortedRecords.getUptoXdigit(curname2, 4))){
			    tpmsNew.set(tpmsNew.size()-1, new Double(tpmsNew.get(tpmsNew.size()-1).doubleValue() + this.tpms.get(j).doubleValue()));
			    visited.put(curname2, curname2);
			}
		    }
		}
	    }
	}
	this.sort(tpmsNew, namesNew);
	this.names = namesNew;
	this.tpms = tpmsNew;
    }

    private void sort(ArrayList<Double> t, ArrayList<String> n){
	for(int i=1; i<t.size(); i++){
	    double tmp = t.get(i).doubleValue();
	    String name = n.get(i);
	    int j;
	    for(j= i - 1; j >= 0 && tmp > t.get(j).doubleValue();j--){
		t.set(j+1, new Double(t.get(j).doubleValue()));
		n.set(j+1, n.get(j));
	    }
	    t.set(j+1, new Double(tmp));
	    n.set(j+1, name);
	}
    }


    public static String getUptoXdigit(String name, int x){
	if(x%2 != 0 && x < 2)
	    return null;
	else{
	    int sPos = name.indexOf("*")+1;
	    String digits = name.substring(sPos);
	    String[] tokens = digits.split(":");
	    StringBuffer buffer = new StringBuffer(tokens[0]);
	    for(int i=1; i<x/2; i++){
		buffer.append(":" + tokens[i]);
	    }
	    return buffer.toString();
	}
    }

    public String extractUptoXdigitFull(int n){
	StringBuffer out = new StringBuffer();
	if(this.tpms.size() > 0){
	    for(int i=0; i<tpms.size() && i<8;i++){
		if(tpms.get(i).doubleValue() > 0){
		    if(i > 0)
			out.append("\t");
		    //int starPos = this.names.get(i).indexOf("*")+1;
		    //out.append(this.names.get(i).substring(0,this.names.get(i).indexOf("*")) + "*" + this.names.get(i).substring(starPos, starPos+n+((n-1)/2)) + "(" + this.tpms.get(i).doubleValue()+")" );
		    out.append(this.names.get(i).substring(0, this.names.get(i).indexOf("*")) + "*" + getNthbestUptoXDigitWTPMS(i, n));
		}
	    }
	}else
	    out.append("NO ALLELE DETECTED");
	return out.toString();
    }
    
    public String getNthbestUptoXDigitWTPMS(int n, int x){
	//int starPos = this.names.get(n).indexOf("*")+1;
	//return this.names.get(n).substring(starPos, starPos+x+((x-1)/2)) + "(" + this.tpms.get(n).doubleValue()+")";
	return SortedRecords.getUptoXdigit(this.names.get(n), x) + "(" + this.tpms.get(n).doubleValue()+")";
    }

    public String getNthBestUptoXDigit(int n, int x){
	//int starPos = this.names.get(n).indexOf("*")+1;
	//return this.names.get(n).substring(starPos, starPos+x+((x-1)/2));
	return SortedRecords.getUptoXdigit(this.names.get(n), x);
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





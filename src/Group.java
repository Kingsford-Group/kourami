/*
Part of Kourami HLA typer/assembler
(c) 2017 by  Heewook Lee, Carl Kingsford, and Carnegie Mellon University.
See LICENSE for licensing.
*/
import java.util.*;

/*
 * G Group
 */
public class Group{

    private String hlaGeneName; //A B C DQA1 DQB1 DRB1 etc

    private String groupname; //representatitve name: ex) 01:01:01G
    
    private HashSet<String> set; 

    public String getGroupString(){
	return hlaGeneName + "*" + groupname;
    }

    public String getHLAGeneName(){
	return hlaGeneName;
    }
    
    public String getGroupName(){
	return this.hlaGeneName + "*" + this.groupname;
    }

    public String getFirstAllele(){
	return this.hlaGeneName + "*" + set.iterator().next();
    }
    
    public Group(String line, NomG nomG){
	this.set = new HashSet<String>();
	this.process(line, nomG);
    }

    public Group(String alleleName){
	this.set = new HashSet<String>();
	this.hlaGeneName = alleleName.substring(0,alleleName.indexOf("*"));
	this.groupname = alleleName.substring(alleleName.indexOf("*")+1);
	this.set.add(groupname);
    }
    
    public void process(String line, NomG nomG){
	String[] tokens = line.split(";");
	this.hlaGeneName = tokens[0].substring(0,tokens[0].indexOf("*"));
	String gName = null;
	if(tokens.length == 2 && line.endsWith(";"))
	    this.groupname = tokens[1];
	else
	    this.groupname = tokens[2];
	
	String[] elements = tokens[1].split("/");
	for(String e : elements){
	    this.set.add(e);
	    nomG.addToAllele2Group(e, this);
	}
    }
    
}

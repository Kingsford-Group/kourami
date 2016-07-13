/*
 * Simple class storing allele name and its typing sequence
 * 
 */
public class HLASequence{
    
    private Group grp;
    private String sequence;
    
    public HLASequence(Group g, Sequence seq){
	this.grp = g;
	this.sequence = seq.getTypingSequence();
    }
    
    public Group getGroup(){
	return this.grp;
    }
    
    public String toString(){
	return ">" + this.grp.getGroupString() + "\n" + this.sequence + "\n";
    }
}

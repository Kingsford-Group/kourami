public class HLASequence{
    
    private Group grp;
    private String sequence;
    
    public HLASequence(Group g, Sequence seq){
	this.grp = g;
	this.sequence = seq.getTypingSequence();
    }
}

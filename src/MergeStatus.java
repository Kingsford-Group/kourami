public class MergeStatus{

    private boolean split;
    private boolean segregating;
    private int lastSegregationColumnIndex;
    
    public MergeStatus(){
	this.split = false;
	this.segregating = false;
	this.lastSegregationColumnIndex = -1;
    }

    public MergeStatus(boolean segregating, int lastSegCI){
	this();
	this.segregating = segregating;
	this.lastSegregationColumnIndex = lastSegCI;
    }

    public int getLastSegregationColumnIndex(){
	return this.lastSegregationColumnIndex;
    }
    
    public boolean isSplit(){
	return this.split;
    }

    public boolean isSegregating(){
	return this.segregating;
    }

    public void setSplit(boolean b){
	this.split = b;
    }
    public void setSegregating(boolean b){
	this.segregating = b;
    }
    
    public void setLastSegregationColumnIndex(int ci){
	this.lastSegregationColumnIndex = ci;
    }
}

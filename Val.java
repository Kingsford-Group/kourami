public class Val{

    private int whichH;
    
    public Val(){
	whichH=-1;
    }

    public Val(int n){
	this.whichH = n;
    }
    
    public void set(int n){
	this.whichH = n;
    }
    
    public int getWhichH(){
	return this.whichH;
    }
}

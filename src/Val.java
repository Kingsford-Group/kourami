/*
Part of Kourami HLA typer/assembler
(c) 2017 by  Heewook Lee, Carl Kingsford, and Carnegie Mellon University.
See LICENSE for licensing.
*/
public class Val{

    private int whichH;//-1 if same score (split), 0 for H1, 1 for H2
    
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

/** All-pass plane-wave destruction filter coefficients
 * Modified from Sergey Fomel's apfilt.c
 * @author Xinming Wu and Dave Hale
 * @version 2015.10.01
 */

package pp;

public class PdFilter {

  public PdFilter (int nw) {
    _n = nw*2;
    _b = new double[_n+1];
    setB(_n,_b);
  }

  public void passfilter (float p, float[] a) {
    for (int k=0; k<=_n; k++) {
      double ak = _b[k];
 	    for (int j=0; j<_n; j++) {
	      if(j<_n-k) {ak *= (_n-j-p);} 
        else       {ak *= (p+j+1);}
	    }
	    a[k] = (float)ak;
    }
  }

  public void aderfilter (float p, float[] a) {
    for (int k=0; k<=_n; k++) {
      double ak = 0.;
	    for (int i=0; i<_n; i++) {
	      double ai = -1.0;
	      for (int j=0; j<_n; j++) {
	      if (j!=i) {			
	        if (j<_n-k){ai *= (_n-j-p);} 
          else       {ai *= (p+j+1);}
	      }else if (j<_n-k) {
	       ai *= -1;
	     }
	    }
	    ak += ai;
	  }
	  a[k] = (float)(ak*_b[k]);
    }
  }

  private int _n; 
  private double[] _b;

  private void setB(int n, double[] b){
    for (int k=0; k<=n; k++) {
	  double bk = 1.0;
	  for (int j=0; j<n; j++) {
	    if (j<n-k) {bk *= (k+j+1.0)/(2*(2*j+1)*(j+1));}
      else       {bk *= 1.0/(2*(2*j+1));}
	  }
	  b[k] = bk;
    }
  }
}


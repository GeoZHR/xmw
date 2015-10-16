/** Matrix definition for plane-wave destruction 
 * Modified from Sergey Fomel's pwd.c
 * @author Xinming Wu and Dave Hale
 * @version 2015.10.01
 */

package pp;

import static edu.mines.jtk.util.ArrayMath.*;

public class PdMatrix {
  int n, na;
  float[] b;
  float[][] a;

  public PdMatrix (int n1, int nw) {
    this.n = n1;
    this.na = 2*nw+1;
    this.b = new float[na];
    this.a = new float[na][n1];
  }

  public static void define (boolean adj, PdMatrix w,
    float[] pp,float[] diag, float[][] offd)
  {
    int n = w.n;
    int na = w.na;
    int nw = (na-1)/2;
    PdFilter pf = new PdFilter(nw);
    for (int i=0; i<n; i++) {
	    pf.passfilter(pp[i],w.b);
	  for (int j=0; j<na; j++) {
	    if(adj) {w.a[j][i] = w.b[na-1-j];} 
      else    {w.a[j][i] = w.b[j];}
	  }}
    for (int i=0; i<n; i++) {
	  for (int j=0; j<na; j++) {
	    int k = i+j-nw;
	    if(k>=nw && k<n-nw){diag[i] += w.a[j][k]*w.a[j][k];}
	  }}
    for (int i=0; i<n; i++) {
	  for (int m=0; m<2*nw; m++) {
	  for (int j=m+1; j<na; j++) {
	   int k = i+j-nw;
	   if (k>=nw && k<n-nw) {
	     float aj = w.a[j][k];
	     float am = w.a[j-m-1][k];
	     offd[m][i] += am*aj;
	   }
	  }}}
  }

  public static void set(boolean adj, PdMatrix w,
    float[] inp, float[] out, float[] tmp)
  {
    zero(tmp);
    int n  = w.n;
    int na = w.na;
    int nw = (na-1)/2;
    float[][] a = w.a;
    if (adj) {
      zero(inp);
      compute1(n,na,nw,a,out,tmp);
      compute2(n,na,nw,a,tmp,inp);
    } else {
      zero(out);
      compute2(n,na,nw,a,inp,tmp);
      compute1(n,na,nw,a,tmp,out);
    }
  }

  private static void compute1(
    int n, int na, int nw, float[][] a, float[] x, float[] y) {
    for(int i=0; i<n;  i++) {
    for(int j=0; j<na; j++) {
      int k = i+j-nw;
      if(k>=nw && k<n-nw){y[k] += a[j][k]*x[i];}
    }}
  }

  private static void compute2(
    int n, int na, int nw, float[][] a, float[] x, float[] y) {
    for(int i=nw; i<n-nw; i++) {
    for(int j=0;  j<na;   j++) {
      int k = i+j-nw;
		  y[k] += a[j][i]*x[i];
	  }}
  }

}


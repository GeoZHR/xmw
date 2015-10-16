/** Banded matrix solver.
 * Modified from Sergey Fomel's banded.c
 * @author Xinming Wu and Dave Hale
 * @version 2015.10.01
 */

package pp;

import static edu.mines.jtk.util.ArrayMath.*;

public class BandedMatrixSolver {
  int n, band;
  float[] d;
  float[][] o;

  public BandedMatrixSolver (int n, int band) {
    this.n = n;
    this.band = band;
    this.d = new float[n];
    this.o = new float[band][];
    for (int i=0; i<band; i++) {
      this.o[i] = new float[n-i-1];
    }
  }

  public static void define(BandedMatrixSolver slv, float[] diag,float[][] offd) {
    int n = slv.n;
    int band = slv.band;
    for(int k=0; k<n; k++) {
	    float tk = diag[k];
	    int m1 = min(k,band);
	    for(int m=0; m<m1; m++) {
        int km = k-m-1;
        float dk = slv.d[km];
        float ok = slv.o[m][km];
        tk -= ok*ok*dk;
      }
      slv.d[k] = tk;
	    int n1 = min(n-k-1,band);
	    for(int i=0; i<n1; i++) {
	      tk = offd[i][k];
	      m1 = min(k,band-i-1);
	      for(int m=0; m<m1; m++) {
          int km = k-m-1;
          float dk  = slv.d[km];
          float omk = slv.o[m    ][km];
          float opk = slv.o[i+m+1][km];
          tk -= omk*opk*dk;
	      }
	      slv.o[i][k] = tk/slv.d[k];
	    }
    }
  }

  public static void solve (BandedMatrixSolver slv, float[] b) {
    int n = slv.n;
    int band = slv.band;
    for(int k=1; k<n; k++) {
	    float tk = b[k];
	    int m1 = min(k,band);
	    for(int m=0; m<m1; m++)
	      tk -= slv.o[m][k-m-1]*b[k-m-1];
	    b[k] = tk;
    }
    for(int k=n-1; k>=0; k--) {
	    float tk = b[k]/slv.d[k];
	    int m1 = min(n-k-1,band);
	    for(int m=0; m<m1; m++)
	      tk -= slv.o[m][k]*b[k+m+1];
	    b[k] = tk;
    }
  }

}




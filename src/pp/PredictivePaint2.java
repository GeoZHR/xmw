package pp;

import edu.mines.jtk.dsp.*;
import static edu.mines.jtk.util.ArrayMath.*;

public class PredictivePaint2 {

  public PredictivePaint2(int n1, int k2, int ord, float eps) {
    _k2 = k2;
    _ord = ord;
    _eps = eps;
    _w1 = new PdMatrix(n1,_ord);
    _bms = new BandedMatrixSolver(n1,2*ord);
  }

  public float[][] paint(Sampling s1, float[][] p) {
    int n2 = p.length;
    int n1 = p[0].length;
    float[] time = new float[n1];
    float f1 = (float)s1.getFirst();
    float d1 = (float)s1.getDelta();
    float[][] u = new float[n2][n1];
    for(int i1=0; i1<n1; ++i1){time[i1]=f1+d1*i1;}
    float[] trace = copy(time);
    u[_k2] = copy(trace);
    for (int i2=_k2-1; i2>=0; i2--) {
      predict(true,false,trace,p[i2]);
      u[i2] = copy(trace);
    }
    trace = copy(u[_k2]);
	  for (int i2=_k2+1; i2<n2; i2++) {
	    predict(true,true,trace,p[i2-1]);
      u[i2] = copy(trace);
	  }
    return u;
  }

  public void predict(
    boolean adj, boolean forw, float[] trace, float[] pp) 
  {
    int nb = 2*_ord;
    int n1 = pp.length;
    float[] diag = new float[n1];
    float[][] offd = new float[nb][n1];
    regularization(diag,offd);
    PdMatrix.define(forw,_w1,pp,diag,offd);
    BandedMatrixSolver.define(_bms,diag,offd);
    if(adj){BandedMatrixSolver.solve(_bms,trace);}
    float t0 = trace[0];
    float t1 = trace[1];
    float t2 = trace[n1-2];
    float t3 = trace[n1-1];
    PdMatrix.set(adj,_w1,trace,trace,diag);
    float eps = _eps*_eps;
    trace[0] += eps*t0;
    trace[1] += eps*t1;
    trace[n1-2] += eps*t2;
    trace[n1-1] += eps*t3;
    if(!adj){BandedMatrixSolver.solve(_bms,trace);}
  }

  private void regularization(float[] diag, float[][] offd) {
    int n1 = diag.length;
    float eps = _eps*_eps;
    for (int i1=0; i1<n1; ++i1) {
      diag[i1]    =  6.0f*eps;
      offd[1][i1] =       eps;
      offd[0][i1] = -4.0f*eps;
    }
    diag[0] = diag[n1-1] = eps+eps;
    diag[1] = diag[n1-2] = eps+5f*eps;
    offd[0][0] = offd[0][n1-2] = -2f*eps;
  }


  ///////////////////////////////////////////////////////////////////////////
  // private
  private PdMatrix _w1;
  private BandedMatrixSolver _bms;
  private int _ord = 1;        // accuracy order
  private float _eps = 0.01f;  // regularization
  private int _k2;             // reference trace

}

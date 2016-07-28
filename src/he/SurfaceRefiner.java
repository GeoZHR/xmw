package he;

import edu.mines.jtk.dsp.*;
import static edu.mines.jtk.util.ArrayMath.*;

import vec.*;
import util.*;

/**
 * This method adjusts a manually or automatically interpreted 
 * seismic horizon so that it follows seismic amplitude throughs, 
 * peaks, or zero crosings.
 * A horizon after the adjustments can reveal more geologic details.
 * @author Xinming
 * @version 2016.07.27
 */

public class SurfaceRefiner {

  public void setSmoothings(float sigma1, float sigma2){
    _sigma1 = sigma1;
    _sigma2 = sigma2;
  }
  
  public void setCG(float small, int niter){
    _small = small;
    _niter = niter;
  }

  public void setExternalIterations(int exniter){
    _exniter = exniter;
  }



  public void surfaceUpdateFromAmplitudes(float[][][] gx, 
     float[][] wp,
     float[] k1, float[] k2, float[] k3,float[][] sf)
  {	
    int n3 = gx.length; 
    int n2 = gx[0].length; 
    int n1 = gx[0][0].length; 
    float lmt = (float)n1-1.f;
    float[][] sft = copy(sf);
    float[][] b   = new float[n3][n2]; 
    float[][] ks = new float[][]{k1,k2,k3};
    float[][][] g1 = new float[n3][n2][n1];
    RecursiveGaussianFilterP rgf = new RecursiveGaussianFilterP(1.0);
    rgf.apply100(gx,g1);
    //float[][] ks = updateConstraints(k1,k2,k3,p,q,ep,surf);
    for (int n=1; n<=_exniter; n++){
      System.out.println(" Iteration "+n+"......");
      VecArrayFloat2 vb    = new VecArrayFloat2(b);
      VecArrayFloat2 vs = new VecArrayFloat2(sf);
      A2 a2 = new A2(wp);
      M2 m2 = new M2(_sigma1,_sigma2,ks[1],ks[2]);
      CgSolver cs = new CgSolver(_small,_niter);
      vb.zero();
      makeRhs(g1,sf,b);
      cs.solve(a2,m2,vb,vs);
      sf = vs.getArray();
      surfaceCorrect(sf,lmt);
      float ad = sum(abs(sub(sft,sf)))/(n3*n2); 
      System.out.println(" Average adjustments per sample = "+ad);
      if (ad<0.002f) break;
      sft = copy(sf);
    }
  }


  private void surfaceCorrect(float[][] surf, float lmt) {
    int n1 = surf[0].length;
    int n2 = surf.length;
    for(int i1=0; i1<n1; ++i1) {
      for (int i2=0; i2<n2; ++i2) {
        if (surf[i2][i1]<0.f) surf[i2][i1]=0.f;
        if (surf[i2][i1]>lmt) surf[i2][i1]=lmt;
      }
    }
  }

  private void makeRhs(float[][][] g1, float[][] x, float[][] b) {
    int n3 = g1.length;
    int n2 = g1[0].length;
    int n1 = g1[0][0].length;
    SincInterpolator si = new SincInterpolator();
    for (int i3=0; i3<n3; i3++){
    for (int i2=0; i2<n2; i2++){
      double x2 = (double)i2;
      double x3 = (double)i3;
      double x1 = (double)x[i3][i2];
	    b[i3][i2] = si.interpolate(n1,1.0,0.0,n2,1.0,0.0,n3,1.0,0.0,g1,x1,x2,x3);
    }}
  }

  private static class A2 implements CgSolver.A{
    A2(float[][] wp){
      _wp = wp;
    }  
    public void apply(Vec vx, Vec vy){
      VecArrayFloat2 v2x = (VecArrayFloat2) vx;
      VecArrayFloat2 v2y = (VecArrayFloat2) vy;
      float[][] x = v2x.getArray();
      float[][] y = v2y.getArray();
      float[][] z = copy(x);
      v2y.zero();
      applyLaplacian(_wp,x,z);
      applyLaplacian(_wp,z,y);
    }
    private float[][] _wp;
  }


  // Preconditioner; includes smoothers and (optional) constraints.
  private static class M2 implements CgSolver.A {
    M2(float sigma1, float sigma2, float[] k2, float[] k3) {
      _sigma1 = sigma1;
      _sigma2 = sigma2;
      if (k2!=null && k3!=null) {
        _k2 = copy(k2);
        _k3 = copy(k3);
      }
    }
    public void apply(Vec vx, Vec vy) {
      VecArrayFloat2 v2x = (VecArrayFloat2)vx;
      VecArrayFloat2 v2y = (VecArrayFloat2)vy;
      float[][] x = v2x.getArray();
      float[][] y = v2y.getArray();
      copy(x,y);
      constrain(_k2,_k3,y);
      smooth2(_sigma2,y);
      smooth1(2.f*_sigma1,y);
      smooth2(_sigma2,y);
      constrain(_k2,_k3,y);
    }
    private float _sigma1,_sigma2;
    private float[] _k2,_k3;
  }

  private static void constrain(float[] k2, float[] k3, float[][] x) {
    if (k2!=null && k3!=null) {
      int np = k2.length;
      for (int ip=0; ip<np; ++ip) {
        int i2 = (int)k2[ip]; 
        int i3 = (int)k3[ip]; 
        x[i3][i2] = 0.f;
      }
    }
  }

  // Smoothing for dimension 2.
  private static void smooth1(float sigma, float[][] x) {
    if (sigma<1.0f)
      return;
    RecursiveExponentialFilter.Edges edges =
      RecursiveExponentialFilter.Edges.OUTPUT_ZERO_SLOPE;
    RecursiveExponentialFilter ref = new RecursiveExponentialFilter(sigma);
    ref.setEdges(edges);
    ref.apply1(x,x);
  }

  // Smoothing for dimension 2.
  private static void smooth2(float sigma, float[][] x) {
    if (sigma<1.0f)
      return;
    RecursiveExponentialFilter.Edges edges =
      RecursiveExponentialFilter.Edges.OUTPUT_ZERO_SLOPE;
    RecursiveExponentialFilter ref = new RecursiveExponentialFilter(sigma);
    ref.setEdges(edges);
    ref.apply2(x,x);
  }


  private static void applyLaplacian(
    float[][] w, float[][] x, float[][] y){
    int n2 = x.length;
    int n1 = x[0].length;
    for (int i2=1; i2<n2; ++i2) {
      for (int i1=1; i1<n1; ++i1) {
        float wi = w[i2][i1]*0.25f;
        float xa = 0.0f;
        float xb = 0.0f;
        xa += x[i2  ][i1  ];
        xb -= x[i2  ][i1-1];
        xb += x[i2-1][i1  ];
        xa -= x[i2-1][i1-1];
        float x1 = (xa+xb);
        float x2 = (xa-xb);
        float y1 = x1*wi;
        float y2 = x2*wi;
        float ya = (y1+y2);
        float yb = (y1-y2);
        y[i2  ][i1  ] += ya;
        y[i2  ][i1-1] -= yb;
        y[i2-1][i1  ] += yb;
        y[i2-1][i1-1] -= ya;
      }
    }
  }


  private float _sigma1 = 10;
  private float _sigma2 = 10;
  private float _small = 0.001f;
  private int _niter = 200;
  private int _exniter = 10; // external iterations of surface updating


}


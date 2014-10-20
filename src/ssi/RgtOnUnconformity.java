package ssi;

import java.util.*;
import edu.mines.jtk.dsp.*;
import static edu.mines.jtk.util.ArrayMath.*;


/**
 * Extract a single seismic horizon surface with control points
 * the constraints derived from control points are incorporated in a 
 * preconditioner in the conjugate gradient method used to solve the 
 * linear system for horizon extracting.
 * @author Xinming Wu
 * @version 2014.08.05
 */

public class RgtOnUnconformity {
 
  public void setSmoothings(float sigma1, float sigma2){
    _sigma1 = sigma1;
    _sigma2 = sigma2;
  }
  
  public void setCG(float small, int niter){
    _small = small;
    _niter = niter;
  }

  // Updates the surface using the seismic normal vectors and control points.
  public void rgtFromNormals
    (float[][] ep,float[][] u2,float[][] u3, float[][] rgt, boolean NaNs)
  {	
    int n3 = u2.length; 
    int n2 = u2[0].length; 
    ArrayList<int[]> bid = new ArrayList<int[]>();
    ArrayList<int[]> nid = new ArrayList<int[]>();
    ArrayList<int[]> pid = new ArrayList<int[]>();
    if(NaNs) {
      boundsAndNeighbors(pid,bid,nid,rgt);
      if(bid.size()>0){fillBoundaries(bid,nid,rgt);}
    }
    float[][] b   = new float[n3][n2]; 
    VecArrayFloat2 vb = new VecArrayFloat2(b);
    VecArrayFloat2 vr = new VecArrayFloat2(rgt);
    setWeightsForSmoothing(bid, ep);
    mute(0.0f,pid,u2,u3,ep);
    Smoother2 smoother2 =  new Smoother2(_sigma1,_sigma2,ep);
    A2 a2 = new A2(smoother2,ep);
    CgSolver cs = new CgSolver(_small, _niter);
    vb.zero();
    makeRhs(ep,u2,u3,b);
    smoother2.applyTranspose(b);
    cs.solve(a2,vb,vr);
    smoother2.apply(rgt);
    restoreNaNs(pid,rgt);
  }

  private static class A2 implements CgSolver.A{
    A2(Smoother2 s2, float[][] wp) {
      _s2 = s2;
      _wp = wp;
    }
    public void apply(Vec vx, Vec vy) {
      VecArrayFloat2 v2x = (VecArrayFloat2)vx;
      VecArrayFloat2 v2y = (VecArrayFloat2)vy;
      float[][] x = v2x.getArray();
      float[][] y = v2y.getArray();
      float[][] z = copy(x);
      zero(y);
      _s2.apply(z);
      applyLhs(_wp,z,y);
      _s2.applyTranspose(y);
    }
    float[][] _wp;
    private Smoother2 _s2;
  }

  private static class Smoother2 {
    public Smoother2(float sigma1, float sigma2, float[][] ep) {
      _ep = ep;
      _sigma1 = sigma1;
      _sigma2 = sigma2;
    }
    public void apply(float[][] x) {
      smooth2(_sigma2,_ep,x);
      smooth1(_sigma1,_ep,x);
    }
    public void applyTranspose(float[][] x) {
      smooth1(_sigma2,_ep,x);
      smooth2(_sigma1,_ep,x);
    }
    private float _sigma1, _sigma2;
    private float[][] _ep;
  }

  // Smoothing for dimension 1
  private static void smooth1(float sigma, float[][] s, float[][] x){
    if (sigma<1.0f)
      return;
    float c = 0.5f*sigma*sigma;
    int n1 = x[0].length;
    int n2 = x.length;
    float[] xt = zerofloat(n1);
    float[] yt = zerofloat(n1);
    float[] st = zerofloat(n1);
    LocalSmoothingFilter lsf = new LocalSmoothingFilter();
    for (int i2=0; i2<n2; ++i2) {
      for (int i1=0; i1<n1; ++i1) {
        xt[i1] = x[i2][i1];
        st[i1] = s[i2][i1];
      }
      lsf.apply(c,st,xt,yt);
      for (int i1=0; i1<n1; ++i1)
        x[i2][i1] = yt[i1];
    }
  }

  // Smoothing for dimension 2
  private static void smooth2(float sigma, float[][] s, float[][] x){
    if (sigma<1.0f)
      return;
    float c = 0.5f*sigma*sigma;
    int n1 = x[0].length;
    int n2 = x.length;
    float[] xt = zerofloat(n2);
    float[] yt = zerofloat(n2);
    float[] st = zerofloat(n2);
    LocalSmoothingFilter lsf = new LocalSmoothingFilter();
    for (int i1=0; i1<n1; ++i1) {
      for (int i2=0; i2<n2; ++i2) {
        xt[i2] = x[i2][i1];
        st[i2] = s[i2][i1];
      }
      lsf.apply(c,st,xt,yt);
      for (int i2=0; i2<n2; ++i2)
        x[i2][i1] = yt[i2];
    }
  }

  private static void applyLhs( float[][] wp, float[][] x, float[][] y) {
    zero(y);
    int n1 = x[0].length;
    int n2 = x.length;
    for (int i2=1; i2<n2; ++i2) {
      for (int i1=1; i1<n1; ++i1) {
        float wpi = (wp!=null)?wp[i2][i1]:1.000f;
        if(wpi<0.05f && wpi!=0.0f){wpi=0.05f;}
        float wps = wpi*wpi;
        float d11 = wps;
        float d22 = wps;
        float xa = 0.0f;
        float xb = 0.0f;
        xa += x[i2  ][i1  ];
        xb -= x[i2  ][i1-1];
        xb += x[i2-1][i1  ];
        xa -= x[i2-1][i1-1];
        float x1 = 0.5f*(xa+xb);
        float x2 = 0.5f*(xa-xb);
        float y1 = d11*x1;
        float y2 = d22*x2;
        float ya = 0.5f*(y1+y2);
        float yb = 0.5f*(y1-y2);
        y[i2  ][i1  ] += ya;
        y[i2  ][i1-1] -= yb;
        y[i2-1][i1  ] += yb;
        y[i2-1][i1-1] -= ya;
      }
    }
  }
  
  //private static void makeRhs
  private static void makeRhs(float[][] wp, float[][] p2, float[][] p3, float[][] y) {
    zero(y);
    int n1 = y[0].length;
    int n2 = y.length;
    for (int i2=1; i2<n2; ++i2) {
      for (int i1=1; i1<n1; ++i1) {
        float wpi = (wp!=null)?wp[i2][i1]:1.000f;
        if(wpi<0.05f && wpi!=0.0f){wpi=0.05f;}
        float p2i = p2[i2][i1];
        float p3i = p3[i2][i1];
        float b11 = wpi;
        float b22 = wpi;
        float x1 = wpi*p2i;
        float x2 = wpi*p3i;
        float y1 = b11*x1;
        float y2 = b22*x2;
        float ya = 0.5f*(y1+y2);
        float yb = 0.5f*(y1-y2);
        y[i2  ][i1  ] += ya;
        y[i2  ][i1-1] -= yb;
        y[i2-1][i1  ] += yb;
        y[i2-1][i1-1] -= ya;
      }
    }
  }

  private void restoreNaNs(ArrayList<int[]> pid, float[][] x) {
    int np = pid.size();
    for (int ip=0; ip<np; ++ip) {
      int i2 = pid.get(ip)[0];
      int i3 = pid.get(ip)[1];
      x[i3][i2] = Float.NaN;
    }
  }

  private void setWeightsForSmoothing(ArrayList<int[]> bid, float[][] ws) {
    int np = bid.size();
    int n3 = ws.length;
    int n2 = ws[0].length;
    for (int ip=0; ip<np; ++ip) {
      int i2 = bid.get(ip)[0];
      int i3 = bid.get(ip)[1];
      ws[i3][i2] = 0.05f;
    }

    for (int i3=0; i3<n3; ++i3) {
      for (int i2=0; i2<n2; ++i2) {
        float wsi = ws[i3][i2];
        if (Float.isNaN(wsi)) {
          continue;
        } else if (wsi<0.05f) {
          ws[i3][i2] = 0.05f;
        }
      }
    }
  }

  private void mute(float v, ArrayList<int[]> pid, 
    float[][] p, float[][] q, float[][] w) {
    int np = pid.size();
    for (int ip=0; ip<np; ++ip) {
      int i2 = pid.get(ip)[0];
      int i3 = pid.get(ip)[1];
      p[i3][i2] = v;
      q[i3][i2] = v;
      w[i3][i2] = v;
    }
  }

  private void boundsAndNeighbors(
    ArrayList<int[]> pid, ArrayList<int[]> bid, ArrayList<int[]> nid, float[][] sf) 
  {
    int n3 = sf.length;
    int n2 = sf[0].length;
    float[][] sft = copy(sf);
    for (int i3=0; i3<n3; ++i3) {
      for (int i2=0; i2<n2; ++i2) {
        float sfi = sft[i3][i2];
        if(Float.isNaN(sfi)) {
          sf[i3][i2] = 0.0f;
          pid.add(new int[]{i2,i3});
          int[] idn = new int[2];
          boolean isBound = checkBound(i2,i3,n2,n3,sft,idn); 
          if(isBound) {
            nid.add(idn);
            bid.add(new int[]{i2,i3});
          }
        }
      }
    }
  }

  private boolean checkBound(int i2, int i3, int n2, int n3, float[][] sf, int[] ind) {
    int i2m = i2-1;
    int i2p = i2+1;
    int i3m = i3-1;
    int i3p = i3+1;
    if(i2m>=0 && !Float.isNaN(sf[i3][i2m])){
      ind[0]=i2m;
      ind[1]=i3;
      return true;
    }else if(i2p<n2 && !Float.isNaN(sf[i3][i2p])){
      ind[0]=i2p;
      ind[1]=i3;
      return true;
    }else if(i3m>=0 && !Float.isNaN(sf[i3m][i2])){
      ind[0]=i2;
      ind[1]=i3m;
      return true;
    }else if(i3p<n3 && !Float.isNaN(sf[i3p][i2])){
      ind[0]=i2;
      ind[1]=i3p;
      return true;
    }else if(i2m>=0 && i3m>=0 && !Float.isNaN(sf[i3m][i2m])){
      ind[0]=i2m;
      ind[1]=i3m;
      return true;
    }else if(i2p<n2 && i3p<n3 && !Float.isNaN(sf[i3p][i2p])){
      ind[0]=i2p;
      ind[1]=i3p;
      return true;
    }else{
      return false;
    }
  }

  private void fillBoundaries(
    ArrayList<int[]> bid, ArrayList<int[]> nid, float[][] sf) 
  {
    int n = bid.size();
    for (int i=0; i<n; ++i) {
      int bi2 = bid.get(i)[0];
      int bi3 = bid.get(i)[1];
      int ni2 = nid.get(i)[0];
      int ni3 = nid.get(i)[1];
      sf[bi3][bi2] = sf[ni3][ni2];
    }
  }
  
  ///////////////////////////////////////////////////////////////////////////
  // private
  private float _sigma1 = 6.0f; // precon smoothing extent for 1st dim
  private float _sigma2 = 6.0f; // precon smoothing extent for 2nd dim
  private float  _small = 0.01f; // stop CG iterations if residuals small
  private int    _niter = 200; // maximum number of CG iterations
}

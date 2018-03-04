package ipsi;

import java.util.*;
import edu.mines.jtk.dsp.*;
import edu.mines.jtk.interp.*;
import static edu.mines.jtk.util.ArrayMath.*;

import vec.*;
import util.*;



/**
 * Extract a single seismic horizon surface with control points
 * the constraints derived from control points are incorporated in a 
 * preconditioner in the conjugate gradient method used to solve the 
 * linear system for horizon extracting.
 * @author Xinming Wu and Dave Hale
 * @version 2014.03.12
 */

public class SurfaceExtractor {
 
  // scale the curvature term 
  public void setWeights(float w){
    _weight = w;
  }
  
  public void setConstraints(float[] k1, float[] k2, float[] k3) {
    _k1 = k1;
    _k2 = k2;
    _k3 = k3;
  }
  public void setSmoothings(float sigma1, float sigma2){
    _sigma1 = sigma1;
    _sigma2 = sigma2;
  }
  
  public void setCG(float small, int niter){
    _small = small;
    _niter = niter;
  }
  
  public float[] getK1(){
    return _k1;
  }

  public float[] getK2(){
    return _k2;
  }
  public float[] getK3(){
    return _k3;
  }

  public ArrayList<float[][]> accumulateConstraints(float[][] sf) {
    _extendC.add(sf);
    return _extendC;
  }
  public ArrayList<int[]> accumulateInd(int[] ind) {
    _extendInd.add(ind);
    return _extendInd;
  }

  // find the peak or trough nearest to each control point
  public void refineConstraints(float[][][] u) {
    int np = _k1.length;
    int n1 = u[0][0].length;
    for (int ip=0; ip<np; ++ip) {
      float k1i = _k1[ip];
      int i1 = (int)k1i;
      int i2 = (int)_k2[ip];
      int i3 = (int)_k3[ip];
      int i1m = i1-1; 
      int i1p = i1+1; 
      if(i1m<0  ){i1m=0;   i1=1;   i1p=2;}
      if(i1p>=n1){i1m=n1-3;i1=n1-2;i1p=n1-1;}
      float um = u[i3][i2][i1m];
      float ui = u[i3][i2][i1 ];
      float up = u[i3][i2][i1p];
      float kp = parabolicPeak(k1i,um,ui,up);
      if(abs(kp-k1i)<2.0){_k1[ip]=kp;}
    }
  }  

  // Interpolate an initial surface passing through control points
  public float[][] surfaceInitialization(int n2, int n3, float[][][] p, boolean NaNs) {
    int n1 = p[0][0].length;
    float lmt = (float)n1-1.0f;
    if (_k1.length==1) {
      float[][] surf = zerofloat(n2,n3);
      add(surf,_k1[0],surf);
      if(NaNs){restoreNaNs(p,surf);}
      return surf; 
    } else {
      Sampling s2 = new Sampling(n2,1.0f,0.0f);
      Sampling s3 = new Sampling(n3,1.0f,0.0f);
      RadialInterpolator2.Biharmonic bs = new RadialInterpolator2.Biharmonic();
      RadialGridder2 rg = new RadialGridder2(bs,_k1,_k2,_k3);    
      float[][] surf = rg.grid(s2,s3);
      surfaceCorrect(surf,lmt);
      if(NaNs){restoreNaNs(p,surf);}
      checkControlPoints(_k2, _k3, surf); 
      return surf;
    }
  }

  // Updates the surface using the seismic normal vectors and control points.
  public float surfaceUpdateFromSlopes
    (float[][][] ep,float[][][] p,float[][][] q,float[][] surf, boolean NaNs)
  {	
    int n3 = p.length; 
    int n2 = p[0].length; 
    int n1 = p[0][0].length; 
    float lmt = (float)n1-1.f;
    float[][] surft = copy(surf);
    ArrayList<int[]> bid = new ArrayList<int[]>();
    ArrayList<int[]> nid = new ArrayList<int[]>();
    ArrayList<int[]> pid = new ArrayList<int[]>();
    if(NaNs) {
      boundsAndNeighbors(pid,bid,nid,surf);
      if(bid.size()>0){fillBoundaries(bid,nid,surf);}
    }
    float[][] b   = new float[n3][n2]; 
    float[][] pi1 = new float[n3][n2]; 
    float[][] qi1 = new float[n3][n2]; 
    float[][] wi1 = new float[n3][n2]; 
    //float[][] ws1 = new float[n3][n2]; 
    VecArrayFloat2 vb    = new VecArrayFloat2(b);
    VecArrayFloat2 vsurf = new VecArrayFloat2(surf);
    if(!NaNs){
      updateSlopesAndWeights(p,q,ep,surf,pi1,qi1,wi1,0.0f);
    }else{
      updateSlopesAndWeights(mute(p),mute(q),mute(ep),surf,pi1,qi1,wi1,0.0f);
      mute(0.0f,pid,pi1,qi1,wi1);
    }
    setWeightsForSmoothing(bid, wi1);
    A2 a2 = new A2(wi1,0.0f);
    M2 m2 = new M2(_sigma1,_sigma2,wi1,_k2,_k3);
    CgSolver cs = new CgSolver(_small, _niter);
    vb.zero();
    makeRhs(wi1,pi1,qi1,b);
    cs.solve(a2,m2,vb,vsurf);
    surf = vsurf.getArray();
    surfaceCorrect(surf,lmt);
    restoreNaNs(pid,surf);
    float ad = adjustPerSample(surft, surf);
    return ad;
  }

  private static void updateSlopesAndWeights (
    float[][][] p, float[][][] q, float[][][] w,
    float[][] surf, float[][] pi1, float[][] qi1,float[][] wi1,float dth)
  {
    int n3 = p.length;
    int n2 = p[0].length;
    int n1 = p[0][0].length;
    SincInterp si = new SincInterp();
    si.setExtrapolation(SincInterp.Extrapolation.CONSTANT);
    for (int i3=0; i3<n3; i3++){
      for (int i2=0; i2<n2; i2++){
        float x2i = (float) i2;
        float x3i = (float) i3;
        float x1i = surf[i3][i2];
	      float wi = si.interpolate(n1,1,0,n2,1,0,n3,1,0,w,x1i,x2i,x3i);
        if(wi<0.0005f && wi!=0.0f){wi1[i3][i2]=0.0005f;} 
        else                      {wi1[i3][i2]=wi;}
        pi1[i3][i2] = si.interpolate(n1,1,0,n2,1,0,n3,1,0,p,x1i,x2i,x3i);
	      qi1[i3][i2] = si.interpolate(n1,1,0,n2,1,0,n3,1,0,q,x1i,x2i,x3i);
      }
    }
  }

  private void surfaceCorrect(float[][] surf, float lmt) {
    int n1 = surf[0].length;
    int n2 = surf.length;
    for(int i1=0; i1<n1; ++i1) {
      for (int i2=0; i2<n2; ++i2) {
        float sfi = surf[i2][i1];
        if (sfi!=Float.NaN) {
          if (sfi<0.f) surf[i2][i1]=0.f;
          if (sfi>lmt) surf[i2][i1]=lmt;
        }
      }
    }
  }

  private static class A2 implements CgSolver.A{
    A2(float[][] wp,float w){
      _w  = w;
      _wp = wp;
    }  
    public void apply(Vec vx, Vec vy){
      VecArrayFloat2 v2x = (VecArrayFloat2) vx;
      VecArrayFloat2 v2y = (VecArrayFloat2) vy;
      float[][] x = v2x.getArray();
      float[][] y = v2y.getArray();
      int n1 = y[0].length; int n2 = y.length;
      float[][] yy = new float[n2][n1];
      float[][] yt = new float[n2][n1];
      VecArrayFloat2 v2yy = new VecArrayFloat2(yy);
      v2y.zero();
      v2yy.zero();
      applyLhs(_wp,x,y);
      if (_w>0.0f) {
        applyLhs(_wp,x,yt);
        applyLhs(_wp,yt,yy);
        v2y.add(1.f,v2yy,_w);
      }
    }
    private float _w;
    private float[][] _wp;
  }

   // Preconditioner; includes smoothers and (optional) constraints.
  private static class M2 implements CgSolver.A {
    M2(float sigma1, float sigma2, float[][] wp, float[] k2, float[] k3) 
    {
      _wp = wp;
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
      smooth2(_sigma2,_wp,y);
      smooth1(2.f*_sigma1,_wp,y);
      smooth2(_sigma2,_wp,y);
      constrain(_k2,_k3,y);
    }
    private float[][] _wp;
    private float[] _k2,_k3;
    private float _sigma1,_sigma2;
  }

  private static void checkControlPoints(float[] k2, float[] k3, float[][] f) {
    if (k2!=null && k3!=null) {
      int np = k2.length;
      for (int ip=0; ip<np; ++ip) {
        int i2 = (int)k2[ip];
        int i3 = (int)k3[ip];
        System.out.println(" i2="+i2+" i3="+i3+" f1="+f[i3][i2]);
      }
    }
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

  private float adjustPerSample(float[][] sfm, float[][] sfi) {
    float cot = 0.0f; 
    float sum = 0.0f;
    int n3 = sfm.length;
    int n2 = sfm[0].length;
    for (int i3=0; i3<n3; ++i3) {
      for (int i2=0; i2<n2; ++i2) {
        float sfmi = sfm[i3][i2];
        float sfii = sfi[i3][i2];
        if(Float.isNaN(sfmi)||Float.isNaN(sfii)) {
          continue;
        } else {
          cot += 1.0f;
          sum += abs(sfii-sfmi);
        }
      }
    }
    return sum/cot;
  }

  private void restoreNaNs(float[][][] f, float[][] x) {
    int n3 = x.length;
    int n2 = x[0].length;
    for (int i3=0; i3<n3; ++i3) {
      for (int i2=0; i2<n2; ++i2) {
        int i1 = round(x[i3][i2]);
        if(Float.isNaN(f[i3][i2][i1])) {
          x[i3][i2] = Float.NaN;
        }
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


  private static float[][][] mute(float[][][] f) {
    int n3 = f.length;
    int n2 = f[0].length;
    int n1 = f[0][0].length;
    float[][][] fm = new float[n3][n2][n1];
    for (int i3=0; i3<n3; ++i3) {
      for (int i2=0; i2<n2; ++i2) {
        for (int i1=0; i1<n1; ++i1) {
          float fi = f[i3][i2][i1];
          if(Float.isNaN(fi)) {
            continue;
          } else {
            fm[i3][i2][i1] = f[i3][i2][i1];
          }
        }
      }
    }
    return fm;
  }

  private void setWeightsForSmoothing(ArrayList<int[]> bid, float[][] ws) {
    int np = bid.size();
    for (int ip=0; ip<np; ++ip) {
      int i2 = bid.get(ip)[0];
      int i3 = bid.get(ip)[1];
      ws[i3][i2] = 0.0005f;
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
  
  // use 3 points to fit a parabolic curve and find its peak
  private float parabolicPeak(float z, float um, float ui, float up) {
    float a = um-up;
    float b = 2.0f*(um+up)-4.0f*ui;
    return (z+a/b);
  }

  ///////////////////////////////////////////////////////////////////////////
  // private
  private float[] _k1,_k2,_k3;
  private ArrayList<int[]> _extendInd = new ArrayList<int[]>();
  private ArrayList<float[][]> _extendC = new ArrayList<float[][]>();
  private float _weight = 0.5f; // -(fxx+fyy)+weight*(fxxxx+2fxxyy+fyyyy)=px+qy
  private float _sigma1 = 6.0f; // precon smoothing extent for 1st dim
  private float _sigma2 = 6.0f; // precon smoothing extent for 2nd dim
  private float _small = 0.01f; // stop CG iterations if residuals small
  private int _niter = 200; // maximum number of CG iterations
}

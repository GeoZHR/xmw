package pik;

import java.util.*;
import edu.mines.jtk.dsp.*;
import edu.mines.jtk.util.*;
import static edu.mines.jtk.util.ArrayMath.*;

import vec.*;
import util.*;



/**
 * Fault enhancing with ray tracing
 * @author Xinming Wu and Sergey Fomel
 * @version 2016.09.29
 */

public class FaultEnhance {
 
  public FaultEnhance(int gate, float an) {
    _gate = gate;
    _an = an;
    _prob = new float[gate*2-1];
  }

  public float[][] applyForEnhanceX(
    int d1, int d2, float sigma, int[][] seeds, float[][] fx) 
  {
    int n2 = fx.length;
    int n1 = fx[0].length;
    int ns = seeds[0].length;
    float[][] gx = new float[n2][n1];
    float[][] wx = applyForWeight(fx);
    RecursiveGaussianFilterP rgf = new RecursiveGaussianFilterP(sigma);
    for (int is=0; is<ns; ++is) {
      int k1 = seeds[0][is];
      int k2 = seeds[1][is];
      int b1 = max(k2-d2,0);
      int e1 = min(k2+d2+1,n2);

      int bb2 = max(k1-d1,0);
      int eb2 = k1+1;
      int bf2 = eb2;
      int ef2 = min(k1+d1+1,n1);
      int m1 = ef2-bb2; 
      float[] fp = new float[m1];
      float[] pb = backwardPick(k2, b1, bb2, e1, eb2, wx);
      float[] pf =  forwardPick(k2, b1, bf2, e1, ef2, wx);
      for (int i1=bb2; i1<ef2; ++i1) {
        int i2=0;
        if (i1<eb2) {i2 = round(pb[i1-bb2]+b1);}
        else        {i2 = round(pf[i1-bf2]+b1);}
        i2 = max(i2,0);
        i2 = min(i2,n2-1);
        fp[i1-bb2] = fx[i2][i1];
      }
      rgf.apply0(fp,fp);
      for (int i1=bb2; i1<ef2; ++i1) {
        int i2=0;
        if (i1<eb2) {i2 = round(pb[i1-bb2]+b1);}
        else        {i2 = round(pf[i1-bf2]+b1);}
        i2 = max(i2,0);
        i2 = min(i2,n2-1);
        if (abs(i1-eb2)>10) {
          float fpi = fp[i1-bb2];
          gx[i2][i1]+=fpi;
        }

      }
    }
    return gx;
  }


  public float[][] applyForEnhance(float sigma, int[][] seeds, float[][] fx) {
    int n2 = fx.length;
    int n1 = fx[0].length;
    int ns = seeds[0].length;
    float[][] gx = new float[n2][n1];
    float[][] wx = applyForWeight(fx);
    RecursiveGaussianFilterP rgf = new RecursiveGaussianFilterP(sigma);
    for (int is=0; is<ns; ++is) {
      int k1 = seeds[0][is];
      int k2 = seeds[1][is];
      if (k1==0) {
        float[] pt = forwardPick(k2,wx);
        float[] fp = new float[n1];
        for (int i1=0; i1<n1; ++i1) {
          int i2 = round(pt[i1]);
          i2 = max(i2,0);
          i2 = min(i2,n2-1);
          fp[i1] = fx[i2][i1];
        }
        rgf.apply0(fp,fp);
        for (int i1=0; i1<n1; ++i1) {
          if (abs(i1-k1)>20) {continue;}
          int i2 = round(pt[i1]);
          i2 = max(i2,0);
          i2 = min(i2,n2-1);
          float fpi = fp[i1];
          float gxi = gx[i2][i1];
          if (fpi>gxi) {gx[i2][i1] = fpi;}
        }
      }else if (k1==n1-1) {
        float[] pt = backwardPick(k2,wx);
        float[] fp = new float[n1];
        for (int i1=0; i1<n1; ++i1) {
          int i2 = round(pt[i1]);
          i2 = max(i2,0);
          i2 = min(i2,n2-1);
          fp[i1] = fx[i2][i1];
        }
        rgf.apply0(fp,fp);
        for (int i1=0; i1<n1; ++i1) {
          if (abs(i1-k1)>20) {continue;}
          int i2 = round(pt[i1]);
          i2 = max(i2,0);
          i2 = min(i2,n2-1);
          float fpi = fp[i1];
          float gxi = gx[i2][i1];
          if (fpi>gxi) {gx[i2][i1] = fpi;}
        }
      } else {
        int d12 = round((k1+1)/2);
        int b12 = k2-d12; b12 = max(0,b12);
        int e12 = k2+d12; e12 = min(n2-1,e12);
        int d22 = round((n1-k1-1)/2);
        int b22 = k2-d22; b22 = max(0,b22);
        int e22 = k2+d22; e22 = min(n2-1,e22);
        float[][] wx1 = copy(e12-b12+1,k1+1,b12,0,wx); 
        float[][] wx2 = copy(e22-b22+1,n1-k1-1,b22,k1,wx); 
        float[] p1 = backwardPick(k2-b12,wx1);
        float[] p2 = forwardPick(k2-b22,wx2);
        float[] fp = new float[n1];
        for (int i1=0; i1<n1; ++i1) {
          if (abs(i1-k1)>20) {continue;}
          int i2=0;
          if (i1<=k1) {
            i2 = round(p1[i1]+b12);
          } else {
            i2 = round(p2[i1-k1-1]+b22);
          }
          i2 = max(i2,0);
          i2 = min(i2,n2-1);
          fp[i1] = fx[i2][i1];
        }
        rgf.apply0(fp,fp);
        for (int i1=0; i1<n1; ++i1) {
          int i2 = 0;
          if (i1<=k1) {
            i2 = round(p1[i1]+b12);
          } else {
            i2 = round(p2[i1-k1-1]+b22);
          }
          i2 = max(i2,0);
          i2 = min(i2,n2-1);
          float fpi = fp[i1];
          float gxi = gx[i2][i1];
          if (fpi>gxi) {gx[i2][i1] = fpi;}
        }
      }
    }
    return gx;
  }

  public float[][] thin(float fm, float[][] fx) {
    int n2 = fx.length;
    int n1 = fx[0].length;
    float[][] ft = new float[n2][n1];
    float[][] fs = new float[n2][n1];
    RecursiveGaussianFilterP rgf = new RecursiveGaussianFilterP(1);
    rgf.apply00(fx,fs);
    for (int i2=1;i2<n2-1;i2++) {
    for (int i1=0;i1<n1  ;i1++) {
      float fsi = fs[i2][i1];
      float fsm = fs[i2-1][i1];
      float fsp = fs[i2+1][i1];
      if (fsm<fsi&&fsp<fsi&&fsi>fm) {
        ft[i2  ][i1] = fsi;
        ft[i2-1][i1] = fsm;
        ft[i2+1][i1] = fsp;
      }
    }}
    return ft;
  }


  public int[][] findSeeds(int d1, int d2, float fm, 
    float[][] fx, float[][] ft) 
  {
    int n2 = fx.length;
    int n1 = fx[0].length;
    float[][] fs = new float[n2][n1];
    ArrayList<Integer> k1 = new ArrayList<Integer>();
    ArrayList<Integer> k2 = new ArrayList<Integer>();
    RecursiveGaussianFilterP rgf = new RecursiveGaussianFilterP(1);
    rgf.apply00(fx,fs);
    for (int i2=5;i2<n2-5;i2+=d2) {
    for (int i1=5;i1<n1-5;i1+=d1) {
      float fsi = fs[i2][i1];
      float fsm = fs[i2-1][i1];
      float fsp = fs[i2+1][i1];
      if (fsm<fsi&&fsp<fsi&&fsi>fm) {
        k1.add(i1);
        k1.add(i1);
        k1.add(i1);
        k2.add(i2);
        k2.add(i2-1);
        k2.add(i2+1);
        ft[i2][i1] = fsi;
      }
    }}
    int ns = k1.size();
    int[][] seeds = new int[2][ns];
    for (int is=0; is<ns; ++is) {
      seeds[0][is] = k1.get(is); 
      seeds[1][is] = k2.get(is); 
    }
    return seeds;
  }

  public float[][] applyForWeight(float[][] vel) {
    int n2 = vel.length;
    int n1 = vel[0].length;
    float[][] w = new float[n1][n2];
    for (int i2=0; i2<n2; ++i2) {
    for (int i1=0; i1<n1; ++i1) {
      w[i1][i2]  = exp(-vel[i2][i1]);
    }}
    return w;
  }

  public float[][] applyForWeightX(float[][] vel) {
    int n2 = vel.length;
    int n1 = vel[0].length;
    float[][] w = new float[n1][n2];
    for (int i2=0; i2<n2; ++i2) {
    for (int i1=0; i1<n1; ++i1) {
      w[i2][i1]  = exp(-vel[i2][i1]);
    }}
    return w;
  }


  public float[] applyForPath(
    int i0, float sig, float[][] vel) {
    int n1 = vel[0].length;
    float[][] wx = applyForWeight(vel);
    float[] p1 = forwardPick(i0,wx);
    float[] p2 = backwardPick(round(p1[n1-1]),wx);
    float[] b = new float[n1];
    float[] r = new float[n1];
    float[] ws = new float[n1];
    makeRhsWeights(p1,p2,vel,b,ws);
    VecArrayFloat1 vb = new VecArrayFloat1(b);
    VecArrayFloat1 vr = new VecArrayFloat1(r);
    Smoother1 smoother1 = new Smoother1(sig);
    A1 a1 = new A1(smoother1,ws);
    CgSolver cs = new CgSolver(0.001,200);
    smoother1.applyTranspose(b);
    cs.solve(a1,vb,vr);
    return r;
  }

  public float[][][] applyForSurfaceInline(
    float sig1, float sig2, final float[][][] vel) 
  {
    int n3 = vel.length;
    int n2 = vel[0].length;
    int n1 = vel[0][0].length;
    float[][] p1 = new float[n2][n3];
    float[][] w1 = new float[n3][n1];
    for (int i2=0; i2<n2; ++i2) {
      for (int i3=0; i3<n3; ++i3) {
      for (int i1=0; i1<n1; ++i1) {
        w1[i3][i1] = exp(-vel[i3][i2][i1]);
      }}
      p1[i2] = forwardPick(n1-1,w1);
    }
    float[][] b = new float[n3][n2];
    float[][] r = new float[n3][n2];
    float[][] ws = new float[n3][n2];
    makeRhsWeightsInline(p1,vel,b,ws);
    VecArrayFloat2 vb = new VecArrayFloat2(b);
    VecArrayFloat2 vr = new VecArrayFloat2(r);
    Smoother2 smoother2 = new Smoother2(sig1,sig2);
    A2 a2 = new A2(smoother2,ws);
    CgSolver cs = new CgSolver(0.001,200);
    smoother2.applyTranspose(b);
    cs.solve(a2,vb,vr);
    float[][] p1t = new float[n3][n2];
    for (int i3=0; i3<n3; ++i3) {
    for (int i2=0; i2<n2; ++i2) {
      p1t[i3][i2] = p1[i2][i3];
    }}
    return new float[][][]{p1t,r};
  }


  public float[][][] applyForSurface(
    float sig1, float sig2, final float[][][] vel) 
  {
    int n3 = vel.length;
    int n2 = vel[0].length;
    int n1 = vel[0][0].length;
    float[][] p1 = new float[n2][n3];
    float[][] w1 = new float[n3][n1];
    for (int i2=0; i2<n2; ++i2) {
      for (int i3=0; i3<n3; ++i3) {
      for (int i1=0; i1<n1; ++i1) {
        w1[i3][i1] = exp(-vel[i3][i2][i1]);
      }}
      p1[i2] = forwardPick(n1-1,w1);
    }
    float[][] p2 = new float[n3][n2];
    float[][] w2 = new float[n2][n1];
    for (int i3=0; i3<n3; ++i3) {
      for (int i2=0; i2<n2; ++i2) {
      for (int i1=0; i1<n1; ++i1) {
        w2[i2][i1] = exp(-vel[i3][i2][i1]);
      }}
      int i0 = round(p1[0][i3]);
      i0 = min(i0,n1-1);
      i0 = max(i0,0);
      p2[i3] = forwardPick(i0,w2);
    }
    float[][] b = new float[n3][n2];
    float[][] r = new float[n3][n2];
    float[][] ws = new float[n3][n2];
    makeRhsWeights(p1,p2,vel,b,ws);
    VecArrayFloat2 vb = new VecArrayFloat2(b);
    VecArrayFloat2 vr = new VecArrayFloat2(r);
    Smoother2 smoother2 = new Smoother2(sig1,sig2);
    A2 a2 = new A2(smoother2,ws);
    CgSolver cs = new CgSolver(0.001,200);
    smoother2.applyTranspose(b);
    cs.solve(a2,vb,vr);
    float[][] p1t = new float[n3][n2];
    for (int i3=0; i3<n3; ++i3) {
    for (int i2=0; i2<n2; ++i2) {
      p1t[i3][i2] = p1[i2][i3];
    }}
    return new float[][][]{p1t,p2,r};
  }

  public float[] backwardPick(int i0, float[][] wx) {
    int n2 = wx.length;
    int n1 = wx[0].length;
    float[] p = new float[n2];
    float[][] tt = new float[n2][n1];
    float[][] what = new float[n2][n1];
    float[] prev = new float[n1];
    float[] next = new float[n1];
    float[] dist = new float[n1];
    for (int i1=0; i1<n1; ++i1) {
      dist[i1] = sqrt(i1*i1+_an*_an); 
    }
	  for (int i1=0; i1<n1; i1++) {
	    float wi = 0.5f*(wx[n2-1][i1]+wx[n2-1][i0]);
	    tt[n2-1][i1] = abs(i1-i0)*wi;
	  }
    
    for (int i1=0; i1<n1; i1++) {
	    float wi = 0.5f*(wx[n2-2][i1]+wx[n2-1][i0]);
	    prev[i1] = dist[abs(i1-i0)]*wi;
	    what[n2-2][i1] = i0;
	    tt[n2-2][i1]=prev[i1];
    }

    for (int i2=n2-3; i2>=0; i2--) {
	    for (int i1=0; i1<n1; i1++) {
	      float wi = wx[i2][i1];
	      int ib = max(i1-_gate,-1);
	      int ie = min(i1+_gate,n1);
	      float c = FLT_MAX;
	      int ic = -1;
	      for (int i=ib+1; i<ie; i++) {
		      float w2 = 0.5f*(wi+wx[i2+1][i]);
		      float d = dist[abs(i1-i)]*w2+prev[i];
		      int it = i-ib-1;
		      if (d < c) {
		        c =	d;
		        ic = it;
		      }
		      _prob[it]=d;
	      }
        float[] vs = find_minimum(ic,ie-ib-1,ib+1,c,what[i2][i1]);
	      next[i1]= vs[0];
        what[i2][i1] = vs[1];
	    }
	    for (int i1=0; i1<n1; i1++) {
	      prev[i1]=next[i1];
	      tt[i2][i1]=prev[i1];
	    }
    }
    forwardTrack(p, next, what);
    return p;
  }

  public float[] backwardPick(int i0, int b1, int b2, int e1, int e2, float[][] wx) {
    int n2 = e2-b2;
    int n1 = e1-b1;
    float[] p = new float[n2];
    float[][] tt = new float[n2][n1];
    float[][] what = new float[n2][n1];
    float[] prev = new float[n1];
    float[] next = new float[n1];
    float[] dist = new float[n1];
    for (int k1=0; k1<n1; ++k1) {
      dist[k1] = sqrt(k1*k1+_an*_an); 
    }
	  for (int i1=b1; i1<e1; i1++) {
      int k1 = i1-b1;
	    float wi = 0.5f*(wx[e2-1][i1]+wx[e2-1][i0]);
	    tt[n2-1][k1] = abs(i1-i0)*wi;
	  }
    for (int i1=b1; i1<e1; i1++) {
      int k1 = i1-b1;
	    float wi = 0.5f*(wx[e2-2][i1]+wx[e2-1][i0]);
	    prev[k1] = dist[abs(i1-i0)]*wi;
	    what[n2-2][k1] = i0-b1;
	    tt[n2-2][k1]=prev[k1];
    }
    for (int i2=e2-3; i2>=b2; i2--) {
      int k2 = i2-b2;
	    for (int i1=b1; i1<e1; i1++) {
        int k1 = i1-b1;
	      float wi = wx[i2][i1];
	      int kb = max(k1-_gate,-1);
	      int ke = min(k1+_gate,n1);
	      float c = FLT_MAX;
	      int kc = -1;
	      for (int k=kb+1; k<ke; k++) {
		      float w2 = 0.5f*(wi+wx[i2+1][k+b1]);
		      float d = dist[abs(k1-k)]*w2+prev[k];
		      int kt = k-kb-1;
		      if (d < c) {
		        c =	d;
		        kc = kt;
		      }
		      _prob[kt]=d;
	      }
        float[] vs = find_minimum(kc,ke-kb-1,kb+1,c,what[k2][k1]);
	      next[k1]= vs[0];
        what[k2][k1] = vs[1];
	    }
	    for (int k1=0; k1<n1; k1++) {
	      prev[k1]=next[k1];
	      tt[k2][k1]=prev[k1];
	    }
    }
    forwardTrack(p, next, what);
    return p;
  }



  public float[] backwardPick(int i0, float[][] wx, float[][] tx) {
    int n2 = wx.length;
    int n1 = wx[0].length;
    float[] p = new float[n2];
    float[][] tt = new float[n2][n1];
    float[][] what = new float[n2][n1];
    float[] prev = new float[n1];
    float[] next = new float[n1];
    float[] dist = new float[n1];
    for (int i1=0; i1<n1; ++i1) {
      dist[i1] = sqrt(i1*i1+_an*_an); 
    }
	  for (int i1=0; i1<n1; i1++) {
	    float wi = 0.5f*(wx[n2-1][i1]+wx[n2-1][i0]);
	    tt[n2-1][i1] = abs(i1-i0)*wi;
	  }
    
    for (int i1=0; i1<n1; i1++) {
	    float wi = 0.5f*(wx[n2-2][i1]+wx[n2-1][i0]);
	    prev[i1] = dist[abs(i1-i0)]*wi;
	    what[n2-2][i1] = i0;
	    tt[n2-2][i1]=prev[i1];
    }

    for (int i2=n2-3; i2>=0; i2--) {
	    for (int i1=0; i1<n1; i1++) {
	      float wi = wx[i2][i1];
	      int ib = max(i1-_gate,-1);
	      int ie = min(i1+_gate,n1);
	      float c = FLT_MAX;
	      int ic = -1;
	      for (int i=ib+1; i<ie; i++) {
		      float w2 = 0.5f*(wi+wx[i2+1][i]);
		      float d = dist[abs(i1-i)]*w2+prev[i];
		      int it = i-ib-1;
		      if (d < c) {
		        c =	d;
		        ic = it;
		      }
		      _prob[it]=d;
	      }
        float[] vs = find_minimum(ic,ie-ib-1,ib+1,c,what[i2][i1]);
	      next[i1]= vs[0];
        what[i2][i1] = vs[1];
	    }
	    for (int i1=0; i1<n1; i1++) {
	      prev[i1]=next[i1];
	      tt[i2][i1]=prev[i1];
	    }
    }
    forwardTrack(p, next, what);
    for (int i2=0; i2<n2;++i2) {
    for (int i1=0; i1<n1;++i1) {
      tx[i1][i2]  = tt[i2][i1];
    }}
    return p;
  }

  public float[] forwardPick(int i0, float[][] wx) {
    int n2 = wx.length;
    int n1 = wx[0].length;
    float[] p = new float[n2];
    float[][] tt = new float[n2][n1];
    float[][] what = new float[n2][n1];
    float[] prev = new float[n1];
    float[] next = new float[n1];
    float[] dist = new float[n1];
    for (int i1=0; i1<n1; ++i1) {
      dist[i1] = sqrt(i1*i1+_an*_an); 
    }
	  for (int i1=0; i1<n1; i1++) {
	    float wi = 0.5f*(wx[0][i1]+wx[0][i0]);
	    tt[0][i1] = abs(i1-i0)*wi;
	  }
    
    for (int i1=0; i1<n1; i1++) {
	    float wi = 0.5f*(wx[1][i1]+wx[0][i0]);
	    prev[i1] = dist[abs(i1-i0)]*wi;
	    what[1][i1] = i0;
	    tt[1][i1]=prev[i1];
    }

    for (int i2=2; i2<n2; i2++) {
	    for (int i1=0; i1<n1; i1++) {
	      float wi = wx[i2][i1];
	      int ib = max(i1-_gate,-1);
	      int ie = min(i1+_gate,n1);
	      int ic = -1;
	      float c = FLT_MAX;
	      for (int i=ib+1; i<ie; i++) {
		      float w2 = 0.5f*(wi+wx[i2-1][i]);
		      float d = dist[abs(i1-i)]*w2+prev[i];
		      int it = i-ib-1;
		      if (d < c) {
		        c =	d;
		        ic = it;
		      }
		      _prob[it]=d;
	      }
        float[] vs = find_minimum(ic,ie-ib-1,ib+1,c,what[i2][i1]);
	      next[i1]= vs[0];
        what[i2][i1] = vs[1];
	    }
	    for (int i1=0; i1<n1; i1++) {
	      prev[i1]=next[i1];
	      tt[i2][i1]=prev[i1];
	    }
    }
    backwardTrack(p, next, what);
    return p;
  }

  public float[] forwardPick(int i0, int b1, int b2, int e1, int e2, float[][] wx) {
    int n2 = e2-b2;
    int n1 = e1-b1;
    float[] p = new float[n2];
    float[][] tt = new float[n2][n1];
    float[][] what = new float[n2][n1];
    float[] prev = new float[n1];
    float[] next = new float[n1];
    float[] dist = new float[n1];
    for (int k1=0; k1<n1; ++k1) {
      dist[k1] = sqrt(k1*k1+_an*_an); 
    }
	  for (int i1=b1; i1<e1; i1++) {
	    float wi = 0.5f*(wx[b2][i1]+wx[b2][i0]);
	    tt[0][i1-b1] = abs(i1-i0)*wi;
	  }
    
    for (int i1=b1; i1<e1; i1++) {
      int k1 = i1-b1;
	    float wi = 0.5f*(wx[b2+1][i1]+wx[b2][i0]);
	    prev[k1] = dist[abs(i1-i0)]*wi;
	    what[1][k1] = i0-b1;
	    tt[1][k1]=prev[k1];
    }

    for (int i2=b2+2; i2<e2; i2++) {
      int k2 = i2-b2; 
	    for (int i1=b1; i1<e1; i1++) {
        int k1 = i1-b1; 
	      float wi = wx[i2][i1];
	      int kb = max(k1-_gate,-1);
	      int ke = min(k1+_gate,n1);
	      int kc = -1;
	      float c = FLT_MAX;
	      for (int k=kb+1; k<ke; k++) {
		      float w2 = 0.5f*(wi+wx[i2-1][k+b1]);
		      float d = dist[abs(k1-k)]*w2+prev[k];
		      int kt = k-kb-1;
		      if (d < c) {
		        c =	d;
		        kc = kt;
		      }
		      _prob[kt]=d;
	      }
        float[] vs = find_minimum(kc,ke-kb-1,kb+1,c,what[k2][k1]);
	      next[k1]= vs[0];
        what[k2][k1] = vs[1];
	    }
	    for (int k1=0; k1<n1; k1++) {
	      prev[k1]=next[k1];
	      tt[k2][k1]=prev[k1];
	    }
    }
    backwardTrack(p, next, what);
    return p;
  }


  public float[] forwardPick(int i0, float[][] wx, float[][] tx) {
    int n2 = wx.length;
    int n1 = wx[0].length;
    float[] p = new float[n2];
    float[][] tt = new float[n2][n1];
    float[][] what = new float[n2][n1];
    float[] prev = new float[n1];
    float[] next = new float[n1];
    float[] dist = new float[n1];
    for (int i1=0; i1<n1; ++i1) {
      dist[i1] = sqrt(i1*i1+_an*_an); 
    }
	  for (int i1=0; i1<n1; i1++) {
	    float wi = 0.5f*(wx[0][i1]+wx[0][i0]);
	    tt[0][i1] = abs(i1-i0)*wi;
	  }
    
    for (int i1=0; i1<n1; i1++) {
	    float wi = 0.5f*(wx[1][i1]+wx[0][i0]);
	    prev[i1] = dist[abs(i1-i0)]*wi;
	    what[1][i1] = i0;
	    tt[1][i1]=prev[i1];
    }

    for (int i2=2; i2<n2; i2++) {
	    for (int i1=0; i1<n1; i1++) {
	      float wi = wx[i2][i1];
	      int ib = max(i1-_gate,-1);
	      int ie = min(i1+_gate,n1);
	      float c = FLT_MAX;
	      int ic = -1;
	      for (int i=ib+1; i<ie; i++) {
		      float w2 = 0.5f*(wi+wx[i2-1][i]);
		      float d = dist[abs(i1-i)]*w2+prev[i];
		      int it = i-ib-1;
		      if (d < c) {
		        c =	d;
		        ic = it;
		      }
		      _prob[it]=d;
	      }
        float[] vs = find_minimum(ic,ie-ib-1,ib+1,c,what[i2][i1]);
	      next[i1]= vs[0];
        what[i2][i1] = vs[1];
	    }
	    for (int i1=0; i1<n1; i1++) {
	      prev[i1]=next[i1];
	      tt[i2][i1]=prev[i1];
	    }
    }
    backwardTrack(p, next, what);
    for (int i2=0; i2<n2;++i2) {
    for (int i1=0; i1<n1;++i1) {
      tx[i1][i2]  = tt[i2][i1];
    }}
    return p;
  }



  private float[] find_minimum(
    int ic, int nc, int jc, float c, float pick)
  {
    float fm, f0, fp, a, b;
    if (0==ic) {
	    ic++;
	    fm=c;
	    f0=_prob[ic];
	    fp=_prob[ic+1];
    } else if (nc-1==ic) {
	    ic--;
	    fm=_prob[ic-1];
	    f0=_prob[ic];
	    fp=c;
    } else {
	    fm=_prob[ic-1];
	    f0=c;
	    fp=_prob[ic+1];
    }

    ic += jc;
    a = fm+fp-2f*f0;
    if (a <= 0.) { /* no minimum */
	    if (fm < f0 && fm < fp) {
	      pick = ic-1;
	      return new float[]{fm,pick};
	    } 
	    if (fp < f0 && fp < fm) {
	      pick = ic+1;
	      return new float[]{fp,pick};
	    } 
	    pick = ic;
	    return new float[]{f0,pick};
    }

    b = 0.5f*(fm-fp);
    a = b/a;
    if (a > 1.) {
	    pick = ic+1;
	    return new float[]{fp,pick};
    }

    if (a < -1.) {
	    pick = ic-1;
	    return new float[]{fm,pick};
    }

    if (f0 < 0.5*b*a) {
	    pick = ic;
	    return new float[]{f0,pick};
    }

    f0 -= 0.5*b*a;
    pick=ic+a;
    return new float[]{f0,pick};
  }

  private void forwardTrack(
    float[] path, float[] next, float[][] what)
  {
    float c, d, fc;
    int n1 = next.length;
    int n2 = path.length;
    c = FLT_MAX;
    fc = 0;
    /* minimum at the bottom */
    for (int i1=0; i1 < n1; i1++) {
	    d = next[i1];
	    if (d < c) {
	      c = d;
	      fc = i1;
	    }
    }
    /* coming up */
    for (int i2=0; i2<n2; i2++) {
	    path[i2]=fc;
	    fc = interpolate(fc,i2,what);
    }
  }


  private void backwardTrack(
    float[] path, float[] next, float[][] what)
  {
    float c, d, fc;
    int n1 = next.length;
    int n2 = path.length;
    c = FLT_MAX;
    fc = 0;
    /* minimum at the bottom */
    for (int i1=0; i1 < n1; i1++) {
	    d = next[i1];
	    if (d < c) {
	      c = d;
	      fc = i1;
	    }
    }
    /* coming up */
    for (int i2=n2-1; i2 >= 0; i2--) {
	    path[i2]=fc;
	    fc = interpolate(fc,i2,what);
    }
  }

  private float interpolate(float fc, int i2, float[][] what) {
    int n1 = what[0].length;
    int ic = round(fc-0.5f);
    fc -= ic;
    if (n1-1 <= ic) return what[i2][n1-1];
    if (0 > ic) return what[i2][0];
    fc = what[i2][ic]*(1f-fc)+what[i2][ic+1]*fc;
    return fc;
  }

  private void makeRhsWeightsInline(
    float[][] p23, float[][][] vel, float[][] b, float[][] ws) 
  {
    int n3 = vel.length;
    int n2 = vel[0].length;
    int n1 = vel[0][0].length;
    for (int i3=0; i3<n3; ++i3) {
    for (int i2=0; i2<n2; ++i2) {
      int k23 = round(p23[i2][i3]);
      k23 = max(0,k23);
      k23 = min(n1-1,k23);
      float w23i = vel[i3][i2][k23];
      w23i *= w23i;
      ws[i3][i2] = w23i;
      b[i3][i2] = p23[i2][i3]*w23i;
    }}
  }


  private void makeRhsWeights(
    float[][] p23, float[][] p32, float[][][] vel, float[][] b, float[][] ws) 
  {
    int n3 = vel.length;
    int n2 = vel[0].length;
    int n1 = vel[0][0].length;
    for (int i3=0; i3<n3; ++i3) {
    for (int i2=0; i2<n2; ++i2) {
      int k23 = round(p23[i2][i3]);
      int k32 = round(p32[i3][i2]);
      k23 = max(0,k23);
      k23 = min(n1-1,k23);
      k32 = max(0,k32);
      k32 = min(n1-1,k32);
      float w23i = vel[i3][i2][k23];
      float w32i = vel[i3][i2][k32];
      w23i *= w23i;
      w32i *= w32i;
      ws[i3][i2] = w23i+w32i;
      b[i3][i2] = p23[i2][i3]*w23i+p32[i3][i2]*w32i;
    }}
  }


  private void makeRhsWeights(
    float[] p1, float[] p2, float[][] vel, float[] b, float[] ws) 
  {
    int n2 = vel.length;
    int n1 = vel[0].length;
    for (int i1=0; i1<n1; ++i1) {
      int k12 = round(p1[i1]);
      int k22 = round(p2[i1]);
      k12 = max(0,k12);
      k12 = min(n2-1,k12);
      k22 = max(0,k22);
      k22 = min(n2-1,k22);
      float w1i = vel[k12][i1];
      float w2i = vel[k22][i1];
      w1i *= w1i;
      w2i *= w2i;
      ws[i1] = w1i+w2i;
      b[i1] = p1[i1]*w1i+p2[i1]*w2i;
    }
  }

  // Conjugate-gradient operators.
  private static class A1 implements CgSolver.A {
    A1(Smoother1 s1, float[] wp) 
    {
      _s1 = s1;
      _wp = wp;
      float n1 = wp.length;
      _sc = sum(wp)/(n1);
    }
    public void apply(Vec vx, Vec vy) {
      VecArrayFloat1 v1x = (VecArrayFloat1)vx;
      VecArrayFloat1 v1y = (VecArrayFloat1)vy;
      float[] x = v1x.getArray();
      float[] y = v1y.getArray();
      float[] z = copy(x);
      v1y.zero();
      _s1.apply(z);
      addAndScale(-_sc,z,y);
      applyLhs(_wp,z,y);
      _s1.applyTranspose(y);
      addAndScale( _sc,x,y);
    }
    private float _sc;
    private Smoother1 _s1;
    private float[] _wp;
  }


  // Conjugate-gradient operators.
  private static class A2 implements CgSolver.A {
    A2(Smoother2 s2, float[][] wp) 
    {
      _s2 = s2;
      _wp = wp;
      float n2 = wp.length;
      float n1 = wp[0].length;
      _sc = 4f*sum(wp)/(n1*n2);
    }
    public void apply(Vec vx, Vec vy) {
      VecArrayFloat2 v2x = (VecArrayFloat2)vx;
      VecArrayFloat2 v2y = (VecArrayFloat2)vy;
      float[][] x = v2x.getArray();
      float[][] y = v2y.getArray();
      float[][] z = copy(x);
      v2y.zero();
      _s2.apply(z);
      addAndScale(-_sc,z,y);
      applyLhs(_wp,z,y);
      _s2.applyTranspose(y);
      addAndScale( _sc,x,y);
    }
    private float _sc;
    private Smoother2 _s2;
    private float[][] _wp;
  }

  private static void applyLhs(float[] wp, float[] x, float[] y) {
    int n1 = wp.length;
    for (int i1=0; i1<n1; ++i1)
      y[i1] += wp[i1]*x[i1];
  }


  private static void applyLhs(float[][] wp, float[][] x, float[][] y) {
    int n2 = wp.length;
    for (int i2=0; i2<n2; ++i2)
      applyLhs(wp[i2],x[i2],y[i2]);
  }

  private static void addAndScale(float sc, float[][] x, float[][] y) {
    int n2 = x.length;
    for (int i2=0; i2<n2; ++i2) {
      addAndScale(sc,x[i2],y[i2]);
    }
  }

  private static void addAndScale(float sc, float[] x, float[] y) {
    int n1 = x.length;
    for (int i1=0; i1<n1; ++i1) {
      y[i1] += sc*x[i1];
    }
  }


  // Smoother used as a preconditioner. After smoothing, enforces zero-shift
  // boundary conditions at top and bottom.
  private static class Smoother2 {
    public Smoother2(float sigma1, float sigma2) {
      _sigma1 = sigma1;
      _sigma2 = sigma2;
    }
    public void apply(float[][] x) {
      smooth2(_sigma2,x);
      smooth1(_sigma1,x);
    }
    public void applyTranspose(float[][] x) {
      smooth1(_sigma1,x);
      smooth2(_sigma2,x);
    }
    private float _sigma1,_sigma2;
  }


  private static class Smoother1 {
    public Smoother1(float sigma) {
      _sigma = sigma;
    }
    public void apply(float[] x) {
      smooth1(_sigma,x);
    }
    public void applyTranspose(float[] x) {
      smooth1(_sigma,x);
    }
    private float _sigma;
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
  private static void smooth1(float sigma, float[] x) {
    if (sigma<1.0f)
      return;
    RecursiveExponentialFilter.Edges edges =
      RecursiveExponentialFilter.Edges.OUTPUT_ZERO_SLOPE;
    RecursiveExponentialFilter ref = new RecursiveExponentialFilter(sigma);
    ref.setEdges(edges);
    ref.apply(x,x);
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



  ///////////////////////////////////////////////////////////////////////////
  // private
  private int _gate;
  private float _an;
  private float[] _prob;


}

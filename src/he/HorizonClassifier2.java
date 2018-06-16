package he;

import java.util.*;

import edu.mines.jtk.dsp.*;
import edu.mines.jtk.interp.*;
import static edu.mines.jtk.util.ArrayMath.*;
import vec.*;
import util.*;

/**
 * Extract a single seismic horizon curve by iterative 
 * seismic wavefome classification.
 * @author Xinming Wu
 * @version 2017.12.22
 */

public class HorizonClassifier2{

  public HorizonClassifier2 (int w1, int w2) {
    _w1 = w1;
    _w2 = w2;
    _s1 = new float[w1*2+1];
    _s2 = new float[w2*2+1];
    computeGaussianWeights(w1,w2);
  }

  public void setGaussianWeights(float sigma1, float sigma2) {
    computeGaussianWeights(sigma1,sigma2);
  }

  public float[][] initialPick(
    int c1, int c2, int[][] pks, float[][] p2, float[][] fx) 
  {
    int n2 = fx.length;
    float[][] hp = new float[2][n2];
    hp[0][c1] = c1;
    hp[1][c1] = 1f;
    // pick right
    for (int i2=c2+1; i2<n2; ++i2) {
      int p1 = 0;
      float dm = 1f;
      int np = pks[i2].length;
      for (int ip=0; ip<np; ip++) {
        int i1 = pks[i2][ip];
        float di = getDistance(c1,c2,i1,i2,fx);
        if(di<dm) {dm=di; p1 = i1;}
      }
      hp[0][i2] = p1;
      hp[1][i2] = 1f-dm;
    }
    // pick left
    for (int i2=c2-1; i2>0; --i2) {
      int p1 = 0;
      float dm = 1f;
      int np = pks[i2].length;
      for (int ip=0; ip<np; ip++) {
        int i1 = pks[i2][ip];
        float di = getDistance(c1,c2,i1,i2,fx);
        if(di<dm) {dm=di; p1=i1;}
      }
      hp[0][i2] = p1;
      hp[1][i2] = 1f-dm;
    }
    for (int iter=0; iter<1; iter++)
      hp = updatePick(c1,c2,pks,hp,fx);
    return hp;
  }

  public float[][] updatePick(
    int c1, int c2, int[][] pks, float[][] hp, float[][] fx) 
  {
    int n2 = fx.length;
    // pick right
    for (int i2=0; i2<n2; ++i2) {
      int p1 = 0;
      float dm = 1f;
      int np = pks[i2].length;
      for (int ip=0; ip<np; ip++) {
        int i1 = pks[i2][ip];
        float si = 5f;
        float di = getDistance(c1,c2,i1,i2,fx)*si;
        for (int k2=0; k2<n2; k2+=8) {
          float sp = hp[1][k2];
          int k1 = round(hp[0][k2]);
          if(sp>0.7f) {
            si += sp;
            di += getDistance(k1,k2,i1,i2,fx)*sp;
          }
        }
        di /= si;
        if(di<dm) {dm=di; p1=i1;}
      }
      hp[0][i2] = p1;
      hp[1][i2] = 1f-dm;
    }
    hp[0][c2] = c1;
    hp[1][c2] = 1f;
    return hp;
  }


  private float getDistance(
    int a1, int a2, int b1, int b2, float[][] fx) {
    DynamicWarping dw = new DynamicWarping(-1,1);
    dw.setStrainMax(1.0);
    dw.setShiftSmoothing(10);
    int n1 = fx[0].length;
    int fa1 = a1-_w1; int dfa1 = max(fa1,0)-fa1; fa1 += dfa1;
    int fb1 = b1-_w1; int dfb1 = max(fb1,0)-fb1; fb1 += dfb1;
    int la1 = a1+_w1; int dla1 = la1-min(la1,n1-1); la1 -= dla1;
    int lb1 = b1+_w1; int dlb1 = lb1-min(lb1,n1-1); lb1 -= dlb1;
    int df1 = max(dfa1,dfb1);
    int dl1 = max(dla1,dlb1);
    fa1 = a1-_w1+df1;
    la1 = a1+_w1-dl1;
    fb1 = b1-_w1+df1;
    lb1 = b1+_w1-dl1;
    int m1 = _w1+_w1-dl1-df1+1;
    float[] fa = new float[m1];
    float[] fb = new float[m1];
    for (int i1=0; i1<m1; i1++) {
      fa[i1] = fx[a2][fa1+i1];
      fb[i1] = fx[b2][fb1+i1];
    }
    float[] fs = dw.findShifts(fa,fb);
    float[] ft = dw.applyShifts(fs,fb);
    float ds = 0f;
    float ss = 0f;
    for (int i1=0; i1<m1; i1++) {
      ds += abs(fa[i1]-ft[i1])*_s1[i1+df1];
      ss += _s1[i1+df1];
    }
    return ds/ss;
  }


  private float[] getTrace(int c1, int c2, float[][] fx) {
    int m1 = _w1*2+1;
    float[] fc = fillfloat(FLT_MAX,m1);
    for (int i1=0; i1<m1; ++i1)
      fc[i1] = fx[c2][c1+i1-_w1];
    return fc;
  }


  public int[][] findTroughs(float[][] fx) {
    int n2 = fx.length;
    int n1 = fx[0].length;
    float[][] gx = new float[n2][n1];
    RecursiveGaussianFilter rgf = new RecursiveGaussianFilter(1);
    rgf.apply0X(fx,gx);
    int[][] pks = new int[n2][];
    for (int i2=0; i2<n2  ; ++i2) {
      float[] gx2 = gx[i2];
      ArrayList<Integer> ps = new ArrayList<Integer>();
      for (int i1=1; i1<n1-1; ++i1) {
        int i1m = i1-1;
        int i1p = i1+1;
        float gxi = gx2[i1 ];
        float gxm = gx2[i1m];
        float gxp = gx2[i1p];
        float gum = gxi-gxm;
        float gup = gxp-gxi;
        if(gup*gum<0f) {
          if(gxi<gxm&&gxi<gxp) ps.add(i1);
        }
      }
      int np = ps.size();
      pks[i2] = new int[np];
      for (int ip=0; ip<np; ++ip)
        pks[i2][ip] = ps.get(ip);
    }
    return pks;
  }

  public float[][] getPmap(int n1, int n2, float[][] hps) {
    float[][] mp = new float[n2][n1];
    for (int i2=0; i2<n2; ++i2) {
        mp[i2][round(hps[0][i2])] = hps[1][i2]; 
    }
    return mp;
  }


  public float[][] getPmap(int n1, int n2, int[][] mks) {
    float[][] mp = new float[n2][n1];
    for (int i2=0; i2<n2; ++i2) {
      int[] ps = mks[i2];
      int np = ps.length;
      for (int ip=0; ip<np; ++ip)
        mp[i2][ps[ip]] = 1; 
    }
    return mp;
  }


  public float[][] predict(int c1, int c2, float[][] p2, float[][] fx) {
    int n2 = fx.length;
    int n1 = fx[0].length;
    int m1 = _w1*2+1;
    int m2 = _w2*2+1;
    float[] x1l = new float[m1];
    float[] p2s = new float[m1];
    float[][] fps = new float[m2][m1];
    SincInterpolator si = new SincInterpolator();
    float ps = parabolicShift(fx[c2][c1-1],fx[c2][c1],fx[c2][c1+1]);
    for (int i1=0; i1<m1; ++i1)
      x1l[i1] = c1+i1-_w1+ps;
    float[] x1r = copy(x1l);
    si.interpolate(n1,1,0,fx[c2],m1,x1l,fps[_w2]);
    int k = 0;
    for (int i2=_w2-1; i2>=0; i2--) {
      if(c2-k-1<0) {
        fps[i2] = fps[i2+1];
      } else {
        si.interpolate(n1,1,0,p2[c2-k],m1,x1l,p2s);
        sub(x1l,p2s,x1l);
        si.interpolate(n1,1,0,fx[c2-k-1],m1,x1l,fps[i2]);
        k++;
      }
    }
    k = 0;
    for (int i2=_w2+1; i2<m2; i2++) {
      if(c2+k+1>=n2) {
        fps[i2]=fps[i2-1];
      }else {
        si.interpolate(n1,1,0,p2[c2+k],m1,x1r,p2s);
        add(x1r,p2s,x1r);
        si.interpolate(n1,1,0,fx[c2+k+1],m1,x1r,fps[i2]);
        k++;
      }
    }
    return fps;
  }

  private float pearsonCorrelation(float[] g, float[] f) {
    int n1 = g.length;
    float gf = 0f;
    float ga = 0f;
    float fa = 0f;
    float gs = 0f;
    float fs = 0f;
    for (int i1=0; i1<n1; ++i1) {
      float gi = g[i1];
      float fi = f[i1];
      ga += gi;
      fa += fi;
      gf += gi*fi;
      gs += gi*gi;
      fs += fi*fi;
    }
    ga /= n1;
    fa /= n1;
    float pcn = gf-ga*fa*n1;
    float pcd = n1*sqrt(gs/n1-ga*ga)*sqrt(fs/n1-fa*fa);
    return pcn/pcd;

  }
 
  private float semblance(float[] g, float[] f) {
    float sn = 0f;
    float sd = 0f;
    float ga = 0f;
    float fa = 0f;
    int n1 = g.length;
    for (int i1=0; i1<n1; ++i1) {
      ga += g[i1];
      fa += f[i1];
    }
    ga /= n1;
    fa /= n1;
    for (int i1=0; i1<n1; ++i1) {
      float gi = g[i1]-ga;
      float fi = f[i1]-fa;
      sn += (fi+gi)*(fi+gi)*_s1[i1];
      sd += (fi*fi+gi*gi)*_s1[i1];
    }
    return sn/sd/2f;
  }

  private float semblanceX(float[] g, float[] f) {
    float sn = 0f;
    float sd = 0f;
    float ga = 0f;
    float fa = 0f;
    int n1 = g.length;
    for (int i1=0; i1<n1; ++i1) {
      ga += g[i1];
      fa += f[i1];
    }
    ga /= n1;
    fa /= n1;
    for (int i1=0; i1<n1; ++i1) {
      float gi = g[i1]-ga;
      float fi = f[i1]-fa;
      sn += (fi+gi)*(fi+gi);
      sd += (fi*fi+gi*gi);
    }
    return sn/sd/2f;
  }

  private float distance(float[] g, float[] f) {
    int n1 = g.length;
    int d1 = 10;
    float ds = 0f;
    float ss = 0f;
    for (int i1=0; i1<n1; ++i1) {
      int b1 = i1-d1; b1 = max(b1,0);
      int e1 = i1+d1; e1 = min(e1,n1-1);
      float di = 0f;
      float sk = 0f;
      for (int k1=b1; k1<=e1; k1++) {
        float gk = g[k1];
        float fk = f[k1];
        if(gk!=FLT_MAX&&fk!=FLT_MAX) {
          di += abs(gk-fk);
          sk += 1f;
        }
      }
      if (sk>0f) {
        di /= sk;
        ds += di*_s1[i1];
        ss += _s1[i1];
      }
    }
    return ds/ss;
  }

  //A Shape-based Similarity Measure for Time Series Data 
  //with Ensemble Learning, 2013, proposed by Tetsuya Nakamura
  private float shapeBasedSimilarity(float[] g, float[] f) {
    float ss = 0;
    int n1 = g.length;
    float d = 2f;
    float ds = d*d;
    for (int i1=0; i1<n1; ++i1) {
      float gi = g[i1];
      float fi = f[i1];
      float ga = sqrt(ds+gi*gi);
      float fa = sqrt(ds+fi*fi);
      float si = (ds+gi*fi)/(ga*fa); 
      //if(si>0) 
      ss += si*_s1[i1];
    }
    return ss/n1;
  }

  private float shapeBasedSimilarity(int d, float[] g, float[] f) {
    float ss = 0;
    int n1 = g.length;
    for (int i1=d; i1<n1; i1+=d) {
      float gi=0f;
      float fi=0f;
      for (int k=i1-d;k<i1;k++) {
        gi += g[k];
        fi += f[k];
      }
      float ga = sqrt(1*1+gi*gi);
      float fa = sqrt(1*1+fi*fi);
      float si = (1*1+gi*fi)/(ga*fa); 
      if(si>0) ss += si*_s1[i1];
    }
    return ss/n1;
  }


  private void computeGaussianWeights(float sigma1, float sigma2) {
    float pi = (float)Math.PI;
    float sc1 = 1f/(sigma1*sqrt(2*pi));
    float sc2 = 1f/(sigma2*sqrt(2*pi));
    for (int i1=0; i1<_w1; i1++) {
      float di = i1-_w1;
      float ds = di*di;
      float ss = 1f/(sigma1*sigma1);
      _s1[i1] = sc1*exp(-0.5f*ds*ss);
      _s1[2*_w1-i1] = _s1[i1];
    }
    _s1[_w1] = sc1;
    for (int i2=0; i2<_w2; i2++) {
      float di = i2-_w2;
      float ds = di*di;
      float ss = 1f/(sigma1*sigma1);
      _s2[i2] = sc2*exp(-0.5f*ds*ss);
      _s2[2*_w2-i2] = _s2[i2];
    }
    _s2[_w2] = sc2;
  }

  private float parabolicShift(float fa, float fb, float fc) {
    return 0f;//0.5f*(fa-fc)/(fa+fc-2*fb);
  }

  private float[] _s1;
  private float[] _s2;
  private int _w1;
  private int _w2;
}

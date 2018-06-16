package tjxd;

import edu.mines.jtk.dsp.*;
import edu.mines.jtk.util.*;
import edu.mines.jtk.interp.*;
import static edu.mines.jtk.util.ArrayMath.*;

public class Helper {

  public void fault(int v, int n1, int n2, 
    float[][] p1s, float[][] p2s, float[][] fs) {
    int np = p1s.length;
    for (int ip=0; ip<np; ip++) {
      int nk = p1s[ip].length;
      for (int ik=0; ik<nk-1; ik++) {
        int p1i = round(p1s[ip][ik]);
        int p2i = round(p2s[ip][ik]);
        int p1p = round(p1s[ip][ik+1]);
        int p2p = round(p2s[ip][ik+1]);
        float dk = (float)(p1p-p1i)/(float)(p2p-p2i);
        int d = 1; if(p2p<p2i) d=-1;
        int i2 = p2i;
        while(i2>=p2p) {
          int i1 = round(p1i+dk*(i2-p2i));
          fs[i2][i1] = v;
          i2 +=d;
        }
      }
    }
  }

  public float[][] faultBlocker(float v, float[][] fm, float[][] fx) {
    int b2 = 0;
    int e2 = 0;
    int n2 = fx.length;
    int n1 = fx[0].length;
    for (int i2=0; i2<n2; ++i2) {
      if(fm[i2][0]==v) {b2=i2;break;}
    }
    for (int i2=0; i2<n2; ++i2) {
      if(fm[i2][n1-1]==v) {e2=i2;break;}
    }
    float[][] fxl = copy(fx);
    float[][] fxr = copy(fx);
    for (int i1=0; i1<n1; ++i1) {
      int k2 = 0;
      for (int i2=0; i2<n2; ++i2) {
        if(fm[i2][i1]==v) {k2=i2;break;}
      }
      for (int i2=k2; i2<n2; i2++) {
        fxl[i2][i1] = fx[k2-1][i1];
      }
      for (int i2=0; i2<=k2; i2++) {
        fxr[i2][i1] = fx[k2+1][i1];
      }
    }
    float[][] fxs = new float[n2-e2+b2][n1];
    for (int i2=0; i2<=b2; i2++) {
      fxs[i2] = fxl[i2];
    }
    for (int i2=n2-e2+b2-1; i2>b2; i2--)
      fxs[i2] = fxr[i2+e2-b2];
    return fxs;
  }

  public float[][] unfaultBlocker(float v, float[][] fm, float[][] fx) {
    int b2 = 0;
    int e2 = 0;
    int n2 = fx.length;
    int n1 = fx[0].length;
    for (int i2=0; i2<n2; ++i2) {
      if(fm[i2][0]==v) {b2=i2;break;}
    }
    for (int i2=0; i2<n2; ++i2) {
      if(fm[i2][n1-1]==v) {e2=i2;break;}
    }
    float[][] fxs = new float[n2+e2-b2][n1];
    for (int i1=0; i1<n1; ++i1) {
      int k2 = 0;
      for (int i2=0; i2<n2; ++i2) {
        if(fm[i2][i1]==v) {k2=i2;break;}
      }
      for (int i2=0; i2<=k2; i2++) {
        fxs[i2][i1] = fx[i2][i1];
      }
      for (int i2=k2+1; i2<n2+e2-b2; i2++) {
        fxs[i2][i1] = fx[i2+b2-e2][i1];
      }
    }
    return fxs;
  }

  public float[][] rgtInterpolate(float[][] fm, float[][] hs) {
    int n2 = fm.length;
    int n1 = fm[0].length;
    int nh = hs.length;
    float[][] tx = new float[n2][n1];
    for (int i2=0; i2<n2; ++i2) {
      float xh1 = hs[0][i2];
      for(int i1=0; i1<xh1; i1++) {
        tx[i2][i1] = hs[0][0]-xh1+i1;
      }
      for (int ih=1; ih<nh; ih++) {
        float xhm = hs[ih-1][i2];
        float xhi = hs[ih  ][i2];
        for (int i1=round(xhm); i1<=xhi; i1++) {
          float d1 = xhi-xhm;
          float dt = hs[ih][0 ]-hs[ih-1][0 ];
          float dx = i1-xhm;
          tx[i2][i1] = dx*dt/d1+hs[ih-1][0];
        }
      }
      float xhn = hs[nh-1][i2];
      for(int i1=round(xhn); i1<n1; i1++) {
        tx[i2][i1] = hs[nh-1][0]+i1-xhn;
      }
    }
    return tx;
  }


  public float[][] rgtInterpolate(int n1, int n2, float[][] hs) {
    int nh = hs.length;
    float[][] tx = new float[n2][n1];
    for (int i2=0; i2<n2; ++i2) {
      float xh1 = hs[0][i2];
      for(int i1=0; i1<xh1; i1++) {
        tx[i2][i1] = hs[0][0]-xh1+i1;
      }
      for (int ih=1; ih<nh; ih++) {
        float xhm = hs[ih-1][i2];
        float xhi = hs[ih  ][i2];
        for (int i1=round(xhm); i1<=xhi; i1++) {
          float d1 = xhi-xhm;
          float dt = hs[ih][0 ]-hs[ih-1][0 ];
          float dx = i1-xhm;
          tx[i2][i1] = dx*dt/d1+hs[ih-1][0];
        }
      }
      float xhn = hs[nh-1][i2];
      for(int i1=round(xhn); i1<n1; i1++) {
        tx[i2][i1] = hs[nh-1][0]+i1-xhn;
      }
    }
    return tx;
  }

  public float[][] interpWithRgt(int c2, float[] wl, float[][] tx) {
    int n2 = tx.length;
    int n1 = tx[0].length;
    float[][] vx = new float[n2][n1];
    float[] tc = zerofloat(n1);
    int k = 1;
    tc[0] = tx[c2][0];
    wl[0] = wl[0];
    for (int i1=1; i1<n1; i1++) {
      if(tx[c2][i1]>tc[k-1]) {
        tc[k] = tx[c2][i1];
        wl[k] = wl[i1];
        k++;
      }

    }
    tc = copy(k,0,tc);
    wl = copy(k,0,wl);
    CubicInterpolator ci = new CubicInterpolator(tc,wl);
    for (int i2=0; i2<n2; ++i2) 
      vx[i2] = ci.interpolate(tx[i2]);
    return vx;

  }

  public float[][] maskout(int c2, int d2, int b1, int e1, float[][] fx, float[][] tx) {
    int n2 = tx.length;
    int n1 = tx[0].length;
    float fxmin = min(fx);
    float[][] fy = fillfloat(fxmin-100,n1,n2);
    float tmin = tx[c2][b1];
    float tmax = tx[c2][e1];
    int b2 = c2-d2; b2 = max(0,b2);
    int e2 = c2+d2; e2 = min(n2,e2);
    float[] w = new float[n2*2];
    int h2 = n2;
    w[h2] = 1f;
    RecursiveGaussianFilter rgf = new RecursiveGaussianFilter(100);
    rgf.apply0(w,w);
    w = div(w,max(w));
    for (int i2=b2; i2<e2; ++i2) {
    for (int i1= 0; i1<n1; ++i1) {
      if(tx[i2][i1]<=tmax&&tx[i2][i1]>=tmin) {
        fy[i2][i1] = fx[i2][i1];
      }
    }}
    int d=0;
    for (int i2=e2; i2<n2; ++i2) {
      d++;
    for (int i1= 0; i1<n1; ++i1) {
      if(tx[i2][i1]<=tmax&&tx[i2][i1]>=tmin) {
        fy[i2][i1] = fx[i2][i1]*w[h2+d];
      }
    }}

    return fy;
  }


}

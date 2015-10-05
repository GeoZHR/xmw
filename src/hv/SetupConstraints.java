package hv;

import java.util.ArrayList;

import edu.mines.jtk.dsp.*;
import edu.mines.jtk.util.*;
import edu.mines.jtk.interp.*;
import static edu.mines.jtk.util.ArrayMath.*;

/**
 * Set up control points for generating a horizon volume
 * Method rearrange:
 *   Compute the average depth of all the control points of each set, and set 
 *   the point with the average depth as the reference point for that set.
 * Method extend:
 *   Optional but usually useful: extend scattered control points to a control 
 *   surfaces (each set produces one surface) by using the "HorizonExtractorC" method.
 *
 * @author Xinming Wu
 * @version 2014.03.13
 */

public class SetupConstraints {
  
  public void setForExtension (float sigma1, float sigma2, float scale) {
    _scale  = scale ;
    _sigma1 = sigma1;
    _sigma2 = sigma2;
  } 

  public void setMask(float[][][] mask) {
    _mask = mask;
  }

  public float[][] constraintsFromSurface(float[][] surf) {
    int n3 = surf.length;
    int n2 = surf[0].length;
    int np = n2*n3;
    float[] k1 = new float[np];
    float[] k2 = new float[np];
    float[] k3 = new float[np];
    int ip = 0;
    for (int i3=0; i3<n3; ++i3) {
      for (int i2=0; i2<n2; ++i2) {
        k3[ip] = i3;
        k2[ip] = i2;
        k1[ip] = surf[i3][i2];
        ip ++;
      }
    }
    return rearrange(k1,k2,k3);
  }

  public float[][] rearrange(float[] k1, float[] k2, float[] k3) {
    int ai = 0;
    int np = k1.length;
    float[][] k = zerofloat(np,4);
    float k1a = sum(k1)/(float)np;
    float min = Float.POSITIVE_INFINITY;
    for (int ip=0; ip<np; ip++) {
      int k1i = round(k1[ip]);
      k[0][ip] = (float)k1i; 
      k[1][ip] = k2[ip];
      k[2][ip] = k3[ip];
      k[3][ip] = k1[ip]-k1i;
      float dk = abs(k1[ip]-k1a);
      if (dk<min){min = dk; ai=ip;}
    }
    float t0 = k[0][0];
    float t1 = k[1][0];
    float t2 = k[2][0];
    float t3 = k[3][0];
    k[0][0] = k[0][ai];  k[0][ai] = t0;
    k[1][0] = k[1][ai];  k[1][ai] = t1;
    k[2][0] = k[2][ai];  k[2][ai] = t2;
    k[3][0] = k[3][ai];  k[3][ai] = t3;
    return k;
  }

  public float[][] firstExtend(float[] k1, float[] k2, float[] k3,
    int w2, int w3, float[][][] p2, float[][][] p3, float[][][] ep, 
    float[][][] u, float[][][] wp) {
    int np = k1.length;
    float[][][] kk = new float[np][4][];
    float[] k1i = new float[1];
    float[] k2i = new float[1];
    float[] k3i = new float[1];
    int npp = 0;
    for (int ip=0; ip<np; ++ip) {
      k1i[0] = k1[ip];
      k2i[0] = k2[ip];
      k3i[0] = k3[ip];
      kk[ip] = extend(k1i,k2i,k3i,w2,w3,p2,p3,ep);
      npp += kk[ip][0].length;
    }
    int k = 0;
    float[][] kkp = new float[3][npp];
    for (int ip=0; ip<np; ++ip) {
      int npi = kk[ip][0].length;
      for (int kp=0; kp<npi; ++kp) {
        kkp[0][k] = kk[ip][0][kp]+kk[ip][3][kp];
        kkp[1][k] = kk[ip][1][kp];
        kkp[2][k] = kk[ip][2][kp];
        k++;
      }
    }
    return kkp;
  }

  public float[][] extend(float[] k1, float[] k2, float[] k3,
    int w2, int w3, float[][][] p2, float[][][] p3, float[][][] ep) 
  {
    int n3  = p2.length;
    int n2  = p2[0].length;
    int n1  = p2[0][0].length;
    int ib2 = (int)min(k2)-w2; if(ib2<0   ) {ib2=0;   }
    int ib3 = (int)min(k3)-w3; if(ib3<0   ) {ib3=0;   }
    int ie2 = (int)max(k2)+w2; if(ie2>n2-1) {ie2=n2-1;} 
    int ie3 = (int)max(k3)+w3; if(ie3>n3-1) {ie3=n3-1;} 
    int n2m = ie2-ib2+1;
    int n3m = ie3-ib3+1;
    float[][][] p2m = copy(n1,n2m,n3m,0,ib2,ib3,p2);
    float[][][] p3m = copy(n1,n2m,n3m,0,ib2,ib3,p3);
    float[][][] epm = copy(n1,n2m,n3m,0,ib2,ib3,ep);
    sub(k2,ib2,k2);
    sub(k3,ib3,k3);
    SurfaceExtractorC se = new SurfaceExtractorC();
    se.setCG(0.01f,200);
    se.setExternalIterations(10);
    se.setSmoothings(_sigma1,_sigma2);
    se.setWeights(_scale);
    float[][] surf = checkPoints(k1,k2,k3,p2m,p3m,epm);
    se.surfaceUpdateFromSlopes(epm,p2m,p3m,k1,k2,k3,surf);
    updateSlopes(ib2,ib3,surf,p2,p3);
    //return surf;
    int ci = 0;
    int ai = 0;
    float[][] k = zerofloat(n2m*n3m,4);
    float min = Float.POSITIVE_INFINITY;
    float avg = heightAvg(n1-1,surf);
    for (int i3=0; i3<n3m; i3++) {
      for (int i2=0; i2<n2m; i2++) {
        int i1m = round(surf[i3][i2]);
        int i2m = i2+ib2;     
        int i3m = i3+ib3;     
        if (i1m<n1-1){
          k[0][ci] = (float)i1m; 
          k[1][ci] = (float)i2m;
          k[2][ci] = (float)i3m;
          k[3][ci] = surf[i3][i2]-(float)i1m;
          float df = abs(surf[i3][i2]-avg);
          if (df<min){min = df; ai=ci;}
          ci++;
        }
      }
    }
    float t0 = k[0][0];
    float t1 = k[1][0];
    float t2 = k[2][0];
    float t3 = k[3][0];
    k[0][0] = k[0][ai]; k[0][ai] = t0;
    k[1][0] = k[1][ai]; k[1][ai] = t1;
    k[2][0] = k[2][ai]; k[2][ai] = t2;
    k[3][0] = k[3][ai]; k[3][ai] = t3;
    return copy(ci,4,0,0,k);
  } 

  // make seismic slopes near the control surface to be 
  // consistent to the surface slopes. 
  private void updateSlopes(int ib2, int ib3, 
    float[][] sf, float[][][] p2, float[][][] p3) 
  {
    int n3 = sf.length;
    int n2 = sf[0].length;
    int n1 = p2[0][0].length;
    float[][] g2 = new float[n3][n2];
    float[][] g3 = new float[n3][n2];
    despike(3,sf);
    RecursiveExponentialFilter ref = new RecursiveExponentialFilter(1.0);
    ref.setEdges(RecursiveExponentialFilter.Edges.OUTPUT_ZERO_SLOPE);
    ref.apply(sf,sf);
    int d = 5;
    float[] ip = new float[101]; ip[50] = 1f;
    RecursiveGaussianFilter rgf1 = new RecursiveGaussianFilter(1);
    RecursiveGaussianFilter rgf2 = new RecursiveGaussianFilter(max(d-2,1));
    rgf1.apply10(sf,g2); 
    rgf1.apply01(sf,g3);
    setBounds(g2);
    setBounds(g3);
    rgf2.apply0(ip,ip);
    float[] sc = copy(d*2+1,50-d,ip);
    sc = mul(sc,0.9f/max(sc));
    for (int i3=0; i3<n3; ++i3) {
    for (int i2=0; i2<n2; ++i2) {
      int i2m = i2+ib2;
      int i3m = i3+ib3;
      int i1m = round(sf[i3][i2]);
      int i1b = i1m-d;
      int i1e = i1m+d;
      i1b=max(i1b,0);
      i1e=min(i1e,n1-1);
      float pf2 = g2[i3][i2];
      float pf3 = g3[i3][i2];
      for (int j1=i1b; j1<=i1e; ++j1){
        int di = j1-i1b;
        float ps2 = p2[i3m][i2m][j1];
        float ps3 = p3[i3m][i2m][j1];
        float[] ps = combineSlopes(5f,5f,sc[di],pf2,pf3,ps2,ps3);
        p2[i3m][i2m][j1] = ps[0];
        p3[i3m][i2m][j1] = ps[1];
      }
    }}
  }

  // use surface slopes to better estimate 
  // seismic slopes near the surface
  float[] combineSlopes(float p2max, float p3max,
    float sc, float pf2, float pf3, float ps2, float ps3) 
  {
    float au = 1.000f;
    float av = 0.600f;
    float aw = 0.200f;
    float auv = au-av;
    float awv = aw-av;
    float uf1 = 1f/sqrt(pf2*pf2+pf3*pf3+1f);
    float uf2 = -pf2*uf1;
    float uf3 = -pf3*uf1;
    float us1 = 1f/sqrt(ps2*ps2+ps3*ps3+1f);
    float us2 = -ps2*us1;
    float us3 = -ps3*us1;
    float ufs =  1f/sqrt(uf1*uf1+uf2*uf2);
    float wf1 = -uf2*ufs;
    float wf2 =  uf1*ufs;
    float wf3 =  0.0f;
    float uss =  1f/sqrt(us1*us1+us2*us2);
    float ws1 = -us2*uss;
    float ws2 =  us1*uss;
    float ws3 =  0.0f;
    float uf11 = uf1*uf1;
    float uf12 = uf1*uf2;
    float uf13 = uf1*uf3;
    float uf22 = uf2*uf2;
    float uf23 = uf2*uf3;
    float uf33 = uf3*uf3;
    float us11 = us1*us1;
    float us12 = us1*us2;
    float us13 = us1*us3;
    float us22 = us2*us2;
    float us23 = us2*us3;
    float us33 = us3*us3;
    float wf11 = wf1*wf1;
    float wf12 = wf1*wf2;
    float wf13 = wf1*wf3;
    float wf22 = wf2*wf2;
    float wf23 = wf2*wf3;
    float wf33 = wf3*wf3;
    float ws11 = ws1*ws1;
    float ws12 = ws1*ws2;
    float ws13 = ws1*ws3;
    float ws22 = ws2*ws2;
    float ws23 = ws2*ws3;
    float ws33 = ws3*ws3;
    double[] e = new double[3];
    double[][] z = new double[3][3];
    double[][] a = new double[3][3];

    a[0][1] += sc*(auv*uf12+awv*wf12);
    a[0][2] += sc*(auv*uf13+awv*wf13);
    a[1][2] += sc*(auv*uf23+awv*wf23);

    a[0][1] += (1f-sc)*(auv*us12+awv*ws12);
    a[0][2] += (1f-sc)*(auv*us13+awv*ws13);
    a[1][2] += (1f-sc)*(auv*us23+awv*ws23);

    a[1][0] = a[0][1];
    a[2][0] = a[0][2];
    a[2][1] = a[1][2];

    a[0][0] += sc*(auv*uf11+awv*wf11+av);
    a[1][1] += sc*(auv*uf22+awv*wf22+av);
    a[2][2] += sc*(auv*uf33+awv*wf33+av);

    a[0][0] += (1f-sc)*(auv*us11+awv*ws11+av);
    a[1][1] += (1f-sc)*(auv*us22+awv*ws22+av);
    a[2][2] += (1f-sc)*(auv*us33+awv*ws33+av);
    Eigen.solveSymmetric33(a,z,e);
    float u1i = (float)z[0][0];
    float u2i = (float)z[0][1];
    float u3i = (float)z[0][2];
    if (u1i<0.0f) {
       u1i = -u1i;
       u2i = -u2i;
       u3i = -u3i;
    }
    float p2=0.0f;
    float p3=0.0f;
    if (-u2i<-p2max*u1i) u2i =  p2max*u1i;
    if (-u2i> p2max*u1i) u2i = -p2max*u1i;
    if (-u3i<-p3max*u1i) u3i =  p3max*u1i;
    if (-u3i> p3max*u1i) u3i = -p3max*u1i;
    if (u1i==0.0f) {
      p2 = (u2i<0.0f)?p2max:-p2max;
      p3 = (u3i<0.0f)?p3max:-p3max;
    } else {
      p2 = -u2i/u1i;
      p3 = -u3i/u1i;
    }
    return new float[]{p2,p3};
  }

  private void setBounds(float[][] g) {
    int n3 = g.length;
    int n2 = g[0].length;
    for (int i3=0; i3<4; ++i3) {
    for (int i2=0; i2<n2; ++i2) {
      g[i3][i2] = g[4][i2];
    }}
    for (int i3=n3-4; i3<n3; ++i3) {
    for (int i2=0; i2<n2; ++i2) {
      g[i3][i2] = g[n3-5][i2];
    }}
    for (int i2=0; i2<4; ++i2) {
    for (int i3=0; i3<n3; ++i3) {
      g[i3][i2] = g[i3][4];
    }}
    for (int i2=n2-4; i2<n2; ++i2) {
    for (int i3=0; i3<n3; ++i3) {
      g[i3][i2] = g[i3][n2-5];
    }}
  }

  
  /**
   * Applies a despking filter to a control surface.
   * @param nmed number of median-of-nine filter passes.
   * @param surf input surface
   */

  private void despike(int nmed, float[][] surf) {
    for (int imed=0; imed<nmed; ++imed) {
      despike(surf);
    }
  }

  private void despike(float[][] f) {
    int n3 = f.length;
    int n2 = f[0].length;
    float[] fs = new float[9];
    MedianFinder mf = new MedianFinder(9);
    for (int i3=1; i3<n3-1; ++i3) {
    for (int i2=1; i2<n2-1; ++i2) {
      fs[0] = f[i3-1][i2-1];
      fs[1] = f[i3-1][i2  ];
      fs[2] = f[i3-1][i2+1];
      fs[3] = f[i3  ][i2-1];
      fs[4] = f[i3  ][i2  ];
      fs[5] = f[i3  ][i2+1];
      fs[6] = f[i3+1][i2-1];
      fs[7] = f[i3+1][i2  ];
      fs[8] = f[i3+1][i2+1];
      f[i3][i2] = mf.findMedian(fs);
    }}
  }

  // This method is applied when many control points are privated and 
  // some of them might not be correct
  private float[][] checkPoints(float[] k1, float[] k2, float[] k3,
    float[][][] p2, float[][][] p3, float[][][] ep) 
  {
    int n3  = p2.length;
    int n2  = p2[0].length;
    int n1  = p2[0][0].length;
    SurfaceExtractorC se = new SurfaceExtractorC();
    se.setCG(0.01f,200);
    se.setExternalIterations(10);
    se.setSmoothings(_sigma1,_sigma2);
    se.setWeights(_scale);
    float lmt = n1-1.f;
    float[][] surf = se.surfaceInitialization(n2,n3,lmt,k1,k2,k3);
    se.surfaceUpdateFromSlopes(ep,p2,p3,k1,k2,k3,surf);
    despike(3,surf);
    return pointCorrection(k1,k2,k3,surf);
  }


  private float[][] pointCorrection(
    float[] k1, float[] k2, float[] k3, float[][] sf) 
  {
    int np = k1.length;
    int n3 = sf.length;
    int n2 = sf[0].length;
    for (int ip=0; ip<np; ++ip) {
      int i2 = (int)k2[ip];
      int i3 = (int)k3[ip];
      if(i2<4||i3<4){continue;}
      if(i2>n2-5||i3>n3-5){continue;}
      int i2m = i2-8;
      int i2p = i2+8;
      int i3m = i3-8;
      int i3p = i3+8;
      i2m = max(i2m,0);
      i3m = max(i3m,0);
      i2p = min(i2p,n2-1);
      i3p = min(i3p,n3-1);
      FloatList x2s = new FloatList();
      FloatList x3s = new FloatList();
      FloatList fxs = new FloatList();
      for (int j3=i3m; j3<i3-2; ++j3) {
      for (int j2=i2m; j2<i2-2; ++j2) {
        x2s.add(j2);
        x3s.add(j3);
        fxs.add(sf[j3][j2]);
      }}
      for (int j3=i3+3; j3<=i3p; ++j3) {
      for (int j2=i2+3; j2<=i2p; ++j2) {
        x2s.add(j2);
        x3s.add(j3);
        fxs.add(sf[j3][j2]);
      }}
      SibsonInterpolator2 si=
        new SibsonInterpolator2(fxs.trim(),x2s.trim(),x3s.trim());    
      k1[ip] = si.interpolate(i2,i3);
      for (int j3=i3-3; j3<=i3+2; ++j3) {
      for (int j2=i2-3; j2<=i2+2; ++j2) {
        sf[j3][j2] = si.interpolate(j2,j3);
      }}
    }
    return sf;
  }

  private float heightAvg(float lmt, float[][] x) {
    int n2=x.length;
    int n1=x[0].length;
    float sum = 0.0f;
    float ci = 0.0f;
    for (int i2=0; i2<n2; ++i2) {
      for (int i1=0; i1<n1; ++i1) {
        float xi = x[i2][i1];
        if(xi<lmt) {sum +=xi; ci +=1.0f;}
      }
    }
    return sum/ci;
  }
  private float _scale  = 0.0f;
  private float _sigma1 = 6.0f;
  private float _sigma2 = 6.0f;
  private float[][][] _mask=null;

  private class FloatList {
    public int n;
    public float[] a = new float[1024];
    public void add(float f) {
      if (n==a.length) {
        float[] t = new float[2*n];
        System.arraycopy(a,0,t,0,n);
        a = t;
      }
      a[n++] = f;
    }
    public float[] trim() {
      if (n==0)
        return null;
      float[] t = new float[n];
      System.arraycopy(a,0,t,0,n);
      return t;
    }
  }

}

package rgi;

import edu.mines.jtk.dsp.*;
import static edu.mines.jtk.util.ArrayMath.*;

import ipfx.*;
import java.util.*;

/**
 * Convert points between input, unfault, and flattened spaces.
 * @author Xinming Wu
 * @version 2015.07.26
 */
public class ConvertPoints {

  public float[][] getSamplesX(Point[][] ps) {
    int ns = ps.length;
    FloatList fvl = new FloatList();
    FloatList x1l = new FloatList();
    FloatList x2l = new FloatList();
    FloatList x3l = new FloatList();
    for (int is=0; is<ns; ++is) {
      int np = ps[is].length;
    for (int ip=0; ip<np; ++ip) {
      Point pi = ps[is][ip];
      fvl.add(pi.fv);
      x1l.add(pi.x1);
      x2l.add(pi.x2);
      x3l.add(pi.x3);
    }}
    float[] fv = fvl.trim();
    float[] x1 = x1l.trim();
    float[] x2 = x2l.trim();
    float[] x3 = x3l.trim();
    return new float[][]{fv,x1,x2,x3};
  }

  public float[][] getSamplesW(Point[][] ps) {
    int ns = ps.length;
    FloatList fvl = new FloatList();
    FloatList x1l = new FloatList();
    FloatList x2l = new FloatList();
    FloatList x3l = new FloatList();
    for (int is=0; is<ns; ++is) {
      int np = ps[is].length;
    for (int ip=0; ip<np; ++ip) {
      Point pi = ps[is][ip];
      fvl.add(pi.fv);
      x1l.add(pi.w1);
      x2l.add(pi.w2);
      x3l.add(pi.w3);
    }}
    float[] fv = fvl.trim();
    float[] x1 = x1l.trim();
    float[] x2 = x2l.trim();
    float[] x3 = x3l.trim();
    return new float[][]{fv,x1,x2,x3};
  }


  public float[][] getSamplesU(Point[][] ps) {
    int ns = ps.length;
    FloatList fvl = new FloatList();
    FloatList x1l = new FloatList();
    FloatList x2l = new FloatList();
    FloatList x3l = new FloatList();
    for (int is=0; is<ns; ++is) {
      int np = ps[is].length;
    for (int ip=0; ip<np; ++ip) {
      Point pi = ps[is][ip];
      fvl.add(pi.fv);
      x1l.add(pi.u1);
      x2l.add(pi.u2);
      x3l.add(pi.u3);
    }}
    float[] fv = fvl.trim();
    float[] x1 = x1l.trim();
    float[] x2 = x2l.trim();
    float[] x3 = x3l.trim();
    return new float[][]{fv,x1,x2,x3};
  }

  /**
   * Set points with coordinates in unfaulted space
   */
  public Point[][] setUnfaultCoord(float[][][] fx, FaultSkin[] skins,
    float[][][] sw1, float[][][] sw2, float[][][] sw3) 
  {
    int n3 = sw1.length;
    int n2 = sw1[0].length;
    int n1 = sw1[0][0].length;
    short[][][] mk = new short[n3][n2][n1];
    faultMark(skins,mk);
    int ns = fx[0].length;
    int np = fx[0][0].length;
    Point[][] ps = new Point[ns][np];
    for (int is=0; is<ns; ++is) {
      boolean onFault = false;
      for (int ip=0; ip<np; ++ip) {
        float fv = fx[0][is][ip];
        float x1 = fx[1][is][ip];
        float x2 = fx[2][is][ip];
        float x3 = fx[3][is][ip];
        int i1 = (int)x1;
        int i2 = (int)x2;
        int i3 = (int)x3;
        float s1 = sw1[i3][i2][i1];
        float s2 = sw2[i3][i2][i1];
        float s3 = sw3[i3][i2][i1];
        float w1 = x1-s1;
        float w2 = x2-s2;
        float w3 = x3-s3;
        if(mk[i3][i2][i1]==1){onFault=true;}
        ps[is][ip] = new Point(fv,x1,x2,x3,w1,w2,w3,onFault);
      }
    }
    return ps;
  }

  public Point[][] setFlattenedCoord(
    Sampling s1, Sampling s2, Sampling s3, float[][][] rgt, Point[][] ps) {
    int ns = ps.length;
    SincInterpolator si = new SincInterpolator();
    si.setExtrapolation(SincInterpolator.Extrapolation.CONSTANT);
    for (int is=0; is<ns; ++is) {
      int np = ps[0].length;
      Point p0 = ps[is][0];
      float w1 = p0.w1;
      float w2 = p0.w2;
      float w3 = p0.w3;
      float u1 = round(si.interpolate(s1,s2,s3,rgt,w1,w2,w3));
      p0.setCoordU(u1,w2,w3);
      Point pl = ps[is][np-1];
      w1 = pl.w1;
      w2 = pl.w2;
      w3 = pl.w3;
      u1 = round(si.interpolate(s1,s2,s3,rgt,w1,w2,w3));
      pl.setCoordU(u1,w2,w3);
    for (int ip=1; ip<np-1; ++ip) {
      Point pm = ps[is][ip-1];
      Point pi = ps[is][ip  ];
      Point pp = ps[is][ip+1];
      float w1m = pm.w1;
      float w2m = pm.w2;
      float w3m = pm.w3;
      float w1i = pi.w1;
      float w2i = pi.w2;
      float w3i = pi.w3;
      float w1p = pp.w1;
      float w2p = pp.w2;
      float w3p = pp.w3;
      float u1mf = si.interpolate(s1,s2,s3,rgt,w1m,w2m,w3m);
      float u1if = si.interpolate(s1,s2,s3,rgt,w1i,w2i,w3i);
      float u1pf = si.interpolate(s1,s2,s3,rgt,w1p,w2p,w3p);
      int u1m = round(u1mf);
      int u1i = round(u1if);
      int u1p = round(u1pf);
      if(u1i==u1m&&(u1i+1.0f)!=u1p){
        pi.setCoordU(u1if+1.0f,w2i,w3i);
      } else {
        pi.setCoordU(u1if,w2i,w3i);
      }
    }}
    /*
    double fu1 = s1.getFirst();
    double du1 = s1.getDelta();
    int nu1 = (int)round((max(rgt)-fu1)/du1); 
    Sampling su1 = new Sampling(nu1,du1,fu1);
    */
    double fu1 = 0.0;
    double lu1 = max(rgt);
    double du1 = s1.getDelta();
    int nu1 = (int)((lu1-fu1)/du1)+1;
    Sampling su1 = new Sampling(nu1,du1,fu1);

    //return ps;
    return checkPoints(su1,ps);
  }

  public Point[][] setFlattenedCoord(float[][][] x1s, Point[][] ps) {
    int ns = ps.length;
    int n3 = x1s.length;
    int n2 = x1s[0].length;
    int n1 = x1s[0][0].length;
    float[][] xk = new float[1][n1];
    for (int is=0; is<ns; ++is) {
      int np = ps[is].length;
      for (int ip=0; ip<np; ++ip) {
        Point pi = ps[is][ip];
        float w1 = pi.w1;
        float w2 = pi.w2;
        float w3 = pi.w3;
        int i2 = round(w2);
        int i3 = round(w3);
        if(i2<0){i2=0;} if(i2>=n2){i2=n2-1;}
        if(i3<0){i3=0;} if(i3>=n3){i3=n3-1;}
        xk[0] = x1s[i3][i2];
        KdTree kt = new KdTree(xk);
        int ik = kt.findNearest(new float[]{w1});
        pi.setCoordU(ik,i2,i3);
      }
    }
    return ps;
  }

  public Point[][] checkPoints(int n1, Point[][] ps) {
    int ns = ps.length;
    Point[][] pp = new Point[ns][];
    for (int is=0; is<ns; ++is) {
      int np = ps[is].length;
      float[][] x = new float[1][np];
      for (int ip=0; ip<np; ++ip) {
        x[0][ip] = ps[is][ip].u1; 
      }
      float xb = min(x);
      float xe = max(x);
      KdTree kt = new KdTree(x);
      ArrayList<Point> pa = new ArrayList<Point>();
      for (int i1=0; i1<n1; ++i1) {
        int xk = kt.findNearest(new float[]{i1});
        float xi = x[0][xk];
        Point pi = ps[is][xk];
        float di = abs(xi-i1);
        if((xi==xb||xi==xe||pi.onFault) && di>0.0f){continue;}
        Point pk = new Point(pi.fv,pi.x1,pi.x2,pi.x3,pi.onFault);
        pk.setCoordW(pi.w1,pi.w2,pi.w3);
        pk.setCoordU(i1,   pi.u2,pi.u3);
        pa.add(pk);
      }
      int ik = 0;
      int npp = pa.size();
      pp[is] = new Point[npp];
      for (Point pi:pa) {
        pp[is][ik] = pi;
        ik++;
      }
    }
    return pp;
  }


  public Point[][] checkPoints(Sampling s1, Point[][] ps) {
    int ns = ps.length;
    int n1 = s1.getCount();
    double f1 = s1.getFirst();
    double d1 = s1.getDelta();
    Point[][] pp = new Point[ns][];
    for (int is=0; is<ns; ++is) {
      int np = ps[is].length;
      float[][] x = new float[1][np];
      for (int ip=0; ip<np; ++ip) {
        x[0][ip] = ps[is][ip].u1; 
      }
      float xb = min(x);
      float xe = max(x);
      KdTree kt = new KdTree(x);
      ArrayList<Point> pa = new ArrayList<Point>();
      for (int i1=0; i1<n1; ++i1) {
        float u1 = (float)(f1+d1*i1);
        int xk = kt.findNearest(new float[]{u1});
        float xi = x[0][xk];
        Point pi = ps[is][xk];
        float di = abs(round(xi-u1));
        if((xi==xb||xi==xe||pi.onFault)&&di>1.0f){continue;}
        if(di>1.0f){continue;}
        Point pk = new Point(pi.fv,pi.x1,pi.x2,pi.x3,pi.onFault);
        pk.setCoordW(pi.w1,pi.w2,pi.w3);
        pk.setCoordU(u1,   pi.u2,pi.u3);
        pa.add(pk);
      }
      int ik = 0;
      int npp = pa.size();
      pp[is] = new Point[npp];
      for (Point pi:pa) {
        pp[is][ik] = pi;
        ik++;
      }
    }
    return pp;
  }

 /*
  public Point[][] checkPoints(Sampling s1, Point[][] ps) {
    int ns = ps.length;
    int n1 = s1.getCount();
    double f1 = s1.getFirst();
    double d1 = s1.getDelta();
    Point[][] pp = new Point[ns][];
    for (int is=0; is<ns; ++is) {
      int np = ps[is].length;
      float[][] x = new float[1][np];
      for (int ip=0; ip<np; ++ip) {
        x[0][ip] = (int)round((ps[is][ip].u1-f1)/d1); 
      }
      float xb = min(x);
      float xe = max(x);
      KdTree kt = new KdTree(x);
      ArrayList<Point> pa = new ArrayList<Point>();
      for (int i1=0; i1<n1; ++i1) {
        int xk = kt.findNearest(new float[]{i1});
        float xi = x[0][xk];
        Point pi = ps[is][xk];
        if(xi==xb||xi==xe||pi.onFault){continue;}
        float u1 = (float)(f1+d1*i1);
        pi.setU1(u1);
        pa.add(pi);
      }
      int ik = 0;
      int npp = pa.size();
      pp[is] = new Point[npp];
      for (Point pi:pa) {
        pp[is][ik] = pi;
        ik++;
      }
    }
    return pp;
  }
 */

  private void faultMark(FaultSkin[] skins, short[][][] mk) {
    int n3 = mk.length;
    int n2 = mk[0].length;
    for (FaultSkin skin:skins) {
    for (FaultCell cell:skin) {
      int i1 = cell.getI1();
      int i2 = cell.getI2();
      int i3 = cell.getI3();
      mk[i3][i2][i1] = 1;
    }}
  }

  private class FloatList {
    public int n = 0;
    public float[] a = new float[1024];
    public void add(float f) {
      if (n==a.length) {
        float[] t = new float[2*n];
        System.arraycopy(a,0,t,0,n);
        a = t;
      }
      a[n++] = f;
    }
    public void add(float[] f) {
      int m = f.length;
      for (int i=0; i<m; ++i)
        add(f[i]);
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
      

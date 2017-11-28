package hp;

import java.util.*;
import edu.mines.jtk.dsp.*;
import static edu.mines.jtk.util.ArrayMath.*;

/**
 * Extract a single seismic horizon curve with control points,
 * the constraints derived from control points are incorporated in a 
 * preconditioner in the conjugate gradient method used to solve the 
 * linear system for horizon extracting.
 * @author Xinming Wu, University of Texas at Austin
 * @version 2017.11.21
 */

public class HorizonPicker2{

  public int[][][] findPTZs(float[][] fx) {
    int n2 = fx.length;
    int n1 = fx[0].length;
    float[][] gx = new float[n2][n1];
    RecursiveGaussianFilter rgf = new RecursiveGaussianFilter(1);
    rgf.apply0X(fx,gx);
    int[][][] mks = new int[3][n2][];
    for (int i2=0; i2<n2  ; ++i2) {
      float[] gx2 = gx[i2];
      ArrayList<Integer> ps = new ArrayList<Integer>();
      ArrayList<Integer> ns = new ArrayList<Integer>();
      ArrayList<Integer> zs = new ArrayList<Integer>();
      for (int i1=1; i1<n1-1; ++i1) {
        int i1m = i1-1;
        int i1p = i1+1;
        float gxi = gx2[i1 ];
        float gxm = gx2[i1m];
        float gxp = gx2[i1p];
        float gum = gxi-gxm;
        float gup = gxp-gxi;
        if(gup*gum<0f) {
          if(gxi>gxm&&gxi>gxp) ps.add(i1);
          if(gxi<gxm&&gxi<gxp) ns.add(i1);
        }
        if(gxp*gxm<0f&&gup*gum>0f) zs.add(i1);
      }
      int np = ps.size();
      int nn = ns.size();
      int nz = zs.size();
      mks[0][i2] = new int[np];
      mks[1][i2] = new int[nn];
      mks[2][i2] = new int[nz];
      for (int ip=0; ip<np; ++ip)
        mks[0][i2][ip] = ps.get(ip);
      for (int in=0; in<nn; ++in)
        mks[1][i2][in] = ns.get(in);
      for (int iz=0; iz<nz; ++iz)
        mks[2][i2][iz] = zs.get(iz);
    }
    return mks;
  }

  public float[][][] getTroughs(int w1, int[][] mps, float[][] gx) {
    int n2 = gx.length;
    ArrayList<float[]> tx = new ArrayList<float[]>();
    ArrayList<float[]> tf = new ArrayList<float[]>();
    for (int i2=0; i2<n2; ++i2) {
      float[] gx2 = gx[i2];
      int[] ps = mps[i2];
      int np = ps.length;
      // add points
      for (int ip=0; ip<np; ip++) {
        int p1 = ps[ip];
        float[] pt = getTrace(w1,p1,gx2);
        tx.add(new float[]{p1,i2});
        tf.add(pt);
      }
    }
    int np = tx.size();
    float[][][] txf = new float[2][np][];
    for (int ip=0; ip<np; ++ip) {
      txf[0][ip] = tx.get(ip);
      txf[1][ip] = tf.get(ip);
    }
    return txf;
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

  public List<Point> setDataPoints(int w1, int[][] mps, float[][] gx) {
    int n2 = gx.length;
    List<Point> points = new LinkedList<Point>();
    for (int i2=0; i2<n2; ++i2) {
      float[] gx2 = gx[i2];
      int[] ps = mps[i2];
      int np = ps.length;
      // add points
      for (int ip=0; ip<np; ip++) {
        int p1 = ps[ip];
        float[] pt = getTrace(w1,p1,gx2);
        Point point = new Point(new int[]{ip,i2},pt);
        points.add(point);
      }
    }
    return points;
  }

  public List<Point> setDataPoints(int w1, int d1, int[][] mps, float[][] gx) {
    int n2 = gx.length;
    List<Point> points = new LinkedList<Point>();
    for (int i2=0; i2<n2; ++i2) {
      float[] gx2 = gx[i2];
      int[] ps = mps[i2];
      int np = ps.length;
      // add points
      for (int ip=0; ip<np; ip++) {
        int p1 = ps[ip];
        float[][] pt = new float[w1*2+1][d1*2+1];
        for (int k1=-w1; k1<=w1; k1++) {
          pt[k1+w1] = copy(getTrace(d1,p1+k1,gx2));
        }
        Point point = new Point(new int[]{ip,i2},pt);
        points.add(point);
      }
    }
    return points;
  }


  public Cluster[] setInitialClusters(int w1, int[][] mps, float[][] gx) {
    int m2 = 0;
    int mp = 0;
    int n2 = gx.length;
    for (int i2=0; i2<n2; ++i2) {
      int np = mps[i2].length;
      if(np>mp) {mp=np; m2=i2;}
    }
    Cluster[] clusters = new Cluster[mp];
    for (int ip=0; ip<mp; ++ip) {
      int p1 = mps[m2][ip];
      float[] pt = getTrace(w1,p1,gx[m2]);
      Point point = new Point(pt);
      clusters[ip] = new Cluster(ip,point);
    }
    return clusters;
  }

  public Cluster[] setInitialClusters(int w1, int d1, int[][] mps, float[][] gx) {
    int m2 = 0;
    int mp = 0;
    int n2 = gx.length;
    for (int i2=0; i2<n2; ++i2) {
      int np = mps[i2].length;
      if(np>mp) {mp=np; m2=i2;}
    }
    Cluster[] clusters = new Cluster[mp];
    for (int ip=0; ip<mp; ++ip) {
      int p1 = mps[m2][ip];
      float[][] pt = new float[w1*2+1][d1*2+1];
      for (int k1=-w1; k1<=w1; k1++) {
        pt[k1+w1] = copy(getTrace(d1,p1+k1,gx[m2]));
      }
      Point point = new Point(pt);
      clusters[ip] = new Cluster(ip,point);
    }
    return clusters;
  }



  public void setLinkConstraints(
    int[][] mps, List<Point> points, float[][] gx) {
    int n2 = gx.length;
    int n1 = gx[0].length;
    int ip = 0;
    int[][] pmap = new int[n2][n1];
    for (int i2=0; i2<n2; ++i2) {
      int[] ps = mps[i2];
      int np = ps.length;
      for (int kp=0; kp<np; kp++) {
        int i1 = ps[kp];
        pmap[i2][i1] = ip;
        ip++;
      }
    }
    for (int i2=1; i2<n2-1; ++i2) {
      int[] ps = mps[i2];
      int np = ps.length;
      for (int kp=0; kp<np; kp++) {
        int i1 = ps[kp];
        if(i1<1||i1>=n1-1) continue;
        float gxi = gx[i2][i1];
        float dgl = FLT_MAX;
        float dgr = FLT_MAX;
        int i1l = 0;
        int i1r = 0;
        int i2m = i2-1;
        int i2p = i2+1;
        int i1m = i1-1;
        int i1p = i1+1;
        float dglm = abs(gx[i2m][i1m]-gxi);
        float dgli = abs(gx[i2m][i1 ]-gxi);
        float dglp = abs(gx[i2m][i1p]-gxi);
        float dgrm = abs(gx[i2p][i1m]-gxi);
        float dgri = abs(gx[i2p][i1 ]-gxi);
        float dgrp = abs(gx[i2p][i1p]-gxi);
        if (pmap[i2m][i1m]>0&&dglm<dgl) {dgl = dglm; i1l = i1m;}
        if (pmap[i2m][i1 ]>0&&dgli<dgl) {dgl = dgli; i1l = i1 ;}
        if (pmap[i2m][i1p]>0&&dglp<dgl) {dgl = dglp; i1l = i1p;}
        if (pmap[i2p][i1m]>0&&dgrm<dgr) {dgr = dgrm; i1r = i1m;}
        if (pmap[i2p][i1 ]>0&&dgri<dgr) {dgr = dgri; i1r = i1 ;}
        if (pmap[i2p][i1p]>0&&dgrp<dgr) {dgr = dgrp; i1r = i1p;}
        int ipi = pmap[i2][i1  ];
        int ipl = pmap[i2m][i1l];
        int ipr = pmap[i2p][i1r];
        Point pi = points.get(ipi);
        Point pl = points.get(ipl);
        Point pr = points.get(ipr);
        if (dgl<0.5f) pi.setLeftPoint(pl);
        if (dgr<0.5f) pi.setRightPoint(pr);
      }
    }
  }


  public float[][] getClusterMap(int n1, int n2, int[][] mps, Cluster[] clusters) {
    float[][] cmap = new float[n2][n1];
    int nc = clusters.length;
    for (int ic=0; ic<nc; ++ic) {
      List<Point> points = clusters[ic].getPoints();
      for (Point point:points) {
        int i2 = point.getI2();
        int p1 = point.getI1();
        int i1 = mps[i2][p1];
        cmap[i2][i1] = ic+1;
      }
    }
    return cmap;
  }

  public Cluster[] applyClustering(int[] w1s, int[][] mps, float[][] gx) {
    int n2 = gx.length;
    int nw = w1s.length;
    CKMeanClusterer kmc = new CKMeanClusterer();
    kmc.setConvergence(100,0.000001f);
    float[][][] mpc = new float[n2][][];
    for (int i2=0; i2<n2; ++i2) {
      int np = mps[i2].length;
      mpc[i2] = new float[np][nw];
    }
    for (int iw=0; iw<nw; ++iw) {
      System.out.println("iw="+iw);
      int w1i = w1s[iw];
      List<Point> points = setDataPoints(w1i,mps,gx);
      Cluster[] clustert = setInitialClusters(w1i,mps,gx);
      Cluster[] clusters = kmc.applyClustering(false,copy(mps),clustert,points);
      for (Cluster cluster:clusters) {
        List<Point> cpoints = cluster.getPoints();
        for (Point point:cpoints) {
          int ip = point.getI1();
          int i2 = point.getI2();
          mpc[i2][ip][iw] = point.getClusterId();
        }
      }
    }

    List<Point> points = new LinkedList<Point>();
    for (int i2=0; i2<n2; ++i2) {
      int np = mpc[i2].length;
      for (int ip=0; ip<np; ++ip) {
        Point point = new Point(new int[]{ip,i2},mpc[i2][ip]);
        points.add(point);
      }
    }
    int m2 = 0, mp = 0;
    for (int i2=0; i2<n2; ++i2) {
      int np = mps[i2].length;
      if(np>mp) {mp=np; m2=i2;}
    }
    Cluster[] clustert = new Cluster[mp];
    for (int ip=0; ip<mp; ++ip) {
      Point point = new Point(mpc[m2][ip]);
      clustert[ip] = new Cluster(ip,point);
    }
    Cluster[] clusters = kmc.applyClustering(true,copy(mps),clustert,points);
    return clusters;
  }

  public float[][] pickPeakHorizon(
    int c1, int c2, int w1, int[][] mps, float[][] gx, float[][] hp) {
    int n2 = gx.length;
    int n1 = gx[0].length;
    int[] hs = new int[n2];
    float[][] ts = new float[n2][w1*2+1];
    ts[c2] = getTrace(w1,mps[c2][c1],gx[c2]);
    hp[c2][mps[c2][c1]] = 1;
    hs[c2] = mps[c2][c1];
    int smax = round(w1*0.5f);
    DynamicWarping dw = new DynamicWarping(-smax,smax);
    dw.setStrainMax(0.1);
    //pick left
    float[] ws = new float[w1*2+1];
    ws[w1] = 1;
    RecursiveGaussianFilter rgf = new RecursiveGaussianFilter(w1);
    rgf.apply0(ws,ws);
    float[] sms = new float[n2];
    sms[c2] = 1f;
    for (int i2=c2-1; i2>=0; i2--) {
      float sm = FLT_MIN;
      int[] mp2 = mps[i2];
      int np = mp2.length;
      int m1 = 0;
      for (int ip=0; ip<np; ++ip) {
        int i1 = mp2[ip];
        float[] t2 = getTrace(w1,i1,gx[i2]);
        float ss = 0f;
        float cs = 0f;
        for (int k2=i2+1; k2<=c2; k2++){
          if(sms[k2]<0.7f) continue;
          cs += sms[k2];
          ss += semblance(ts[k2],t2,ws)*sms[k2];
        }
        ss /= cs;
        if(ss>sm) {m1 = i1;sm=ss;} 
      }
      hs[i2] = m1;
      sms[i2] = sm;
      ts[i2] = getTrace(w1,m1,gx[i2]);
      if(sm>0.5f) {hp[i2][m1] = c1;}
    }
    //pick right

    for (int i2=c2+1; i2<n2; i2++) {
      float sm = FLT_MIN;
      int[] mp2 = mps[i2];
      int np = mp2.length;
      int m1 = 0;
      for (int ip=0; ip<np; ++ip) {
        int i1 = mp2[ip];
        float[] t2 = getTrace(w1,i1,gx[i2]);
        float ss = 0f;
        float cs = 0f;
        for (int k2=i2-1; k2>=0; k2--) {
          if(sms[k2]<0.7f) continue;
          cs += sms[k2];
          ss += semblance(ts[k2],t2,ws)*sms[k2];
        }
        ss /= cs;
        if(ss>sm) {m1 = i1;sm=ss;} 
      }
      sms[i2] = sm;
      hs[i2] = m1;
      ts[i2] = getTrace(w1,m1,gx[i2]);
      if(sm>0.5f) {hp[i2][m1] = c1;}
    }
    return ts;
  }
  
  public float distance(float[] f, float[] g) {
    int n1 = f.length;
    float ds = 0f;
    for (int i1=0; i1<n1; ++i1) {
      float fi = f[i1];
      float gi = g[i1];
      float di = fi-gi;
      ds += di*di;
    }
    return sqrt(ds);
  }


  public float semblance(float[] f, float[] g, float[] ws) {
    int n1 = f.length;
    float ds = 0f;
    float ns = 0f;
    for (int i1=0; i1<n1; ++i1) {
      float fi = f[i1];
      float gi = g[i1];
      float si = (fi+gi)*0.5f;
      ds += si*si*ws[i1];
      float fs = fi*fi;
      float gs = gi*gi;
      ns += (fs+gs)*0.5f*ws[i1];
    }
    return ds/ns;
  }

  public float[] getTrace(int w1, int c1, float[] gx) {
    int n1 = gx.length;
    float[] tc = new float[w1*2+1];
    Random rd = new Random();
    for (int k1=-w1; k1<=w1; k1++) {
      int i1 = k1+c1;
      if(i1<0||i1>=n1) i1 = rd.nextInt(n1);
      tc[k1+w1] = gx[i1];
    }
    return tc;
  }
 
 
  ///////////////////////////////////////////////////////////////////////////
  // private
  private float _weight = 0.5f; // -(fxx+fyy)+weight*(fxxxx+2fxxyy+fyyyy)=px+qy
  private float _sigma1 = 6.0f; // precon smoothing extent for 1st dim
  private float _small = 0.01f; // stop CG iterations if residuals small
  private int _niter = 200; // maximum number of CG iterations
  private int _exniter = 10; // external iterations of surface updating
  private float _wc = 1f;
  private int _d2 = 2;
  private int _dm = 20;


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

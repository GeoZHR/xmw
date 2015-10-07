package swt;

import edu.mines.jtk.dsp.*;
import edu.mines.jtk.util.*;
import static edu.mines.jtk.util.ArrayMath.*;

/**
 * Seismic depth conversion.
 * 
 * @author Xinming Wu, Colorado School of Mines
 * @version 2015.09.29
 */

public class SeismicDepthConversion {

  public float[][][] convert(
    Sampling sz, Sampling st, Sampling s1, Sampling s2, Sampling s3, 
    final float[][][] gx, final float[][][] zt) 
  {
    final int n3 = gx.length;
    final int n2 = gx[0].length;
    final int n1 = gx[0][0].length;
    final int nz = sz.getCount();
    final InverseInterpolator ii = new InverseInterpolator(st,sz);
    final float[][][] tz = new float[n3][n2][nz];
    Parallel.loop(n3,new Parallel.LoopInt() {
    public void compute(int i3) {
    for (int i2=0; i2<n2; ++i2) {
      ii.invert(zt[i3][i2],tz[i3][i2]);
    }}});
    final double d1 = s1.getDelta();
    final double f1 = s1.getFirst();
    final SincInterpolator si = new SincInterpolator();
    final float[][][] gz = new float[n3][n2][nz];
    Parallel.loop(n3,new Parallel.LoopInt() {
    public void compute(int i3) {
      for (int i2=0; i2<n2; ++i2) {
        si.interpolate(n1,d1,f1,gx[i3][i2],nz,tz[i3][i2],gz[i3][i2]);
    }}});
    return gz;
  }
  
  public float[][][] depthTime(
    Sampling s1, Sampling s2, Sampling s3, 
    float[][][] vt, float[][][] gt) 
  {
    int n1 = s1.getCount();
    int n2 = s2.getCount();
    int n3 = s3.getCount();
    float[][][] zt = new float[n3][n2][n1];
    for (int i3=0; i3<n3; ++i3) {
    for (int i2=0; i2<n2; ++i2) {
    for (int i1=0; i1<n1; ++i1) {
      float ti = (float)s1.getValue(i1); 
      zt[i3][i2][i1] = 0.5f*vt[i3][i2][i1]*ti;
      if(i1>0 && zt[i3][i2][i1]<=zt[i3][i2][i1-1])
        zt[i3][i2][i1] = zt[i3][i2][i1-1]+0.000001f;
    }}}
    return zt;
  }

  public float[][][] averageVelInterp(
    Sampling s1, Sampling s2, Sampling s3, SynSeis.Model[] mds, float[][][] gt) 
  {
    float[][] fvx = averageVelAtWells(s1,mds);
    RgtInterp3 rgi = new RgtInterp3(fvx[0],fvx[1],fvx[2],fvx[3]);
    rgi.setRgt(gt);
    rgi.setScales(0.2f,1.0f);
    float[][][][] ftpq = rgi.gridBlend(s1,s2,s3);
    return ftpq[2];
  }

  public float[][] averageVelAtWellsA(Sampling s1, SynSeis.Model[] mds) {
    float ftm =  FLT_MAX;
    float ltm = -FLT_MAX;
    int nl = mds.length;
    float f1 = (float)s1.getFirst();
    float d1 = (float)s1.getDelta();
    for (int il=0; il<nl; ++il) {
      SynSeis.Model md = mds[il];
      float ft = f1+round((min(md.t)-f1)/d1)*d1+d1;
      int nt = round((max(md.t)-ft)/d1);
      float lt = ft+nt*d1;
      if(ft<ftm){ftm=ft;}
      if(lt>ltm){ltm=lt;}
    }
    int ntm = round((ltm-ftm)/d1)+1;
    float[][] fvx = new float[nl][ntm];
    float[][] ft1 = new float[nl][ntm];
    float[][] fx2 = new float[nl][ntm];
    float[][] fx3 = new float[nl][ntm];
    for (int il=0; il<nl; ++il) {
      SynSeis.Model md = mds[il];
      Sampling sz = md.sz;
      float ft = f1+round((min(md.t)-f1)/d1)*d1+d1;
      int nt = round((max(md.t)-ft)/d1);
      int jt = round((ft-ftm)/d1);
      Sampling st = new Sampling(nt,d1,ft);
      InverseInterpolator ii = new InverseInterpolator(sz,st);
      float[] tz = md.t;
      float[] zt = new float[nt];
      float x2 = (float)md.x2;
      float x3 = (float)md.x3;
      for (int i=0; i<md.nz; ++i)
        if(i>0 && tz[i]<=tz[i-1])
          tz[i] = tz[i-1]+0.000001f;
      ii.invert(tz,zt);
      float dv = 0.0f;
      for (int it=0; it<nt; ++it) {
        float zi = zt[it];
        float ti = (float)st.getValue(it);
        float vi = 2f*zi/ti;
        fvx[il][it+jt] = vi;
        ft1[il][it+jt] = ti;
        fx2[il][it+jt] = x2;
        fx3[il][it+jt] = x3;
        if(it>0)
          dv += fvx[il][it+jt]-fvx[il][it+jt-1];
      }
      dv /= (nt-1f);
      System.out.println("dv="+dv);
      for (int j=jt-1; j>=0; --j) {
        fvx[il][j] = fvx[il][j+1]-dv;
        ft1[il][j] = ft1[il][j+1]-d1;
        fx2[il][j] = x2;
        fx3[il][j] = x3;
      }
      for (int j=nt+jt; j<ntm; ++j) {
        fvx[il][j] = fvx[il][j-1]+dv;
        ft1[il][j] = ft1[il][j-1]+d1;
        fx2[il][j] = x2;
        fx3[il][j] = x3;
      }
    }
    FloatList fvxl = new FloatList();
    FloatList ft1l = new FloatList();
    FloatList fx2l = new FloatList();
    FloatList fx3l = new FloatList();
    for (int il=0; il<nl; ++il) {
    for (int it=0; it<ntm; ++it) {
      fvxl.add(fvx[il][it]);
      ft1l.add(ft1[il][it]);
      fx2l.add(fx2[il][it]);
      fx3l.add(fx3[il][it]);
    }}
    return new float[][]{fvxl.trim(),ft1l.trim(),fx2l.trim(),fx3l.trim()};
  }

  public float[][][] averageVelAtWellsC(Sampling s1, SynSeis.Model[] mds) {
    float ftm =  FLT_MAX;
    float ltm = -FLT_MAX;
    int nl = mds.length;
    float f1 = (float)s1.getFirst();
    float d1 = (float)s1.getDelta();
    for (int il=0; il<nl; ++il) {
      SynSeis.Model md = mds[il];
      float ft = f1+round((min(md.t)-f1)/d1)*d1+d1;
      int nt = round((max(md.t)-ft)/d1);
      float lt = ft+nt*d1;
      if(ft<ftm){ftm=ft;}
      if(lt>ltm){ltm=lt;}
    }
    int ntm = round((ltm-ftm)/d1)+1;
    float[][] fvx = new float[nl][ntm];
    float[][] ft1 = new float[nl][ntm];
    float[][] fx2 = new float[nl][ntm];
    float[][] fx3 = new float[nl][ntm];
    for (int il=0; il<nl; ++il) {
      SynSeis.Model md = mds[il];
      Sampling sz = md.sz;
      float ft = f1+round((min(md.t)-f1)/d1)*d1+d1;
      int nt = round((max(md.t)-ft)/d1);
      int jt = round((ft-ftm)/d1);
      Sampling st = new Sampling(nt,d1,ft);
      InverseInterpolator ii = new InverseInterpolator(sz,st);
      float[] tz = md.t;
      float[] zt = new float[nt];
      float x2 = (float)md.x2;
      float x3 = (float)md.x3;
      for (int i=0; i<md.nz; ++i)
        if(i>0 && tz[i]<=tz[i-1])
          tz[i] = tz[i-1]+0.000001f;
      ii.invert(tz,zt);
      for (int it=0; it<nt; ++it) {
        float zi = zt[it];
        float ti = (float)st.getValue(it);
        float vi = 2f*zi/ti;
        fvx[il][it+jt] = vi;
        ft1[il][it+jt] = ti;
        fx2[il][it+jt] = x2;
        fx3[il][it+jt] = x3;
      }
      float dv1 = fvx[il][jt+5]-fvx[il][jt+4];
      float dv2 = fvx[il][nt-5+jt]-fvx[il][nt+jt-6];
      for (int j=jt-1; j>=0; --j) {
        fvx[il][j] = fvx[il][j+1]-dv1;
        ft1[il][j] = ft1[il][j+1]-d1;
        fx2[il][j] = x2;
        fx3[il][j] = x3;
      }
      for (int j=nt+jt; j<ntm; ++j) {
        fvx[il][j] = fvx[il][j-1]+dv2;
        ft1[il][j] = ft1[il][j-1]+d1;
        fx2[il][j] = x2;
        fx3[il][j] = x3;
      }
    }
    return new float[][][]{fvx,ft1,fx2,fx3};
  }



  public float[][] averageVelAtWells(Sampling s1, SynSeis.Model[] mds) {
    int nl = mds.length;
    float f1 = (float)s1.getFirst();
    float d1 = (float)s1.getDelta();
    FloatList fvx = new FloatList();
    FloatList ft1 = new FloatList();
    FloatList fx2 = new FloatList();
    FloatList fx3 = new FloatList();
    for (int il=0; il<nl; ++il) {
      SynSeis.Model md = mds[il];
      Sampling sz = md.sz;
      float ft = f1+round((min(md.t)-f1)/d1)*d1+d1;
      int nt = round((max(md.t)-ft)/d1);
      Sampling st = new Sampling(nt,d1,ft);
      InverseInterpolator ii = new InverseInterpolator(sz,st);
      float[] tz = md.t;
      float[] zt = new float[nt];
      float x2 = (float)md.x2;
      float x3 = (float)md.x3;
      for (int i=0; i<md.nz; ++i)
        if(i>0 && tz[i]<=tz[i-1])
          tz[i] = tz[i-1]+0.000001f;
      ii.invert(tz,zt);
      for (int it=0; it<nt; ++it) {
        float zi = zt[it];
        float ti = (float)st.getValue(it);
        float vi = 2f*zi/ti;
        fvx.add(vi);
        ft1.add(ti);
        fx2.add(x2);
        fx3.add(x3);
      }
    }
    return new float[][]{fvx.trim(),ft1.trim(),fx2.trim(),fx3.trim()};
  }

  public float[][][] averageVelAtWellsX(Sampling s1, SynSeis.Model[] mds) {
    int nl = mds.length;
    float[][][] sps = new float[4][nl][];
    float f1 = (float)s1.getFirst();
    float d1 = (float)s1.getDelta();
    for (int il=0; il<nl; ++il) {
      SynSeis.Model md = mds[il];
      Sampling sz = md.sz;
      System.out.println("fz="+sz.getFirst());
      float ft = f1+round((min(md.t)-f1)/d1)*d1+d1;
      int nt = round((max(md.t)-ft)/d1);
      Sampling st = new Sampling(nt,d1,ft);
      InverseInterpolator ii = new InverseInterpolator(sz,st);
      float[] ts = md.t;
      float[] zs = new float[nt];
      float x2 = (float)md.x2;
      float x3 = (float)md.x3;
      sps[0][il] = new float[nt];
      sps[1][il] = new float[nt];
      sps[2][il] = new float[nt];
      sps[3][il] = new float[nt];
      for (int i=0; i<md.nz; ++i)
        if(i>0 && ts[i]<=ts[i-1])
          ts[i] = ts[i-1]+0.000001f;
      ii.invert(ts,zs);
      for (int it=0; it<nt; ++it) {
        float zi = zs[it];
        float ti = (float)st.getValue(it);
        float vi = 2.0f*zi/ti;
        sps[0][il][it] = vi;
        sps[1][il][it] = ti;
        sps[2][il][it] = x2;
        sps[3][il][it] = x3;
      }
    }
    return sps;
  }


  public float[][][][] imageCut(
    Sampling s1, SynSeis.Model[] mds, 
    float[][][] gx, float[][][] gt, double[] ndf) 
  {
    int n3 = gx.length;
    int n2 = gx[0].length;
    double tmin =  FLT_MAX;
    double tmax = -FLT_MAX;
    for (SynSeis.Model md:mds) {
      int nt = md.t.length;
      float tb = md.t[0];
      float te = md.t[nt-1];
      if (tb<tmin) tmin=tb;
      if (te>tmax) tmax=te;
    }
    double f1 = s1.getFirst();
    double d1 = s1.getDelta();
    int jt = (int)((tmin-f1)/d1);
    double ft = f1+jt*d1;
    int nt = (int)((tmax-ft)/d1+1);
    ndf[0] = nt; ndf[1] = d1; ndf[2] = ft;
    float[][][] gxc = copy(nt,n2,n3,jt,0,0,gx);
    float[][][] gtc = copy(nt,n2,n3,jt,0,0,gt);
    return new float[][][][]{gxc,gtc};
  }

}


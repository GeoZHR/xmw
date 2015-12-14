package slt;

import ipf.*;
import java.util.*;
import edu.mines.jtk.dsp.*;
import edu.mines.jtk.util.Parallel;
import static edu.mines.jtk.util.ArrayMath.*;

/**
 * Estimate normals of points computed from salt likelihoods, these points 
 * might be at the boundaries of salts.
 * <em>EXPERIMENTAL</em>
 *
 * @author Xinming Wu, Colorado School of Mines
 * @version 2015.12.13
 */

public class SaltNormal {

  public float[][] applyForNormals(float[][][] fx) {
    float[][] ps = setPoints(fx);
    return applyForNormals(ps);
  }

  public float[][] applyForNormals( float[][] ps){
    return applyForNormals(ps[0],ps[1],ps[2],ps[3]);
  }

  public float[][] applyForNormals(FaultCell[] cells) {
    int nc = cells.length;
    float[] x1 = new float[nc];
    float[] x2 = new float[nc];
    float[] x3 = new float[nc];
    float[] fx = new float[nc];
    for (int ic=0; ic<nc; ++ic) {
      x1[ic] = cells[ic].getX1();
      x2[ic] = cells[ic].getX2();
      x3[ic] = cells[ic].getX3();
      fx[ic] = cells[ic].getFl();
    }
    return applyForNormals(x1,x2,x3,fx);
  }


  public float[][] applyForNormals(
    float[] x1, float[] x2, float[] x3, float[] fx) 
  {
    int np = x1.length;
    float[][] x = new float[][]{x1,x2,x3};
    KdTree kt = new KdTree(x);
    float[] xms = new float[3];
    float[] xps = new float[3];
    float[][] ps = new float[7][np];
    int c = 0;
    for (int ip=0; ip<np; ++ip) {
      int d = 6;
      float x1i = x1[ip];
      float x2i = x2[ip];
      float x3i = x3[ip];
      xms[0] = x1i-d; xps[0] = x1i+d;
      xms[1] = x2i-d; xps[1] = x2i+d;
      xms[2] = x3i-d; xps[2] = x3i+d;
      int[] ids = kt.findInRange(xms,xps);
      int nd = ids.length;
      while(nd<20 && d<20) {
        d++;
        xms[0]--; xps[0]++;
        xms[1]--; xps[1]++;
        xms[2]--; xps[2]++;
        ids = kt.findInRange(xms,xps);
        nd = ids.length;
      }
      if(nd<10){continue;}
      float[] x1s = new float[nd];
      float[] x2s = new float[nd];
      float[] x3s = new float[nd];
      float[] fxs = new float[nd];
      for (int ik=0; ik<nd; ++ik) {
        int ic = ids[ik];
        x1s[ik] = x1[ic];
        x2s[ik] = x2[ic];
        x3s[ik] = x3[ic];
        fxs[ik] = fx[ic];
      }
      float[] ups = findNormals(d,x1s,x2s,x3s,fxs);
      ps[0][c] = x1i;
      ps[1][c] = x2i;
      ps[2][c] = x3i;
      ps[3][c] = ups[0];
      ps[4][c] = ups[1];
      ps[5][c] = ups[2];
      ps[6][c] = ups[3];
      c++;
    }
    return ps;
  }

  public float[][] setPoints(float[][][] f) {
    int n3 = f.length;
    int n2 = f[0].length;
    int n1 = f[0][0].length;
    ArrayList<Float> x1l = new ArrayList<Float>();
    ArrayList<Float> x2l = new ArrayList<Float>();
    ArrayList<Float> x3l = new ArrayList<Float>();
    ArrayList<Float> fxl = new ArrayList<Float>();
    for (int i3=0; i3<n3; ++i3) {
    for (int i2=0; i2<n2; ++i2) {
    for (int i1=0; i1<n1; ++i1) {
      float fx= f[i3][i2][i1];
      if(fx>0.0f) {
        fxl.add(fx);
        x1l.add((float)i1);
        x2l.add((float)i2);
        x3l.add((float)i3);
      }
    }}}
    int ns = fxl.size();
    float[][] ps = new float[4][ns];
    for (int is=0; is<ns; ++is) {
      ps[0][is] = x1l.get(is);
      ps[1][is] = x2l.get(is);
      ps[2][is] = x3l.get(is);
      ps[3][is] = fxl.get(is);
    }
    return ps;
  }

  private float[] findNormals(
    float sigma, float[] x1, float[] x2, float[] x3, float[] fx) 
  {
    int np = x1.length;
    float c1 = 0.0f;
    float c2 = 0.0f;
    float c3 = 0.0f;
    float cs = 0.0f;
    for (int ip=0; ip<np; ++ip) {
      float fxi = fx[ip];
      c1 += fxi*x1[ip];
      c2 += fxi*x2[ip];
      c3 += fxi*x3[ip];
      cs += fxi;
    }
    c1 /= cs;
    c2 /= cs;
    c3 /= cs;
    float sigs = 0.5f/(sigma*sigma);
    float gsci = (float)(1.0/(sigma*sqrt(2*Math.PI)));
    double[][] a = new double[3][3];
    double[][] z = new double[3][3];
    double[] e = new double[3];
    for (int ip=0; ip<np; ++ip) {
      float fxi = fx[ip];
      float x1i = x1[ip]-c1;
      float x2i = x2[ip]-c2;
      float x3i = x3[ip]-c3;
      float x1s = x1i*x1i;
      float x2s = x2i*x2i;
      float x3s = x3i*x3i;
      float xsi = x1s+x2s+x3s;
      float sci = fxi*gsci*exp(-xsi*sigs);
      a[0][0] += x1s*sci;
      a[1][1] += x2s*sci;
      a[2][2] += x3s*sci;
      a[0][1] += x1i*x2i*sci;
      a[0][2] += x1i*x3i*sci;
      a[1][2] += x2i*x3i*sci;
    }
    a[1][0] = a[0][1];
    a[2][0] = a[0][2];
    a[2][1] = a[1][2];
    Eigen.solveSymmetric33(a,z,e);
    float w1i = (float)z[2][0];
    float w2i = (float)z[2][1];
    float w3i = (float)z[2][2];
    if (w1i>0.0f) {
      w1i = -w1i;
      w2i = -w2i;
      w3i = -w3i;
    }
    float eui = (float)e[0];
    float evi = (float)e[1];
    float ewi = (float)e[2];
    if (ewi<0.0f) ewi = 0.0f;
    if (evi<ewi) evi = ewi;
    if (eui<evi) eui = evi;
    float esi = (eui>0.0f)?1.0f/eui:1.0f;
    float epi = (eui-evi)*esi;
    return new float[]{w1i,w2i,w3i,epi};
  } 

}

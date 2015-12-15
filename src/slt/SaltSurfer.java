package slt;

import java.util.*;

import edu.mines.jtk.dsp.*;
import edu.mines.jtk.util.*;
import static edu.mines.jtk.util.ArrayMath.*;


/**
 * Extract salt boundary surfaces from oriented salt points. 
 * <em>EXPERIMENTAL</em>
 *
 * @author Xinming Wu, Colorado School of Mines
 * @version 2015.12.10
 */

public class SaltSurfer {

  public float[][][] findSurfacesX(
    float[][][] fx, float[][][] u1, float[][][] u2, float[][][] u3) 
  {
    int n3 = fx.length;
    int n2 = fx[0].length;
    int n1 = fx[0][0].length;
    float[][][] fp = copy(fx);
    float[][][] g1 = new float[n3][n2][n1];
    float[][][] g2 = new float[n3][n2][n1];
    float[][][] g3 = new float[n3][n2][n1];
    RecursiveGaussianFilter rgf = new RecursiveGaussianFilter(4.0);
    rgf.apply100(fp,g1);
    rgf.apply010(fp,g2);
    rgf.apply001(fp,g3);
    mul(g1,u1,g1);
    mul(g2,u2,g2);
    mul(g3,u3,g3);
    add(g1,g2,g2);
    add(g2,g3,g3);
    for (int i3=0; i3<n3; ++i3) {
    for (int i2=0; i2<n2; ++i2) {
    for (int i1=0; i1<n1; ++i1) {
      if (fx[i3][i2][i1]<0.3f) {
        g3[i3][i2][i1] = 0f;
      }
    }}}
    System.out.println("test1");
    float[][][] g11 = new float[n3][n2][n1];
    float[][][] g12 = new float[n3][n2][n1];
    float[][][] g13 = new float[n3][n2][n1];
    float[][][] g22 = new float[n3][n2][n1];
    float[][][] g23 = new float[n3][n2][n1];
    float[][][] g33 = new float[n3][n2][n1];
    float[][] xs = findPoints(fx,u1,u2,u3);
    int np = xs[0].length;
    for (int ip=0; ip<np; ++ip) {
      int i1 = round(xs[0][ip]);
      int i2 = round(xs[1][ip]);
      int i3 = round(xs[2][ip]);
      float sc = xs[3][ip];
      float u1i = u1[i3][i2][i1];
      float u2i = u2[i3][i2][i1];
      float u3i = u3[i3][i2][i1];
      g11[i3][i2][i1] = u1i*u1i*sc;
      g12[i3][i2][i1] = u1i*u2i*sc;
      g13[i3][i2][i1] = u1i*u3i*sc;
      g22[i3][i2][i1] = u2i*u2i*sc;
      g23[i3][i2][i1] = u2i*u3i*sc;
      g33[i3][i2][i1] = u3i*u3i*sc;
    }
    System.out.println("test2");
    RecursiveGaussianFilter rgf1 = new RecursiveGaussianFilter(8.0);
    rgf1.apply000(g11,g11);
    rgf1.apply000(g12,g12);
    rgf1.apply000(g13,g13);
    rgf1.apply000(g22,g22);
    rgf1.apply000(g23,g23);
    rgf1.apply000(g33,g33);
    System.out.println("test3");
    EigenTensors3 ets = applyForTensors(g11,g12,g13,g22,g23,g33);
    LocalSmoothingFilter lsf = new LocalSmoothingFilter();
    float[][][] gs = new float[n3][n2][n1];
    lsf.apply(ets,20,g3,gs);
    System.out.println("test4");
    return gs;
  }

  public EigenTensors3 applyForTensors(
    float[][][] g11, float[][][] g12, float[][][] g13,
    float[][][] g22, float[][][] g23, float[][][] g33) 
  {
    int n3 = g11.length;
    int n2 = g11[0].length;
    int n1 = g11[0][0].length;
    float[][][] u2 = new float[n3][n2][n1];
    float[][][] u3 = new float[n3][n2][n1];
    float[][][] w1 = new float[n3][n2][n1];
    float[][][] w2 = new float[n3][n2][n1];
    float[][][] eu = fillfloat(0.001f,n1,n2,n3);
    float[][][] ev = fillfloat(1.000f,n1,n2,n3);
    float[][][] ew = fillfloat(1.000f,n1,n2,n3);
    solveEigenproblems(g11,g12,g13,g22,g23,g33,u2,u3,w1,w2);
    // Compute u1 such that u3 > 0.
    float[][][] u1 = u3;
    for (int i3=0; i3<n3; ++i3) {
    for (int i2=0; i2<n2; ++i2) {
    for (int i1=0; i1<n1; ++i1) {
      float u2i = u2[i3][i2][i1];
      float u3i = u3[i3][i2][i1];
      float u1s = 1.0f-u2i*u2i-u3i*u3i;
      float u1i = (u1s>0.0f)?sqrt(u1s):0.0f;
      if (u3i<0.0f) {
        u1i = -u1i;
        u2i = -u2i;
      }
      u1[i3][i2][i1] = u1i;
      u2[i3][i2][i1] = u2i;
    }}}
    return new EigenTensors3(u1,u2,w1,w2,eu,ev,ew,true);
  }


  private void solveEigenproblems(
    final float[][][] g11, final float[][][] g12, final float[][][] g13,
    final float[][][] g22, final float[][][] g23, final float[][][] g33,
    final float[][][] u2, final float[][][] u3, final float[][][] w1, 
    final float[][][] w2) 
  {
    final int n1 = g11[0][0].length;
    final int n2 = g11[0].length;
    final int n3 = g11.length;
    Parallel.loop(n3,new Parallel.LoopInt() {
      public void compute(int i3) {
        double[][] a = new double[3][3];
        double[][] z = new double[3][3];
        double[] e = new double[3];
        for (int i2=0; i2<n2; ++i2) {
          for (int i1=0; i1<n1; ++i1) {
            a[0][0] = g11[i3][i2][i1];
            a[0][1] = g12[i3][i2][i1];
            a[0][2] = g13[i3][i2][i1];
            a[1][0] = g12[i3][i2][i1];
            a[1][1] = g22[i3][i2][i1];
            a[1][2] = g23[i3][i2][i1];
            a[2][0] = g13[i3][i2][i1];
            a[2][1] = g23[i3][i2][i1];
            a[2][2] = g33[i3][i2][i1];
            Eigen.solveSymmetric33(a,z,e);
            u2[i3][i2][i1] = (float)z[0][1];
            u3[i3][i2][i1] = (float)z[0][2];
            w1[i3][i2][i1] = (float)z[2][0];
            w2[i3][i2][i1] = (float)z[2][1];
          }
        }
      }
    });
  }


   

  public float[][][] findSurfaces(
    float[][][] fx, float[][][] u1, float[][][] u2, float[][][] u3) 
  {
    int n3 = fx.length;
    int n2 = fx[0].length;
    int n1 = fx[0][0].length;
    float[][] xs = findPoints(fx,u1,u2,u3);
    SaltNormal sn = new SaltNormal();
    float[][] xus = sn.applyForNormals(xs);
    PointSetSurface pss = new PointSetSurface();
    return pss.findScalarField(n1,n2,n3,xus);
  }

  public void findSurfaces(float[][] xus) {

  }

  public float[][] findPoints (
    final float[][][] ss, 
    final float[][][] u1, final float[][][] u2, final float[][][] u3) 
  {
    final int n3 = ss.length;
    final int n2 = ss[0].length;
    final int n1 = ss[0][0].length;
    final Sampling s1 = new Sampling(n1);
    final Sampling s2 = new Sampling(n2);
    final Sampling s3 = new Sampling(n3);
    final float[][] xps = new float[4][n2*n3*20];
    final SincInterpolator si = new SincInterpolator();
    int k = 0;
    //Parallel.loop(n3,new Parallel.LoopInt() {
    //public void compute(int i3) {
      for (int i3=0; i3<n3 ;++i3) {
      for (int i2=0; i2<n2 ;++i2) {
      for (int i1=0; i1<n1 ;++i1) {
        float u1i = u1[i3][i2][i1];
        float u2i = u2[i3][i2][i1];
        float u3i = u3[i3][i2][i1];
        float x1m = i1-u1i;
        float x2m = i2-u2i;
        float x3m = i3-u3i;
        float x1p = i1+u1i;
        float x2p = i2+u2i;
        float x3p = i3+u3i;
        float sxi = ss[i3][i2][i1];
        float sxm = si.interpolate(s1,s2,s3,ss,x1m,x2m,x3m);
        float sxp = si.interpolate(s1,s2,s3,ss,x1p,x2p,x3p);
        if (sxi>sxm && sxi>sxp && sxi>0.3f) {
          xps[0][k] = i1;
          xps[1][k] = i2;
          xps[2][k] = i3;
          xps[3][k] = sxi;
          k++;
        }
      }}}
    //}});
    return copy(k,4,xps);
  }
  
}

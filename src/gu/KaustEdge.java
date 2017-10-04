package gu;

import edu.mines.jtk.dsp.*;
import edu.mines.jtk.util.*;
import static edu.mines.jtk.util.ArrayMath.*;

/**
 * Kaust edge detection
 * @author Xinming Wu, University of Texas at Austin
 * @version 2017.08.11
 */
public class KaustEdge {

  public float[][][][] findSlopes(float sig1, float sig2, float[][][] gx) {
    final int n3 = gx.length;
    final int n2 = gx[0].length;
    final int n1 = gx[0][0].length;
    final float[][][] p2 = new float[n3][n2][n1];
    final float[][][] p3 = new float[n3][n2][n1];
    final PlaneWaveDestructor pd = new PlaneWaveDestructor(-3,3);
    pd.setSmoothness(sig1,sig2);
    Parallel.loop(n3,new Parallel.LoopInt() {
    public void compute(int i3) {
      p3[i3] = pd.findSlopes(gx[i3]);
    }});
    Parallel.loop(n2,new Parallel.LoopInt() {
    public void compute(int i2) {
      float[][] g2i = new float[n3][n1];
      for (int i3=0; i3<n3; ++i3)
        g2i[i3] = gx[i3][i2];
      float[][] p2i = pd.findSlopes(g2i);
      for (int i3=0; i3<n3; ++i3)
        p2[i3][i2] = p2i[i3];
    }});
    return new float[][][][]{p2,p3};
  }

  public float[][][] applyDestruction2(
    float[][][] gx, float[][][] p2) {
    final int n3 = gx.length;
    final int n2 = gx[0].length;
    final int n1 = gx[0][0].length;
    final float[][][] d2 = new float[n3][n2][n1];
    final PlaneWaveDestructor pd = new PlaneWaveDestructor(-3,3);
    Parallel.loop(n2,new Parallel.LoopInt() {
    public void compute(int i2) {
      float[][] g2i = new float[n3][n1];
      float[][] p2i = new float[n3][n1];
      for (int i3=0; i3<n3; ++i3) {
        g2i[i3] = gx[i3][i2];
        p2i[i3] = p2[i3][i2];
      }
      float[][] d2i = pd.applyFilter(p2i,g2i);
      for (int i3=0; i3<n3; ++i3)
        d2[i3][i2] = d2i[i3];
    }});
    return d2;
  }

  public float[][][] applyDestruction3(
    float[][][] gx, float[][][] p3) {
    final int n3 = gx.length;
    final int n2 = gx[0].length;
    final int n1 = gx[0][0].length;
    final float[][][] d3 = new float[n3][n2][n1];
    final PlaneWaveDestructor pd = new PlaneWaveDestructor(-3,3);
    Parallel.loop(n3,new Parallel.LoopInt() {
    public void compute(int i3) {
      d3[i3] = pd.applyFilter(p3[i3],gx[i3]);
    }});
    return d3;
  }

  public float[][][][] amplitudeCurvature(
    float[][][] gx, float[][][] p2, float[][][] p3) {
    final int n3 = gx.length;
    final int n2 = gx[0].length;
    final int n1 = gx[0][0].length;
    final float[][][] d2  = applyDestruction2(gx,p2);
    final float[][][] d3  = applyDestruction3(gx,p3);
    final float[][][] s2  = applyDestruction2(d2,p2);
    final float[][][] s3  = applyDestruction3(d3,p3);
    final float[][][] s23 = applyDestruction3(d2,p3);
    final float[][][] s32 = applyDestruction2(d3,p2);
    final float[][][] pc  = new float[n3][n2][n1];
    final float[][][] nc  = new float[n3][n2][n1];
    Parallel.loop(n3,new Parallel.LoopInt() {
    public void compute(int i3) {
    for (int i2=0; i2<n2; ++i2) {
    for (int i1=0; i1<n1; ++i1) {
      float s2i = s2[i3][i2][i1]*0.5f;
      float s3i = s3[i3][i2][i1]*0.5f;
      float d23 = (s23[i3][i2][i1]+s32[i3][i2][i1])*0.5f;
      d23 *= d23;
      float m23 = s2i-s3i;
      float p23 = s2i+s3i;
      pc[i3][i2][i1] = p23+sqrt(m23*m23+d23);
      nc[i3][i2][i1] = p23-sqrt(m23*m23+d23);
    }}
    }});
    return new float[][][][]{pc,nc};
  }

  public float[][][] planeWaveDestruction(float sig1, float sig2, float[][][] gx) {
    final int n3 = gx.length;
    final int n2 = gx[0].length;
    final int n1 = gx[0][0].length;
    final float[][][] d2 = new float[n3][n2][n1];
    final float[][][] d3 = new float[n3][n2][n1];
    final float[][][] ds = new float[n3][n2][n1];
    final PlaneWaveDestructor pd = new PlaneWaveDestructor(-3,3);
    pd.setSmoothness(sig1,sig2);
    Parallel.loop(n3,new Parallel.LoopInt() {
    public void compute(int i3) {
      float[][] p3 = pd.findSlopes(gx[i3]);
      d3[i3] = pd.applyFilter(p3,gx[i3]);
    }});
    Parallel.loop(n2,new Parallel.LoopInt() {
    public void compute(int i2) {
      float[][] g2i = new float[n3][n1];
      for (int i3=0; i3<n3; ++i3)
        g2i[i3] = gx[i3][i2];
      float[][] p2i = pd.findSlopes(g2i);
      float[][] d2i = pd.applyFilter(p2i,g2i);
      for (int i3=0; i3<n3; ++i3)
        d2[i3][i2] = d2i[i3];
    }});
    Parallel.loop(n3,new Parallel.LoopInt() {
    public void compute(int i3) {
    for (int i2=0; i2<n2; ++i2) {
    for (int i1=0; i1<n1; ++i1) {
      float d2i = d2[i3][i2][i1];
      float d3i = d3[i3][i2][i1];
      ds[i3][i2][i1] = sqrt(d2i*d2i+d3i*d3i);
    }}
    }});
    return ds;
  }

  public float[][][] directionalDifference(
    float[][][] gx, float[][][] v1, float[][][] v2, float[][][] v3,
    float[][][] w1, float[][][] w2, float[][][] w3) {
    final int n3 = gx.length;
    final int n2 = gx[0].length;
    final int n1 = gx[0][0].length;
    final float[][][] dd = new float[n3][n2][n1];
    RecursiveGaussianFilter rgf = new RecursiveGaussianFilter(1);
    final float[][][] g1 = new float[n3][n2][n1];
    final float[][][] g2 = new float[n3][n2][n1];
    final float[][][] g3 = new float[n3][n2][n1];
    rgf.apply100(gx,g1);
    rgf.apply010(gx,g2);
    rgf.apply001(gx,g3);
    Parallel.loop(n3,new Parallel.LoopInt() {
    public void compute(int i3) {
    for (int i2=0; i2<n2; ++i2) {
    for (int i1=0; i1<n1; ++i1) {
      float v1i = v1[i3][i2][i1];
      float v2i = v2[i3][i2][i1];
      float v3i = v3[i3][i2][i1];
      float w1i = w1[i3][i2][i1];
      float w2i = w2[i3][i2][i1];
      float w3i = w3[i3][i2][i1];
      float g1i = g1[i3][i2][i1];
      float g2i = g2[i3][i2][i1];
      float g3i = g3[i3][i2][i1];
      float gvi = v1i*g1i+v2i*g2i+v3i*g3i;
      float gwi = w1i*g1i+w2i*g2i+w3i*g3i;
      dd[i3][i2][i1] = sqrt(gvi*gvi+gwi*gwi);
    }}
    }});
    return dd;
  }
}

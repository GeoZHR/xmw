package unct;

import static edu.mines.jtk.util.ArrayMath.*;

import java.util.ArrayList;

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
  
  public float[][][] constraintsFromSurfaces(float[][][] surf) {
    int ns = surf.length;
    float[][][] cs = new float[4][ns][];
    for(int is=0; is<ns; ++is) {
      float[][] ks = constraintsFromSurface(surf[is]);
      cs[0][is] = ks[0];
      cs[1][is] = ks[1];
      cs[2][is] = ks[2];
      cs[3][is] = ks[3];
    }
    return cs;
  }

  public int[][][] uncConstraints(float[][][] surf) {
    int ns = surf.length;
    int n3 = surf[0].length;
    int n2 = surf[0][0].length;
    int[][][] unc = new int[3][ns][];
    for (int is=0; is<ns; ++is) {
      ArrayList<Integer> k1s = new ArrayList<Integer>();
      ArrayList<Integer> k2s = new ArrayList<Integer>();
      ArrayList<Integer> k3s = new ArrayList<Integer>();
      for (int i3=0; i3<n3; ++i3) {
      for (int i2=0; i2<n2; ++i2) {
        int i1 = round(surf[is][i3][i2]);
        if(i1<0){continue;}
        k1s.add(i1);
        k2s.add(i2);
        k3s.add(i3);
        int np = k1s.size();
        unc[0][is] = new int[np];
        unc[1][is] = new int[np];
        unc[2][is] = new int[np];
      }}
      int ik = 0;
      for (int k1:k1s) {unc[0][is][ik] = k1;ik++;}
      ik = 0;
      for (int k2:k2s) {unc[1][is][ik] = k2;ik++;}
      ik = 0;
      for (int k3:k3s) {unc[2][is][ik] = k3;ik++;}
    }
    return unc;
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
}

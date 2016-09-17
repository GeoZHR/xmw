package aii;


import edu.mines.jtk.dsp.*;
import static edu.mines.jtk.util.ArrayMath.*;
/**
 * Method to display 3D tensors
 * @author Xinming Wu
 * @version 2016.07.30
 */



public class TensorView {

  public EigenTensors3 tensorsFromNormal(
    float[][][] u1, float[][][] u2, float[][][] u3, 
    float[][][] w1, float[][][] w2, float[][][] w3, 
    float[][][] au, float[][][] av, float[][][] aw) 
  {
    int n3 = u1.length;
    int n2 = u1[0].length;
    int n1 = u1[0][0].length;
    EigenTensors3 ets = new EigenTensors3(n1,n2,n3,true);
    for (int i3=0; i3<n3; ++i3) {
    for (int i2=0; i2<n2; ++i2) {
    for (int i1=0; i1<n1; ++i1) {
      float aui = au[i3][i2][i1];
      float avi = av[i3][i2][i1];
      float awi = aw[i3][i2][i1];
      float u1i = u1[i3][i2][i1];
      float u2i = u2[i3][i2][i1];
      float u3i = u3[i3][i2][i1];
      float w1i = w1[i3][i2][i1];
      float w2i = w2[i3][i2][i1];
      float w3i = w3[i3][i2][i1];
      ets.setEigenvalues(i1,i2,i3,aui,avi,awi);
      ets.setEigenvectorU(i1,i2,i3,u1i,u2i,u3i);
      ets.setEigenvectorW(i1,i2,i3,w1i,w2i,w3i);
    }}}
    return ets;
  }

  public float[][] applyForSegmentsW(
    float scale, EigenTensors3 et, float[][] hz) 
  {
    int n3 = et.getN3();
    int n2 = et.getN2();
    int n1 = et.getN1();
    float[][][] w1 = new float[n3][n2][n1];
    float[][][] w2 = new float[n3][n2][n1];
    float[][][] w3 = new float[n3][n2][n1];
    for (int i3=0; i3<n3; i3++) {
    for (int i2=0; i2<n2; i2++) {
    for (int i1=0; i1<n1; i1++) {
      float[] w = et.getEigenvectorW(i1,i2,i3);
      w1[i3][i2][i1] = w[0];
      w2[i3][i2][i1] = w[1];
      w3[i3][i2][i1] = w[2];
    }}}
    SincInterpolator si = new SincInterpolator();
    Sampling s1 = new Sampling(n1);
    Sampling s2 = new Sampling(n2);
    Sampling s3 = new Sampling(n3);
    float[][] ss = new float[n2*n3][6];
    int k = 0;
    int d = round(scale+1);
    for (int i3=d; i3<n3-d; i3+=d) {
    for (int i2=d; i2<n2-d; i2+=d) {
      int p = 0;
      float x1 = hz[i3][i2];
      float w1i = si.interpolate(s1,s2,s3,w1,x1,i2,i3);
      float w2i = si.interpolate(s1,s2,s3,w2,x1,i2,i3);
      float w3i = si.interpolate(s1,s2,s3,w3,x1,i2,i3);
      w1i *= scale;
      w2i *= scale;
      w3i *= scale;
      ss[k][p++] = i3-w3i;
      ss[k][p++] = i2-w2i;
      //ss[k][p++] = x1-w1i;
      ss[k][p++] = si.interpolate(s2,s3,hz,i2-w2i,i3-w3i);
      ss[k][p++] = i3+w3i;
      ss[k][p++] = i2+w2i;
      //ss[k][p++] = x1+w1i;
      ss[k][p++] = si.interpolate(s2,s3,hz,i2+w2i,i3+w3i);
      k++;
    }}
    return copy(6,k,0,0,ss);
  }

  public float[][] applyForSegmentsV(
    float scale, EigenTensors3 et, float[][] hz) 
  {
    int n3 = et.getN3();
    int n2 = et.getN2();
    int n1 = et.getN1();
    float[][][] v1 = new float[n3][n2][n1];
    float[][][] v2 = new float[n3][n2][n1];
    float[][][] v3 = new float[n3][n2][n1];
    for (int i3=0; i3<n3; i3++) {
    for (int i2=0; i2<n2; i2++) {
    for (int i1=0; i1<n1; i1++) {
      float[] v = et.getEigenvectorV(i1,i2,i3);
      v1[i3][i2][i1] = v[0];
      v2[i3][i2][i1] = v[1];
      v3[i3][i2][i1] = v[2];
    }}}
    SincInterpolator si = new SincInterpolator();
    Sampling s1 = new Sampling(n1);
    Sampling s2 = new Sampling(n2);
    Sampling s3 = new Sampling(n3);
    float[][] ss = new float[n2*n3][6];
    int k = 0;
    int d = round(scale+1);
    for (int i3=d; i3<n3-d; i3+=d) {
    for (int i2=d; i2<n2-d; i2+=d) {
      int p = 0;
      float x1 = hz[i3][i2];
      float v1i = si.interpolate(s1,s2,s3,v1,x1,i2,i3);
      float v2i = si.interpolate(s1,s2,s3,v2,x1,i2,i3);
      float v3i = si.interpolate(s1,s2,s3,v3,x1,i2,i3);
      v1i *= scale;
      v2i *= scale;
      v3i *= scale;
      ss[k][p++] = i3-v3i;
      ss[k][p++] = i2-v2i;
      //ss[k][p++] = x1-v1i;
      ss[k][p++] = si.interpolate(s2,s3,hz,i2-v2i,i3-v3i);
      ss[k][p++] = i3+v3i;
      ss[k][p++] = i2+v2i;
      ss[k][p++] = si.interpolate(s2,s3,hz,i2+v2i,i3+v3i);
      //ss[k][p++] = x1+v1i;
      k++;
    }}
    return copy(6,k,0,0,ss);
  }


   
}

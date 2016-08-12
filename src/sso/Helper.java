package sso;


import edu.mines.jtk.dsp.*;
import static edu.mines.jtk.util.ArrayMath.*;

/**
 * Method to display 3D tensors
 * @author Xinming Wu
 * @version 2016.07.30
 */



public class Helper {

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

  public float[][] channelAzimuth(
    float[][][] w2, float[][][] w3, float[][] hz) {
    int n3 = w2.length;
    int n2 = w3[0].length;
    float[][] ha = fillfloat(100f,n2,n3);
    float pi = (float)Math.PI;
    float dp = 2f*pi/(float)n2;
    for (int i2=0; i2<n2; ++i2) {
      float x3 = 0.6f*n3+20f*sin(dp*i2)+30;
      int j3 = (int)(x3+0.5);
    for (int i3=j3-2; i3<j3+2; ++i3) {
      float hi = hz[i3][i2];
      int i1 = round(hi);
      float w2i = w2[i3][i2][i1];
      float w3i = w3[i3][i2][i1];
      if(w2i<0f) {w2i = -w2i;w3i = -w3i;}
      ha[i3][i2] = atan(w3i/w2i);
    }}
    return ha;
  }

  public float[][][] channelAzimuth(
    float[][][] w2, float[][][] w3) {
    int n3 = w2.length;
    int n2 = w3[0].length;
    int n1 = w3[0][0].length;
    float[][][] az = new float[n3][n2][n1];
    for (int i3=0; i3<n3; ++i3) {
    for (int i2=0; i2<n2; ++i2) {
    for (int i1=0; i1<n1; ++i1) {
      float w2i = w2[i3][i2][i1];
      float w3i = w3[i3][i2][i1];
      if(w2i<0f) {w2i = -w2i;w3i = -w3i;}
      az[i3][i2][i1] = atan(w3i/w2i);
    }}}
    return az;
  }

  public float[][] channelAzimuth(
    float[][] w2, float[][] w3, float[][] hz) {
    int n3 = w2.length;
    int n2 = w3[0].length;
    float[][] ha = fillfloat(100f,n2,n3);
    float pi = (float)Math.PI;
    float dp = 2f*pi/(float)n2;
    for (int i2=0; i2<n2; ++i2) {
      float x3 = 0.6f*n3+20f*sin(dp*i2)+30;
      int j3 = (int)(x3+0.5);
    for (int i3=j3-2; i3<j3+2; ++i3) {
      float w2i = w2[i3][i2];
      float w3i = w3[i3][i2];
      if(w2i<0f) {w2i = -w2i;w3i = -w3i;}
      ha[i3][i2] = atan(w3i/w2i);
    }}
    return ha;
  }


  public float[][] channelAzimuth(float[][] hz) {
    int n3 = hz.length;
    int n2 = hz[0].length;
    float[][] ha = fillfloat(100f,n2,n3);
    float pi = (float)Math.PI;
    float dp = 2f*pi/(float)n2;
    for (int i2=0; i2<n2; ++i2) {
      float x3 = 0.6f*n3+20f*sin(dp*i2)+30;
      int j3 = (int)(x3+0.5);
    for (int i3=j3-2; i3<j3+2; ++i3) {
      ha[i3][i2] = atan(20*2f*pi*cos(dp*i2)/(float)n2);
    }}
    return ha;
  }


   
}

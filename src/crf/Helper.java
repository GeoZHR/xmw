package crf;

import util.*;
import ipfx.*;
import static edu.mines.jtk.util.ArrayMath.*;

public class Helper {

  public void setWeights(FaultSkin[] skins, float[][][] wp) {
    int n3 = wp.length;
    int n2 = wp[0].length;
    int n1 = wp[0][0].length;
    float[][][] fl = new float[n3][n2][n1];
    for (FaultSkin skin:skins) {
    for (FaultCell cell:skin) {
      int i1 = cell.getI1();
      int i2 = cell.getI2();
      int i3 = cell.getI3();
      fl[i3][i2][i1] = cell.getFl();
    }}
    RecursiveGaussianFilterP rgf = new RecursiveGaussianFilterP(2.0);
    rgf.apply000(fl,fl);
    for (int i3=0; i3<n3; ++i3) {
    for (int i2=0; i2<n2; ++i2) {
    for (int i1=0; i1<n1; ++i1) {
      float fli = 1f-fl[i3][i2][i1];
      fli *= fli;
      fli *= fli;
      fli *= fli;
      wp[i3][i2][i1] *= fli;
    }}}
  }

  public float[][][] resample(float[][][] gx) {
    int n3 = gx.length;
    int n2 = gx[0].length;
    int n1 = gx[0][0].length;
    int m2 = round(n2/2f);
    float[][][] fx = new float[n3][m2][n1];
    for (int i3=0; i3<n3; i3++) {
      int k2=-1;
      for (int i2=0; i2<n2; i2+=2) {
        k2++;
      for (int i1=0; i1<n1; ++i1) {
        fx[i3][k2][i1] = gx[i3][i2][i1];
      }}
    }
    return fx;
  }

  public void rotate(float phi, float[][][] fps) {
    int n3 = fps.length;
    int n2 = fps[0].length;
    int n1 = fps[0][0].length;
    for (int i3=0; i3<n3; ++i3) {
    for (int i2=0; i2<n2; ++i2) {
    for (int i1=0; i1<n1; ++i1) {
      float fpi = fps[i3][i2][i1];
      if(fpi>=0.0f) {
        fpi += phi;
        if (fpi>=360f) fpi-=360f;
        fps[i3][i2][i1] = fpi;
      }
    }}}
  }


  public void convert(float[][][] fps) {
    int n3 = fps.length;
    int n2 = fps[0].length;
    int n1 = fps[0][0].length;
    for (int i3=0; i3<n3; ++i3) {
    for (int i2=0; i2<n2; ++i2) {
    for (int i1=0; i1<n1; ++i1) {
      float fpi = fps[i3][i2][i1];
      if(fpi>180.0f) {
        fpi -= 180f;
        fps[i3][i2][i1] = fpi;
      }
    }}}
  }


  public float[][] getOceanBottom(float dv, float[][][] gx) {
    int n3 = gx.length;
    int n2 = gx[0].length;
    int n1 = gx[0][0].length;
    float[][] ob = new float[n3][n2];
    for (int i3=0; i3<n3; ++i3) {
    for (int i2=0; i2<n2; ++i2) {
    for (int i1=0; i1<n1-1; ++i1) {
      float dx = gx[i3][i2][i1+1]-gx[i3][i2][i1];
      if(abs(dx)>dv) {
        ob[i3][i2] = i1;
        i1 = n1;
      }
    }}}
    RecursiveGaussianFilterP rgf = new RecursiveGaussianFilterP(3.0);
    rgf.apply00(ob,ob);
    return ob;
  }

}

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

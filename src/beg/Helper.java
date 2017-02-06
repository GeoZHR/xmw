package beg;

import edu.mines.jtk.dsp.*;
import util.*;
import ipfx.*;
import java.io.*;

public class Helper {

  public void writeAsciiHorizon(String name, float[][] hz) throws IOException  {
    int n3 = hz.length;
    int n2 = hz[0].length;
    PrintWriter writer = new PrintWriter(name,"UTF-8");
    for (int i2= 0; i2<n2; i2++) {
      String row = String.valueOf(hz[0][i2])+"  ";
      for (int i3= 1; i3<n3; i3++) {
        row += String.valueOf(hz[i3][i2])+"  ";
      }
      writer.println(row);
    }
    writer.close();
  }

  public float[][][] convertToAsciiHorizons(
    Sampling s1, Sampling s2, float f3, float d3, float[][][] hs) 
  {
    int ns = hs.length;
    int n3 = hs[0].length;
    int n2 = hs[0][0].length;
    float d1 = (float)s1.getDelta();
    float f1 = (float)s1.getFirst();
    float[][][] ha = new float[ns][3][n2*n3];
    for (int is=0; is<ns; ++is) {
      int ik = 0;
    for (int i3=0; i3<n3; ++i3) {
    for (int i2=0; i2<n2; ++i2) {
      ha[is][0][ik] = f3+i3*d3;
      ha[is][1][ik] = (float)s2.getValue(i2);
      ha[is][2][ik] = hs[is][i3][i2]*d1+f1;
      ik++;
    }}}
    return ha;
  }
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

  public float[][][][] slopesFromNormals(float pmax, 
    float[][][] u1, float[][][] u2, float[][][] u3) {
    int n3 = u1.length;
    int n2 = u1[0].length;
    int n1 = u1[0][0].length;
    float p2min = -pmax;
    float p3min = -pmax;
    float p2max =  pmax;
    float p3max =  pmax;
    float[][][] p2 = new float[n3][n2][n1];
    float[][][] p3 = new float[n3][n2][n1];
    for (int i3=0; i3<n3; ++i3) {
    for (int i2=0; i2<n2; ++i2) {
    for (int i1=0; i1<n1; ++i1) {
      float u1i = u1[i3][i2][i1];
      float u2i = u2[i3][i2][i1];
      float u3i = u3[i3][i2][i1];
      if (u1i<0f) {
        u1i = -u1i;
        u2i = -u2i;
        u3i = -u3i;
      }
      if (-u2i<p2min*u1i) u2i = -p2min*u1i;
      if (-u2i>p2max*u1i) u2i = -p2max*u1i;
      if (-u3i<p3min*u1i) u3i = -p3min*u1i;
      if (-u3i>p3max*u1i) u3i = -p3max*u1i;
      if (u1i==0.0f) {
        p2[i3][i2][i1] = (u2i<0.0f)?p2max:p2min;
        p3[i3][i2][i1] = (u3i<0.0f)?p3max:p3min;
      } else {
        p2[i3][i2][i1] = -u2i/u1i;
        p3[i3][i2][i1] = -u3i/u1i;
      }
    }}}
    return new float[][][][]{p2,p3};
  }

}

package cfd;

import util.*;
import edu.mines.jtk.dsp.*;
import static edu.mines.jtk.util.ArrayMath.*;

public class TensorMaker {


  public EigenTensors3 applyForTensors(
    float sigma1, float sigma2, float sigma3, float[][][] mk, float[][][] gx) 
  {
    int n3 = mk.length;
    int n2 = mk[0].length;
    int n1 = mk[0][0].length;
    LocalOrientFilterP lof = new LocalOrientFilterP(sigma1,sigma2,sigma3);
    EigenTensors3 ets = lof.applyForTensors(gx);
    for (int i3=0; i3<n3; ++i3) {
    for (int i2=0; i2<n2; ++i2) {
    for (int i1=0; i1<n1; ++i1) {
      if (mk[i3][i2][i1]==0f) {
        ets.setEigenvectorU(i1,i2,i3,1f,0f,0f);
        ets.setEigenvectorW(i1,i2,i3,0f,0f,1f);
      }
    }}}
    return ets;
  }

  public float[][][] mask(
    double small, double sigma1, double sigma2, double sigma3,
    float[][][] x) 
  {
    int n3 = x.length;
    int n2 = x[0].length;
    int n1 = x[0][0].length;
    float[][][] t = abs(x);
    float a = ((sum(t)/n1)/n2)/n3; // global mean absolute amplitude
    RecursiveGaussianFilterP rgf1 = new RecursiveGaussianFilterP(sigma1);
    RecursiveGaussianFilterP rgf2 = new RecursiveGaussianFilterP(sigma2);
    RecursiveGaussianFilterP rgf3 = new RecursiveGaussianFilterP(sigma3);
    float[][][] b = zerofloat(n1,n2,n3);
    rgf1.apply0XX(t,b);
    rgf2.applyX0X(b,t);
    rgf3.applyXX0(t,b); // local mean absolute amplitude
    float[][][] mask = new float[n3][n2][n1];
    for (int i3=0; i3<n3; ++i3) {
    for (int i2=0; i2<n2; ++i2) {
    for (int i1=0; i1<n1; ++i1) {
      if (b[i3][i2][i1]>small*a) {
        mask[i3][i2][i1] = 1;
      }
    }}}
    return mask;
  }

}

package figi;


import edu.mines.jtk.dsp.*;
import edu.mines.jtk.interp.*;
import static edu.mines.jtk.util.ArrayMath.*;


public class Helper {

  public float[][][] sliceExtraction(
    float[] c2, float[] c3, float[][][] gx, float[][][] wx) {
    int n3 = gx.length;
    int n2 = gx[0].length;
    int n1 = gx[0][0].length;
    float[][] gs = new float[n2][n1];
    float[][] ws = new float[n2][n1];
    Sampling s1 = new Sampling(n1);
    Sampling s2 = new Sampling(n2);
    Sampling s3 = new Sampling(n3);
    SincInterpolator si = new SincInterpolator();
    CubicInterpolator ci = new CubicInterpolator(c2,c3);
    for (int i2=0; i2<n2; ++i2) {
      float x3 = ci.interpolate(i2);
      for (int i1=0; i1<n1; ++i1) {
        gs[i2][i1] = si.interpolate(s1,s2,s3,gx,i1,i2,x3);
        int k3 = round(x3);
        k3 = min(k3,n3-1);
        k3 = max(k3,0);
        ws[i2][i1] = wx[k3][i2][i1];
      }
    }
    for (int i1=0; i1<n1; ++i1) {
      ws[208][i1] = 0f;
      //ws[285][i1] = 0f;
    }
    return new float[][][]{gs,ws};
  }
}

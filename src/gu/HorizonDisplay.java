package gu;

import edu.mines.jtk.awt.*;
import edu.mines.jtk.dsp.*;
import edu.mines.jtk.interp.*;
import static edu.mines.jtk.util.ArrayMath.*;



/**
 * Horizon display
 * @author Xinming Wu, University of Texas at Austin
 * @version 2017.08.11
 */
public class HorizonDisplay {
 

  public float[][][] heightRgb(
    ColorMap mp, float[][] sf) {
    int n3 = sf.length;
    int n2 = sf[0].length;
    float[] sa = new float[n3*n2];
    int k = 0;
    for (int i3=0; i3<n3; ++i3) 
    for (int i2=0; i2<n2; ++i2) 
      sa[k++] = sf[i3][i2];
    float[] rgb = mp.getRgbFloats(sa);
    float[][] r = new float[n3][n2];
    float[][] g = new float[n3][n2];
    float[][] b = new float[n3][n2];
    k = 0;
    for (int i3=0; i3<n3; ++i3) {
    for (int i2=0; i2<n2; ++i2) {
      r[i3][i2] = rgb[k++];
      g[i3][i2] = rgb[k++];
      b[i3][i2] = rgb[k++];
    }}
    return new float[][][]{r,g,b};
  }

  public float[][][] amplitudeRgb(
    ColorMap mp, float[][][] fx, float[][] sf) {
    int n3 = fx.length;
    int n2 = fx[0].length;
    int n1 = fx[0][0].length;
    float[] sa = new float[n3*n2];
    int k = 0;
    SincInterpolator si = new SincInterpolator();
    for (int i3=0; i3<n3; ++i3) 
    for (int i2=0; i2<n2; ++i2) 
      sa[k++] = si.interpolate(n1,1,0,fx[i3][i2],sf[i3][i2]);
    float[] rgb = mp.getRgbFloats(sa);
    float[][] r = new float[n3][n2];
    float[][] g = new float[n3][n2];
    float[][] b = new float[n3][n2];
    k = 0;
    for (int i3=0; i3<n3; ++i3) {
    for (int i2=0; i2<n2; ++i2) {
      r[i3][i2] = rgb[k++];
      g[i3][i2] = rgb[k++];
      b[i3][i2] = rgb[k++];
    }}
    return new float[][][]{r,g,b};
  }

  public void fillHoles(float[][] sf) {
    int n2 = sf.length;
    int n1 = sf[0].length;
    float[] fx = new float[n1*n2];
    float[] x1 = new float[n1*n2];
    float[] x2 = new float[n1*n2];
    int ip = 0;
    for (int i2=0; i2<n2; ++i2) {
    for (int i1=0; i1<n1; ++i1) {
      float sfi = sf[i2][i1];
      if(sfi>0) {
        x1[ip] = i1;
        x2[ip] = i2;
        fx[ip] = sfi;
        ip++;
      }
    }}
    int np = ip;
    if (np==n1*n2) return;
    fx = copy(np,0,fx);
    x1 = copy(np,0,x1);
    x2 = copy(np,0,x2);
     RadialInterpolator2.Biharmonic basis = new RadialInterpolator2.Biharmonic();
    RadialInterpolator2 si = new RadialInterpolator2(basis,fx,x1,x2);
    for (int i2=0; i2<n2; ++i2) {
    for (int i1=0; i1<n1; ++i1) {
      float sfi = sf[i2][i1];
      if(sfi==0) {
        sf[i2][i1] = si.interpolate(i1,i2);
      }
    }}
  }

}

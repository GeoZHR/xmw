package gu;

import edu.mines.jtk.awt.*;
import edu.mines.jtk.dsp.*;
import edu.mines.jtk.interp.*;
import static edu.mines.jtk.util.ArrayMath.*;

import java.util.*;



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

  public float[][] amplitudeOnHorizon(float[][][] fx, float[][] sf) {
    int n3 = fx.length;
    int n2 = fx[0].length;
    int n1 = fx[0][0].length;
    float[][] sa = new float[n3][n2];
    SincInterpolator si = new SincInterpolator();
    for (int i3=0; i3<n3; ++i3) { 
    for (int i2=0; i2<n2; ++i2) {
      sa[i3][i2] = si.interpolate(n1,1,0,fx[i3][i2],sf[i3][i2]);
      if(sf[i3][i2]==0f) sa[i3][i2] = Float.NaN;
    }}
    return sa;
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
    RadialInterpolator2.Biharmonic basis = new RadialInterpolator2.Biharmonic();
    for (int i2=0; i2<n2; ++i2) {
    for (int i1=0; i1<n1; ++i1) {
      float sfi = sf[i2][i1];
      if(sfi==0) {
        ArrayList<Float> fxa = new ArrayList<Float>();
        ArrayList<Float> x1a = new ArrayList<Float>();
        ArrayList<Float> x2a = new ArrayList<Float>();
        int m2 = max(i2-20,0);
        int m1 = max(i1-20,0);
        int p2 = min(i2+20,n2-1);
        int p1 = min(i1+20,n1-1);
        for (int k2=m2; k2<=p2; k2+=2) {
        for (int k1=m1; k1<=p1; k1+=2) {
          float sfk = sf[k2][k1];
          if(sfk>0f) {
            fxa.add(sfk);
            x1a.add((float)k1);
            x2a.add((float)k2);
          }
        }}
        int np = fxa.size();
        float[] fx = new float[np];
        float[] x1 = new float[np];
        float[] x2 = new float[np];
        for (int ip=0; ip<np; ++ip) {
          fx[ip] = fxa.get(ip);
          x1[ip] = x1a.get(ip);
          x2[ip] = x2a.get(ip);
        }
        RadialInterpolator2 si = new RadialInterpolator2(basis,fx,x1,x2);
        //SibsonInterpolator2 si = new SibsonInterpolator2(fx,x1,x2);
        sf[i2][i1] = si.interpolate(i1,i2);
      }
    }}
  }

}

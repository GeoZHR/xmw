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
 

  public float[][] getTimeSlice(int k1, float[][][] gx) {
    int n3 = gx.length;
    int n2 = gx[0].length;
    float[][] g1 = new float[n3][n2];
    for (int i3=0; i3<n3; ++i3) {
    for (int i2=0; i2<n2; ++i2) {
      g1[i3][i2] = gx[i3][i2][k1];
    }}
    return g1;
  }

  public float[][] getTimeSliceW(int k1, int d1, float[][][] gx) {
    int n3 = gx.length;
    int n2 = gx[0].length;
    float[][] g1 = new float[n3][n2];
    for (int i3=0; i3<n3; ++i3) {
    for (int i2=0; i2<n2; ++i2) {
      float c = 0f;
      for (int i1=k1-d1; i1<=k1+d1; i1++) {
        g1[i3][i2] = gx[i3][i2][k1];
      }
    }}
    return g1;
  }


  public float[][] flip1(float[][] gx) {
    int n2 = gx.length;
    int n1 = gx[0].length;
    float[][] gf = new float[n2][n1];
    for (int i2=0; i2<n2; ++i2) {
    for (int i1=0; i1<n1; ++i1) {
      gf[i2][i1] = gx[i2][n1-i1-1];
    }}
    return gf;
  }

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

  public float[][] buildTrigs(
    Sampling sx, Sampling sy, int nz,
    float fmin, float fmax, float[][] z, float[][] f) {
    int i = 0;
    int k = 0;
    int c = 0;
    int nx = sx.getCount();
    int ny = sy.getCount();
    float[] zas = new float[nx*ny*6];
    float[] zfs = new float[nx*ny*6];
    float[] xyz = new float[nx*ny*6*3];
    for (int ix=0;ix<nx-1; ++ix) {
      float x0 = (float)sx.getValue(ix  );
      float x1 = (float)sx.getValue(ix+1);
      for (int iy=0; iy<ny-1; ++iy) {
        float y0 = (float)sy.getValue(iy  );
        float y1 = (float)sy.getValue(iy+1);
        float z1 = z[ix  ][iy  ];
        float z2 = z[ix  ][iy+1];
        float z3 = z[ix+1][iy  ];
        float z4 = z[ix+1][iy  ];
        float z5 = z[ix  ][iy+1];
        float z6 = z[ix+1][iy+1];
        if(z1<=0||z2<=0||z3<=0){continue;}
        if(z4<=0||z5<=0||z6<=0){continue;}
        if(z1>nz||z2>nz||z3>nz){continue;}
        if(z4>nz||z5>nz||z6>nz){continue;}
        zas[k++] = z1;  zas[k++] = z2;  zas[k++] =z3;
        zas[k++] = z4;  zas[k++] = z5;  zas[k++] =z6;
        if(f!=null) {
          zfs[c++] = f[ix  ][iy  ];
          zfs[c++] = f[ix  ][iy+1];
          zfs[c++] = f[ix+1][iy  ];
          zfs[c++] = f[ix+1][iy  ];
          zfs[c++] = f[ix  ][iy+1];
          zfs[c++] = f[ix+1][iy+1];
        }
        xyz[i++] = x0;  xyz[i++] = y0;  xyz[i++] = z[ix  ][iy  ];
        xyz[i++] = x0;  xyz[i++] = y1;  xyz[i++] = z[ix  ][iy+1];
        xyz[i++] = x1;  xyz[i++] = y0;  xyz[i++] = z[ix+1][iy  ];
        xyz[i++] = x1;  xyz[i++] = y0;  xyz[i++] = z[ix+1][iy  ];
        xyz[i++] = x0;  xyz[i++] = y1;  xyz[i++] = z[ix  ][iy+1];
        xyz[i++] = x1;  xyz[i++] = y1;  xyz[i++] = z[ix+1][iy+1];
      }
    }
    float[] rgb;
    zas = copy(k,0,zas);
    xyz = copy(i,0,xyz);
    float zmin = Float.MAX_VALUE;
    float zmax = -Float.MAX_VALUE;
    for (int ix=0; ix<nx; ++ix) {
    for (int iy=0; iy<ny; ++iy) {
      float zi = z[ix][iy];
      if (Float.isNaN(zi)) {continue;}
      if (zi<zmin) {zmin=zi;}
      if (zi>zmax) {zmax=zi;}
    }}
    ColorMap cp = new ColorMap(fmin,fmax,ColorMap.GRAY);
    rgb = cp.getRgbFloats(zfs);
    return new float[][]{xyz,rgb};
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

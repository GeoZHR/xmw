package crf;

import edu.mines.jtk.dsp.*;
import edu.mines.jtk.util.*;
import edu.mines.jtk.awt.ColorMap;
import static edu.mines.jtk.util.ArrayMath.*;

import ipfx.*;
import util.*;

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

  public void resample(
    final Sampling s1, final Sampling s2, final Sampling s3, 
    final float d3i, final float[][][] fx, final float[][][] fi) 
  {
    final int n3 = fi.length;
    final int n2 = fi[0].length;
    final int n1 = fi[0][0].length;
    final SincInterpolator si = new SincInterpolator();
    si.setExtrapolation(SincInterpolator.Extrapolation.CONSTANT);
    Parallel.loop(n3,new Parallel.LoopInt() {
    public void compute(int i3) {
      double x3i = i3*d3i;
      for (int i2=0; i2<n2; ++i2) {
        double x2i = s2.getValue(i2);
        for (int i1=0; i1<n1; ++i1) {
          double x1i = s1.getValue(i1);
          fi[i3][i2][i1] = si.interpolate(s1,s2,s3,fx,x1i,x2i,x3i);
        }
      }
    }});
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

  public void rotateX(float phi, float[][][] fps) {
    int n3 = fps.length;
    int n2 = fps[0].length;
    int n1 = fps[0][0].length;
    for (int i3=0; i3<n3; ++i3) {
    for (int i2=0; i2<n2; ++i2) {
    for (int i1=0; i1<n1; ++i1) {
      float fpi = fps[i3][i2][i1];
      if(fpi>=0.0f) {
        fpi = phi-fpi;
        if (fpi<0f) fpi+=360f;
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

  public void mask(float[][] ob, float[][][] fx) {
    int n3 = fx.length;
    int n2 = ob[0].length;
    for (int i3=0; i3<n3; ++i3) {
    for (int i2=0; i2<n2; ++i2) {
    for (int i1=0; i1<round(ob[i3][i2]); ++i1) {
      fx[i3][i2][i1] = -0.001f;
    }}}
  }

  public void setStrikes(float[][] ps, float[][][] fp) {
    int np = ps[0].length;
    for (int ip=0; ip<np; ++ip) {
      int k1 = (int)ps[0][ip];
      int k2 = (int)ps[1][ip];
      int k3 = (int)ps[2][ip];
      if (k1==99) {
        fp[k3][k2][k1] = ps[3][ip];
      }
    }

  }

  public void horizonToImage(float[][] hz, float[][][] fx) {
    int n3 = fx.length;
    int n2 = fx[0].length;
    int n1 = fx[0][0].length;
    for (int i3=0; i3<n3; ++i3) {
    for (int i2=0; i2<n2; ++i2) {
      int i1 = round(hz[i3][i2]);
      i1 = min(i1,n1-1);
      int im = max(i1-1,0);
      int ip = min(i1+1,n1-1);
      fx[i3][i2][i1] = 1f;
      fx[i3][i2][im] = 1f;
      fx[i3][i2][ip] = 1f;
    }}
  }

  public float[][] faultDensity(float[][] ob, float[][] sf, float[][][] fp) {
    int n3 = fp.length;
    int n2 = fp[0].length;
    float[][] fd = new float[n3][n2];
    for (int i3=0; i3<n3; ++i3) {
    for (int i2=0; i2<n2; ++i2) {
      int nc = 0;
      int b1 = round(ob[i3][i2]);
      int e1 = round(sf[i3][i2]);
      for (int i1=b1; i1<=e1; ++i1) {
        float fpi = fp[i3][i2][i1];
        if (fpi>0f) {nc++;}
      }
      fd[i3][i2] = (float)nc/(float)(e1-b1);
    }}
    return fd;
  }

  public float[][] horizonWithFaultDensity(
    int nz, float[] mfs, float[][] hz, float[][] fd) 
  {
    int nx = fd.length;
    int ny = fd[0].length;
    Sampling sx = new Sampling(nx);
    Sampling sy = new Sampling(ny);
    return buildTrigs(nz-2,sx,sy,hz,-1,mfs,fd); 
  }

  public float[][] buildTrigs(
    int nz, Sampling sx, Sampling sy, float[][] z,  
    float color, float[] mfs, float[][] f) 
  {
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
        if(Float.isNaN(z1)){continue;}
        if(Float.isNaN(z2)){continue;}
        if(Float.isNaN(z3)){continue;}
        if(Float.isNaN(z4)){continue;}
        if(Float.isNaN(z5)){continue;}
        if(Float.isNaN(z6)){continue;}
        if(z1<0||z2<0||z3<0){continue;}
        if(z4<0||z5<0||z6<0){continue;}
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
    if(color>0.0f) {
      zero(zas);
      add(zas,color,zas);
      ColorMap cp = new ColorMap(0.0f,1.0f,ColorMap.JET);
      rgb = cp.getRgbFloats(zas);
    } else if (f==null) {
      ColorMap cp = new ColorMap(-zmax,-zmin,ColorMap.JET);
      rgb = cp.getRgbFloats(mul(zas,-1f));
    } else {
      //ColorMap cp = new ColorMap(mfs[0],mfs[1],ColorMap.GRAY);
      ColorMap cp = new ColorMap(mfs[0],mfs[1],ColorMap.JET);
      rgb = cp.getRgbFloats(zfs);
    }
    return new float[][]{xyz,rgb};
  }



}

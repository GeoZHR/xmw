package uff;

import ipfx.*;
import edu.mines.jtk.dsp.*;
import edu.mines.jtk.interp.*;
import edu.mines.jtk.sgl.*;
import edu.mines.jtk.awt.ColorMap;
import static edu.mines.jtk.util.ArrayMath.*;

public class HorizonExtractor {

  public float[][][][] computeCompositeShifts(
    float[][][][] ts, float[][][][] rs) 
  {
    int n3 = ts[0].length;
    int n2 = ts[0][0].length;
    int n1 = ts[0][0][0].length;
    float[][][] r1 = rs[0];
    float[][][] r2 = rs[1];
    float[][][] r3 = rs[2];
    float[][][] t1 = ts[0];
    float[][][] t2 = ts[1];
    float[][][] t3 = ts[2];
    float[][][][] cs = new float[3][n3][n2][n1];
    for (int i3=0; i3<n3; ++i3) {
    for (int i2=0; i2<n2; ++i2) {
    for (int i1=0; i1<n1; ++i1) {
      float r1i = r1[i3][i2][i1];
      float r2i = r2[i3][i2][i1];
      float r3i = r3[i3][i2][i1];
      float w1i = i1-r1i;
      float w2i = i2-r2i;
      float w3i = i3-r3i;
      float t1i = nearestInterp(w1i,w2i,w3i,t1);
      float t2i = nearestInterp(w1i,w2i,w3i,t2);
      float t3i = nearestInterp(w1i,w2i,w3i,t3);
      cs[0][i3][i2][i1] = r1i-t1i;
      cs[1][i3][i2][i1] = r2i-t2i;
      cs[2][i3][i2][i1] = r3i-t3i;
    }}}
    return cs;
  }

  public float[][][] findInUnfaultSpace(int i1, float[][][][] rs) {
    float[][][] r1 = rs[0];
    float[][][] r2 = rs[1];
    float[][][] r3 = rs[2];
    int n3 = r1.length;
    int n2 = r1[0].length;
    float[][] w1 = new float[n3][n2];
    float[][] w2 = new float[n3][n2];
    float[][] w3 = new float[n3][n2];
    for (int i3=0; i3<n3; ++i3) {
    for (int i2=0; i2<n2; ++i2) {
      w1[i3][i2] = i1-r1[i3][i2][i1];
      w2[i3][i2] = i2-r2[i3][i2][i1];
      w3[i3][i2] = i3-r3[i3][i2][i1];

    }}
    return new float[][][]{w1,w2,w3};
  }

  public float[][][] findInCurrentSpace(int i1,
    float[][][] ws, float[][][][] ts, float[][][][] rs)
  {
    float[][] w1 = ws[0];
    float[][] w2 = ws[1];
    float[][] w3 = ws[2];
    float[][][] t1 = ts[0];
    float[][][] t2 = ts[1];
    float[][][] t3 = ts[2];
    float[][][] r1 = rs[0];
    float[][][] r2 = rs[1];
    float[][][] r3 = rs[2];
    int n3 = w1.length;
    int n2 = w1[0].length;
    float[][] x1 = new float[n3][n2];
    float[][] x2 = new float[n3][n2];
    float[][] x3 = new float[n3][n2];
    for (int i3=0; i3<n3; ++i3) {
    for (int i2=0; i2<n2; ++i2) {
      float w1i = w1[i3][i2];
      float w2i = w2[i3][i2];
      float w3i = w3[i3][i2];
      float t1i = nearestInterp(w1i,w2i,w3i,t1);
      float t2i = nearestInterp(w1i,w2i,w3i,t2);
      float t3i = nearestInterp(w1i,w2i,w3i,t3);
      float r1i = r1[i3][i2][i1];
      float r2i = r2[i3][i2][i1];
      float r3i = r3[i3][i2][i1];
      x1[i3][i2] = i1+t1i-r1i;
      x2[i3][i2] = i2+t2i-r2i;
      x3[i3][i2] = i3+t3i-r3i;
    }}
    return new float[][][]{x1,x2,x3};
  }

  public void setHorizon(float[][][] xs, float[][][] gx) {
    int n3 = gx.length;
    int n2 = gx[0].length;
    int n1 = gx[0][0].length;
    int ns = n2*n3;
    float[] k1 = new float[ns];
    float[] k2 = new float[ns];
    float[] k3 = new float[ns];
    int ik = 0;
    for (int i3=0; i3<n3; ++i3) {
    for (int i2=0; i2<n2; ++i2) {
      k1[ik] = xs[0][i3][i2];
      k2[ik] = xs[1][i3][i2];
      k3[ik] = xs[2][i3][i2];
      ik ++;
    }}
    Sampling s2 = new Sampling(n2);
    Sampling s3 = new Sampling(n3);
    DiscreteSibsonGridder2 ds = new DiscreteSibsonGridder2(k1,k2,k3);
    float[][] x1 = ds.grid(s2,s3);
    for (int i3=0; i3<n3; ++i3) {
    for (int i2=0; i2<n2; ++i2) {
      int i1 = round(x1[i3][i2]);
      if(i1<0){i1=0;}
      if(i1>=n1){i1=n1-1;}
      gx[i3][i2][i1] = 10f;
    }}
  }


  public TriangleGroup applyForTg( float color,
    float[][][] xs, float[][][] gx, FaultSkin[] sks) {
    _color = color;
    int n3 = gx.length;
    int n2 = gx[0].length;
    int n1 = gx[0][0].length;
    short[][][] mk = zeroOnFaults(n1,n2,n3,sks);
    Sampling s2 = new Sampling(n2);
    Sampling s3 = new Sampling(n3);
    float[][] hz = horizonInterp(xs,mk);
    float[][] tg = buildTrigs(s3,s2,_color,hz,gx,mk);
    //float[][][] rgb = rgbFromHeight(hz);
    return new TriangleGroup(true,tg[0],tg[1]);
  }

  /*
  public QuadGroup applyForQg(float[][][] xs, float[][][] gx) {
    int n3 = gx.length;
    int n2 = gx[0].length;
    Sampling s2 = new Sampling(n2);
    Sampling s3 = new Sampling(n3);
    float[][] hz = horizonInterp(xs);
    float[][][] rgb = rgbFromHeight(hz);
    return new QuadGroup(false,s2,s3,hz,rgb[0],rgb[1],rgb[2]);
  }
  */

  private short[][][] zeroOnFaults(int n1, int n2, int n3, FaultSkin[] sks) {
    short[][][] mk = new short[n3][n2][n1];
    for (FaultSkin sk:sks) {
      for (FaultCell fc:sk) {
        int[] is = fc.getI();
        int i1i = is[0];
        int i2i = is[1];
        int i3i = is[2];
        mk[i3i][i2i][i1i] = 1;
      }
    }
    return mk;
  }

  private float[][] buildTrigs(Sampling sx, Sampling sy, 
    float color, float[][] z, float[][][] gx, short[][][] mk) {
    int i = 0;
    int k = 0;
    int nz = mk[0][0].length;
    int nx = sx.getCount()-1;
    int ny = sy.getCount()-1;
    float[] zas = new float[nx*ny*6];
    float[] xyz = new float[nx*ny*6*3];
    for (int ix=0;ix<nx; ++ix) {
      float x0 = (float)sx.getValue(ix  );
      float x1 = (float)sx.getValue(ix+1);
      for (int iy=0; iy<ny; ++iy) {
        float y0 = (float)sy.getValue(iy  );
        float y1 = (float)sy.getValue(iy+1);
        float z1 = z[ix  ][iy  ];
        float z2 = z[ix  ][iy+1];
        float z3 = z[ix+1][iy  ];
        float z4 = z[ix+1][iy  ];
        float z5 = z[ix  ][iy+1];
        float z6 = z[ix+1][iy+1];
        if(abs(z1-z2)>2f){continue;}
        if(abs(z1-z3)>2f){continue;}
        if(abs(z2-z3)>2f){continue;}
        if(abs(z4-z5)>2f){continue;}
        if(abs(z4-z6)>2f){continue;}
        if(abs(z5-z6)>2f){continue;}
        if(onFault(x0,y0,z1,mk)){continue;}
        if(onFault(x0,y1,z2,mk)){continue;}
        if(onFault(x1,y0,z3,mk)){continue;}
        if(onFault(x1,y0,z4,mk)){continue;}
        if(onFault(x0,y1,z5,mk)){continue;}
        if(onFault(x1,y1,z6,mk)){continue;}
        if(z1<0||z2<0||z3<0){continue;}
        if(z4<0||z5<0||z6<0){continue;}
        if(z1>nz||z2>nz||z3>nz){continue;}
        if(z4>nz||z5>nz||z6>nz){continue;}
        //float v1 = sincInterp(z1,y0,x0,gx);
        //float v2 = sincInterp(z2,y1,x0,gx);
        //float v3 = sincInterp(z3,y0,x1,gx);
        //float v4 = sincInterp(z4,y0,x1,gx);
        //float v5 = sincInterp(z5,y1,x0,gx);
        //float v6 = sincInterp(z6,y1,x1,gx);
        //zas[k++] = v1;  zas[k++] = v2;  zas[k++] =v3;
        //zas[k++] = v4;  zas[k++] = v5;  zas[k++] =v6;
        zas[k++] = z1;  zas[k++] = z2;  zas[k++] =z3;
        zas[k++] = z4;  zas[k++] = z5;  zas[k++] =z6;

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
    float zmin = -max(zas);
    float zmax = -min(zas);
    if(color>0.0f) {
      zero(zas);
      add(zas,color,zas);
      ColorMap cp = new ColorMap(0.0f,1.0f,ColorMap.JET);
      rgb = cp.getRgbFloats(zas);
    }else {
      ColorMap cp = new ColorMap(zmin,zmax,ColorMap.JET);
      rgb = cp.getRgbFloats(mul(zas,-1f));
    }
    //ColorMap cp = new ColorMap(min(zas),max(zas),ColorMap.RED_WHITE_BLUE);
    //float[] rgb = cp.getRgbFloats(zas);
    return new float[][]{xyz,rgb};
  }

  private boolean onFault(float x3, float x2, float x1, short[][][] mk) {
    int n3 = mk.length;
    int n2 = mk[0].length;
    int n1 = mk[0][0].length;
    int i1 = round(x1);
    int i2 = round(x2);
    int i3 = round(x3);
    if(i1<0){i1=0;}if(i1>=n1){i1=n1-1;}
    if(i2<0){i2=0;}if(i2>=n2){i2=n2-1;}
    if(i3<0){i3=0;}if(i3>=n3){i3=n3-1;}
    if(mk[i3][i2][i1]==1){return true;}
    else                {return false;}
  }

  private float[][][] rgbFromHeight(float[][] hz) {
    int n3 = hz.length;
    int n2 = hz[0].length;
    float[][] r = new float[n3][n2];
    float[][] g = new float[n3][n2];
    float[][] b = new float[n3][n2];
    float[] ha = new float[n2*n3];
    float hmin = -max(hz);
    float hmax = -min(hz);
    ColorMap cp = new ColorMap(hmin,hmax,ColorMap.JET);
    int ik = 0;
    for (int i3=0; i3<n3; ++i3) {
    for (int i2=0; i2<n2; ++i2) {
      ha[ik] = -hz[i3][i2];
      ik++;
    }}
    ik = 0;
    float[] htRGB = cp.getRgbFloats(ha);
    for (int i3=0; i3<n3; ++i3) {
    for (int i2=0; i2<n2; ++i2) {
      r[i3][i2] = htRGB[ik];
      g[i3][i2] = htRGB[ik+1];
      b[i3][i2] = htRGB[ik+2];
      ik += 3;
    }}
    return new float[][][]{r,g,b};
  }


  private float[][] horizonInterp(float[][][] xs, short[][][] mk) {
    int n3 = mk.length;
    int n2 = mk[0].length;
    int n1 = mk[0][0].length;
    int ns = n2*n3;
    float[] k1 = new float[ns];
    float[] k2 = new float[ns];
    float[] k3 = new float[ns];
    int ik = 0;
    for (int i3=0; i3<n3; ++i3) {
    for (int i2=0; i2<n2; ++i2) {
      k1[ik] = xs[0][i3][i2];
      k2[ik] = xs[1][i3][i2];
      k3[ik] = xs[2][i3][i2];
      ik ++;
    }}
    Sampling s2 = new Sampling(n2);
    Sampling s3 = new Sampling(n3);
    NearestGridder2 ng = new NearestGridder2(k1,k2,k3);
    float[][] hz = ng.grid(s2,s3);
    float[][] sc = fillfloat(1.0f,n2,n3);
    for (int i3=0; i3<n3; ++i3) {
    for (int i2=0; i2<n2; ++i2) {
      int i1 = round(hz[i3][i2]);
      if(i1<0){i1=0;}if(i1>=n1){i1=n1-1;}
      if(mk[i3][i2][i1]==1){sc[i3][i2]=0f;}
    }}
    LocalSmoothingFilter lsf = new LocalSmoothingFilter();
    lsf.apply(2.0f,sc,hz,hz);
    return hz;
  }

  private float sincInterp(
    float x1, float x2, float x3, float[][][] gx) 
  {
    int n3 = gx.length;
    int n2 = gx[0].length;
    int n1 = gx[0][0].length;
    Sampling s1 = new Sampling(n1);
    Sampling s2 = new Sampling(n2);
    Sampling s3 = new Sampling(n3);
    SincInterpolator si = new SincInterpolator();
    return si.interpolate(s1,s2,s3,gx,x1,x2,x3);
  }

  private float nearestInterp(
    float x1, float x2, float x3, float[][][] gx) {
    int n3 = gx.length;
    int n2 = gx[0].length;
    int n1 = gx[0][0].length;
    int i1 = round(x1); 
    int i2 = round(x2); 
    int i3 = round(x3); 
    if(i1<0){i1=0;}if(i1>=n1){i1=n1-1;}
    if(i2<0){i2=0;}if(i2>=n2){i2=n2-1;}
    if(i3<0){i3=0;}if(i3>=n3){i3=n3-1;}
    return gx[i3][i2][i1];
  }

  private float _color = -1f;
}

/****************************************************************************
Copyright (c) 2013, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is 
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/

package hv;

import edu.mines.jtk.util.*;
import edu.mines.jtk.dsp.*;
import edu.mines.jtk.awt.ColorMap;
import edu.mines.jtk.sgl.TriangleGroup;
import static edu.mines.jtk.util.ArrayMath.*;

import ipfx.*;

/**
 * From a precomputed horizon volume x1(u1,x2,x3) or 
 * a relative geologic time volume u1(x1,x2,x3), seismic horizons can be easily 
 * extracted by simply horizontal slicing in the horizon volume at each relative 
 * geologic time u1
 * @author Xinming Wu, Colorado School of Mines
 * @version 2015.05.27
 */

public class HorizonExtraction {
  /**
   * Sets horizon volume and relative geologc time volume for horizon extraction.
   * @param u1 an array of relative geologic time volume u1(x1,x2,x3).
   * @param x1 an array of horizon volume x1(u1,x2,x3).
   */
  public HorizonExtraction(
    Sampling s1, Sampling s2, Sampling s3, float[][][] u1, float[][][] x1) 
  {
    _s1 = s1;
    _s2 = s2;
    _s3 = s3;
    _u1 = u1;
    if(x1!=null){_x1=x1;}
    else {_x1=horizonVolumeFromRgt(s1,u1);}
  }

  /**
   * Extract a single horizon corresponding to relative geologic time u1i
   * defined in the sampling st
   * @param u1i relative geologic time 
   */
  public float[][] singleHorizon(float u1i) {
    int n3 = _s3.getCount();
    int n2 = _s2.getCount();
    double d1 = _s1.getDelta();
    double f1 = _s1.getFirst();
    float[][] surf = new float[n3][n2];
    int i1 = (int)round((u1i-f1)/d1);
      for (int i3=0; i3<n3; ++i3) {
      for (int i2=0; i2<n2; ++i2) {
        surf[i3][i2] = _x1[i3][i2][i1];
      }}
    return surf;
  }


  /**
   * Extract horizons corresponding relative geologic time 
   * defined in the sampling st
   * @param st sampling in relative geolgoic time for horizon extraction
   */
  public float[][][] multipleHorizons(Sampling st) {
    int n3 = _s3.getCount();
    int n2 = _s2.getCount();
    double[] u1s = st.getValues();
    double d1 = _s1.getDelta();
    double f1 = _s1.getFirst();
    int n1s = u1s.length;
    float[][][] surfs = new float[n1s][n3][n2];
    for (int i1s=0; i1s<n1s; ++i1s) {
      double u1i = u1s[i1s];
      int i1 = (int)round((u1i-f1)/d1);
      for (int i3=0; i3<n3; ++i3) {
      for (int i2=0; i2<n2; ++i2) {
        surfs[i1s][i3][i2] = _x1[i3][i2][i1];
      }}
    }
    return surfs;
  }

  /**
   * Build triangle groups for a single horizon surface 
   * corresponging to the relative geologic time u1i
   * @param u1i relative geologic time
   * @param sf array of the corresponding horizon surface 
   */

  public TriangleGroup applyForTgs(float u1i, float[][] sf) {
    float umin = min(_u1);
    float umax = max(_u1);
    ColorMap cp = new ColorMap(umin,umax,ColorMap.JET);
    float[] xyz = getVertices(sf);
    float[] vs  = fillfloat(u1i,xyz.length);
    float[] rgb = cp.getRgbFloats(vs);
    return new TriangleGroup(true,xyz,rgb);
  }

  /**
   * Build triangle groups for horizon surfaces 
   * defined in the sampling st
   * @param st sampling in relative geolgoic time for horizon extraction
   * @param sfs array of corresponding horizon surfaces 
   */

  public TriangleGroup[] applyForTgs(Sampling st, float[][][] sfs) {
    int ns = sfs.length;
    double[] u1s = st.getValues();
    TriangleGroup[] tgs = new TriangleGroup[ns];
    for (int is=0; is<ns; ++is) {
      float[][] sfi = sfs[is];
      float u1i = (float)u1s[is];
      tgs[is] = applyForTgs(u1i,sfi);
    }
    return tgs;
  }

  /**
   * Get triangle vertices on the horizon surface
   * @param z horizon surface z(x,y)
   */
  private float[] getVertices(float[][] z) {
    int i = 0;
    int nx = _s3.getCount();
    int ny = _s2.getCount();
    int nz = _s1.getCount();
    float dz = (float)_s1.getDelta();
    float[] xyz = new float[nx*ny*6*3];
    for (int ix=4;ix<nx-4; ++ix) {
      float x0 = (float)_s3.getValue(ix  );
      float x1 = (float)_s3.getValue(ix+1);
      for (int iy=4; iy<ny-4; ++iy) {
        float y0 = (float)_s2.getValue(iy  );
        float y1 = (float)_s2.getValue(iy+1);
        float z1 = z[ix  ][iy  ];
        float z2 = z[ix  ][iy+1];
        float z3 = z[ix+1][iy  ];
        float z4 = z[ix+1][iy  ];
        float z5 = z[ix  ][iy+1];
        float z6 = z[ix+1][iy+1];
        if(abs(z1-z2)/dz>2f){continue;}
        if(abs(z1-z3)/dz>2f){continue;}
        if(abs(z2-z3)/dz>2f){continue;}
        if(abs(z4-z5)/dz>2f){continue;}
        if(abs(z4-z6)/dz>2f){continue;}
        if(abs(z5-z6)/dz>2f){continue;}
        /*
        if(onFault(x0,y0,z1,mk)){continue;}
        if(onFault(x0,y1,z2,mk)){continue;}
        if(onFault(x1,y0,z3,mk)){continue;}
        if(onFault(x1,y0,z4,mk)){continue;}
        if(onFault(x0,y1,z5,mk)){continue;}
        if(onFault(x1,y1,z6,mk)){continue;}
        */
        if(z1<0||z2<0||z3<0){continue;}
        if(z4<0||z5<0||z6<0){continue;}
        if(z1>nz||z2>nz||z3>nz){continue;}
        if(z4>nz||z5>nz||z6>nz){continue;}

        xyz[i++] = x0;  xyz[i++] = y0;  xyz[i++] = z[ix  ][iy  ];
        xyz[i++] = x0;  xyz[i++] = y1;  xyz[i++] = z[ix  ][iy+1];
        xyz[i++] = x1;  xyz[i++] = y0;  xyz[i++] = z[ix+1][iy  ];
        xyz[i++] = x1;  xyz[i++] = y0;  xyz[i++] = z[ix+1][iy  ];
        xyz[i++] = x0;  xyz[i++] = y1;  xyz[i++] = z[ix  ][iy+1];
        xyz[i++] = x1;  xyz[i++] = y1;  xyz[i++] = z[ix+1][iy+1];
      }
    }
    return copy(i,0,xyz);
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
    else                 {return false;}
  }


  private short[][][] markOnFaults(
    int n1, int n2, int n3, FaultSkin[] sks) 
  {
    float d1 = (float)_s1.getDelta();
    short[][][] mk = new short[n3][n2][n1];
    for (FaultSkin sk:sks) {
      for (FaultCell fc:sk) {
        float s1 = fc.getS1();
        if(s1/d1<0.8f){continue;}
        int[] is = fc.getI();
        int i1i = is[0];
        int i2i = is[1];
        int i3i = is[2];
        mk[i3i][i2i][i1i] = 1;
      }
    }
    return mk;
  }


  private float[][][] horizonVolumeFromRgt(
    final Sampling s1, final float[][][] u1) 
  {
    cleanRGT(u1);
    final int n3 = u1.length;
    final int n2 = u1[0].length;
    final int n1 = u1[0][0].length;
    final float[][][] x1 = new float[n3][n2][n1];
    final InverseInterpolator ii = new InverseInterpolator(s1,s1);
    Parallel.loop(n3,new Parallel.LoopInt() {
    public void compute(int i3) {
      for (int i2=0; i2<n2; ++i2) 
        ii.invert(u1[i3][i2],x1[i3][i2]);
    }});
    return x1;
  }

  private void cleanRGT(float[][][] u) {
    int n3 = u.length;
    int n2 = u[0].length;
    int n1 = u[0][0].length;
    float dMin = 0.0001f*(float)_s1.getDelta();
    for (int i3=0; i3<n3; ++i3) {
      for (int i2=0; i2<n2; ++i2) {
        for (int i1=1; i1<n1; ++i1) {
          float ui = u[i3][i2][i1  ];
          float um = u[i3][i2][i1-1];
          float ut = um+dMin;
          if (ui<ut) {u[i3][i2][i1]=ut;}
        }
      }
    }
  }

  private Sampling _s1;
  private Sampling _s2;
  private Sampling _s3;
  private float[][][] _x1;
  private float[][][] _u1;
}

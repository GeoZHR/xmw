/****************************************************************************
Copyright 2007, Colorado School of Mines and others.
Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
****************************************************************************/
package sso;

import edu.mines.jtk.dsp.*;
import edu.mines.jtk.util.Parallel;
import static edu.mines.jtk.util.ArrayMath.*;

import ad.FastExplicitDiffusion;

/**
 * Waveform shape based coherence
 * @author Xinming Wu, Colorado School of Mines
 * @version 2016.08.18
 */
public class ShapeCoherence {

  public float[][] applyForCoherence(int d, float[][] p2, float[][] gx) {
    int n2 = gx.length;
    int n1 = gx[0].length;
    float[][] ch = new float[n2][n1];
    SincInterpolator si = new SincInterpolator();
    Sampling s1 = new Sampling(n1);
    for (int i2=1; i2<n2; i2++) {
      for (int i1=0; i1<n1; i1++) {
        float i1m = i1-p2[i2][i1];
        float chi = 0f;
        float gms = 0f;
        float gis = 0f;
        for (int k=-d; k<=d; k++) {
          int   k1i = i1 +k;
          float k1m = i1m+k;
          if(k1i<0||k1i>n1-1) continue;
          if(k1m<0||k1m>n1-1) continue;
          float ks = k*k*0.001f;
          float gi = gx[i2][k1i];
          //float gm = gx[i2-1][k1i];
          float gm = si.interpolate(s1,gx[i2-1],k1m);
          chi += gi*gm;
          gms += gm*gm;
          gis += gi*gi;
        }
        gms = sqrt(gms);
        gis = sqrt(gis);
        ch[i2][i1] = max(chi/(gms*gis),0);
      }
    }
    return ch;
  }

  public float[][] applyForCoherenceD(int dv, int dh, float[][] p2, float[][] gx) {
    int n2 = gx.length;
    int n1 = gx[0].length;
    float[][] ch = new float[n2][n1];
    SincInterpolator si = new SincInterpolator();
    RecursiveGaussianFilter rgf = new RecursiveGaussianFilter(dv);
    for (int i2=1; i2<n2-1; i2++) {
      float[] gms = new float[n1];
      float[] gps = new float[n1];
      float[] x1m = rampfloat(0,1,n1);
      float[] x1p = rampfloat(0,1,n1);
      for (int k=1; k<=dh; k++) {
        if(i2-k<0||i2+k>=n2) continue;
        float[] gmk = new float[n1];
        float[] gpk = new float[n1];
        float[] pmk = new float[n1];
        float[] ppk = new float[n1];
        si.interpolate(n1,1,0,p2[i2-k+1],n1,x1m,pmk);
        si.interpolate(n1,1,0,p2[i2+k-1],n1,x1p,ppk);
        sub(x1m,pmk,x1m);
        add(x1p,ppk,x1p);
        si.interpolate(n1,1,0,gx[i2-k],n1,x1m,gmk);
        si.interpolate(n1,1,0,gx[i2+k],n1,x1p,gpk);
        add(gms,gmk,gms);
        add(gps,gpk,gps);
      }
      for (int i1=0; i1<n1; i1++) {
        float chm = 0f;
        float chp = 0f;
        float gia = 0f;
        float gma = 0f;
        float gpa = 0f;
        for (int k=-dv; k<=dv; k++) {
          int   k1i = i1+k;
          if(k1i<0||k1i>n1-1) continue;
          float gi = gx[i2][k1i];
          float gm = gms[k1i];
          float gp = gps[k1i];
          chm += gi*gm;
          chp += gi*gp;
          gma += gm*gm;
          gpa += gp*gp;
          gia += gi*gi;
        }
        gma = sqrt(gma);
        gpa = sqrt(gpa);
        gia = sqrt(gia);
        chm /= (gma*gia);
        chp /= (gpa*gia);
        chm = max(chm,0);
        chp = max(chp,0);
        ch[i2][i1] = chm;//*chp;//max(chi/(gma*gpa),0);
      }

      /*
      float[] gn = abs(sub(gps,gms));
      float[] gd = abs(add(gps,gms));
      rgf.apply0(gn,gn);
      rgf.apply0(gd,gd);
      for (int i1=0; i1<n1; i1++) {
        float chi = gn[i1]/gd[i1];
        chi = min(chi,1);
        ch[i2][i1] = chi;
      }
      */
    }
    return ch;
  }


}

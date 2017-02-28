/****************************************************************************
Copyright 2012, Colorado School of Mines and others.
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
package hdw;

import java.util.Random;

import edu.mines.jtk.util.*;

import util.SmoothWithShaping;

import static edu.mines.jtk.util.ArrayMath.*;

import swt.*;

/**
 * Seismic flattening with dynamic warping.
 * @author Xinming Wu, Colorado School of Mines
 * @version 2017.02.16
 */
public class WellHelper {

  public float[][][] toArray(WellLog[] logs, double[] ndf) {
    int nc = 4;
    int nw = logs.length;
    float dz = 0f;
    float zmin =  FLT_MAX;
    float zmax = -FLT_MAX;
    for (int iw=0; iw<nw; ++iw) {
      float[] z = logs[iw].z;
      int nz = z.length;
      float zf = z[0];
      float zl = z[nz-1];
      if(zf<zmin) zmin=zf;
      if(zl>zmax) zmax=zl;
      dz = z[1]-z[0];
    }
    int nz = round((zmax-zmin)/dz)+1;
    float[][][] wd = fillfloat(_nullValue,nz,nw,nc);
    for (int iw=0; iw<nw; ++iw) {
      WellLog log = logs[iw];
      int nzi = log.z.length;
      for (int iz=0; iz<nzi; iz++) {
        float[] v = log.getCurve("vel");
        float[] d = log.getCurve("den");
        float[] g = log.getCurve("gar");
        float[] p = log.getCurve("por");
        int kz = round((log.z[iz]-zmin/dz));
        if(v!=null) wd[0][iw][kz] = v[iz];
        if(d!=null) wd[1][iw][kz] = d[iz];
        if(g!=null) wd[2][iw][kz] = g[iz];
        if(p!=null) wd[3][iw][kz] = p[iz];
      }
    }
    int izf = nz-1;
    int izl = 0;
    for (int ic=0; ic<nc; ++ic) {
    for (int iw=0; iw<nw; ++iw) {
      int wzf = 0;
      int wzl = nz-1;
      while(wzf<nz&&wd[ic][iw][wzf]==_nullValue)
        wzf++;
      while(wzl>=0&&wd[ic][iw][wzl]==_nullValue)
        wzl--;
      if(wzf<izf) izf=wzf;
      if(wzl>izl) izl=wzl;
    }}
    int mz = izl-izf+1;
    ndf[0] = mz;
    ndf[1] = dz;
    ndf[2] = izf*dz+zmin;
    return copy(mz,nw,nc,izf,0,0,wd);
  }

  public float[][][] resample(int d, float[][][] wd) {
    int nc = wd.length;
    int nw = wd[0].length;
    int nz = wd[0][0].length;
    float[][][] wc = new float[nc][nw][nz];
    for (int ic=0; ic<nc; ++ic) {
    for (int iw=0; iw<nw; ++iw) {
      int ik=0;
      for (int iz=0; iz<nz; iz+=d)
        wc[ic][iw][ik++] = wd[ic][iw][iz];
    }}
    int mz=0;
    for (int iz=0; iz<nz; iz+=d)
      mz++;
    return copy(mz,nw,nc,0,0,0,wc);
  }
  
  public float[][] trim(double[] ndf, float[][] wd) {
    int nw = wd.length;
    int nz = wd[0].length;
    int izf = nz-1;
    int izl = 0;
    for (int iw=0; iw<nw; ++iw) {
      int wzf = 0;
      int wzl = nz-1;
      while(wzf<nz&&wd[iw][wzf]==_nullValue)
        wzf++;
      while(wzl>=0&&wd[iw][wzl]==_nullValue)
        wzl--;
      if(wzf<izf) izf=wzf;
      if(wzl>izl) izl=wzl;
    }
    ndf[2] = ndf[2]+ndf[1]*izf;
    return copy(izl-izf+1,nw,izf,0,wd);

  }

  public float[][] resample(int d, float[][] wd) {
    int nw = wd.length;
    int nz = wd[0].length;
    float[][] wc = new float[nw][nz];
    for (int iw=0; iw<nw; ++iw) {
      int ik=0;
      for (int iz=0; iz<nz; iz+=d)
        wc[iw][ik++] = wd[iw][iz];
    }
    int mz=0;
    for (int iz=0; iz<nz; iz+=d)
      mz++;
    return copy(mz,nw,0,0,wc);
  }

  public float[][] sortLogs(float[][] wd) {
    int nw = wd.length;
    int nz = wd[0].length;
    int[] np = new int[nw];
    int[] ids = new int[nw];
    for (int iw=0; iw<nw; ++iw) {
      ids[iw] = iw;
    for (int iz=0; iz<nz; ++iz) {
      if(wd[iw][iz]!=_nullValue)
        np[iw]+=1;
    }}
    quickIndexSort(np,ids);
    float[][] wr = new float[nw][nz];
    for (int iw=0; iw<nw; ++iw) {
      wr[iw] = wd[ids[nw-iw-1]];
    }
    return wr;
  }


  private int _nl;
  private int _lmin;
  private int _lmax;
  private float _nullValue=-999.25f;
  private int _bstrain1;
}


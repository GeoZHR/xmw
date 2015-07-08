package hv;

import java.util.ArrayList;

import edu.mines.jtk.dsp.*;
import edu.mines.jtk.util.*;
import edu.mines.jtk.interp.*;
import static edu.mines.jtk.util.ArrayMath.*;

/**
 * Set up control points for generating a horizon volume
 * Method rearrange:
 *   Compute the average depth of all the control points of each set, and set 
 *   the point with the average depth as the reference point for that set.
 * Method extend:
 *   Optional but usually useful: extend scattered control points to a control 
 *   surfaces (each set produces one surface) by using the "HorizonExtractorC" method.
 *
 * @author Xinming Wu
 * @version 2014.03.13
 */

public class SetupConstraints {
  
  public void setForExtension (float sigma1, float sigma2, float scale) {
    _scale  = scale ;
    _sigma1 = sigma1;
    _sigma2 = sigma2;
  } 

  public void setMask(float[][][] mask) {
    _mask = mask;
  }

  public float[][] constraintsFromSurface(float[][] surf) {
    int n3 = surf.length;
    int n2 = surf[0].length;
    int np = n2*n3;
    float[] k1 = new float[np];
    float[] k2 = new float[np];
    float[] k3 = new float[np];
    int ip = 0;
    for (int i3=0; i3<n3; ++i3) {
      for (int i2=0; i2<n2; ++i2) {
        k3[ip] = i3;
        k2[ip] = i2;
        k1[ip] = surf[i3][i2];
        ip ++;
      }
    }
    return rearrange(k1,k2,k3);
  }

  public float[][] rearrange(float[] k1, float[] k2, float[] k3) {
    int ai = 0;
    int np = k1.length;
    float[][] k = zerofloat(np,4);
    float k1a = sum(k1)/(float)np;
    float min = Float.POSITIVE_INFINITY;
    for (int ip=0; ip<np; ip++) {
      int k1i = round(k1[ip]);
      k[0][ip] = (float)k1i; 
      k[1][ip] = k2[ip];
      k[2][ip] = k3[ip];
      k[3][ip] = k1[ip]-k1i;
      float dk = abs(k1[ip]-k1a);
      if (dk<min){min = dk; ai=ip;}
    }
    float t0 = k[0][0];
    float t1 = k[1][0];
    float t2 = k[2][0];
    float t3 = k[3][0];
    k[0][0] = k[0][ai];  k[0][ai] = t0;
    k[1][0] = k[1][ai];  k[1][ai] = t1;
    k[2][0] = k[2][ai];  k[2][ai] = t2;
    k[3][0] = k[3][ai];  k[3][ai] = t3;
    return k;
  }

  public float[][] firstExtend(float[] k1, float[] k2, float[] k3,
    int w2, int w3, float[][][] p2, float[][][] p3, float[][][] ep, 
    float[][][] u, float[][][] wp) {
    int np = k1.length;
    float[][][] kk = new float[np][4][];
    float[] k1i = new float[1];
    float[] k2i = new float[1];
    float[] k3i = new float[1];
    int npp = 0;
    for (int ip=0; ip<np; ++ip) {
      k1i[0] = k1[ip];
      k2i[0] = k2[ip];
      k3i[0] = k3[ip];
      kk[ip] = extend(k1i,k2i,k3i,w2,w3,p2,p3,ep,u);
      npp += kk[ip][0].length;
    }
    int k = 0;
    float[][] kkp = new float[3][npp];
    for (int ip=0; ip<np; ++ip) {
      int npi = kk[ip][0].length;
      for (int kp=0; kp<npi; ++kp) {
        kkp[0][k] = kk[ip][0][kp]+kk[ip][3][kp];
        kkp[1][k] = kk[ip][1][kp];
        kkp[2][k] = kk[ip][2][kp];
        k++;
      }
    }
    return kkp;
  }

  public float[][] extend(float[] k1, float[] k2, float[] k3,
    int w2, int w3, float[][][] p2, float[][][] p3, float[][][] ep, 
    float[][][] u) {
    int np = k1.length;
    int n3  = p2.length;
    int n2  = p2[0].length;
    int n1  = p2[0][0].length;
    int ib2 = (int)min(k2)-w2; if(ib2<0   ) {ib2=0;   }
    int ib3 = (int)min(k3)-w3; if(ib3<0   ) {ib3=0;   }
    int ie2 = (int)max(k2)+w2; if(ie2>n2-1) {ie2=n2-1;} 
    int ie3 = (int)max(k3)+w3; if(ie3>n3-1) {ie3=n3-1;} 
    int n2m = ie2-ib2+1;
    int n3m = ie3-ib3+1;
    float[][][] um  = copy(n1,n2m,n3m,0,ib2,ib3,u );
    float[][][] p2m = copy(n1,n2m,n3m,0,ib2,ib3,p2);
    float[][][] p3m = copy(n1,n2m,n3m,0,ib2,ib3,p3);
    float[][][] epm = copy(n1,n2m,n3m,0,ib2,ib3,ep);
    sub(k2,ib2,k2);
    sub(k3,ib3,k3);
    SurfaceExtractorC se = new SurfaceExtractorC();
    se.setCG(0.01f,200);
    se.setExternalIterations(20);
    se.setSmoothings(_sigma1,_sigma2);
    se.setWeights(_scale);
    float lmt = n1-1.f;
    k1 = se.refineConstraints(k1,k2,k3,u);
    float[][] surf = se.surfaceInitialization(n2m,n3m,lmt,k1,k2,k3);
    se.surfaceUpdateFromSlopes(epm,p2m,p3m,k1,k2,k3,surf);
    /*
    return surf;
    SurfaceRefine sr = new SurfaceRefine();
    sr.setExternalIterations(10,100);
    sr.setWeights(0.5,0.5);
    float[][] scale = fillfloat(1.0f,n2,n3);
    float[][] amp = se.updateWeights(surf,um);
    float sign = -1.0f*sum(amp)/abs(sum(amp));
    sr.surfaceRefine(scale,mul(um,sign),surf,0.02f);
    */
    int ci = 0;
    int ai = 0;
    float[][] k = zerofloat(n2m*n3m,4);
    float min = Float.POSITIVE_INFINITY;
    float avg = heightAvg(lmt,surf);
    for (int i3=0; i3<n3m; i3++) {
      for (int i2=0; i2<n2m; i2++) {
        int i1m = round(surf[i3][i2]);
        int i2m = i2+ib2;     
        int i3m = i3+ib3;     
        if (i1m<n1-1){
          k[0][ci] = (float)i1m; 
          k[1][ci] = (float)i2m;
          k[2][ci] = (float)i3m;
          k[3][ci] = surf[i3][i2]-(float)i1m;
          float df = abs(surf[i3][i2]-avg);
          if (df<min){min = df; ai=ci;}
          ci++;
        }
      }
    }
    float t0 = k[0][0];
    float t1 = k[1][0];
    float t2 = k[2][0];
    float t3 = k[3][0];
    k[0][0] = k[0][ai]; k[0][ai] = t0;
    k[1][0] = k[1][ai]; k[1][ai] = t1;
    k[2][0] = k[2][ai]; k[2][ai] = t2;
    k[3][0] = k[3][ai]; k[3][ai] = t3;

    return copy(ci,4,0,0,k);
  } 

  private float heightAvg(float lmt, float[][] x) {
    int n2=x.length;
    int n1=x[0].length;
    float sum = 0.0f;
    float ci = 0.0f;
    for (int i2=0; i2<n2; ++i2) {
      for (int i1=0; i1<n1; ++i1) {
        float xi = x[i2][i1];
        if(xi<lmt) {sum +=xi; ci +=1.0f;}
      }
    }
    return sum/ci;
  }
  private float _scale  = 0.0f;
  private float _sigma1 = 6.0f;
  private float _sigma2 = 6.0f;
  private float[][][] _mask=null;
}

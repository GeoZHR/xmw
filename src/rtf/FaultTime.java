package rtf;

import edu.mines.jtk.dsp.*;
import edu.mines.jtk.util.*;
import edu.mines.jtk.interp.*;
import static edu.mines.jtk.util.ArrayMath.*;

public class FaultTime {

  public EigenTensors2 makeImageTensors(float[][]fl, float[][] ft) {
    int n2 = fl.length;
    int n1 = fl[0].length;
    EigenTensors2 et = new EigenTensors2(n1,n2);
    for (int i2=0; i2<n2; ++i2) {
    for (int i1=0; i1<n1; ++i1) {
      float fli = fl[i2][i1];
      float fti = ft[i2][i1];
      fti = (float)Math.toRadians(fti);
      float u1i = 1.0f;
      float u2i = 0.0f;
      float eui = fli*0.01f;
      float evi = fli*0.01f;
      if (fli<0.30f) {
        //eui = fli;
        //evi = fli;
        et.setEigenvalues(i1,i2,eui,evi);
        et.setEigenvectorU(i1,i2,u1i,u2i);
      } else {
        u1i = -sin(fti);
        u2i =  cos(fti);
        evi = fli;
        et.setEigenvalues(i1,i2,eui,evi);
        et.setEigenvectorU(i1,i2,u1i,u2i);
      }
    }}
    return et;
  }

  public float[][] timeImage(EigenTensors2 et, float[][] fl) {
    int n2 = fl.length;
    int n1 = fl[0].length;
    float[][] fp = new float[n2][n1];
    FloatList fls = new FloatList();
    FloatList x1s = new FloatList();
    FloatList x2s = new FloatList();
    for (int i2=0; i2<n2; ++i2) {
    for (int i1=0; i1<n1; ++i1) {
      float fli = fl[i2][i1];
      if(fli>0.6f) {
        x1s.add(i1);
        x2s.add(i2);
        fls.add(fli);
        fp[i2][i1] = fli;
      }
    }}
    float[] fv = fls.trim();
    float[] x1 = x1s.trim();
    float[] x2 = x2s.trim();
    BlendedGridder2 bd2 = new BlendedGridder2(fv,x1,x2);
    bd2.setTensors(et);
    return bd2.gridNearest(0.0f,fp);
  }

  private class FloatList {
    public int n;
    public float[] a = new float[1024];
    public void add(float f) {
      if (n==a.length) {
        float[] t = new float[2*n];
        System.arraycopy(a,0,t,0,n);
        a = t;
      }
      a[n++] = f;
    }
    public float[] trim() {
      if (n==0)
        return null;
      float[] t = new float[n];
      System.arraycopy(a,0,t,0,n);
      return t;
    }
  }

}

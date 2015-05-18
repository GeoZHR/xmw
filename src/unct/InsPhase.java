package unct;

import edu.mines.jtk.dsp.*;
import static edu.mines.jtk.util.ArrayMath.*;

// compute seismic instantaneous phasae from 
// seismic amplitude using Hilbert Transform.
//
public class InsPhase {
  public void applyForPhase1(float[][] u, float[][] ph){
    int n2 = u.length;
    int n1 = u[0].length; 
    float[] ui = new float[n1];
    HilbertTransformFilter hbt = new HilbertTransformFilter();
    for (int i2=0; i2<n2; i2++){
      hbt.apply(n1,u[i2],ui);
      for (int i1=0; i1<n1; i1++){
        float uri =   u[i2][i1];
        float uii =  ui[i1];
        ph[i2][i1] = -atan2(uii,uri);
      }
    } 
  }
  public void applyForPhase2(float[][] u, float[][] ph){
    int n2 = u.length;
    int n1 = u[0].length; 
    float[][] u1i = new float[n2][n1];
    float[][] u2i = new float[n2][n1];
    HilbertTransformFilter hbt = new HilbertTransformFilter();
    for (int i2=0; i2<n2; i2++)
      hbt.apply(n1,u[i2],u1i[i2]);
    float[] ui1  = new float[n2];
    float[] ui1i = new float[n2];
    for (int i1=0; i1<n1; i1++){
      for (int i2=0; i2<n2; i2++)
        ui1[i2] = u[i2][i1];	
      hbt.apply(n2,ui1,ui1i);
      for (int i2=0; i2<n2; i2++)
        u2i[i2][i1] = ui1i[i2];
    }
    float[][] u1 = new float[n2][n1];
    float[][] u2 = new float[n2][n1];
    LocalOrientFilter lof = new LocalOrientFilter(2.0f,4.0f);
    lof.applyForNormal(u,u1,u2);
    for (int i2=0; i2<n2; i2++){
      for (int i1=0; i1<n1; i1++){
        float uri =  u[i2][i1]*u1[i2][i1] + u[i2][i1]*u2[i2][i1];
        float uii =  u1i[i2][i1]*u1[i2][i1] + u2i[i2][i1]*u2[i2][i1];
        ph[i2][i1] = -atan2(uii,uri);
      }
    } 
  }

  public void applyForPhase(float[][][] u, float[][][] ph){
    int n1 = u[0][0].length; 
    int n2 = u[0].length; 
    int n3 = u.length;
    float[] ui = new float[n1];
    HilbertTransformFilter hbt = new HilbertTransformFilter();
    for (int i3=0; i3<n3; i3++){
      for (int i2=0; i2<n2; i2++){
        hbt.apply(n1,u[i3][i2],ui);
        for (int i1=0; i1<n1; i1++){
          float uri =   u[i3][i2][i1];
          float uii =  ui[i1];
          ph[i3][i2][i1] = -atan2(uii,uri);
        }
      }
    } 
  }
  
  public void applyForCosine(float[][][] u, float[][][] cs){
    int n1 = u[0][0].length; int n2 = u[0].length; int n3 = u.length;
    float[] traceR = new float[n1]; float[] traceI = new float[n1];
    HilbertTransformFilter hbt = new HilbertTransformFilter();
    for (int i3=0; i3<n3; i3++){
      for (int i2=0; i2<n2; i2++){
        columnCopy(n1,i2,i3,u,traceR);
        hbt.apply(n1,traceR,traceI);
        for (int i1=0; i1<n1; i1++){
          cs[i3][i2][i1] = cos(atan2(traceI[i1],traceR[i1]));
        }
      }
    } 
  }

  public void applyForAmplitude(float[][][] u, float[][][] a){
    int n1 = u[0][0].length; int n2 = u[0].length; int n3 = u.length;
    float[] traceR = new float[n1]; float[] traceI = new float[n1];
    HilbertTransformFilter hbt = new HilbertTransformFilter();
    for (int i3=0; i3<n3; i3++){
      for (int i2=0; i2<n2; i2++){
        columnCopy(n1,i2,i3,u,traceR);
        hbt.apply(n1,traceR,traceI);
        for (int i1=0; i1<n1; i1++){
          a[i3][i2][i1] = sqrt(traceI[i1]*traceI[i1]+traceR[i1]*traceR[i1]);
        }
      }
    } 
  }
  
  public void applyForImg(
    float[][][] u, float[][][] uh1, float[][][] uh2, float[][][] uh3)
  {
    int n1 = u[0][0].length; int n2 = u[0].length; int n3 = u.length;
    float[] traceR2 = new float[n2];
    float[] traceR3 = new float[n3];
    float[] traceI2 = new float[n2];
    float[] traceI3 = new float[n3];
    HilbertTransformFilter hbt = new HilbertTransformFilter();
    for (int i3=0; i3<n3; i3++){
      for (int i2=0; i2<n2; i2++){
        hbt.apply(n1,u[i3][i2],uh1[i3][i2]);
      }
    }

    for (int i1=0; i1<n1; i1++){
      for (int i3=0; i3<n3; i3++){
        for (int i2=0; i2<n2; i2++){
          traceR2[i2] = u[i3][i2][i1];
        }
        hbt.apply(n2,traceR2,traceI2);
        for (int i2=0; i2<n2; i2++){
          uh2[i3][i2][i1] = traceI2[i2];
        }
      }
    }


    for (int i1=0; i1<n1; i1++){
      for (int i2=0; i2<n2; i2++){
        for (int i3=0; i3<n3; i3++){
          traceR3[i3] = u[i3][i2][i1];
        }
        hbt.apply(n3,traceR3,traceI3);
        for (int i3=0; i3<n3; i3++){
          uh3[i3][i2][i1] = traceI3[i3];
        }
      }
    }
  }
 

  private static void columnCopy(
      int n1, int i2, int i3, float[][][] u, float[] u3)
  {
    for (int i1=0; i1<n1; i1++){u3[i1] = u[i3][i2][i1];}
  } 
  
  private double _small = 0.01;
  private int _niter = 500;

}

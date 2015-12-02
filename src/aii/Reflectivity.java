package aii;

import edu.mines.jtk.util.*;


public class Reflectivity {

  public void setParams (float fmax, float spike) {
    _fmax = fmax;
    _spike = spike;
  }

  public float[][][] apply3D(
    final float dt, final float[][][] gx) 
  {
    final int n3 = gx.length;
    final int n2 = gx[0].length;
    final int n1 = gx[0][0].length;
    final float[][][] rx = new float[n3][n2][n1];
    Parallel.loop(n3,new Parallel.LoopInt() {
    public void compute(int i3) {
      for (int i2=0; i2<n2; ++i2)
        rx[i3][i2]=apply1D(dt,gx[i3][i2]);
    }});

    return rx;
  }

  public float[] apply1D(float dt, float[] gx) {
    ReflectivityEstimation re = new ReflectivityEstimation(gx,_fmax,_spike,dt);
    return re.compute();
  }

  private float _fmax = 55f;
  private float _spike = 1f;
} 

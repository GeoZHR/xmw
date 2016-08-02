package sso;


import edu.mines.jtk.dsp.*;
/**
 * Method to display 3D tensors
 * @author Xinming Wu
 * @version 2016.07.30
 */



public class TensorView {

  public EigenTensors3 tensorsFromNormal(
    float[][][] u1, float[][][] u2, float[][][] u3, 
    float[][][] w1, float[][][] w2, float[][][] w3, 
    float[][][] au, float[][][] av, float[][][] aw) 
  {
    int n3 = u1.length;
    int n2 = u1[0].length;
    int n1 = u1[0][0].length;
    EigenTensors3 ets = new EigenTensors3(n1,n2,n3,true);
    for (int i3=0; i3<n3; ++i3) {
    for (int i2=0; i2<n2; ++i2) {
    for (int i1=0; i1<n1; ++i1) {
      float aui = au[i3][i2][i1];
      float avi = av[i3][i2][i1];
      float awi = aw[i3][i2][i1];
      float u1i = u1[i3][i2][i1];
      float u2i = u2[i3][i2][i1];
      float u3i = u3[i3][i2][i1];
      float w1i = w1[i3][i2][i1];
      float w2i = w2[i3][i2][i1];
      float w3i = w3[i3][i2][i1];
      ets.setEigenvalues(i1,i2,i3,aui,avi,awi);
      ets.setEigenvectorU(i1,i2,i3,u1i,u2i,u3i);
      ets.setEigenvectorW(i1,i2,i3,w1i,w2i,w3i);
    }}}
    return ets;
  }

   
}

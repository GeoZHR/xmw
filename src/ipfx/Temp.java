package ipfx;

import java.util.*;

public class Temp {

  public void getFlFromCells(float[][][] fl, FaultCell[] cells) {
    int n3 = fl.length;
    int n2 = fl[0].length;
    int n1 = fl[0][0].length;
    float[][][] ss = new float[n3][n2][n1];
    for (FaultCell cell:cells) {
        // Get sample indices for the minus and plus sides of the cell.
        int i1 = cell.i1;
        int i2m = cell.i2m;
        int i3m = cell.i3m;
        int i2p = cell.i2p;
        int i3p = cell.i3p;
        if(i2m<0||i2p<0){continue;}
        if(i3m<0||i3p<0){continue;}
        if(i2m>=n2||i2p>=n2){continue;}
        if(i3m>=n3||i3p>=n3){continue;}
        fl[i3m][i2m][i1] += cell.fl;
        ss[i3m][i2m][i1] += 1.0f;
        fl[i3p][i2p][i1]  = cell.fl;
       ss[i3p][i2p][i1]  = 1.0f;
    }

    // Where more than one slip was accumulated, compute the average.
    for (int i3=0; i3<n3; ++i3) {
    for (int i2=0; i2<n2; ++i2) {
    for (int i1=0; i1<n1; ++i1) {
    if (ss[i3][i2][i1]>1.0f) {
      fl[i3][i2][i1] /= ss[i3][i2][i1];
    }}}}
  }

}

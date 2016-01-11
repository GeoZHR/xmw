package f3d;

import java.util.*;
import static edu.mines.jtk.util.ArrayMath.*;

import ipfx.*;

public class FaultSelect {

  public void setValueOnFaults(float v, FaultSkin[] skins, float[][][] x) {
    for (FaultSkin skin:skins) {
      for (FaultCell cell:skin) {
        int i1 = cell.getI1();
        int i2 = cell.getI2();
        int i3 = cell.getI3();
        x[i3][i2][i1] = v;
      }
    }
  }

  public void setValueOnFaults(float v, FaultCell[] cells, float[][][] x) {
    for (FaultCell cell:cells) {
      int i1 = cell.getI1();
      int i2 = cell.getI2();
      int i3 = cell.getI3();
      x[i3][i2][i1] = v;
    }
  }


  public FaultSkin[] getSubSkin(FaultSkin[] sks) {
    ArrayList<FaultSkin> sl = new ArrayList<FaultSkin>();
    for (FaultSkin sk:sks){
      int i1m = 0;
      for (FaultCell fc:sk) {
        int i1 = fc.getI1();
        if(i1>i1m) {i1m=i1;}
      }
      if(i1m>200) sl.add(sk);
    }
    return sl.toArray(new FaultSkin[0]);
  }

}

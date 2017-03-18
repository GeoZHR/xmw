package mef;

import edu.mines.jtk.dsp.*;
import edu.mines.jtk.util.*;
import static edu.mines.jtk.util.ArrayMath.*;

import java.util.*;

public class FaultCleaner {

  public FaultSkin[] filterBySize(int minSize,FaultSkin[] skins) {
    ArrayList<FaultSkin> fsl = new ArrayList<FaultSkin>();
    for (FaultSkin skin:skins) {
      if(skin.size()>minSize) {
        fsl.add(skin);
      }
    }
    return fsl.toArray(new FaultSkin[0]);
  }

  public FaultSkin[] filterBySlip(float minSlip, FaultSkin[] skins) {
    ArrayList<FaultSkin> fsl = new ArrayList<FaultSkin>();
    for (FaultSkin skin:skins) {
      float fs = 0.0f;
      float sc = skin.size();
      for (FaultCell cell:skin)
        fs += cell.getS1(); 
      if (fs/sc>minSlip)
        fsl.add(skin);
    }
    return fsl.toArray(new FaultSkin[0]);
  }
}

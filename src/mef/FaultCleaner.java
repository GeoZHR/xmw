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
}

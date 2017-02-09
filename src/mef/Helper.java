package mef;

import edu.mines.jtk.dsp.*;
import edu.mines.jtk.util.*;
import static edu.mines.jtk.util.ArrayMath.*;

public class Helper {

  public void applyTF(
    final int nf, final float fmin, final float fmax, final int[] ks,
    final float[][] fx, final float[][] pr){
    final int n2 = fx.length;
    final int n1 = fx[0].length; 
    final Sampling st = new Sampling(n1,0.004,0.0);
    final Sampling sf = MorletTransform.frequencySampling(nf,fmin,fmax);
    final MorletTransform mt = new MorletTransform(st,sf);
    Parallel.loop(n2,new Parallel.LoopInt() {
      public void compute(int i2) {
        float[][][] fi = mt.apply(fx[i2]);
        for (int ik:ks) {
          for (int i1=0; i1<n1; ++i1)
            pr[i2][i1] += fi[0][ik][i1];
        }
      }
    }); 
    float nk = ks.length;
    div(pr,nk,pr);
  }

}

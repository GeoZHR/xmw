package he;

import java.awt.*;
import java.awt.event.*;
import java.io.File;
import java.util.ArrayList;

import javax.swing.*;

import edu.mines.jtk.awt.*;
import edu.mines.jtk.dsp.*;
import edu.mines.jtk.mosaic.*;
import edu.mines.jtk.util.Check;
import edu.mines.jtk.util.Cdouble;
import static edu.mines.jtk.util.ArrayMath.*;

public class InteractiveHorizonPicker {

  /**
   * Runs the program.
   * @param args arguments (ignored).
   */
  public static void main(String[] args) {
    Check.argument(args.length>0,"type of input specified");
    String what = args[0];
    SwingUtilities.invokeLater(new Runnable() {
      public void run() {
        new InteractiveHorizonPicker();
      }
    });
  }


  private InteractiveHorizonPicker() {
  }


}

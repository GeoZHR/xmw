/****************************************************************************
Copyright (c) 2009, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is 
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package util;

import java.io.*;
import java.util.*;

public class TxtReader {

  public float[][] readLogs(String file) throws IOException {
    FileInputStream fis = new FileInputStream(file);
    InputStreamReader isr = new InputStreamReader(fis);
    BufferedReader reader = new BufferedReader(isr);
    String str=reader.readLine();
    ArrayList<Float> rha = new ArrayList<Float>();
    ArrayList<Float> vpa = new ArrayList<Float>();
    ArrayList<Float> vsa = new ArrayList<Float>();
    while(str!=null && str.length()>0) {
      str = str.replaceFirst(" ","");
      str = str.replaceAll("          "," ");
      str = str.replaceAll("         "," ");
      str = str.replaceAll("        "," ");
      String[] sts = str.split(" ");
      if(sts.length>1) {
        rha.add(new Float(sts[1]));
        vpa.add(new Float(sts[2]));
        vsa.add(new Float(sts[3]));
        System.out.println("Rho="+sts[1]+"  Vp="+sts[2]+"  Vs="+sts[3]);
      }
      str = reader.readLine();
    }
    int np = rha.size();
    float[][] las = new float[3][np];
    for (int ip=0; ip<np; ++ip) {
      las[0][ip] = rha.get(ip);
      las[1][ip] = vpa.get(ip);
      las[2][ip] = vsa.get(ip);
    }
    return las;
  }
}

/****************************************************************************
Copyright (c) 2009, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is 
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package csu;

import java.io.*;
import java.util.*;

public class TxtReader {

  public float[][] readIps(String file) throws IOException {
    FileInputStream fis = new FileInputStream(file);
    InputStreamReader isr = new InputStreamReader(fis);
    BufferedReader reader = new BufferedReader(isr);
    String str=reader.readLine();
    ArrayList<Float> ip1 = new ArrayList<Float>();
    ArrayList<Float> ip2 = new ArrayList<Float>();
    ArrayList<Float> ip3 = new ArrayList<Float>();
    ArrayList<Float> ip4 = new ArrayList<Float>();
    int k=0;
    while(str!=null && str.length()>0) {
      str = str.replaceAll("        "," ");
      str = str.replaceAll("         "," ");
      str = str.replaceAll("  "," ");
      str = str.replaceAll("	"," ");
      str = str.replaceAll("  "," ");
      str = str.replaceAll("  "," ");
      str = str.replaceAll("  "," ");
      String[] sts = str.split(" ");
      if(sts.length>1&&k>0) {
        ip1.add(new Float(sts[0]));
        ip2.add(new Float(sts[1]));
        ip3.add(new Float(sts[2]));
        ip4.add(new Float(sts[3]));
      }
      str = reader.readLine();
      k++;
    }
    int np = ip1.size();
    float[][] ips = new float[4][np];
    for (int ip=0; ip<np; ++ip) {
      ips[0][ip] = ip1.get(ip);
      ips[1][ip] = ip2.get(ip);
      ips[2][ip] = ip3.get(ip);
      ips[3][ip] = ip4.get(ip);
    }
    return ips;
  }
}

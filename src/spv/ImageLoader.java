package spv;

import java.io.File;
import java.awt.image.*;
import java.io.IOException;
import javax.imageio.ImageIO;

/**
 * Description:
 * Read and write image.
 * @author Xinming Wu
 * @version 2017.07.15
 */

public class ImageLoader{

  public float[][][] readThreeChannels(String sid) throws IOException {
    File fid = new File(sid);
    BufferedImage image = null;
    image = ImageIO.read(fid);
    Object dataElements = null;
    Raster raster = image.getRaster();
    ColorModel colorModel = image.getColorModel();
    int n1 = image.getHeight();
    int n2 = image.getWidth();
    float[][][] fx = new float[3][n2][n1];
    for (int i2=0; i2<n2; ++i2) {
    for (int i1=0; i1<n1; ++i1) {
      dataElements = raster.getDataElements(i2,i1,dataElements);
      fx[0][i2][i1] = colorModel.getRed(dataElements);
      fx[1][i2][i1] = colorModel.getBlue(dataElements);
      fx[2][i2][i1] = colorModel.getGreen(dataElements);
    }}
    return fx;
  }

  public float[][] readColorImage(String sid) throws IOException {
    File fid = new File(sid);
    BufferedImage image = null;
    image = ImageIO.read(fid);
    int n1 = image.getHeight();
    int n2 = image.getWidth();
    float[][] fx = new float[n2][n1];
    for (int i2=0; i2<n2; ++i2) {
    for (int i1=0; i1<n1; ++i1) {
      fx[i2][i1] = image.getRGB(i2,i1);
    }}
    return fx;
  }

}//class ends here

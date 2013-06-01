package gvorn;

import java.awt.*;
import testbed.*;
import base.*;

public class VornPointSite extends VornSite {

  public VornPointSite(FPoint2 loc) {
    this.loc = new FPoint2(loc);
  }
  private FPoint2 loc;

  public FPoint2 location() {
    return loc;
  }

  @Override
  public DArray buildBisector(VornSite otherSite) {
    return DArray.build(new VornPointBisector(this, (VornPointSite) otherSite));
  }

  public void render(Color c, int stroke, int markType) {
    V.pushColor(c, MyColor.cBLUE);
    V.mark(location(), markType);
    V.pop();
  }

  @Override
  public double distanceFrom(FPoint2 pt) {
    return location().distance(pt);
  }

}

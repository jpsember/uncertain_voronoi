package gvorn;

import java.awt.*;
import testbed.*;
import base.*;

public class VornSegSite extends VornSite {

  public VornSegSite(FPoint2 p0, FPoint2 p1, String label) {
    this.s0 = p0;
    this.s1 = p1;
    this.label = label;
  }

  @Override
  public DArray buildBisector(VornSite otherSite) {
    
    VornBisector bs = new VornSegBisector(this, (VornSegSite)otherSite);
    return DArray.build(bs);
//    return VornSegBisector.build(this, (VornSegSite) otherSite);
  }

  @Override
  public double distanceFrom(FPoint2 pt) {
    return MyMath.ptDistanceToSegment(pt, s0, s1, null);
  }

  public String toString() {
    return label;
  }
  public void render(Color c, int stroke, int markType) {
    V.pushColor(c, MyColor.cPURPLE);
    V.pushStroke(stroke, Globals.STRK_NORMAL);
    V.drawLine(s0, s1);
    V.pop(2);
  }
  public FPoint2 pt(int index) {
    return index == 0 ? s0 : s1;
  }

  private FPoint2 s0;
  private FPoint2 s1;
  private String label;

}

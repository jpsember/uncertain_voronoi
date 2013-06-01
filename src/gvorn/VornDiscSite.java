package gvorn;

import java.awt.*;
import testbed.*;
import base.*;

public class VornDiscSite extends VornSite {

  public VornDiscSite(EdDisc d) {
    this(d.getOrigin(), d.getRadius(), d.getLabel());
  }
  
  public VornDiscSite(FPoint2 origin, double radius, String label) {
    if (radius < 0)
      throw new IllegalArgumentException();
    this.radius = radius;
    this.origin = new FPoint2(origin);
    this.label = label;
  }
  private FPoint2 origin;
  private double radius;

  public FPoint2 origin() {
    return origin;
  }

  public double radius() {
    return radius;
  }

  //  @Override
  //  public VornBisector buildBisector(VornSite otherSite) {
  //    return new VornPointBisector(this, (VornPointSite) otherSite);
  //  }

  public void render(Color c, int stroke, int markType) {
    V.pushColor(c, MyColor.cBLUE);
    V.drawCircle(origin(), radius());
    V.pop();
  }

  @Override
  public double distanceFrom(FPoint2 pt) {
    return pt.distance(origin) - radius;
  }

  @Override
  public DArray buildBisector(VornSite otherSite) {
    return DArray.build(new VornDiscBisector(this, (VornDiscSite) otherSite));
  }

  public String toString() {
    return label;
  //  return "Disc(" + origin + " r:" + Tools.f(radius) + ")";
  }
  private String label;
}

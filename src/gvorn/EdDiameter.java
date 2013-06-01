package gvorn;

import base.*;
import testbed.*;
import java.awt.*;

public class EdDiameter extends EdObject implements Globals, Renderable {

  static final boolean db = false;

  private EdDiameter() {
  }

  public void setPoint(int ptIndex, FPoint2 pt, boolean useGrid, TBAction action) {

    switch (ptIndex) {
    case 0:
      super.setPoint(ptIndex, pt, useGrid, action);
      calcPoint2(radius);
      break;
    default:
      {
        FPoint2 or = origin();
        double dist = or.distance(pt);

        theta = 0;
        if (dist > 0) {
          theta = MyMath.polarAngle(or, pt);
          if (ptIndex == 2)
            theta += Math.PI;
        }
        if (ptIndex == 2)
          dist = radius;
        calcPoint2(dist);
      }
      break;
    }
  }

  private static final double MINRADIUS = 2.2;

  private void calcPoint2(double radius) {
    radius = Math.max(radius, MINRADIUS);
    FPoint2 or = origin();
    super.setPoint(1, MyMath.ptOnCircle(or, theta, radius), true, null);
    super.setPoint(2, MyMath.ptOnCircle(or, theta + Math.PI, radius), true,
        null);
    this.radius = radius;
  }

  public boolean complete() {
    return nVert() >= 3;
  }

  public FPoint2 origin() {
    return getPoint(0);
  }

  private double theta;

  private double radius;

  public double distFrom(FPoint2 pt) {
    FPoint2 p1 = getPoint(1);
    FPoint2 p2 = getPoint(2);
    return MyMath.ptDistanceToSegment(pt, p1, p2, null);
  }

  public EdObjectFactory getFactory() {
    return FACTORY;
  }

  /**
   * Get bounding rectangle of object
   * @return FRect
   */
  public FRect getBounds() {
    FRect r = null;
    do {
      for (int i = 0; i < nVert(); i++)
        r = FRect.add(r, getPoint(i));
    } while (false);
    return r;
  }

  /** 
   * @return
   */
  public int nVert() {
    return nPoints();
  }

  /**
   * Move entire object by a displacement
   * Default implementation just adjusts each point.
   * @param orig : a copy of the original object
   * @param delta : amount to move by
   */
  public void moveBy(EdObject orig, FPoint2 delta) {
    setPoint(0, FPoint2.add(orig.getPoint(0), delta, null));
  }

  public void render(Color color, int stroke, int markType) {
    if (color == null)
      color = isActive() ? MyColor.cPURPLE : Color.DARK_GRAY;

    do {
      if (!complete())
        break;
      FPoint2 p0 = origin(), p1 = getPoint(1), p2 = getPoint(2); //ptFwd(), p2 = ptBwd();

      {
        V.pushScale(.8);
        V.pushStroke(STRK_RUBBERBAND);
        V.pushColor(Color.black);
        V.drawCircle(p0, radius - 2.0);
        V.popColor();
        V.popStroke();
        V.popScale();
      }
      if (color != null)
        V.pushColor(color);
      if (stroke >= 0)
        V.pushStroke(stroke);

      EdSegment.plotDirectedLine(p2, p1, true, true);
      if (markType >= 0) {
        V.mark(p0, markType);
        V.mark(p1, markType);
        V.mark(p2, markType);
      }
      if (stroke >= 0)
        V.popStroke();
      if (color != null)
        V.popColor();
    } while (false);
  }

  public static EdObjectFactory FACTORY = new EdObjectFactory() {

    public EdObject construct() {
      return new EdDiameter();
    }

    public String getTag() {
      return "diam";
    }

    public EdObject parse(Tokenizer s, int flags) {
      final boolean db = false;
      if (db)
        Streams.out.println("EdDiameter, parse, next=" + s.peek().debug());

      EdDiameter obj = new EdDiameter();
      obj.setFlags(flags);
      obj.theta = s.extractDouble();
      obj.radius = s.extractDouble();
      obj.addPoint(s.extractFPoint2());
      return obj;
    }

//    public void write(StringBuilder sb, EdObject obj) {
//      EdDiameter seg = (EdDiameter) obj;
//      sb.append(seg.toString());
//    }
    public void write(StringBuilder sb, EdObject obj) {
      EdDiameter d = (EdDiameter) obj;
      if (!d.isActive())
        toString(sb, true);
      toString(sb, d.theta);
      toString(sb, d.radius);
      toString(sb, d.getPoint(0));
    }

    public String getMenuLabel() {
      return "Add diameter";
    }
    public String getKeyEquivalent() {return "m";}

  };

}

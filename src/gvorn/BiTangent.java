package gvorn;

import java.awt.*;
import testbed.*;
import base.*;

public class BiTangent implements Renderable, Globals { //extends DirectedLineSegment {
  //implements NewTraceable, Globals {

  private static final boolean TEST = false;
  static {
    if (TEST)
      Tools.warn("TEST mode is on");
  }

  public BiTangent(EdDisc a, EdDisc b) {

    final boolean db = false;

    if (db)
      Streams.out.println("BiTangent for " + a + ", " + b);

    this.discA = a;
    this.discB = b;

    if (TEST || EdDisc.partiallyDisjoint(a, b)) {
      double th = calcTheta(a, b);
      if (th != BADTHETA) {
        LineEqn eqn = new LineEqn(a.polarPoint(th + Math.PI / 2), th);
        double ta = eqn.parameterFor(a.getOrigin());
        double tb = eqn.parameterFor(b.getOrigin());
        seg = new DirSeg(eqn.pt(ta), eqn.pt(tb));
        //        setLineEqn( new LineEqn(a.polarPoint(th + Math.PI / 2), th));
        //        setPoints(lineEqn.calcClosestPointTo(a.getOrigin()),
        //          lineEqn.calcClosestPointTo(b.getOrigin()));
      }
    }
  }

  private DirSeg seg;

  public boolean defined() {
    return seg != null;
  }

  //  private double ta, tb;
  public FPoint2 tangentPt(int index) {
    return seg.endPoint(index); // == 0 ? ta : tb);
  }
  public String getLabel(int pt) {
    return disc(pt).getLabel();
  }

  public LineEqn lineEqn() {
    return seg.lineEqn();
  }

  //  private LineEqn lineEqn;
  public void render(Color c, int stroke, int markType) {
    if (c == null)
      c = Color.RED;
    V.pushColor(c);

    if (defined()) {
      seg.render(c, stroke, markType);
      //      plotDirectedHalfPlane(vp, lineEqn.pt(ta),lineEqn.pt(tb), markType);
    } else {
      if (false) {

        if (stroke < 0)
          stroke = STRK_THICK;
        discA.render(c, stroke, -1);
        discB.render(c, stroke, -1);
      }

      for (int i = 0; i < 2; i++) {
        if (i == 1 && discA == discB)
          continue;
        EdDisc d = (i == 0) ? discA : discB;
        if (stroke < 0)
          stroke = STRK_RUBBERBAND;
        V.pushStroke(stroke);
        double r = Math.max(d.getRadius() - 4, 2.0);
        V.drawCircle(d.getOrigin(), r);
        V.popStroke();
      }
    }
    V.popColor();

  }
  /**
   * Determine if a disc lies strictly to the right of this bitangent,
   * thus invalidating it as a hull edge
   * @param di
   * @return
   */
  public boolean invalidatedBy(EdDisc di) {
    return seg.lineEqn().signedDistanceFrom(di.getOrigin()) + di.getRadius() < 0;
  }
/**
 * @deprecated see EdSegment
 * @param p0
 * @param p1
 */
  public static void plotDirectedLine(FPoint2 p0, FPoint2 p1) {
    plotDirectedLine(p0, p1, false, true);
  }
  
  /**
   * @deprecated see EdSegment
   * @param p0
   * @param p1
   * @param p0Head
   * @param p1Head
   */
  public static void plotDirectedLine(FPoint2 p0, FPoint2 p1, boolean p0Head,
      boolean p1Head) {
    V.drawLine(p0, p1);
    double len = p0.distance(p1);
    // draw arrowheads
    if (len > 0) {
      final double AH_LEN = 1.2;
      final double AH_ANG = Math.PI * .85;
      double theta = MyMath.polarAngle(p0, p1);

      for (int h = 0; h < 2; h++) {
        FPoint2 ep = h == 0 ? p0 : p1;
        if ((h == 0 ? p0Head : p1Head)) {
          double th = theta; //(h == 0) ? theta + Math.PI : theta;

          FPoint2 a0 = MyMath.ptOnCircle(ep, th + AH_ANG, AH_LEN);
          V.drawLine(ep, a0);
          FPoint2 a1 = MyMath.ptOnCircle(ep, th - AH_ANG, AH_LEN);
          V.drawLine(ep, a1);
        }
      }
    }
  }

  public double theta() {
    if (!defined())
      return 0;
    return seg.lineEqn().polarAngle();
  }

  public String toString() {
    StringBuilder sb = new StringBuilder();
    sb.append("b<");
    sb.append(getLabel(0));
    sb.append(getLabel(1));
    sb.append('>');
    return sb.toString();
  }

  private static final double BADTHETA = -999;
  private static final double TOLERANCE = 1e-12;

  /**
   * Calculate angle of bitangent to two discs.
   * Assumes discs are partially disjoint.
   * Uses Newton's method.
   * @param a, b : discs
   * @return polar angle of bitangent, or BADTHETA if problem
   */
  public static double calcTheta(EdDisc a, EdDisc b) {
    final boolean db = TEST;

    double ret = BADTHETA;

    double deltaR = b.getRadius() - a.getRadius();
    double fa = b.getOrigin().x - a.getOrigin().x;
    double fb = b.getOrigin().y - a.getOrigin().y;
    if (db)
      Streams.out.println("calcTheta\n a=" + a + "\n b=" + b + "\n fa="
          + Tools.f(fa) + " fb=" + Tools.f(fb));

    double p0 = 0;

    // use angle of segment connecting origins as initial approximation
    p0 = MyMath.polarAngle(a.getOrigin(), b.getOrigin());

    if (db)
      Streams.out.println("initial approximation is " + Tools.fa(p0));

    for (int step = 0; step < 15; step++) {
      double c = Math.cos(p0), s = Math.sin(p0);
      double fval = deltaR + fb * c - fa * s;
      double fprime = -(fb * s + fa * c);
      double p = p0;
      if (Math.abs(fprime) >= TOLERANCE) {
        p = p0 - fval / fprime;
      }

      if (db)
        Streams.out.println(" " + Tools.f(step) + "  f=" + Tools.f(fval, 3, 12)
            + " f'=" + Tools.f(fprime, 3, 12) + " p0(deg)=" + Tools.fa(p0)
            + " chg=" + Tools.f(Math.abs(p - p0), 3, 12));

      if (Math.abs(p - p0) < TOLERANCE) {
        if (Math.abs(fval) < 1e-5)
          ret = p;
        break;
      }
      p0 = p;
    }

    if (ret == BADTHETA)
      Tools.warn("problem calculating theta for:\n" + a + "\n" + b);

    return ret;
  }

  /**
   * Test program for bitangent calculations.
   * Set TEST = true to see the output!
   * @param args
   */
  public static void main(String[] args) {

    if (true) {
      double a = 0;
      int n = 1;
      while (true) {
        a += 1.0 / n - 1.0 / (n + 2);
        Streams.out.println("a=" + a + " a*4=" + a * 4 + " " + Math.PI
            / (a * 4));
        n += 4;
      }
    }

    double[] testPairs = { //
    0, 0, 0, 8.66, 5, 0,//
        47.7, 51.45, 5.26, 78.74, 69.91, 4.03,//
        0, 0, 2, 5, 0, 2, //
        0, 0, 5, 3, 0, 7.9, //
        0, 0, 5, 3, 0, 7.9999, //
        0, 0, 3, 2, 1, .5, //
    };
    for (int i = 0; i < testPairs.length; i += 6) {
      EdDisc a = new EdDisc(testPairs, i + 0);
      EdDisc b = new EdDisc(testPairs, i + 3);
      new BiTangent(a, b);
    }
  }

  public static void plotDirectedHalfPlane(FPoint2 p0, FPoint2 p1, int markType) {
    double SEP = .4 * V.getScale();
    double ang = MyMath.polarAngle(p0, p1);
    FPoint2 d0 = MyMath.ptOnCircle(p0, ang + Math.PI / 2, SEP);
    FPoint2 d1 = MyMath.ptOnCircle(p1, ang + Math.PI / 2, SEP);

    plotDirectedLine(p0, p1);
    V.pushStroke(STRK_RUBBERBAND);
    V.drawLine(d0, d1);
    V.popStroke();

    if (markType >= 0) {
      V.mark(d0, markType);
      V.mark(d1, markType);
    }
  }

  private EdDisc discA, discB;

  public EdDisc disc(int i) {
    return (i == 0) ? discA : discB;
  }
}

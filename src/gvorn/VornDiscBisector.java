package gvorn;

import java.awt.*;
import testbed.*;
import base.*;

public class VornDiscBisector extends VornBisector {
  public  IPlaneCurve curve(int i){return curve;}

  public VornDiscBisector(VornDiscSite sa, VornDiscSite sb) {
    super(sa, sb);

    final boolean db = true;

    FPoint2 f1 = sa.origin();
    FPoint2 f2 = sb.origin();
    double r1 = sa.radius();
    double r2 = sb.radius();

    if (EdDisc.contains(f1, r1, f2, r2) || EdDisc.contains(f2, r2, f1, r1)) {
      throw new IllegalArgumentException("discs are nested");
    }

    double dist = FPoint2.distance(f1, f2);
    double t = .5 * (r1 + dist - r2);

    if (db && T.update())
      T.msg("r1=" + sa.radius() + " r2=" + sb.radius() + " d=" + dist + " t="
          + Tools.f(t));

    if (Math.abs(t - dist * .5) < 1e-3) {
      curve = new LineCurve(FPoint2.midPoint(f1, f2), MyMath.polarAngle(f1, f2)
          + Math.PI / 2);
    } else

      curve = new Hyperb(f1, f2, t);

    //    FPoint2 pt = FPoint2.interpolate(f1, f2, t);
    //
    //    // if point on arm is closer to f2 than f1, swap f1 & f2.
    //    boolean flipped = (t > .5);
    //
    ////    FPoint2 focus1 = new FPoint2(f1);
    ////    FPoint2 focus2 = new FPoint2(f2);
    ////    
    ////    this.foci[flipped ? LEFT : RIGHT] = new FPoint2(f1);
    ////    this.foci[flipped ? RIGHT : LEFT] = new FPoint2(f2);
    //
    //    double fociDist = FPoint2.distance(f1, f2);
    //    if (fociDist == 0) {
    //      throw new FPError("Hyperbola foci are same point");
    //    }
    //
    //    c = fociDist * .5;
    //
    //    // calculate the translation of the hyperbola away from
    //    // standard position.
    //
    //    FPoint2 rFocus = flipped ? f1 : f2,
    //        lFocus = flipped ? f2 : f1;
    //    //foci[0], lFocus = foci[1];
    //
    //    FPoint2 origin = new FPoint2(.5 * (rFocus.x + lFocus.x),
    //        .5 * (rFocus.y + lFocus.y));
    //
    //    // calculate the angle of rotation of the hyperbola away
    //    // from the standard position.
    //
    //    double theta = Math.atan2(rFocus.y - lFocus.y, rFocus.x - lFocus.x);
    //
    //    Matrix fromCenterInW = Matrix.getTranslate(origin, true);
    //    Matrix rotToE = Matrix.getRotate(-theta);
    //
    //    toE2 = rotToE;
    //    Matrix.mult(toE2, fromCenterInW, toE2);
    //
    //    // calculate inverse
    //
    //    toW2 = toE2.invert(null);
    //
    //    {
    //      // get the arm point in hyperbola space.
    //      FPoint2 workPt = toE2.apply(pt, null);
    //
    //      double xs = workPt.x * workPt.x;
    //      double cs = c * c;
    //
    //      Polyn q = new Polyn(1, -(cs + xs + workPt.y * workPt.y), cs * xs);
    //      DArray qsoln = new DArray();
    //      q.solve(qsoln);
    //      double val = q.c(1) * -.5;
    //      int ql = qsoln.size();
    //      if (ql >= 1) {
    //        val = qsoln.getDouble(0);
    //      }
    //
    //      // choose the root that is less than c*c.
    //
    //      if (ql == 2) {
    //        if (val > qsoln.getDouble(1)) {
    //          val = qsoln.getDouble(1);
    //        }
    //      }
    //
    //      a = Polyn.sqrt(val);
    //      A = a * a;
    //      B = A / (c * c - A);
    //    }
    //
    //    if (false) {
    //      Tools.warn("Always getting");
    //      getPlaneCurve();
    //    }
  }

  //  @Override
  //  public DArray intersectWith(VornBisector b) {
  //    VornDiscBisector b2 = (VornDiscBisector) b;
  //
  //    VornDiscBisector b1 = this;
  //
  //    DArray ip = findIntersections(b1, b2);
  //
  //    DArray ret = new DArray();
  //
  //    for (int i = 0; i < ip.size(); i++) {
  //      FPoint2 pt = ip.getFPoint2(i);
  //
  //      FPoint2 data = b1.calcParameterAndDistance(pt);
  //      ret.addDouble(data.x);
  //
  //      data = b2.calcParameterAndDistance(pt);
  //      ret.addDouble(data.x);
  //    }
  //    return ret;
  //  }

  //  /**
  //   * Calculate parameter for a point, and indicate how close the
  //   * point is to the curve
  //   * @param pt : point to find parameter for
  //   * @return FPoint2; x contains the parameter, y is its distance
  //   *  from the curve
  //   */
  //  private FPoint2 calcParameterAndDistance(FPoint2 pt) {
  //    final FPoint2 ept = new FPoint2();
  //    toCurveSpace(pt, ept);
  //    double t = ept.y;
  //    FPoint2 as = new FPoint2(Polyn.sqrt(A + t * t * B), t);
  //
  //    return new FPoint2(t, Math.abs(ept.x - as.x));
  //  }

  //  /**
  //   * Find intersections between two hyperbolas 
  //   * @param a Hyperbola
  //   * @param b Hyperbola
  //   * @param iPts where to store intersection points; null to construct
  //   */
  //  private static DArray findIntersections(VornDiscBisector a, VornDiscBisector b) {
  //
  //    DArray iPts = new DArray();
  //    DArray jPts = new DArray();
  //
  //    PlaneCurve.findIntersect(a.getCurve(), b.getCurve(), jPts);
  //
  //    // filter out false intersections
  //    for (int i = 0; i < jPts.size(); i++) {
  //      FPoint2 pt = jPts.getFPoint2(i);
  //      if (a.onBisector(pt) && b.onBisector(pt))
  //        iPts.add(pt);
  //    }
  //    return iPts;
  //  }

  @Override
  public FPoint2 pointAt(double t) {
    return curve.pt(t); // pointAt(t, null);
  }

  @Override
  public void render(Color color, int stroke, int markType) {
    // final boolean db = false;

    V.pushStroke(stroke, Globals.STRK_NORMAL);
    V.pushColor(color, MyColor.cBLUE);

    //    double step = Math.max(1, ((c - a) * (c - a) / c) * .3);
    //
    //    if (db)
    //      Streams.out.println(" step=" + step);
    //
    // plot each visible segment

    for (int seg = 0; seg < nComponents(); seg++) {
      VUtil2.render(curve,cStart(seg),cEnd(seg));
      
      //curve.render(cStart(seg), cEnd(seg));
      //      double t0 = this.cStart(seg), t1 = this.cEnd(seg);
      //      t0 = MyMath.clamp(t0, -500.0, 500.0);
      //      t1 = MyMath.clamp(t1, -500.0, 500.0);
      //
      //      double s0, s1;
      //
      //      s0 = t0;
      //      s1 = t1;
      //
      //      if (s0 > s1) {
      //        double tmp = s0;
      //        s0 = s1;
      //        s1 = tmp;
      //      }
      //      FPoint2 p0 = this.pointAt(s0), p1 = pointAt(s1);
      //      if (db)
      //        Streams.out.println(" p0=" + p0 + ", p1=" + p1);
      //      {
      //        int count = 0;
      //        boolean first = true;
      //        for (double t = t0;; t += step, count++) {
      //          boolean last = (t >= t1);
      //          if (last)
      //            t = t1;
      //          pointAt(t, p1);
      //
      //          if (!p1.isValid()) {
      //            if (last) {
      //              break;
      //            }
      //            continue;
      //          }
      //
      //          if (db) {
      //            System.out.println(" calcPt " + Tools.f(t) + " = " + p1.x + ","
      //                + p1.y);
      //          }
      //          if (!first)
      //            V.drawLine(p0, p1);
      //          if (last)
      //            break;
      //          p0.setLocation(p1);
      //          first = false;
      //        }
    }
    V.pop(2);

  }
  // private static final int RIGHT = 0, LEFT = 1;

  /**
   * Transform a point from world space to curve space
   * @param pt : point in world space
   * @param out : where to store curve space point
   */
  public double parameterFor(FPoint2 pt) {
    return curve.parameterFor(pt);
    //    FPoint2 p2 = toE2.apply(pt);
    //    return p2.y;
  }

  //  /**
  //   * Calculate a point on the hyperbola
  //   * @param t : parameter in internal space (after flipping has occurred)
  //   * @param dest : where to store the calculated point
  //   */
  //  private FPoint2 pointAt(double t, FPoint2 dest) {
  //    return curve.pt(t);
  //    dest = toW2.apply(Polyn.sqrt(A + t * t * B), t, dest);
  //    return dest;
  //  }

  //  /**
  //   * Get the IPlaneCurve associated with this hyperbola
  //   * @return IPlaneCurve
  //   */
  //  private IPlaneCurve getCurve() {
  //    if (hPoly == null) {
  //      {
  //        // We note that x', y' are 'curve space' coordinates, which
  //        // are derived by transforming x,y by the toE matrix;
  //        // that is,
  //        //
  //        //  x' = x * A + y * B + C
  //        //  y' = x * D + y * E + F
  //        //
  //        //
  //        double A, B, C, D, E, F, G, H;
  //
  //        A = toE2.get(0, 0);
  //        B = toE2.get(0, 1);
  //        C = toE2.get(0, 2);
  //        D = toE2.get(1, 0);
  //        E = toE2.get(1, 1);
  //        F = toE2.get(1, 2);
  //        G = a * a - c * c;
  //        H = this.A;
  //        if (Double.isNaN(G) || Double.isNaN(H)) {
  //          System.out.println("not a number: " + G + " or H=" + H);
  //        }
  //
  //        hPoly = Conic.construct(A * A * G + D * D * H, 2 * (A * B * G + D * E
  //            * H), B * B * G + E * E * H, 2 * (A * C * G + D * F * H), 2 * (B
  //            * C * G + E * F * H), C * C * G + H * (F * F - G));
  //
  //      }
  //
  //    }
  //    return hPoly;
  //  }
  //  private static void test(PlaneCurve p) {
  //    System.out.println(p);
  //
  //    double A = p.c(5);
  //    double B = p.c(4);
  //    double C = p.c(3), D = p.c(2), E = p.c(1), F = p.c(0);
  //
  //    System.out.println("A=" + A + "\nB=" + B + "\nC=" + C + "\nD=" + D + "\nE="
  //        + E + "\nF=" + F);
  //
  //    if (false) {
  //      for (double y = 86.6 - 10; y < 86.6 + 10; y += 1) {
  //        double x = 52.6;
  //        double ev = p.eval(x, y);
  //        System.out.println("eval " + new FPoint2(x, y) + " =" + ev);
  //      }
  //    }
  //
  //    double a = C, b = B, c = A;
  //
  //    double disc = b * b - 4 * a * c;
  //
  //    double q = Math.sqrt(disc);
  //
  //    double j1 = (-b + q) / (2 * a);
  //    double j2 = (-b - q) / (2 * a);
  //
  //    System.out.println("a=" + a + " b=" + b + " c=" + c + " disc=" + disc
  //        + " q=" + q + "\nj1=" + j1 + "\nj2=" + j2);
  //
  //    System.out.println("D-BE/2C=" + (D - B * E / (2 * C)));
  //    a = C;
  //    b = E;
  //    c = F;
  //
  //    q = Math.sqrt(b * b - 4 * a * c);
  //    double k1 = (-b + q) / (2 * a);
  //    double k2 = (-b - q) / (2 * a);
  //
  //    for (int i = 0; i < 4; i++) {
  //      double j = (i % 2 == 0 ? j1 : j2);
  //      double k = (i / 2 == 0 ? k1 : k2);
  //      double sum = D + 2 * C * j * k + B * k + E * j;
  //      System.out.println("j=" + j + " k=" + k + " sum=" + sum);
  //    }
  //  }
  //  private double a, c;
  //
  //  // position of foci in world space
  //  //private FPoint2[] foci = new FPoint2[2];
  //
  //  // a*a / (c*c - a*a)
  //  private double B;
  //  private double A;
  //
  //  private Matrix toE2, toW2;

  private IPlaneCurve curve;

//  @Override
//  public IPlaneCurve getPlaneCurve() {
//    return curve;
//    //    if (hPoly == null) 
//    //    {
//    //      {
//    //        // We note that x', y' are 'curve space' coordinates, which
//    //        // are derived by transforming x,y by the toE matrix;
//    //        // that is,
//    //        //
//    //        //  x' = x * A + y * B + C
//    //        //  y' = x * D + y * E + F
//    //        //
//    //        //
//    //        double A, B, C, D, E, F, G, H;
//    //
//    //        A = toE2.get(0, 0);
//    //        B = toE2.get(0, 1);
//    //        C = toE2.get(0, 2);
//    //        D = toE2.get(1, 0);
//    //        E = toE2.get(1, 1);
//    //        F = toE2.get(1, 2);
//    //        G = a * a - c * c;
//    //        H = this.A;
//    //        if (Double.isNaN(G) || Double.isNaN(H)) {
//    //          System.out.println("not a number: " + G + " or H=" + H);
//    //        }
//    //
//    //        hPoly = Conic.construct(A * A * G + D * D * H, 2 * (A * B * G + D * E
//    //            * H), B * B * G + E * E * H, 2 * (A * C * G + D * F * H), 2 * (B
//    //            * C * G + E * F * H), C * C * G + H * (F * F - G));
//    //
//    //      }
//    //
//    //    }
//    //    return hPoly;
//  }
}

package gvorn;

import java.awt.*;
import base.*;
import testbed.*;

/**
 */
public class Hyperbola implements GeomPrimitive, Globals, Renderable {

  public static final int RIGHT = 0, LEFT = 1, MINCLIPPED = 2, MAXCLIPPED = 3;

  public static final double CLIP_MIN = -100000, CLIP_MAX = -CLIP_MIN;

  /**
   * Find intersections of hyperbolic arm with a line.
   *
   * @param s0 : first point on line
   * @param s1 : second point on line
   * @param ipts : intersection points returned here
   */
  public DArray findLineIntersect(FPoint2 s0, FPoint2 s1, DArray ipts) {
    if (ipts == null)
      ipts = new DArray();
    ipts.clear();

    // transform both line points to curve space

    FPoint2 c0 = toCurveSpace(s0, null), c1 = toCurveSpace(s1, null);

    double a = c0.x, b = c0.y, c = c1.x, d = c1.y;

    double e = d - b;
    double f = c - a;

    double D = f * f - B * e * e;
    double E = 2 * a * f - 2 * b * B * e;
    double F = a * a - A - B * b * b;
    Polyn q = new Polyn(D, E, F);

    DArray roots = new DArray();
    q.solve(roots);

    for (int i = 0; i < roots.size(); i++) {
      double k = roots.getDouble(i);
      FPoint2 pt =
      //      ipts.add(
      new FPoint2(s0.x + (s1.x - s0.x) * k, s0.y + (s1.y - s0.y) * k);
      //      );

      // Make sure this point is actually on the arm.
      FPoint2 cpt = toCurveSpace(pt, null);
      if (cpt.x < 0) {
        continue;
      }
      ipts.add(pt);
    }
    return ipts;
  }

  /**
   * Construct a hyperbola
   *
   * @param f1 FPoint2 : first focus
   * @param f2 FPoint2 : second focus
   * @param interceptDistance : closest distance of point on arm
   *   to f1
   */
  public Hyperbola(FPoint2 f1, FPoint2 f2, double interceptDistance) {
    double fDist = f2.distance(f1);
    if (!(interceptDistance >= 0 && interceptDistance <= fDist))
      throw new FPError("Hyperbola construction: icept="
          + Tools.f(interceptDistance) + " of max " + Tools.f(fDist) + "\n f1="
          + f1 + " f2=" + f2);

    double ratio = 0;
    if (fDist > 0) {
      ratio = interceptDistance / fDist;
    }
    FPoint2 pt = FPoint2.interpolate(f1, f2, ratio);
    construct(f1, f2, pt);
  }

  /**
   * Construct a hyperbola that is a bisector of two points
   * @param f1 FPoint2
   * @param f2 FPoint2
   */
  public Hyperbola(FPoint2 f1, FPoint2 f2) {
    construct(f1, f2, null);
  }

  public Hyperbola(Hyperbola src) {
    this.a = src.a;
    this.A = src.A;
    this.c = src.c;
    this.B = src.B;
    this.flipped = src.flipped;
    this.foci = src.foci;
    this.label = src.label;
    this.origin = src.origin;
    this.poly = src.poly;
    this.pt = src.pt;
    this.toE2 = src.toE2;
    this.toW2 = src.toW2;
    this.userData = new DArray();
    userData.addAll(src.userData);
    this.valid = src.valid;
    this.visSeg = (DArray) src.visSeg.clone();
  }

  /**
   * Constructor
   * @param f1 FPoint2
   * @param f2 FPoint2
   * @param pt FPoint2, or null for bisector
   */
  private void construct(FPoint2 f1, FPoint2 f2, FPoint2 pt) {

    //    userData[LEFT] = new DArray();
    //    userData[RIGHT]  =new DArray();

    final boolean db = false;
    if (db) {
      System.out.println("Hyperbola constructor\n f1=" + f1 + "\n f2=" + f2
          + "\n pt=" + pt);
    }
    boolean bisector = (pt == null);
    initializeVisibleSegments();

    // if point on arm is closer to f2 than f1, swap f1 & f2.

    if (!bisector
        && FPoint2.distanceSquared(f1, pt) > FPoint2.distanceSquared(f2, pt)) {
      flipped = true;
    }

    this.foci[RIGHT] = new FPoint2(f1);
    this.foci[LEFT] = new FPoint2(f2);
    if (!bisector) {
      this.pt = new FPoint2(pt);
    }

    double fociDist = FPoint2.distance(f1, f2);
    if (fociDist == 0) {
      throw new FPError("Hyperbola foci are same point");
    }

    c = fociDist * .5;

    // calculate the translation of the hyperbola away from
    // standard position.

    FPoint2 rFocus = getFocus(0), lFocus = getFocus(1);

    origin = new FPoint2(.5 * (rFocus.x + lFocus.x), .5 * (rFocus.y + lFocus.y));

    // calculate the angle of rotation of the hyperbola away
    // from the standard position.

    double theta = Math.atan2(rFocus.y - lFocus.y, rFocus.x - lFocus.x);

    Matrix fromCenterInW = Matrix.getTranslate(origin, true);
    Matrix rotToE = Matrix.getRotate(-theta);

    toE2 = rotToE;
    Matrix.mult(toE2, fromCenterInW, toE2);
    // calculate inverse

    toW2 = toE2.invert(null);

    //      Matrix toCenterInW = Matrix.translationMatrix(origin, false);
    //      Matrix rotToW = Matrix.getRotate2D(theta);
    //
    //      toW2 = toCenterInW;
    //      Matrix.mult(toW2, rotToW, toW2);
    //      Tools.warn("just invert matrix here");
    //

    if (bisector) {
      valid = true;
    } else {
      // get the arm point in hyperbola space.
      FPoint2 workPt = toE2.apply(pt, null);

      double xs = workPt.x * workPt.x;
      double cs = c * c;

      Polyn q = new Polyn(1, -(cs + xs + workPt.y * workPt.y), cs * xs);
      if (db) {
        System.out.println("a2 quadratic:\n" + q);
      }
      final DArray qsoln = new DArray();
      q.solve(qsoln);
      if (db) {
        Streams.out.println(qsoln);
      }
      double val = q.c(1) * -.5;
      int ql = qsoln.size();
      if (ql >= 1) {
        val = qsoln.getDouble(0);
      }

      // choose the root that is less than c*c.

      if (ql == 2) {
        if (val > qsoln.getDouble(1)) {
          val = qsoln.getDouble(1);
          if (db) {
            System.out.println(" two roots, choosing smaller.");
          }
        }
      }
      if (db) {
        System.out.println(" root chosen=" + val);
      }

      a = Polyn.sqrt(val);
      A = a * a;
      B = A / (c * c - A);
    }
    valid = true;
    if (db) {
      System.out.println(" ==> " + this);
    }
  }

  /**
   * Construct a hyperbolic arm by specifying the two focii and a point on
   * the arm.
   * @param f1 : location of focus to right of arm
   * @param f2 : location of focus to left of arm
   * @param pt : location of a point on the hyperbola
   */
  public Hyperbola(FPoint2 f1, FPoint2 f2, FPoint2 pt) {
    construct(f1, f2, pt);
  }

  /**
   * Find intersections of hyperbolic arm with an axes-aligned
   * line segment.  Arm must intersect segment, not its
   * underlying (infinite) line.  The intersection points are
   * sorted into increasing parameter values w.r.t. the arm.
   *
   * @param s0 : start of line segment
   * @param s1 : end of line segment
   * @param ipts : intersection points returned here
   * @param reverseOrder : if true, points are sorted by
   *  decreasing parameter values
   * @param dbFlag : true to print debug information
   */
  public void findOrthogonalIntersect(FPoint2 s0, FPoint2 s1, DArray ipts,
      boolean reverseOrder, boolean dbFlag) {

    final boolean db = true && dbFlag;

    if (db) {
      System.out.println("findOrthogonalIntersect " + s0 + " -> " + s1
          + " rev:" + Tools.f(reverseOrder));
    }

    // Determine whether this is a vertical or horizontal line segment.

    boolean vert = (Math.abs(s1.y - s0.y) > Math.abs(s1.x - s0.x));
    if (db) {
      System.out.println(" vert=" + Tools.f(vert));
    }

    // Find the quadratic to solve.
    Polyn p;
    PlaneCurve cv = getCurve();
    if (vert) {
      p = cv.solveForX(s0.x);
    } else {
      p = cv.solveForY(s0.y);
    }

    DArray lst = new DArray();
    p.solve(lst);

    if (db) {
      System.out.println(" curve=" + cv.toString(true));
      System.out.println(" polyn=" + p.toString(true));
      System.out.println(" roots=" + Tools.d(lst));
    }

    // Sort points, discarding those not on line segment,
    // and those not on the correct arm.
    ipts.clear();

    for (int i = 0; i < lst.size(); i++) {
      double ta = lst.getDouble(i);
      FPoint2 pt;
      if (vert) {
        pt = new FPoint2(s0.x, ta);
      } else {
        pt = new FPoint2(ta, s0.y);
      }

      if (db) {
        System.out.println(" position on arm for ta=" + ta + " is " + pt);
      }
      double t = MyMath.positionOnSegment(pt, s0, s1);
      if (db) {
        System.out.println("  pt=" + pt + " t=" + t);
      }
      if (t < 0 || t > 1) {
        if (db) {
          System.out.println("   not on segment, skipping");
        }
        continue;
      }

      FPoint2 cpt = toCurveSpace(pt, null);
      if (db) {
        System.out.println("  curveSpace=" + cpt);
      }
      if (cpt.x < 0) {
        if (db) {
          System.out.println("   skipping...");
        }
        continue;
      }

      double t1 = calcParameter(pt);
      int j = 0;
      while (true) {
        if (j == ipts.size()) {
          break;
        }
        double t2 = calcParameter(ipts.getFPoint2(j));
        if (!reverseOrder) {
          if (t1 < t2) {
            break;
          }
        } else if (t1 > t2) {
          break;
        }
        j++;
      }
      ipts.add(j, pt);
    }
    if (db) {
      System.out.println(" ipts=" + Tools.d(ipts));
    }
  }

  //  public Object getData(int side) {
  //    return getData(0, side);
  //  }

  //  public void setData(Object dRight, Object dLeft) {
  //    setData(RIGHT, dRight);
  //    setData(LEFT,dLeft);
  //  }

  public Object getData(int field) {
    Object ret = null;
    if (userData.exists(field))
      ret = userData.get(field);
    return ret;
  }

  public Object getOtherData(Object first) {
    int t = 0;
    if (first == userData.get(t))
      t++;
    return userData.get(t);
  }

  public void setData(int field, Object data) {
    userData.growSet(field, data);
  }

  /**
   * @param dRight
   * @param dLeft
   */
  public void setData(Object dRight, Object dLeft) {
    setData(RIGHT, dRight);
    setData(LEFT, dLeft);
  }

  /**
   * Get number of segments in clipped arm
   * @return # segments; 0 if clipped out
   */
  private int nClip() {
    int len = visSeg.size() / 2;
    return len;
  }

  public boolean isEmpty() {
    return nClip() == 0;
  }

  private boolean flipped;

  /**
   * Get a particular focus
   * @param side : RIGHT or LEFT, after flipping; thus RIGHT focus
   *  is on the INSIDE of the arm
   * @return FPoint2
   */
  private FPoint2 getFocus(int side) {
    return foci[flipped ? side ^ 1 : side];
  }

  /**
   * Determine if hyperbola has been flipped (rotated might be more
   * accurate) from its constructed form to its internal representation.
   * @return true if arm curves to the left as t increases
   */
  public boolean flipped() {
    return flipped;
  }

  /**
   * Convert parameter to internal value (i.e. flip if necessary)
   * @param t double
   * @return double
   */
  private double toInt(double t) {
    if (flipped()) {
      t = -t;
    }
    return t;
  }

  /**
   * Convert parameter to external value (i.e. flip if necessary)
   * @param t double
   * @return double
   */
  private double toExt(double t) {
    if (flipped()) {
      t = -t;
    }
    return t;
  }

  public String visString() {
    StringBuilder sb = new StringBuilder();
    if (visSeg != null) {
      for (int i = 0; i < visSeg.size(); i += 2) {
        double t0 = 0, t1 = 0;
        if (flipped()) {
          int j = visSeg.size() - i - 2;
          t0 = toExt(visSeg.getDouble(j + 1));
          t1 = toExt(visSeg.getDouble(j));
        } else {
          t0 = visSeg.getDouble(i);
          t1 = visSeg.getDouble(i + 1);
        }

        sb.append("<");
        sb.append(Tools.f(t0));
        sb.append("..");

        sb.append(Tools.f(t1));
        sb.append("> ");
      }
    }
    return sb.toString();
  }

  /**
   * Construct a string
   * @param javaMode : if true, generates java code to construct it
   * @return String
   */
  private String toString(boolean javaMode) {

    StringBuilder sb = new StringBuilder();
    //    sb.setLength(0);

    if (javaMode) {
      sb.append("new Hyperbola(new FPoint2(" + foci[RIGHT].x + ","
          + foci[RIGHT].y + "), new FPoint2(" + foci[LEFT].x + ","
          + foci[LEFT].y + "), new FPoint2(" + pt.x + "," + pt.y + "));");
    } else {
      sb.append("hyperbola: ");
      if (!valid) {
        sb.append("***INVALID***");
      } else {
        sb.append("f(r)=" + foci[RIGHT] + " f(l)=" + foci[LEFT] + " a="
            + Tools.f(a) + " c=" + Tools.f(c));
        sb.append(" flip=" + Tools.f(flipped()));
        if (visSeg != null) {
          sb.append("\nVisible segments: ");
          for (int i = 0; i < visSeg.size(); i += 2) {
            double t0 = 0, t1 = 0;
            if (flipped()) {
              int j = visSeg.size() - i - 2;
              t0 = toExt(visSeg.getDouble(j + 1));
              t1 = toExt(visSeg.getDouble(j));
            } else {
              t0 = visSeg.getDouble(i);
              t1 = visSeg.getDouble(i + 1);
            }

            sb.append("<");
            sb.append(Tools.f(t0));
            sb.append("..");

            sb.append(Tools.f(t1));
            sb.append("> ");
          }
          sb.append("\n");
        }
      }
    }
    return sb.toString();
  }

  /**
   * Construct description
   * @return String
   */
  public String toString() {
    if (label != null)
      return label;
    return toString(false);
  }

  /**
   * Calculate a point on the hyperbola
   * @param t : parameter in internal space (after flipping has occurred)
   * @param dest : where to store the calculated point
   */
  private FPoint2 calcPointInternal(double t, FPoint2 dest) {
    final boolean db = false;
    if (!valid)
      throw new FPError("calcPoint of invalid hyperbola");

    //    Tools.ASSERT(valid, "calcPoint of invalid hyperbola");
    //    final FPoint2 work = new FPoint2();
    //
    //    work.setLocation(Polyn.sqrt(A + t * t * B), t);
    dest = toW2.apply(Polyn.sqrt(A + t * t * B), t, dest);
    //
    //    if (dest == null)
    //      dest = new FPoint2();
    //    Matrix.apply(toW2, work, dest);
    if (db) {
      System.out.println("calcPoint t=" + Tools.f(toExt(t)) //+ " armsp=" + work
          + " world=" + dest);
    }
    return dest;
  }

  /**
   * Calculate a point on the hyperbola
   * @param t : parameter: positive for point above x-axis, negative for below
   * @param dest : where to store the calculated point
   */
  public FPoint2 calcPoint(double t, FPoint2 dest) {
    return calcPointInternal(toInt(t), dest);
  }

  /**
   * Calculate a point on the hyperbola, leave in arm space
   * @param t : parameter: positive for point above x-axis, negative for below
   * @return point in arm space
   */
  private FPoint2 calcPointInArmSpace(double t) {
    return calcPointInArmSpace0(toInt(t));
  }

  /**
   * Calculate a point on the hyperbola, leave in arm space
   * @param t : parameter, after conversion to 'internal' value
   * @return point in arm space
   */
  private FPoint2 calcPointInArmSpace0(double t) {
    return new FPoint2(Polyn.sqrt(A + t * t * B), t);
  }

  /**
   * Calculate point on curve
   * @param t : parameter: angle of rotation around unit circle
   * @return point on curve
   */
  public FPoint2 calcPoint(double t) {
    FPoint2 dest = new FPoint2();
    calcPoint(t, dest);
    return dest;
  }

  /**
   * Transform a point from world space to curve space
   * @param pt : point in world space
   * @param out : where to store curve space point
   */
  public FPoint2 toCurveSpace(FPoint2 pt, FPoint2 out) {
    return toE2.apply(pt, out);
  }

  /**
   * Transform a point from curve space to world space
   * @param pt : point in curve space
   * @return point in world space
   */
  public FPoint2 toWorldSpace(FPoint2 pt) {
    return toWorldSpace(pt.x, pt.y, null);
  }

  /**
   * Transform a point from curve space to world space
   * @param pt : point in curve space
   * @param out : where to store world space point
   */
  public FPoint2 toWorldSpace(FPoint2 pt, FPoint2 out) {
    return toWorldSpace(pt.x, pt.y, out);
  }

  public FPoint2 toWorldSpace(double d, double e, FPoint2 out) {
    return toW2.apply(d, e, out);
  }

  /**
   * Get world space coordinates of curve origin.
   * Note: this is NOT the intersection of the arm with the x-axis!
   * For that, use calcPoint(0).
   *
   * @return FPoint2
   */
  public FPoint2 origin() {
    return origin;
  }

  public double slopeOfAsymptote() {
    return Math.sqrt(c * c - a * a) / a;
  }

  /**
   * Determine if hyperbola is valid (radius is large enough for focal
   * points)
   * @return boolean
   */
  public boolean valid() {
    return valid;
  }

  /**
   * Calculate the parameter associated with a point
   * @param pt : point to find parameter for
   * @return best-fit parameter for this point
   */
  public double calcParameter(FPoint2 pt) {
    final FPoint2 ept = new FPoint2();
    toCurveSpace(pt, ept);
    return toExt(ept.y);
  }

  /**
   * Snap point to a point on the hyperbola
   * @param pt point
   * @return point on hyperbola
   */
  public FPoint2 snap(FPoint2 pt) {
    double t = calcParameter(pt);
    return calcPoint(t);
  }

  /**
   * Calculate parameter for a point, and indicate how close the
   * point is to the curve
   * @param pt : point to find parameter for
   * @return FPoint2; x contains the parameter, y is its distance
   *  from the curve
   */
  private FPoint2 calcParameterAndDistance(FPoint2 pt) {
    final FPoint2 ept = new FPoint2();
    toCurveSpace(pt, ept);
    double t = toExt(ept.y);
    FPoint2 as = calcPointInArmSpace(t);
    return new FPoint2(t, Math.abs(ept.x - as.x));
  }

  /**
   * Test program for Hyperbola class
   * @param args String[]
   */
  public static void main(String[] args) {

    final double[] pts = { //
    100, 0, -100, 0, 75, 20,

    //  120, 30, -100, -10, 70, 50,
    };

    for (int i = 0; i < pts.length; i += 6) {
      try {
        Hyperbola h = new Hyperbola(new FPoint2(pts[i + 0], pts[i + 1]),
            new FPoint2(pts[i + 2], pts[i + 3]), new FPoint2(pts[i + 4],
                pts[i + 5]));
        System.out.println("Constructed:\n" + h);

        for (double t = -50; t <= 50; t += 10) {

          FPoint2 pt = h.calcPoint(t);

          FPoint2 pt2 = new FPoint2(pt.x, pt.y + 5);

          double tClosest = h.closestPointTo(pt2);

          System.out.println("t=" + Tools.f(t) + " pt=" + pt + " closest="
              + tClosest);

          if (t == -20) {
            for (double t2 = tClosest - .1; t2 <= tClosest + .1; t2 += .01) {
              FPoint2 pt3 = h.calcPoint(t2);
              Streams.out.println("t2=" + t2 + " dist=" + pt3.distance(pt2));
            }

          }

        }

      } catch (TBError e) {
        System.out.println(e.toString());
      }

    }
  }

  /**
   * Determine where point is relative to arm
   * @param x :
   * @param y : point to test
   * @return an integer, which is 0 if it's on the arm,
   *   1 if it's to the right of the arm,
   *   -1 if it's to the left of the arm
   */
  public int testPoint(double x, double y) {
    final FPoint2 ept = new FPoint2(), eps = new FPoint2();
    eps.setLocation(x, y);
    toCurveSpace(eps, ept);
    FPoint2 as = calcPointInArmSpace(toExt(ept.y));
    int out = 0;
    double diff = ept.x - as.x;
    if (diff > 0) {
      out = 1;
    } else if (diff < 0) {
      out = -1;
    }
    if (flipped()) {
      out = -out;
    }
    return out;
  }

  /**
   * Determine where point is relative to arm
   * @param pt point to test
   * @return an integer, which is 0 if it's on the arm,
   *   1 if it's to the right of the arm,
   *   -1 if it's to the left of the arm
   */
  public int testPoint(FPoint2 pt) {
    return testPoint(pt.x, pt.y);
  }

  /**
   * Determine if a parameter value's point has not been clipped away
   * @param t : parameter
   * @return true if it hasn't been clipped
   */
  public boolean containsPoint(double t) {
    return containsPoint(t, false);
    /*
         t = toInt(t);
         for (int i = 0; i < visSeg.size(); i += 2) {
      if (t >= visSeg.getDouble(i) && t <= visSeg.getDouble(i + 1)) {
        return true;
      }
         }
         return false;
     */
  }

  /**
   * Determine if a parameter value's point has not been clipped away
   * @param t : parameter
   * @param withFuzzFactor : if true, allows t to be slightly outside
   *  of clipped endpoints
   * @return true if it hasn't been clipped
   */
  public boolean containsPoint(double t, boolean withFuzzFactor) {
    final boolean db = false;
    if (db) {
      System.out.println("contains point " + t);
    }
    t = toInt(t);
    for (int i = 0; i < visSeg.size(); i += 2) {
      if (withFuzzFactor) {
        final double FUZZ = .0001;
        if (db) {
          System.out.println(" bounds=" + visSeg.getDouble(i) + " and "
              + visSeg.getDouble(i + 1));
        }
        if (t + FUZZ >= visSeg.getDouble(i)
            && t - FUZZ <= visSeg.getDouble(i + 1)) {
          return true;
        }
      } else {
        if (t >= visSeg.getDouble(i) && t <= visSeg.getDouble(i + 1)) {
          return true;
        }
      }
    }
    return false;
  }

  public double getClipParameter(int clipSeg, boolean max) {

    int i = clipSeg * 2 + (max ? 1 : 0);
    if (flipped()) {
      i = nClip() - 1 - i;
    }

    return toExt(visSeg.getDouble(i));

  }

  /**
   * Return minimum existing parameter value
   * @return double
   */
  public double minParameter() {
    if (visSeg.isEmpty())
      throw new FPError("minParameter with empty arm");
    return toExt(visSeg.getDouble(flipped() ? visSeg.size() - 1 : 0));
  }

  /**
   * Return maximum existing parameter value
   * @return double
   */
  public double maxParameter() {
    if (visSeg.isEmpty())
      throw new FPError("maxParameter with empty arm");
    return toExt(visSeg.getDouble(flipped() ? 0 : visSeg.size() - 1));
  }

  /**
   * Clip out the entire arm
   */
  public void clipAll() {
    //  initClipIfNec();
    visSeg.clear();
  }

  /**
   * Clip out a section of the arm
   * @param c0 : start of section to clip out
   * @param c1 : end of section to clip out; if < c0, c0 & c1 are swapped
   */
  public void clip(double c0, double c1) {
    clip0(toInt(c0), toInt(c1));
  }

  public void clipAllBut(double c0, double c1) {
    if (c0 > c1) {
      double tmp = c1;
      c1 = c0;
      c0 = tmp;
    }

    clip(CLIP_MIN, c0);
    clip(c1, CLIP_MAX);
  }

  public void setVisibleRange(double c0, double c1) {
    visSeg = new DArray();
    c0 = toInt(c0);
    c1 = toInt(c1);
    if (c0 > c1) {
      double t = c0;
      c0 = c1;
      c1 = t;
    }
    visSeg.addDouble(c0);
    visSeg.addDouble(c1);
  }

  public boolean containsRange(double c0, double c1) {
    final boolean db = false;
    if (c0 > c1) {
      double t = c0;
      c0 = c1;
      c1 = t;
    }
    boolean ret = false;
    final double EPS = 1e-4;
    for (int i = 0; i < visSeg.size(); i += 2) {

      double t0 = visSeg.getDouble(i), t1 = visSeg.getDouble(i + 1);
      if (t0 <= c0 + EPS && t1 >= c1 - EPS) {
        ret = true;
        break;
      }
    }
    if (db)
      Streams.out.println("containsRange,\n " + this + " " + this.visString()
          + "\n" + ", " + c0 + " ... " + c1 + " ret " + ret);
    return ret;
  }

  /**
   * Clip out a section of the arm, using internal parameter values
   * @param c0 : start of section to clip out
   * @param c1 : end of section to clip out; if < c0, c0 & c1 are swapped
   */
  private void clip0(double c0, double c1) {
    final boolean db = false;

    if (c0 > c1) {
      double t = c0;
      c0 = c1;
      c1 = t;
    }

    // array to store new clip elements in
    DArray nc = new DArray();

    if (db) {
      System.out.println("clip " + Tools.f(c0) + "..." + Tools.f(c1));
    }

    for (int i = 0; i < visSeg.size(); i += 2) {

      double t0 = visSeg.getDouble(i), t1 = visSeg.getDouble(i + 1);
      if (db) {

        System.out.println(" seg " + i + " " + Tools.f(t0) + "..."
            + Tools.f(t1));
      }

      // if this seg is outside clip region, leave intact (add it)
      if (t0 >= c1 || t1 <= c0) {
        nc.addDouble(t0);
        nc.addDouble(t1);
        continue;
      }

      if (t0 >= c0 && t1 <= c1) {
        if (db) {
          System.out.println("  clipping entire segment");
        }
        // remove entire segment, by skipping it
        continue;
      }

      // add unclipped portion
      if (t0 < c0) {
        nc.addDouble(t0);
        nc.addDouble(c0);
      }
      if (t1 > c1) {
        nc.addDouble(c1);
        nc.addDouble(t1);
      }

    }
    visSeg = nc;
    if (db) {
      System.out.println(" after clipping:" + visSeg);
    }
  }

  private void initializeVisibleSegments() {
    visSeg = new DArray();
    visSeg.addDouble(CLIP_MIN);
    visSeg.addDouble(CLIP_MAX);
  }

  //  public void render() {
  //    render(0, false);
  //  }

  //  private void render(double step, boolean dashed ) {
  //    final boolean db = false;
  //    // Get array of visible segments.  If no such array exists,
  //    // use default.
  //    DArray vseg = visSeg;
  //
  //    if (step == 0) {
  //      step = renderStep();
  //    }
  ////    vp V = TestBed.view();
  //
  //    if (db)
  //      Streams.out.println(" step=" + step);
  //
  //    // plot each visible segment
  //
  //    for (int seg = 0; seg < vseg.size(); seg += 2) {
  //      double t0 = vseg.getDouble(seg + 0), t1 = vseg.getDouble(seg + 1);
  //      t0 = MyMath.clamp(t0, -500.0, 500.0);
  //      t1 = MyMath.clamp(t1, -500.0, 500.0);
  //
  //      // render() expects external parameters.
  //
  //      double s0 = toExt(t0), s1 = toExt(t1);
  //      if (s0 > s1) {
  //        double tmp = s0;
  //        s0 = s1;
  //        s1 = tmp;
  //      }
  //      FPoint2 p0 = calcPoint(s0), p1 = calcPoint(s1);
  //      if (db)
  //        Streams.out.println(" p0=" + p0 + ", p1=" + p1);
  //      if (isLine() && !dashed) {
  //        V.drawLine(p0, p1);
  //      } else {
  //
  //        /*
  //                   if (Math.abs(s0) >= 500
  //            ||Math.abs(s1) >= 500)
  //          System.out.println("Rendering "+t0+" to "+t1+" step "+step);
  //         */
  //        if (dashed)
  //          V.pushStroke(Globals.STRK_RUBBERBAND);
  //        {
  //          int count = 0;
  //          boolean first = true;
  //          for (double t = t0;; t += step, count++) {
  //            boolean last = (t >= t1);
  //            if (last) {
  //              t = t1;
  //            }
  //            calcPointInternal(t, p1);
  //
  //            if (!p1.isValid()) {
  //              if (last) {
  //                break;
  //              }
  //              continue;
  //            }
  //
  //            if (db) {
  //              System.out.println(" calcPt " + Tools.f(toExt(t)) + " = " + p1.x
  //                  + "," + p1.y);
  //            }
  //            if (!first) {
  //              V.drawLine(p0, p1);
  //              if (false) {
  //                Tools.warn("highlighting int");
  //                V.mark(p0);
  //              }
  //            }
  //            if (last) {
  //              break;
  //            }
  //            p0.setLocation(p1);
  //            first = false;
  //          }
  //        }
  //        if (dashed)
  //          V.popStroke();
  //      }
  //    }
  //
  //  }

  /**
   * Calculate dx, dy values for a point
   * @param t : parameter
   * @param delta : where to store dx, dy values
   */
  private void calcTangentAt(double t, FPoint2 delta) {
    final boolean db = false;
    if (db) {
      System.out.println("calcTangentAt " + t);
    }
    t = toInt(t);
    if (db) {
      System.out.println(" flipped=" + flipped() + " ti=" + t);
    }
    double a = toW2.get(0, 0), b = toW2.get(0, 1), c = toW2.get(0, 2);
    double d = toW2.get(1, 0), e = toW2.get(1, 1), f = toW2.get(1, 2);
    if (db) {
      System.out.println(" a=" + a + " b=" + b + " c=" + c + "\n d=" + d
          + " e=" + e + " f=" + f);
    }
    double rt = Polyn.sqrt(A + B * t * t);
    double dx = a * B * t / rt + b;
    double dy = d * B * t / rt + e;

    if (flipped()) {
      dx = -dx;
      dy = -dy;
    }
    if (db) {
      System.out.println(" dx=" + dx + "\n dy=" + dy);
      System.out.println(" ratio=" + (dy / dx));
    }
    delta.setLocation(dx, dy);
  }

  private static final FPoint2 workPt = new FPoint2();

  /**
   * Calculate angle of tangent line at a particular point
   * @param t : parameter
   * @return angle
   */
  public double calcDirectionAt(double t) {
    final boolean db = false;
    calcTangentAt(t, workPt);
    double a = Math.atan2(workPt.y, workPt.x);
    if (db) {
      System.out.println("calcDirectionAt " + Tools.f(t) + " a=" + Tools.fa(a));
    }
    return a;
  }

  /**
   * Get the PlaneCurve associated with this hyperbola
   * @return PlaneCurve
   */
  public PlaneCurve getCurve() {
    //Tools.warn("testing planecurve cache problem?");
    if (poly == null) {
      if (isLine()) {

        double A, B, C;

        A = toE2.get(0, 0);
        B = toE2.get(0, 1);
        C = toE2.get(0, 2);

        //        double D = toE.getShearY(),  E = toE.getScaleY(), F = toE.getTranslateY();
        poly = new PlaneCurve(0, // x^2
            0, // xy
            0, // y^2
            A, // x
            B, // y
            C);
        /*
            0,0,0,
            E,
            -B,
            (F-C)*E);
         */
      } else {
        //
        // We note that x', y' are 'curve space' coordinates, which
        // are derived by transforming x,y by the toE matrix;
        // that is,
        //
        //  x' = x * A + y * B + C
        //  y' = x * D + y * E + F
        //
        //
        double A, B, C, D, E, F, G, H;

        A = toE2.get(0, 0);
        B = toE2.get(0, 1);
        C = toE2.get(0, 2);
        D = toE2.get(1, 0);
        E = toE2.get(1, 1);
        F = toE2.get(1, 2);
        G = a * a - c * c;
        H = this.A;
        if (Double.isNaN(G) || Double.isNaN(H)) {
          System.out.println("not a number: " + G + " or H=" + H);
        }

        poly = new PlaneCurve(A * A * G + D * D * H,
            2 * (A * B * G + D * E * H), B * B * G + E * E * H,
            2 * (A * C * G + D * F * H), 2 * (B * C * G + E * F * H), C * C * G
                + H * (F * F - G));

      }
    }
    return poly;
  }

  /**
   * Calculate standard step value for rendering the arm.
   * @return amount to add to parameter for each segment
   *   approximating curve
   */
  private double renderStep() {
    double step = Math.max(1, ((c - a) * (c - a) / c) * .3);
    return step;
  }

  /**
   * Find intersections between two hyperbolas
   * @param a Hyperbola
   * @param b Hyperbola
   * @param iPts where to store intersection points; null to construct
   */
  public static DArray findIntersections(Hyperbola a, Hyperbola b, DArray iPts) {
    final boolean db = false;

    if (iPts == null)
      iPts = new DArray();
    iPts.clear();
    final DArray jPts = new DArray();
    //a.initClipIfNec();
    //    b.initClipIfNec();
    if (db) {
      System.out.println("finding intersections between:\n" + a.toString(true)
          + "\n" + b.toString(true));
    }

    PlaneCurve.findIntersect(a.getCurve(), b.getCurve(), jPts);

    // filter out false intersections
    if (db) {
      System.out.println("prefilter # intersections= " + jPts.size());
    }

    FPoint2 ptc = new FPoint2();

    for (int i = 0; i < jPts.size(); i++) {
      FPoint2 pt = jPts.getFPoint2(i);

      // make sure calculated intersect point is actually on both arms
      {
        FPoint2 data = a.calcParameterAndDistance(pt);
        if (data.y > .001)
          continue;
        data = b.calcParameterAndDistance(pt);
        if (data.y > .001)
          continue;
      }

      if (!a.isLine()) {
        // put in curve space to verify it's to the right of the y axis
        a.toCurveSpace(pt, ptc);
        if (db) {
          System.out.println("Filter " + i + " in curveA= " + ptc);
        }
        if (ptc.x <= 0) {
          if (db) {
            System.out.println(" < 0");
          }
          continue;
        }
      }
      if (!b.isLine()) {
        b.toCurveSpace(pt, ptc);
        if (db) {
          System.out.println("        in curveB= " + ptc);
        }
        if (ptc.x <= 0) {
          if (db) {
            System.out.println(" < 0");
          }
          continue;
        }
      }
      if (db) {
        System.out.println(" adding " + pt);
      }
      iPts.add(pt);
    }
    return iPts;
  }

  //  /**
  //   * Find t value such that the disc, Ct, centered at p(t), which has
  //   * inside-tangent a disc centered at f1 with radius f1Radius, is outside-tangent
  //   * to the disc centered at qPt with radius qRad
  //   * @param f1Radius
  //   * @param qPt
  //   * @param qRad
  //   * @return
  //   */
  //  public double findTangentCircle(double f1Radius, FPoint2 qPt, double qRad) {
  //    // get qPt in arm space
  //    qPt = toCurveSpace(qPt,null);
  //
  //
  //  }

  public boolean isLine() {
    return A == 0;
  }

  public void setLabel(String s) {
    this.label = s;
  }

  // list of visible segments, stored as ti...tj pairs of doubles
  private DArray visSeg;

  private double a, c;

  // position of foci in world space
  private FPoint2[] foci = new FPoint2[2];

  // original point on arm (only needed for generating Java code
  // in toString(...))
  private FPoint2 pt;

  private FPoint2 origin;

  // a*a / (c*c - a*a)
  private double B;
  private double A;

  // true if hyperbola construction occurred without problems
  private boolean valid;

  private Matrix toE2, toW2;

  // polynomial representing the hyperbola equation
  // (null if it hasn't been required yet)
  private PlaneCurve poly;

  private DArray userData = new DArray();

  private String label;

  public FPoint2 getSamplePoint() {
    double t = 0;
    if (!containsPoint(t)) {
      double t0 = minParameter();
      double t1 = maxParameter();
      if (t0 > 0)
        t = t0;
      else
        t = t1;
    }
    return calcPoint(t);
  }

  /**
   * Determine closest point on hyperbolic arm to a point pt
   * @param pt : point
   * @return t value for closest point; not necessarily within clipped range
   */
  public double closestPointTo(FPoint2 pt) {

    final boolean db = false;

    pt = toCurveSpace(pt, null);

    double U = B + 1;
    double V = -pt.y;
    double W = B * pt.x;

    Polyn p = new Polyn(//
        B * U * U, //
        2 * B * U * V, //
        B * V * V + A * U * U - W * W, //
        2 * A * U * V, //
        A * V * V//
    );

    if (db)
      Streams.out.println("closestPointTo, pt=" + pt + "\n" + p);
    double ret = 0;

    try {
      DArray r = new DArray();

      if (Math.abs(p.c(0)) < 1e-5)
        r.addDouble(0);
      else
        p.solve(r);

      if (r.isEmpty()) {
        throw new FPError("can't find closest point, poly=\n" + p);
      }

      double bestDist = 0;

      for (int i = 0; i < r.size(); i++) {
        double t = r.getDouble(i);
        FPoint2 apt = calcPoint(t);
        double dist = apt.distance(pt);
        if (i == 0 || dist < bestDist) {
          bestDist = dist;
          ret = t;
        }
      }
    } catch (FPError e) {
      Tools.warn("caught FPError");
      //      Streams.out.println("caught:\n" + e);
      ret = (this.minParameter() + this.maxParameter()) * .5;
    }

    return ret;
  }

  /**
   * Convert array of hyperbolas to hyperbolas with single connected range
   * @param h1
   * @return
   */
  public static Hyperbola[] singlify(Hyperbola[] h1) {
    DArray e = new DArray();
    for (int i = 0; i < h1.length; i++) {
      Hyperbola hi = h1[i];

      if (hi.visSeg.size() == 1) {
        e.add(hi);
      } else {
        for (int j = 0; j < hi.visSeg.size(); j += 2) {
          double t0 = hi.visSeg.getDouble(j + 0);
          double t1 = hi.visSeg.getDouble(j + 1);
          Hyperbola h = new Hyperbola(hi);
          h.visSeg.clear();
          h.visSeg.addDouble(t0);
          h.visSeg.addDouble(t1);
          e.add(h);
        }
      }
    }
    return (Hyperbola[]) e.toArray(Hyperbola.class);
  }

  public void render(Color c, int stroke, int markType) {
    final boolean db = false;
    // Get array of visible segments.  If no such array exists,
    // use default.
    DArray vseg = visSeg;
    boolean dashed = false;

    //    if (step == 0) {
    double step = renderStep();
    //    }
    //    vp V = TestBed.view();

    if (db)
      Streams.out.println(" step=" + step);

    // plot each visible segment

    for (int seg = 0; seg < vseg.size(); seg += 2) {
      double t0 = vseg.getDouble(seg + 0), t1 = vseg.getDouble(seg + 1);
      t0 = MyMath.clamp(t0, -500.0, 500.0);
      t1 = MyMath.clamp(t1, -500.0, 500.0);

      // render() expects external parameters.

      double s0 = toExt(t0), s1 = toExt(t1);
      if (s0 > s1) {
        double tmp = s0;
        s0 = s1;
        s1 = tmp;
      }
      FPoint2 p0 = calcPoint(s0), p1 = calcPoint(s1);
      if (db)
        Streams.out.println(" p0=" + p0 + ", p1=" + p1);
      if (isLine() && !dashed) {
        V.drawLine(p0, p1);
      } else {

        /*
                   if (Math.abs(s0) >= 500
            ||Math.abs(s1) >= 500)
          System.out.println("Rendering "+t0+" to "+t1+" step "+step);
         */
        if (dashed)
          V.pushStroke(Globals.STRK_RUBBERBAND);
        {
//          int count = 0;
          boolean first = true;
          for (double t = t0;; t += step) { //, count++) {
            boolean last = (t >= t1);
            if (last) {
              t = t1;
            }
            calcPointInternal(t, p1);

            if (!p1.isValid()) {
              if (last) {
                break;
              }
              continue;
            }

            if (db) {
              System.out.println(" calcPt " + Tools.f(toExt(t)) + " = " + p1.x
                  + "," + p1.y);
            }
            if (!first) {
              V.drawLine(p0, p1);
              if (false) {
                Tools.warn("highlighting int");
                V.mark(p0);
              }
            }
            if (last) {
              break;
            }
            p0.setLocation(p1);
            first = false;
          }
        }
        if (dashed)
          V.popStroke();
      }
    }

  }

}

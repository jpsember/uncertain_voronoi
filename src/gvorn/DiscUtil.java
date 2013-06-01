package gvorn;

import testbed.*;
import base.*;

public class DiscUtil {
  public static final int DISC_CONTAINED = 1 << (EdDisc.USER_FLAG_BITS_DISC - 1);
  public static final int DISC_OVERLAPPING = 1 << (EdDisc.USER_FLAG_BITS_DISC - 2);

  public static boolean contained(EdDisc d) {
    return d.hasFlags(DISC_CONTAINED);
  }
  public static boolean overlapping(EdDisc d) {
    return d.hasFlags(DISC_OVERLAPPING);
  }
  /**
   * Calculate amount disc c2 must be moved so it is inside tangent disc c
   * @param c
   * @param c2
   * @return
   */
  public static double itanDist(EdDisc c, EdDisc c2) {
    return Math.abs(FPoint2.distance(c2.getOrigin(), c.getOrigin())
        + c2.getRadius() - c.getRadius());
  }
  /**
   * Calculate amount disc c2 must be moved so it is outside tangent disc c
   * @param c
   * @param c2
   * @return
   */
  public static double otanDist(EdDisc c, EdDisc c2) {
    return Math.abs(FPoint2.distance(c2.getOrigin(), c.getOrigin())
        - c2.getRadius() - c.getRadius());
  }


  /**
   * Determine if ca is closer to being itangent to c than otangent to c
   * @param c
   * @param ca
   * @return true if origin of ca is interior to c
   */
  public static boolean itan(EdDisc c, EdDisc ca) {
    return FPoint2.distance(ca.getOrigin(), c.getOrigin()) - c.getRadius() < 0;
  }

  /**
     * Construct hyperbola containing centers of discs that have i/otan ca, cb
     * @param c disc that approximates a disc on the hyperbola, that is used
     *  to determine if ca, cb are itan or otan
     * @param ca
     * @param cb
     * @return
     */
  public static Hyperbola supportingHyperbola(EdDisc c, EdDisc ca, EdDisc cb) {

    boolean ia = itan(c, ca);
    boolean ib = itan(c, cb);

    double idist = FPoint2.distance(ca.getOrigin(), cb.getOrigin());

    if (ia)
      idist -= ca.getRadius();
    else
      idist += ca.getRadius();

    if (ib)
      idist += cb.getRadius();
    else
      idist -= cb.getRadius();

    Hyperbola h = new Hyperbola(ca.getOrigin(), cb.getOrigin(), idist * .5);
    return h;
  }

  /**
   * Construct hyperbola containing centers of discs that have itan ca,cb
   * @param ca
   * @param cb
   * @return
   */
  public static Hyperbola supportingHyperbola(EdDisc ca, EdDisc cb) {

    double idist = FPoint2.distance(ca.getOrigin(), cb.getOrigin());
    idist -= ca.getRadius();
    idist += cb.getRadius();
    Hyperbola h = new Hyperbola(ca.getOrigin(), cb.getOrigin(), idist * .5);
    return h;
  }
  public static EdDisc smallestBoundingDisc(EdDisc a, EdDisc b) {
    EdDisc ret = null;
    if (EdDisc.contains(a, b)) {
      ret = new EdDisc(a);
    } else if (EdDisc.contains(b, a)) {
      ret = new EdDisc(b);
    } else {
      double theta = MyMath.polarAngle(a.getOrigin(), b.getOrigin());
      FPoint2 pa = a.polarPoint(theta + Math.PI);
      FPoint2 pb = b.polarPoint(theta);
      ret = new EdDisc(FPoint2.midPoint(pa, pb), (pa.distance(pb)) * .5);
    }
    return ret;
  }


  public static EdDisc smallestBoundingDisc(EdDisc a, EdDisc b, EdDisc c) {
    EdDisc ret = null;


    final boolean db = false;
    Inf inf = new Inf("smallestBoundingDisc", 100);

    if (db && T.update())
      T.msg("smallestBoundingDisc" + T.show(a) + T.show(b) + T.show(c));
    do {
      EdDisc tmp = DiscUtil.smallestBoundingDisc(a, b);
      if (EdDisc.contains(tmp, c)) {
        ret = tmp;
        break;
      }
      tmp = DiscUtil.smallestBoundingDisc(a, c);
      if (EdDisc.contains(tmp, b)) {
        ret = tmp;
        break;
      }
      tmp = DiscUtil.smallestBoundingDisc(b, c);
      if (EdDisc.contains(tmp, a)) {
        ret = tmp;
        break;
      }

      Hyperbola h = MinMaxDiscBisector.S.getBisector(a, b);
      double t0 = -1000;
      double t1 = 1000;
      while (t0 < t1 - 1e-6) {
        inf.update();
        double tm = .67 * t0 + .33 * t1;
        double tn = .33 * t0 + .67 * t1;
        double dm = distFrom(h, tm, c, a);
        double dn = distFrom(h, tn, c, a);
        if (dm < dn) {
          if (db && T.update())
            T.msg("moving t1 left to " + Tools.f(tn));
          t1 = tn;
        } else {
          if (db && T.update())
            T.msg("moving t0 right to " + Tools.f(tm));
          t0 = tm;
        }

        if (db && T.update()) {
          FPoint2 origin = h.calcPoint(t0);
          double radius = distFrom(h, t0, c, a);
          ret = new EdDisc(origin, radius);
          T.msg("min bound disc" + T.show(ret) + T.show(h));
        }
      }
      FPoint2 origin = h.calcPoint(t0);
      double radius = distFrom(h, t0, c, a);
      ret = new EdDisc(origin, radius);
      if (db && T.update())
        T.msg("min bound disc" + T.show(ret));
    } while (false);
    return ret;
  }

  private static double distFrom(Hyperbola h, double t, EdDisc c, EdDisc a) {
    FPoint2 hpt = h.calcPoint(t);
    double rad = c.getOrigin().distance(hpt) + c.getRadius();
    double rad2 = a.getOrigin().distance(hpt) + a.getRadius();
    return Math.max(rad2, rad);
  }


}

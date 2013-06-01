package gvorn;

import java.awt.*;
import testbed.*;
import base.*;

public class Hyperb implements IPlaneCurve {

  /**
   * Construct hyperbolic arm
   * @param f1 first focus
   * @param f2 second focus
   * @param fd distance of nearest point on arm from first focus
   */
  public Hyperb(FPoint2 f1, FPoint2 f2, double fd) {

    final boolean db = true;

    c = FPoint2.distance(f2, f1) * .5;

    if (fd > c) {
      fd = 2 * c - fd;
      FPoint2 tmp = f1;
      f1 = f2;
      f2 = tmp;
    }
    if (fd <= 0 || fd > c)
      throw new IllegalArgumentException();

    a = c - fd;
    A = a * a;
    B = A / (c * c - A);

    FPoint2 origin = FPoint2.midPoint(f1, f2);
    double theta = MyMath.polarAngle(f1, f2);

    if (db && T.update())
      T.msg("Hyperb fd=" + Tools.f(fd) + " c=" + Tools.f(c) + " a="
          + Tools.f(a) + " theta=" + Tools.fa(theta));

    Matrix fromCenterInW = Matrix.getTranslate(origin, true);
    toE2 = Matrix.getRotate(-theta);
    Matrix.mult(toE2, fromCenterInW, toE2);
    toW2 = toE2.invert(null);

    /* Maple code:
    eq1 := s^2 - B*t^2 - A = 0;
    eq2 := s = c0*x + c1*y + c2;
    eq3 := t = c3*x + c4*y + c5;
    eq4 := subs({eq2,eq3},eq1);
    eq5:= expand(eq4);
     */
    coeff = new double[6];
    double c0 = toE2.get(0, 0);
    double c1 = toE2.get(0, 1);
    double c2 = toE2.get(0, 2);
    double c3 = toE2.get(1, 0);
    double c4 = toE2.get(1, 1);
    double c5 = toE2.get(1, 2);
    coeff[5] = c0 * c0 - B * c3 * c3;
    coeff[4] = 2 * (c0 * c1 - B * c3 * c4);
    coeff[3] = c1 * c1 - B * c4 * c4;
    coeff[2] = 2 * (c0 * c2 - B * c3 * c5);
    coeff[1] = 2 * (c1 * c2 - B * c4 * c5);
    coeff[0] = c2 * c2 - B * c5 * c5 - A;

    if (db) {
      pt(1.0);
    }

  }

  /**
  * Calculate a point on the parabola
  * @param t   parameter
  * @return point on parabola
  */
  public FPoint2 pt(double t) {
    double s2 = A + B * t * t;
    double s = -Math.sqrt(s2);
    return toW2.apply(s, t, null);
  }

  /**
   * Get parameter corresponding to point
   * @param pt  point
   * @return parameter corresponding to that point (if not on curve,
   *   a reasonable approximation)
   */
  public double parameterFor(FPoint2 pt) {
    FPoint2 p2 = toE2.apply(pt);
    return p2.y;
  }

  public double coeff(int n) {
    return coeff[n];
  }

  public int degree() {
    return 2;
  }

  public void render(double t0, double t1) {
    double step = Math.max(1, ((c - a) * (c - a) / c) * .3);
    t0 = MyMath.clamp(t0, -500.0, 500.0);
    t1 = MyMath.clamp(t1, -500.0, 500.0);

    double s0, s1;

    s0 = t0;
    s1 = t1;

    if (s0 > s1) {
      double tmp = s0;
      s0 = s1;
      s1 = tmp;
    }
    FPoint2 p0 = pt(s0), p1 = null;
    {
//      int count = 0;
      boolean first = true;
      for (double t = t0;; t += step) { //, count++) {
        boolean last = (t >= t1);
        if (last)
          t = t1;
        p1 = pt(t);

        if (!p1.isValid()) {
          if (last) {
            break;
          }
          continue;
        }

        if (!first)
          V.drawLine(p0, p1);
        if (last)
          break;
        p0 = p1;
        first = false;
      }
    }
  }

  /**
   * Evaluate point on hyperbola
   * @param x  x coordinate
   * @param y  y coordinate
   * @return double; should be zero for a point on the curve
   */
  public double eval(double x, double y) {
    return x * (coeff(5) * x + coeff(4) * y + coeff(2)) + y
        * (coeff(3) * y + coeff(1)) + coeff(0);
  }

  private double a, c;
  // a*a
  private double A;
  // a*a / (c*c - a*a)
  private double B;

  private Matrix toE2, toW2;
  private double[] coeff;
}

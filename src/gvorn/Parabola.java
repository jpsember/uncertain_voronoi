package gvorn;

import java.awt.*;
import testbed.*;
import base.*;

public class Parabola implements IPlaneCurve {

  public void render(double t0, double t1) {
    final boolean db = false;

    double step = Math.max(.1, f * .08);
    // step *= 4;
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
    FPoint2 p0 = this.pt(s0), p1 = null; //pt(s1);
    if (db)
      Streams.out.println(" p0=" + p0 + ", p1=" + p1);
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

        if (db) {
          System.out.println(" calcPt " + Tools.f(t) + " = " + p1.x + ","
              + p1.y);
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
  public static void main(String[] args) {
    db = true;

    FPoint2 s0 = new FPoint2(0, 1);
    FPoint2 s1 = new FPoint2(5, 8);
    FPoint2 pt = new FPoint2(7, 2);
    Parabola p = new Parabola(s0, s1, pt);
    System.out.println("f=" + p.f);

  }
  private static boolean db;

  /**
   * Constructor
   * @param s0 point on line
   * @param s1 point on line
   * @param pt focus point
   */
  public Parabola(FPoint2 s0, FPoint2 s1, FPoint2 pt) {

    LineEqn line = new LineEqn(s0, s1);
    double tOrigin = line.parameterFor(pt);
    FPoint2 orig = line.pt(tOrigin);

    f = FPoint2.distance(pt, orig);
    a = 1 / (2 * f);
    b = f / 2;

    if (db)
      System.out.println("constructing Parabola, s0=" + s0 + " s1=" + s1
          + " pt=" + pt);
    // calculate the angle of rotation of the parabola away
    // from the standard position.

    double theta = MyMath.polarAngle(orig, pt);

    Matrix fromCenterInW = Matrix.getTranslate(orig, true);
    Matrix rotToE = Matrix.getRotate(-theta);

    toCurve = rotToE;
    Matrix.mult(toCurve, fromCenterInW, toCurve);

    // calculate inverse
    toWorld = toCurve.invert(null);

    if (db)
      Streams.out.println("theta=" + Tools.fa(theta));
    if (db)
      Streams.out.println("toCurve=\n" + toCurve);

    coeff = new double[6];
    {
      // we have an implicit equation for the curve in 'curve space',
      // points defined by (s,t):  t^2 + f^2 - 2ft = 0;
      //
      // and we have a transform matrix for world->curve space; so
      // substitute for curve space parameters the expressions for the world
      // coordinates passed through this matrix.
      double A, B, C, D, E, F;

      A = toCurve.get(0, 0);
      B = toCurve.get(0, 1);
      C = toCurve.get(0, 2);
      D = toCurve.get(1, 0);
      E = toCurve.get(1, 1);
      F = toCurve.get(1, 2);

      // from Maple:
      /*
      eq1 := t^2 + f^2 - 2*f*s = 0;
      eq2 := s = A*x + B*y + C;
      eq3 := t = D*x + E*y + F;
      eq4 := subs({eq2,eq3},eq1);
      eq5:= expand(eq4);
       */
      coeff[5] = D * D; //U * A * A + D * D;
      coeff[4] = 2 * D * E; //2 * (U * A * B + D * E);
      coeff[3] = E * E; //U * B * B + E * E;
      coeff[2] = 2 * D * F - 2 * f * A; //2 * (U * A * C + D * F);
      coeff[1] = 2 * E * F - 2 * B * f; //2 * (U * B * C + E * F);
      coeff[0] = F * F + f * f - 2 * f * C; //U * C * C + F * F + V;

      //      hPoly = Conic.construct(U * A * A + D * D, 2 * (U * A * B + D * E), U * B
      //          * B + E * E, 2 * (U * A * C + D * F), 2 * (U * B * C + E * F), U * C
      //          * C + F * F + V);
    }

  }
  private double[] coeff;
  /**
  * Calculate a point on the parabola
  * @param t   parameter
  * @return point on parabola
  */
  public FPoint2 pt(double t) {
    return toWorld.apply(a * t * t + b, t, null);
  }

  /**
   * Get parameter corresponding to point
   * @param pt  point
   * @return parameter corresponding to that point (if not on curve,
   *   a reasonable approximation)
   */
  public double parameterFor(FPoint2 pt) {
    FPoint2 p2 = toCurve.apply(pt);
    return p2.y;
  }

  //  /**
  //   * Get IPlaneCurve representing parabola
  //   * @return
  //   */
  //  public IPlaneCurve getPlaneCurve() {
  //    if (hPoly == null) {
  //      {
  //        // we have an implicit equation for the curve in 'curve space';
  //        // and we have a transform matrix for world->curve space; so
  //        // substitute for curve space parameters the expressions for the world
  //        // coordinates passed through this matrix.
  //        double A, B, C, D, E, F;
  //
  //        A = toCurve.get(0, 0);
  //        B = toCurve.get(0, 1);
  //        C = toCurve.get(0, 2);
  //        D = toCurve.get(1, 0);
  //        E = toCurve.get(1, 1);
  //        F = toCurve.get(1, 2);
  //
  //        double U = -2 * f;
  //        double V = f * f;
  //
  //        // from Maple:
  //        /*
  //         eq1 := U*s^2 + t^2 + V = 0;
  //         eq2 := s = A*xp + B*yp + C;
  //         eq3 := t = D*xp + E*yp + F;
  //         eq4 := subs({eq2,eq3},eq1);
  //         eq5 := expand(eq4);
  //         */
  //
  //        hPoly = Conic.construct(U * A * A + D * D, 2 * (U * A * B + D * E), U
  //            * B * B + E * E, 2 * (U * A * C + D * F), 2 * (U * B * C + E * F),
  //            U * C * C + F * F + V);
  //      }
  //
  //    }
  //    return hPoly;
  //  }

  public double coeff(int n) {
    return coeff[n];
  }

  public int degree() {
    return 2;
  }

  // lazy-initialized polynomial representing the hyperbola equation
  //private IPlaneCurve hPoly;
  private double f;
  private double a, b;
  private Matrix toCurve, toWorld;
}

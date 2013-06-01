package gvorn;

import gvorn.GabrielOper.*;
import java.awt.*;
import testbed.*;
import base.*;

public class GabrielZoneFunction extends NewtonSinCos implements Renderable,
    Globals {

  public GabrielZoneFunction(EdDisc ca, EdDisc cb) {
    FPoint2 aOrig = ca.getOrigin();
    FPoint2 bOrig = cb.getOrigin();

    midpoint = FPoint2.midPoint(aOrig, bOrig);
    e = ca.getRadius();
    f = cb.getRadius();

    a = FPoint2.distance(aOrig, bOrig);
    abAngle = MyMath.polarAngle(aOrig, bOrig);

    // construct matrix to convert from pair space (disc a at origin, disc b on positive x-axis)
    // to view space

    Matrix m2 = Matrix.getRotate(abAngle);
    Matrix m3 = Matrix.getTranslate(aOrig, false);

    tfmDiscToView = Matrix.mult(m3, m2, null);
    tfmViewToDisc = tfmDiscToView.invert(null);

    double thetaMin = MyMath.normalizeAngle(findRoot(0));
    double thetaMax = MyMath.normalizeAngle(findRoot(Math.PI / 2));
    setDomain(thetaMin, thetaMax);
  }

  /**
    * Determine if a disc intersects this Gabriel zone
    * @param d disc
    * @param ptOnBoundary if not null, nearest boundary point stored here
    * @param lineEqn if not null, and length >= 1, tangent line stored here
    * 
    * @return true if intersection occurs
    */
  public boolean discIntersectsZone(EdDisc d, FPoint2 ptOnBoundary,
      LineEqn[] lineEqn) {

    double theta = nearestPointTo(d.getOrigin());

    FPoint2 bndPt0 = new FPoint2(calcXAt(theta), calcYAt(theta));

    FPoint2 bndPt = new FPoint2(bndPt0);
    if (queryBelowXAxis()) {
      bndPt.y = -bndPt.y;
    }
    bndPt = tfmDiscToView.apply(bndPt);

    //    FPoint2 bndPt = pointAt(theta, queryBelowXAxis());

    if (ptOnBoundary != null)
      ptOnBoundary.setLocation(bndPt);

    // calculate tangent line
    calcSinCos(theta);

    double dx = -a * cos * sin - bndPt0.y;
    double dy = -a * sin * sin + bndPt0.x;

    double viewAngle = MyMath.polarAngle(new FPoint2(dx, dy));

    if (queryBelowXAxis())
      viewAngle = Math.PI - viewAngle;
    viewAngle += abAngle;

    LineEqn ln = new LineEqn(bndPt, viewAngle);
    if (lineEqn != null && lineEqn.length > 0)
      lineEqn[0] = ln;

    if (true) {
    double distFromLine = -ln.signedDistanceFrom(d.getOrigin());
    return distFromLine <= d.getRadius();
    } else {
    // calculate distance to midpoint of disc pair
    double qDist = d.getOrigin().distance(midpoint) - d.getRadius();
    double bDist = bndPt.distance(midpoint);
    return qDist <= bDist;
    }
  }

  /**
   * Calculate nearest point on boundary of zone from a query point
   * @param pt location of point, in view space
   * @return theta of nearest point
   */
  private double nearestPointTo(FPoint2 pt) {
    FPoint2 pt0 = tfmViewToDisc.apply(pt);
    flipped = pt0.y < 0;
    pt0.y = Math.abs(pt0.y);

    Newton fn = new GabrielFn2(this, pt0);

    double th0;
    double bestTheta = 0;
    double bestDist = -1;

    for (th0 = 0; th0 <= Math.PI / 2; th0 += Math.PI / 8) {

      double theta = fn.findRoot(th0);
      FPoint2 pt2 = pointAt(theta, flipped);
      double dist = FPoint2.distance(pt, pt2);
      if (bestDist < 0 || bestDist > dist) {
        bestDist = dist;
        bestTheta = theta;
      }
    }
    return bestTheta;
  }

  /**
   * Get angle that axis from a to b makes with x-axis
   * @return
   */
  public double getAxisAngle() {
    return abAngle;
  }

  @Override
  public double eval(double theta) {
    return calcYAt(theta);
  }

  private double calcYAt(double theta) {
    calcSinCos(theta);

    return a * cos * sin + f * sin + e * cos;
  }

  private double calcXAt(double theta) {
    calcSinCos(theta);
    return (a * cos + f) * cos - e * sin;
  }

  /**
   * Calculate boundary point for a parameter value
   * @param theta parameter
   * @param flip true to return lower boundary point
   * @return boundary point, in view space
   */
  public FPoint2 pointAt(double theta, boolean flip) {
    FPoint2 pt = new FPoint2(calcXAt(theta), eval(theta));
    if (flip)
      pt.y = -pt.y;

    return tfmDiscToView.apply(pt);
  }

  @Override
  public double evalDerivative(double theta) {
    calcSinCos(theta);

    return a * (cos * cos - sin * sin) + f * cos - e * sin;
  }

  public void sample(int res) {
    double step = (Math.PI ) / res;
    double r[] = new double[Math.min(res * 4, 720)];
    double step2 = Math.PI / r.length;

    FPoint2 ao = tfmDiscToView.apply(0, 0, null);
    FPoint2 bo = tfmDiscToView.apply(a, 0, null);

    FPoint2 aMin = tfmDiscToView.apply(-e, 0, null);
    FPoint2 bMax = tfmDiscToView.apply(a + f, 0, null);

    FPoint2 boundingOrigin = FPoint2.midPoint(aMin, bMax);
    //)tfmDiscToView.apply(a/2,0,null);
    double xAxis = MyMath.polarAngle(ao, bo);

    for (double a = 0; a < 2 * Math.PI; a += step) {
      FPoint2 apt = MyMath.ptOnCircle(ao, a, e);

      for (double b = 0; b < 2 * Math.PI; b += step) {
        FPoint2 bpt = MyMath.ptOnCircle(bo, b, f);

        FPoint2 gabrielOrigin = FPoint2.midPoint(apt, bpt);

        double gabrielRadius = FPoint2.distance(gabrielOrigin, apt);

        for (double c = 0; c < 2 * Math.PI; c += step) {

          FPoint2 ptOnGabriel = MyMath.ptOnCircle(gabrielOrigin, c,
              gabrielRadius);

          double distFromOrigin = FPoint2.distanceSquared(ptOnGabriel,
              boundingOrigin);
          double ang = MyMath.polarAngle(boundingOrigin, ptOnGabriel);
          ang = MyMath.normalizeAnglePositive(ang - xAxis);
          int j = (int) Math.floor(ang / step2);
          if (j < 0 || j >= r.length)
            continue;
          if (distFromOrigin > r[j])
            r[j] = distFromOrigin;
        }
      }
    }
    int i = 0;
    FPoint2 prev = null;
    FPoint2 pt;
    for (double c = xAxis; i < r.length; c += step2, i++) {
      if (r[i] == 0)
        continue;
      pt = MyMath.ptOnCircle(boundingOrigin, c, Math.sqrt(r[i]));
      if (prev != null)
        V.drawLine(prev, pt);
      prev = pt;
    }
    prev = null;
    i = 0;
    for (double c = 0; i < r.length; c += step2, i++) {
      if (r[i] == 0)
        continue;
      pt = MyMath.ptOnCircle(boundingOrigin, xAxis - c, Math.sqrt(r[i]));
      if (prev != null)
        V.drawLine(prev, pt);
      prev = pt;
    }
  }

  /**
   * Determine if last argument for nearestPointTo() was below
   * the disc-space x-axis
   * 
   * @return true if so
   */
  private boolean queryBelowXAxis() {
    return flipped;
  }

  public void setUpperOnly(boolean b) {
    upperOnly = b;
  }

  public void render(Color c, int stroke, int markType) {

    V.pushColor(c, MyColor.cLIGHTGRAY);
    V.pushStroke(stroke, STRK_NORMAL);

    double startTheta = xMin();
    double endTheta = xMax();

    for (int pass = 0; pass < 2; pass++) {
      if (upperOnly) {
        if (pass != 0)
          break;
        startTheta = 0;
        endTheta = Math.PI * 2;
      }

      FPoint2 pt0 = null;
      for (double th = startTheta; th < endTheta; th += MyMath.radians(1)) {

        double cos = Math.cos(th);
        double sin = Math.sin(th);

        double u = a * cos * cos;
        double v = a * cos * sin;
        double x = u + f * cos - e * sin;
        double y = v + f * sin + e * cos;
        if (pass == 1)
          y = -y;

        FPoint2 npt = new FPoint2(x, y);
        FPoint2 pt = tfmDiscToView.apply(npt, null);

        if (pt0 != null)
          V.drawLine(pt0, pt);

        pt0 = pt;
      }
    }
    V.pop(2);

  }

  private Matrix tfmDiscToView;
  private Matrix tfmViewToDisc;
  private double a, e, f;
  private double abAngle;

  // midpoint of the two disc origins
  private FPoint2 midpoint;

  private boolean flipped;
  private boolean upperOnly;

  /**
   * Newton function for finding nearest point on boundary
   * of Gabriel zone to a point
   */
  private static class GabrielFn2 extends NewtonSinCos {
    public GabrielFn2(GabrielZoneFunction fn, FPoint2 queryPt) {
      this.w = queryPt.x;
      this.z = queryPt.y;
      double a = fn.a;
      double e = fn.e;
      double f = fn.f;
      A = a * (2 * w - a);
      B = 2 * a * z;
      C = e * w - f * z;
      D = f * w + e * z - a * f;
      E = -a * z;
      this.setDomain(fn.xMin(), fn.xMax());
    }

    @Override
    public double evalDerivative(double theta) {
      calcSinCos(theta);
      return A * (2 * cos * cos - 1) + 2 * B * cos * sin - C * sin + D * cos;
    }
    private double w, z;
    @Override
    public double eval(double theta) {
      calcSinCos(theta);
      return cos * sin * A + sin * sin * B + cos * C + sin * D + E;
    }
    private double A, B, C, D, E;
  }

}

abstract class NewtonSinCos extends Newton {
  private double lastTheta = 0;
  protected void calcSinCos(double theta) {
    if (theta != lastTheta) {
      lastTheta = theta;
      sin = Math.sin(theta);
      cos = Math.cos(theta);
    }
  }
  protected double sin = 0, cos = 1.0;
}

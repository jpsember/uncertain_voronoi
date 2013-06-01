package gvorn;

import base.*;

public abstract class Newton {

  /**
   * Evaluate function 
   * @param x function parameter
   * @return value of function at x
   */
  public abstract double eval(double x);
  /**
   * Evaluate derivative of function 
   * @param x function parameter
   * @return value of derivative of function at x
   */
  public abstract double evalDerivative(double x);

  public void setDomain(double xMin, double xMax) {
    this.xMin = xMin;
    this.xMax = xMax;
    domainBounded = true;
  }

  public double clamp(double x) {
    if (domainBounded)
      x = MyMath.clamp(x, xMin, xMax);
    return x;
  }
  public double findRoot(double x0) {
    return findRoot(x0, false);
  }
  public double findRoot(double x0, boolean db) {

    x0 = clamp(x0);
    int iter = 0;
    final double EPS = 1e-5;
    final double EPS2 = 1e-10;
    double x = x0;
    while (true) {
      if (iter++ > 20)
        break;

      double f = eval(x);
      double fPrime = evalDerivative(x);
      if (Math.abs(fPrime) < EPS)
        break;

      double xNext = x - f / fPrime;
      if (domainBounded)
        xNext = MyMath.clamp(xNext,xMin,xMax);
      
      
      if (db)
        Streams.out.println("Iter=" + Tools.f(iter, 2) + " x="+Tools.fa(x)+" x2="+Tools.fa(xNext)
            + " diff=" + Tools.fa(Math.abs(xNext - x)));

      if (Math.abs(xNext - x) < EPS2)
        break;
      x = xNext;
    }
    if (db)
      Streams.out.println();

    return x;
  }
public double xMin() {return xMin;}public double xMax() {return xMax;}

  private double xMin, xMax;
  private boolean domainBounded;
}

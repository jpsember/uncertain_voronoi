package gvorn;

import base.*;

/**
 * Polynomial of polynomials class
 *
 * Represents a polynomial in x, where coefficients are polynomials in y
 */
public class PolynOfPolyn {


  public PolynOfPolyn() {
    polynomials = new DArray();
    polynomials.add(Polyn.ZERO);
  }

  /**
   * Initialize a coefficient of the polynomial
   * @param coeff : index of coefficient (power of x)
   * @param p Polyn : polynomial in y
   */
  public void set(int coeff, Polyn p) {
    // add zeros if necessary
    while (degree() + 1 < coeff) {
      polynomials.growSet(degree() + 1, Polyn.ZERO);
    }
    polynomials.growSet(coeff, p);

    // trim zeros from upper
    while (degree() > 0) {
      Polyn p2 = poly(degree());
      if (!p2.isZero()) {
        break;
      }
      polynomials.pop();
    }

    // recalculate magnitude
    calcMagnitude();
  }

  private void calcMagnitude() {
    double m = 0;
    for (int i = 0; i < degree() + 1; i++) {
      m = Math.max(m, Math.abs(poly(i).magnitude()));
    }
    magnitude = m;
  }

  private Polyn poly(int power) {
    return (Polyn) polynomials.get(power);
  }

  private int degree() {
    return polynomials.size() - 1;
  }

  /**
   * Get string describing object
   * @return String
   */
  public String toString() {
    StringBuilder sb = new StringBuilder();
    //sb.append("PPoly: ");

    boolean any = false;

    int d = degree();
    for (int i = d; i >= 0; i--) {
      Polyn p = poly(i);
      if (p.isZero()) {
        if (any || i > 0) {
          continue;
        }
      }
      any = true;

      String s = toString(p, "y");
      if (i < d) {
        sb.append(" +");
      }
      sb.append('(');
      sb.append(s);
      sb.append(')');
      if (i > 0) {
        sb.append("x");
      }
      if (i > 1) {
        sb.append(i);
      }
    }
    return sb.toString();
  }

  public static String toString(Polyn p) {
    return toString(p, "y");
  }

  public static String toString(Polyn p, String var) {
    StringBuilder sb = new StringBuilder();
    int d = p.degree();
    boolean any = false;
    for (int i = d; i >= 0; i--) {
      double v = p.c(i);
      if (Polyn.isZero(v)) {
        if (any || i > 0) {
          continue;
        }
      }
      any = true;
      String s = Tools.f(v, 4,3).trim(); // Tools.dblStr(v, false, 4, 3).trim();
      if (v >= 0 && i < d) {
        sb.append('+');
      }
      sb.append(s);
      if (i > 0) {
        //sb.append(' ');
        sb.append(var);
        if (i > 1) {
          sb.append(i);
        }
      }
    }
    return sb.toString();
  }

  public static void main(String[] args) {

    Polyn q2 = new Polyn(2, -3, 1), q1 = new Polyn(3, 4), q0 = new Polyn(2, -1);

    if (false) {
      PolynOfPolyn p = new PolynOfPolyn();
      p.set(0, q0);
      p.set(1, q1);
      p.set(2, q2);

      System.out.println("PPolyn p = \n" + p);

      Polyn prod = Polyn.product(q0, q1), prod2 = Polyn.product(q0, q2);
      System.out.println("q0xq1= " + toString(prod));
      System.out.println("q0xq2= " + toString(prod2));
    }

    if (false) {
      PolynOfPolyn p1 = new PolynOfPolyn();
      p1.set(1, new Polyn(3, 0));
      p1.set(0, new Polyn(4));
      PolynOfPolyn p2 = new PolynOfPolyn();
      p2.set(1, new Polyn(2));
      p2.set(0, new Polyn(-9, 0));

      System.out.println("p1=" + p1 + "\np2=" + p2);
      System.out.println("reduce(p1,p2) = " + PolynOfPolyn.reduce(p1, p2));
    }

    {
      Polyn a2 = new Polyn(1), a1 = new Polyn(0), a0 = new Polyn(1, -4), b2 = new Polyn(
          1, 1), b1 = new Polyn(-1), b0 = new Polyn(1, 0, -12);

      PolynOfPolyn A = new PolynOfPolyn();
      A.set(2, a2);
      A.set(1, a1);
      A.set(0, a0);
      PolynOfPolyn B = new PolynOfPolyn();
      B.set(2, b2);
      B.set(1, b1);
      B.set(0, b0);

      System.out.println("A= " + A + "\nB= " + B);
      Polyn s = reduceToSingleVar(A, B);
      System.out.println("Reduced to single var:\n" + toString(s));
      Polyn px = A.solveForY(3);
      Polyn qx = B.solveForY(3);
      System.out.println("solve for y=3: \n A=" + px + "\n B=" + qx);

      DArray ans = new DArray();
      solve(A, B, ans);
      System.out.println("solve= " + ans);
    }

  }

  private void shiftLeft(int n) {
    int d = degree();
    for (int i = d; i >= 0; i--)
      polynomials.growSet(i + n, polynomials.get(i));

    for (int i = 0; i < n; i++)
      polynomials.growSet(i,Polyn.ZERO);
  }

  public PolynOfPolyn(PolynOfPolyn src) {
    polynomials = (DArray) src.polynomials.clone();
    magnitude = src.magnitude;
  }

  private static PolynOfPolyn reduce(PolynOfPolyn a, PolynOfPolyn b) {
    boolean dbg = false;
    if (dbg) {
      System.out.println("reduce \n a=" + a + "\n b=" + b);
    }
    if (a.degree() < b.degree()) {
      return reduce(b, a);
    }
    int da = a.degree(), db = b.degree();
    Tools.ASSERT(da > 0, "Attempt to reduce degree 0 polynomials");

    if (dbg) {
      System.out.println(" da=" + da + " db=" + db);
    }
    if (db < da) {
      b = new PolynOfPolyn(b);
      b.shiftLeft(da - db);
      if (dbg) {
        System.out.println(" shifted b left, now\n b=" + b);
      }
      db = da;
    }

    Polyn ai = a.poly(da), bi = b.poly(da);

    // construct a * bi - b * ai
    PolynOfPolyn sum = new PolynOfPolyn();
    for (int i = 0; i < da; i++) {
      Polyn ap = Polyn.product(a.poly(i), bi);
      Polyn bp = Polyn.product(b.poly(i), ai);
      Polyn cp = Polyn.add(ap, 1, bp, -1);
      sum.set(i, cp);
    }

    // scale this polynomial so maximum magnitude is reasonable.
    sum = sum.normalize(1000);

    if (dbg) {
      System.out.println(" returning sum\n " + sum);
    }
    return sum;
  }

  private PolynOfPolyn normalize(double mag) {
    if (magnitude == 0) {
      return new PolynOfPolyn(this);
    }
    return scale(mag / magnitude);
  }
  private PolynOfPolyn scale(double s) {
    PolynOfPolyn p = new PolynOfPolyn();
    for (int i = 0; i <= degree(); i++) {
      p.set(i, poly(i).scaleBy(s));
    }
    return p;
  }

  public static void solve(PolynOfPolyn a, PolynOfPolyn b, DArray solutionPts) {
    boolean db = false;

    solutionPts.clear();

    if (db) {
      System.out.println("solve\na=" + a + "\nb=" + b);
    }
    if (db) {
      System.out.println(" normalized\na=" + a + "\nb=" + b);
    }
    Polyn r = reduceToSingleVar(a, b);

    if (db) {
      System.out.println(" reduced to single variable:\n" + r);
    }

    DArray roots = new DArray();
    DArray roots2 = new DArray();

    r.solve(roots);
    if (db) {
      System.out.println(" verifying solutions= " + roots);
    }

    for (int i = 0; i < roots.size(); i++) {
      double y = roots.getDouble(i);
      if (db) {
        System.out.println(" candidate = " + y);
      }
      Polyn ax = a.solveForY(y);
      Polyn bx = b.solveForY(y);
      if (db) {
        System.out.println(" ax= " + ax + "\n bx= " + bx);
      }

      ax.solve(roots2);
      if (db) {
        System.out.println(" ax solved produced roots= " + roots2);
      }
      for (int j = 0; j < roots2.size(); j++) {
        double x = roots2.getDouble(j);

        // is x,y a solution to both original polynomials?

        double v = ax.eval(x);
        if (db) {
          System.out.println(" ax evaluated at x=" + x + " = " + v);
        }

        if (!nearZero(v)) {
          continue;
        }
        double v2 = bx.eval(x);
        if (db) {
          System.out.println(" bx evaluated at x=" + x + " = " + v2);
        }
        if (!nearZero(v2)) {
          continue;
        }

        // is this point already in the solution array?
        // I don't know if this is necessary...
        boolean skip = false;
        for (int k = 0; k < solutionPts.size(); k++) {
          FPoint2 t = solutionPts.getFPoint2(k);
          if (nearZero(t.x - x) && nearZero(t.y - y)) {
            skip = true;
            break;
          }
        }
        if (!skip) {
          solutionPts.add(new FPoint2(x, y));
        }
      }
    }
    if (db) {
      System.out.println(" solutionPts= " + solutionPts + "\n\n");
    }
  }

  private static Polyn reduceToSingleVar(PolynOfPolyn a, PolynOfPolyn b) {
    boolean db = false;

    if (db)
      System.out.println("reduceToSingleVar\n a=" + a + "\n b=" + b);

    Polyn ret = null;
    do {

      if (a.degree() < b.degree()) {
        ret = reduceToSingleVar(b, a);
        break;
      }

      // if degree of both is zero, we can return the degree-zero polynomial
      // for one of them, which is a polynomial in y
      int d = Math.max(a.degree(), b.degree());
      if (d == 0) {
        ret = a.poly(0);
        break;
      }

      if (db)
        System.out.println(" maximum degree = " + d);
      PolynOfPolyn C = reduce(a, b);
      if (db)
        System.out.println(" reduced to\n C=" + C);
      PolynOfPolyn D = reduce(a, C);
      if (db)
        System.out.println(" reduced to\n D=" + D);
      ret = reduceToSingleVar(C, D);
    } while (false);
    return ret;
  }

  /**
   * Evaluate for a particular value of y
   * @param y
   * @return polynomial in x
   */
  public Polyn solveForY(double y) {
    double[] dwork = new double[1+degree()];
    for (int i = 0; i <= degree(); i++) {
      Polyn p = poly(i);
      dwork[i] = p.eval(y);
    }
    return new Polyn(dwork);
  }

  private static boolean nearZero(double v) {
    return Math.abs(v) < 1e-5;
  }

  // polynomials (Polyn) for each entry
  private DArray polynomials;

  // largest magnitude of any component Polyn
  private double magnitude;
}

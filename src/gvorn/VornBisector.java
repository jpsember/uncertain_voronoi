package gvorn;

import java.awt.*;
import testbed.*;
import base.*;

public abstract class VornBisector implements Renderable {
public int nPieces() {return 1;}
public abstract IPlaneCurve curve(int i);
  //  /**
  //   * Determine if bisector is valid.  An invalid bisector is
  //   * one that does not exist or is undefined; for instance, between two
  //   * regions, one of which is interior to the other
  //   * @return true if valid
  //   * @deprecated
  //   */
  //  public boolean isValid() {
  //    return true;
  //  }

  /**
   * Min and max parameter values for bisector
   */
  public static final double TMAX = 1e4, TMIN = -TMAX;

  /**
   * Get one of the sites associated with the bisector
   * @param index site number (0..1)
   * @return site
   */
  public VornSite getSite(int index) {
    return index == 0 ? sa : sb;
  }

  /**
   * Determine if a point is on the bisector (ignoring clip regions),
   * by seeing if distance from point to each site is the same
   * @param pt query point
   * @return true if point lies on the bisector
   */
  public boolean onBisector(FPoint2 pt) {
    final boolean db = false;
    //final double tolerance = 1e-4;

    double saDist = sa.distanceFrom(pt);
    double sbDist = sb.distanceFrom(pt);

    //double scale = Math.max(1.0, Math.max(saDist, sbDist));

    double d1 = Math.abs(saDist - sbDist);
    boolean ret;

    //    if (true) {
    ret = d1 < 1e-2;
    //    } else
    //      ret = (d1 / scale) < tolerance;

    if (db && T.update())
      T.msg("onBisector " + T.show(pt) + T.show(this) + " distA="
          + Tools.f(saDist) + " distB=" + Tools.f(sbDist) + " diff="
          + Tools.f(d1) + " returning " + ret);
    return ret;
  }

  /**
   * Constructor
   * @param sa first site
   * @param sb second site
   */
  public VornBisector(VornSite sa, VornSite sb) {
    this.sa = sa;
    this.sb = sb;
    comp = new DArray();
    comp.add(new Comp(TMIN, TMAX));
  }

  /**
   * Get the number of connected components.
   * It may get fragmented into more than one component by clipping operations.
   * @return
   */
  public int nComponents() {
    return comp.size();
  }
  private Comp comp(int cNum) {
    return (Comp) comp.get(cNum);
  }
  public void setClipRange(double tMin, double tMax) {
    comp.clear();
    comp.add(new Comp(tMin, tMax));
  }
  public double cStart(int cNum) {
    return comp(cNum).start;
  }
  public double cEnd(int cNum) {
    return comp(cNum).end;
  }

  public void clipAll() {
    comp.clear();
  }

  //  /**
  //   * Get the IPlaneCurve representing this bisector
  //   * @return IPlaneCurve
  //   * @deprecated
  //   * 
  //   */
  //  public final  IPlaneCurve getPlaneCurve( ){return null;}
  public abstract FPoint2 pointAt(double t);

  public abstract double parameterFor(FPoint2 pt);

  public abstract void render(Color c, int stroke, int markType);

  /**
   * Render bisector given additional IPlaneCurve parameter
   * @param curve
   * @param color
   * @param stroke
   * @param markType
   */
  protected void render(IPlaneCurve curve, Color color, int stroke, int markType) {
    V.pushStroke(stroke, Globals.STRK_NORMAL);
    V.pushColor(color, MyColor.cBLUE);

    for (int seg = 0; seg < nComponents(); seg++) {
      curve.render(cStart(seg), cEnd(seg));
    }
    V.pop(2);
  }

  /**
   * Class representing connected component of bisector
   */
  private static class Comp {
    public Comp(double s, double e) {
      this.start = s;
      this.end = e;
    }
    double start, end;
    public String toString() {
      StringBuilder sb = new StringBuilder();
      sb.append("[");
      sb.append(Tools.f(start));
      sb.append(" ... ");
      sb.append(Tools.f(end));
      sb.append("]");
      return sb.toString();
    }
  }

  //  /**
  //   * Calculate all points of intersection of this
  //   * bisector with another
  //   * @param b other bisector
  //   * @return array of parameters representing intersection points; (tThis, tOther) pairs
  //   */
  //  public abstract DArray intersectWith(VornBisector b);

  /**
   * Clip out portion of bisector 
   * @param tj0
   * @param tj
   */
  public void clipOut(double c0, double c1) {
    final boolean db = false;

    if (db && T.update())
      T.msg("clipOut " + T.show(this) + T.show(pointAt(c0))
          + T.show(pointAt(c1)) + Tools.f(c0) + " ... " + Tools.f(c1));

    DArray nc = new DArray();
    for (int i = 0; i < nComponents(); i++) {
      Comp c = comp(i);
      if (false && db && T.update())
        T.msg(" component: " + c);
      if (c.start > c1 || c.end < c0)
        nc.add(c);
      else if (c0 <= c.start) {
        if (c1 < c.end)
          nc.add(new Comp(c1, c.end));
      } else {
        nc.add(new Comp(c.start, c0));
        if (c1 < c.end)
          nc.add(new Comp(c1, c.end));
      }
    }
    comp = nc;
    if (false && db && T.update())
      T.msg(" new components:\n" + comp.toString(true));

  }
  public String componentsStr() {
    StringBuilder sb = new StringBuilder();
    sb.append("components: [\n");
    for (int i = 0; i < nComponents(); i++) {
      sb.append(" " + cStart(i) + " ... " + cEnd(i) + "\n");
    }
    sb.append("]");
    return sb.toString();
  }

  private VornSite sa, sb;
  private DArray comp;

  public VornSite otherSite(VornSite site) {
    return (site == sa) ? sb : sa;
  }

}

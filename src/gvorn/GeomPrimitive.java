package gvorn;

import base.*;
public interface GeomPrimitive {

  /**
   * Calculate a point on the primitive
   * @param t : parameter
   * @param dest : where to store the calculated point
   */
  public FPoint2 calcPoint(double t, FPoint2 dest);

  /**
   * Calculate point on the primitive
   * @param t : parameter
   * @return point on curve
   */
  public FPoint2 calcPoint(double t);

  /**
   * Transform a point from world space to curve space
   * @param pt : point in world space
   * @return point in curve space
   */
  public FPoint2 toCurveSpace(FPoint2 pt, FPoint2 out);

  /**
   * Transform a point from world space to curve space
   * @param in : point in curve space
   * @param out : where to store point in world space
   */
  public FPoint2 toWorldSpace(FPoint2 in, FPoint2 out);

  /**
   * Given a point in world space, calculate curve parameter most
   * closely associated with that point
   *
   * @param pt : point in world space
   * @return parameter
   */
  public double calcParameter(FPoint2 pt);

}

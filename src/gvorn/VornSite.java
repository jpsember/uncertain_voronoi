package gvorn;

import base.*;
import testbed.*;

/**
 * Abstract class representing Voronoi sites
 */
public abstract class VornSite implements Renderable {

  public void setId(int id) {
    this.id = id;
  }
  public int id() {
    return id;
  }
  private int id;

  /**
   * Calculate distance of a point from the site.
   * If point is inside site, behaviour is undetermined.
   * @param pt
   * @return distace of point from site
   */
  public abstract double distanceFrom(FPoint2 pt);

  /**
   * Construct Voronoi bisector between this and another site
   * @param otherSite the other site
   * @return array of bisectors associated with the two sites
   *   (for point and disc sites, this array contains a single bisector;
   *    for segments and polygons, this may contain more than one bisector)
   */
  public abstract DArray buildBisector(VornSite otherSite);

}

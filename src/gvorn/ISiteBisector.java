package gvorn;

import testbed.*;
import base.*;

public interface ISiteBisector {
  /**
   * Construct bisector for two sites
   * @param a
   * @param b
   * @return
   */
  public Hyperbola getBisector(EdPoint a, EdPoint b);

  /**
   * Clip two hyperbolas based on # intersections
   * @param hij : first hyperbola
   * @param hik : second hyperbola 
   * @param commonSite : the site associated with their common focus
   * @param iPoints : the points of intersection of the sites
   */
  public void clipPair(Hyperbola hij, Hyperbola hik, EdPoint commonSite, DArray iPoints);

  /**
   * Determine if cell is trivially empty
   */
  public boolean cellIsEmpty(EdPoint[] sites, EdPoint site);
}

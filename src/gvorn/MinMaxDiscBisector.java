package gvorn;

import testbed.*;
import base.*;

public class MinMaxDiscBisector extends StandardDiscBisector {

  public Hyperbola getBisector(EdPoint a, EdPoint b) {
    Hyperbola h = null;
    try {
      EdDisc c1 = (EdDisc) a, c2 = (EdDisc) b;
      if (c1.getRadius() == c2.getRadius()) {
        h = PointBisector.S.getBisector(c1, c2);
      } else {
        double dist = FPoint2.distance(c1.getOrigin(), c2.getOrigin());
        double adj = .5 * dist - (-c2.getRadius() + c1.getRadius()) * .5;
        if (adj > 0 && adj < dist) {
          h = new Hyperbola(c1.getOrigin(), c2.getOrigin(), adj);
          h.setLabel(c1.getLabel() + "/" + c2.getLabel());
          h.setData(c1, c2);
        }
      }
    } catch (FPError e) {
      Tools.warn("FPError caught");
    }
    return h;
  }

  public boolean cellIsEmpty(EdPoint[] sites, EdPoint site) {
    boolean ret = false;

    // if this site contains others, it is empty
    for (int i = 0; i < sites.length; i++) {
      EdDisc ci = (EdDisc) sites[i];
      if (ci == site)
        continue;
      if (EdDisc.contains((EdDisc) site, ci)) {
        ret = true;
        break;
      }
    }

    return ret;
  }

  public static final ISiteBisector S = new MinMaxDiscBisector();
}

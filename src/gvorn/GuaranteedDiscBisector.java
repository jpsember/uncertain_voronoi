package gvorn;

import testbed.*;
import base.*;

public class GuaranteedDiscBisector extends StandardDiscBisector {

  public Hyperbola getBisector(EdPoint a, EdPoint b) {
    EdDisc c1 = (EdDisc) a, c2 = (EdDisc) b;
    Hyperbola h = null;
    do {
      try {
        do {
          if (c1.getRadius() == 0 && c2.getRadius() == 0) {
            h = PointBisector.S.getBisector(c1, c2);
            break;
          }
          double dist = FPoint2.distance(c1.getOrigin(), c2.getOrigin());

          double adj = .5 * ( -c2.getRadius() - c1.getRadius() + dist);
          if (adj >= dist)
            break;
          h = new Hyperbola(c1.getOrigin(), c2.getOrigin(), //
              adj);
          h.setLabel(c1.getLabel() + "/" + c2.getLabel());
          h.setData(c1, c2);
        } while (false);
      } catch (FPError e) {
        Tools.warn("FPError caught");
        h = null;
      }

    } while (false);
    return h;
  }

  public boolean cellIsEmpty(EdPoint[] sites, EdPoint index) {
    boolean ret = false;
    EdDisc di = (EdDisc) index;

    for (int j = 0; j < sites.length; j++) {
      EdDisc dj = (EdDisc) sites[j];
      if (dj != di && EdDisc.overlap(di, dj)) {
        ret = true;
        break;
      }
    }
    return ret;
  }
  public static final ISiteBisector S = new GuaranteedDiscBisector();
}

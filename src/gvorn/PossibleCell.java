package gvorn;

import testbed.*;
import base.*;

public class PossibleCell {

  private static class PossibleDiscBisector2 extends GuaranteedDiscBisector {
    public Hyperbola getBisector(EdPoint a, EdPoint b) {
      EdDisc c1 = (EdDisc) a, c2 = (EdDisc) b;
      Hyperbola h = null;
      do {
        if (EdDisc.overlap(c1, c2))
          break;

        try {
          do {
            if (c1.getRadius() == 0 && c2.getRadius() == 0) {
              h = PointBisector.S.getBisector(c1, c2);
              break;
            }

            double dist = FPoint2.distance(c1.getOrigin(), c2.getOrigin());

            double adj = .5 * (dist + (c2.getRadius() + c1.getRadius()));
            if (adj < 0)
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

    public void clipPair(Hyperbola hij, Hyperbola hik, EdPoint commonSite,
        DArray iPts) {

      final boolean db = false;
      if (db)
        Streams.out.println("PossibleDiscBisector2, clip " + hij + "," + hik
            + " iPts=" + iPts);

      switch (iPts.size()) {
      default:
        super.clipPair(hij, hik, commonSite, iPts);
        break;
      case 0:
        {
          // If arms are nested, clip inner one
          if (!hij.isEmpty() && hik.testPoint(hij.getSamplePoint()) < 0) {
            hij.clipAll();
            break;
          }
          if (!hik.isEmpty() && hij.testPoint(hik.getSamplePoint()) < 0) {
            hik.clipAll();
            break;
          }
        }

        // If one disc contains the other, clip larger disc's site.
        EdDisc dj = (EdDisc) hij.getOtherData(commonSite);
        EdDisc dk = (EdDisc) hik.getOtherData(commonSite);
        if (!hij.isEmpty() && EdDisc.contains(dj, dk)) {
          hij.clipAll();
        } else if (!hik.isEmpty() && EdDisc.contains(dk, dj)) {
          hik.clipAll();
        }
        break;
      }
    }
  }

  private static final ISiteBisector possibleDiscBisector2 = new PossibleDiscBisector2();

  /**
   * @param siteI
   * @return
   */
  public static VornGraph buildPossibleCell(int siteI) {
    Hyperbola[] h = VornUtil.buildCell(Main.getDiscs(), siteI,
        possibleDiscBisector2);
    return new VornGraph(h);
  }
}

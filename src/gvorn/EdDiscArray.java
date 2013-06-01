package gvorn;

import java.awt.*;
import testbed.*;
import base.*;

public class EdDiscArray extends DArray {
  public EdDiscArray(EdDisc[] discs) {
    super(discs);
  }
  public EdDiscArray() {
  }

  public EdDiscArray(EdDisc[] discs, int start, int count) {
    super();
    for (int i = start; i < start + count; i++)
      add(discs[i]);
  }

  public EdDisc disc(int i) {
    return (EdDisc) get (i);
  }

  /**
   * Find the disc to begin the Jarvis march with.  This
   * is the disc with the lowest uppermost point.
   * @return EdDisc, or null if no discs found
   */
  public EdDisc findLowestDisc() {
    final boolean db = false;
    EdDisc lowestDisc = null;
    double lwst = 0;
    for (int i = 0; i < size(); i++) {
      EdDisc d = disc(i);
      double y = d.getOrigin().y + d.getRadius();
      if (db && T.update())
        T.msg("finding lowest disc " + T.show(d, Color.RED));
      if (lowestDisc == null || y < lwst) {
        lowestDisc = d;
        lwst = y;
      }
    }
    return lowestDisc;
  }
  public EdDisc peekDisc(int distFromLast) {
    return (EdDisc) peek(distFromLast);
  }

}

package gvorn;

import java.awt.*;
import java.util.*;
import testbed.*;
import base.*;

public class VornUtil {

  private static String cellLabel(EdPoint[] sites, DArray cell) {
    StringBuilder sb = new StringBuilder();
    for (int i = 0; i < cell.size(); i++) {
      EdPoint v = sites[cell.getInt(i)];
      //      sb.append("#"+cell.getInt(i)+":");
      sb.append(v.getLabel());
    }
    return sb.toString();
  }

  /**
   * Get list of cells
   * @param n : number of regions
   * @param k : number of regions per cell
   * @return DArray of DArray of ints, labels for each cell
   */
  private static DArray getSites(int n, int k) {

    final boolean db = false;
    if (db)
      Streams.out.println("getSites n=" + n + " k=" + k);

    DArray s = new DArray();

    if (n >= k) {
      DArray c = new DArray();

      outer: while (true) {
        if (s.size() == 0) {
          for (int j = 0; j < k; j++)
            c.addInt(j);
        } else {
          int d = k - 1;
          while (true) {
            if (d < 0)
              break outer;
            if (db)
              Streams.out.println(" d=" + d);

            int h = c.getInt(d);
            if (h >= n - (k - d)) {
              d--;
            } else {
              if (db)
                Streams.out.println(" setting to " + (h + 1));

              c.setInt(d, ++h);
              while (d + 1 < k) {
                c.setInt(++d, ++h);
              }
              break;
            }
          }
        }
        if (db)
          Streams.out.println(" adding " + c);

        s.add(c.clone());
      }
    }

    return s;
  }

  /**
   * Build a single cell of a Voronoi diagram
   * @param sites : sites
   * @param siteI : site to build for (index into sites)
   * @param bisector : generates bisectors
   * @return array of Hyperbolas
   */
  public static Hyperbola[] buildCell(EdPoint[] sites, int siteI,
      ISiteBisector bisector) {

    final boolean db = false;
    if (db && T.update())
      T.msg("VornUtil.buildCell for site " + siteI);
    DArray ret = new DArray();

    EdPoint vs = sites[siteI];
    DArray cjk = new DArray();

    DArray cell = new DArray();
    cell.addInt(siteI);

    for (int j = 0; j < sites.length; j++) {
      if (j != siteI)
        cjk.addInt(j);
    }

    Hyperbola[] hy;
    hy = getHyp(bisector, sites, cell, cjk);
    for (int i = 0; i < hy.length; i++) {
      Hyperbola hi = hy[i];
      if (hi == null)
        continue;
      for (int j = i + 1; j < hy.length; j++) {
        Hyperbola hj = hy[j];
        if (hj == null)
          continue;
        if (false) {
          Tools.warn("is this correct?");
          if (hi.isEmpty() || hj.isEmpty())
            continue;
        }
        bisector
            .clipPair(hi, hj, vs, Hyperbola.findIntersections(hi, hj, null));
      }
    }

    for (int i = 0; i < hy.length; i++) {
      Hyperbola h = hy[i];
      if (h == null)
        continue;
      ret.add(h);
    }

    return (Hyperbola[]) ret.toArray(Hyperbola.class);
  }

  public static VornGraph vornForPointSites(EdPoint[] sites) {
    EdPoint[] s = new EdPoint[sites.length];
    for (int i = 0; i < sites.length; i++)
      s[i] = new EdPoint(sites[i].getOrigin());
    Hyperbola[] h = build(s, 1, PointBisector.S);
    return new VornGraph(h);
  }

  /**
   * Build order-k Voronoi diagram
   * @param sites
   * @param k
   */
  public static Hyperbola[] build(EdPoint[] sites, int k, ISiteBisector bisector) {

    boolean db2 = false;

    if (db2)
      Streams.out.println("\n\n\nVornUtil.build for n=" + sites.length
          + " sites, k=" + k + ", bisector=" + bisector);

    k = Math.max(k, 1);

    // build list of cells
    DArray cl = getSites(sites.length, k);
    DArray ret = new DArray();

    for (int ii = 0; ii < cl.size(); ii++) {
      if (false && ii != 0)
        continue;

      boolean db = db2 && (ii == 0);

      //    for (Iterator it = cl.iterator(); it.hasNext();) {
      DArray cell = (DArray) cl.get(ii); //)it.next();
      if (db)
        Streams.out.println(" cell=" + cell);

      EdPoint vs = sites[cell.getInt(0)];
      DArray cjk = new DArray();
      {
        int j = 0;
        for (int i = 0; i < k; i++) {
          int h = cell.getInt(i);
          while (j < h) {
            cjk.addInt(j++);
          }
          j = h + 1;
        }
        while (j < sites.length)
          cjk.addInt(j++);
      }
      Hyperbola[] hy;
      if (bisector.cellIsEmpty(sites, vs))
        hy = new Hyperbola[0];
      else
        hy = getHyp(bisector, sites, cell, cjk);

      for (int i = 0; i < hy.length; i++) {
        Hyperbola hi = hy[i];

        if (db)
          Streams.out.println(" hi=" + hi);

        if (hi == null
        //|| hi.isEmpty()
        ) {
          if (db)
            Streams.out.println("  null or empty");
          continue;
        }
        for (int j = i + 1; j < hy.length; j++) {
          Hyperbola hj = hy[j];
          if (db)
            Streams.out.println(" hj=" + hj);

          if (hj == null
          //|| hj.isEmpty()
          ) {
            if (db)
              Streams.out.println("  null or empty");
            continue;
          }
          if (db)
            Streams.out.println(" clipping");
          bisector.clipPair(hi, hj, vs, Hyperbola.findIntersections(hi, hj,
              null));
        }
      }

      if (db) {
        Streams.out.println("Cell " + cellLabel(sites, cell) + ":");
        for (int i = 0; i < hy.length; i++) {
          Hyperbola h = hy[i];
          if (h == null || h.isEmpty())
            continue;
          Streams.out.println(h + " " + h.minParameter() + "..."
              + h.maxParameter());

          Streams.out.println(" minPt=" + h.calcPoint(h.minParameter())
              + " maxPt=" + h.calcPoint(h.maxParameter()));
        }
      }

      for (int i = 0; i < hy.length; i++) {
        Hyperbola h = hy[i];
        if (h == null)
          continue;
        ret.add(h);
      }
    }

    return (Hyperbola[]) ret.toArray(Hyperbola.class);

  }

  private static Hyperbola[] getHyp(ISiteBisector bi, EdPoint[] sites,
      DArray cell, DArray comp) {

    final boolean db = false;

    if (db && T.update())
      T.msg("VornUtil.getHyp, cell=" + cell + ", comp=" + comp);
    Hyperbola[] h = new Hyperbola[cell.size() * comp.size()];
    int v = 0;
    for (int i = 0; i < cell.size(); i++) {
      int c = cell.getInt(i);
      for (int j = 0; j < comp.size(); j++) {
        int d = comp.getInt(j);

        Hyperbola hy = null;

        if (db && T.update())
          T.msg(" adding bisector for " + sites[c] + "," + sites[d]);
        hy = bi.getBisector(sites[c], sites[d]);
        h[v++] = hy;
      }
    }
    return h;
  }

  /**
   * Construct a dual graph from a Voronoi diagram.
   * Each node's data is the VornSite, and each edge's data is
   * the Hyperbola
   * @param edges : Hyperbolas
   * @return Graph
   */
  public static Graph buildDualGraph(Hyperbola[] edges) {
    Graph g = new Graph();

    // determine ids to associate with sites
    Map siteToIdMap = new HashMap();
    for (int i = 0; i < edges.length; i++) {
      Hyperbola h = edges[i];

      int ida = 0, idb = 0;
      for (int side = 0; side < 2; side++) {
        EdPoint site = (EdPoint) h.getData(side);
        Integer id = (Integer) siteToIdMap.get(site);
        if (id == null) {
          id = new Integer(g.newNode());
          siteToIdMap.put(site, id);
          g.setNodeData(id.intValue(), site);
        }
        if (side == 0)
          ida = id.intValue();
        else
          idb = id.intValue();
      }
      g.addEdgesBetween(ida, idb, h, h);
    }
    return g;
  }

  /**
  * Construct a Voronoi diagram using the 'contests' algorithm.
  * 
  * @param sites : DArray of EdCircles
  * @param clipType: CLIP_xxx
  * @return DArray of Hyperbolas, the edges
  */
  public static Hyperbola[] buildSubset(EdDisc[] sites, boolean withClipping,
      ISiteBisector bisector) {

    final boolean db = false;

    if (db)
      Streams.out.println("\n\nbuildSubset");

    DArray edgeList = new DArray();

    Hyperbola[] hy = constructArmMatrix(sites, withClipping);
    if (withClipping) {
      int ns = sites.length;
      for (int ii = 0; ii < ns; ii++) {

        for (int jj = 0; jj < ns; jj++) {
          if (jj == ii)
            continue;

          Hyperbola hij;
          hij = hy[(ii * ns) + jj];

          if (hij == null || hij.isEmpty())
            continue;

          for (int kk = 0; kk < ns; kk++) {
            if (kk == ii || kk == jj)
              continue;
            Hyperbola hik = hy[ii * ns + kk];
            if (hik == null || hik.isEmpty())
              continue;

            DArray iPts = Hyperbola.findIntersections(hij, hik, null);

            bisector.clipPair(hij, hik, sites[ii], iPts);
            if (hij.isEmpty())
              break;
          }
        }
      }
    }
    for (int i = 0; i < hy.length; i++) {
      Hyperbola hi = hy[i];
      if (hi == null || hi.isEmpty())
        continue;
      edgeList.add(hi);
    }

    return (Hyperbola[]) edgeList.toArray(Hyperbola.class);
  }

  private static Hyperbola[] constructArmMatrix(EdDisc[] sites,
      boolean withClipping) {
    int ns = sites.length;

    Hyperbola[] hy = new Hyperbola[ns * ns];

    for (int s1 = 0; s1 < ns; s1++) {
      EdDisc c1 = sites[s1];
      if (DiscUtil.contained(c1))
        continue;

      if (withClipping && DiscUtil.contained(c1))
        continue;

      for (int s2 = 0; s2 < ns; s2++) {
        if (s2 == s1)
          continue;
        EdDisc c2 = sites[s2];
        if (DiscUtil.contained(c2))
          continue;
        try {
          Hyperbola h = null;
          h = GuaranteedDiscBisector.S.getBisector(c1, c2);
          //constructWCBisector(c1, c2);
          hy[s1 + s2 * ns] = h;
        } catch (FPError e) {
          Tools.warn("FPError caught");
        }
      }
    }
    return hy;
  }

  private static Color[] dc = { MyColor.get(MyColor.PURPLE, .4),
      MyColor.get(MyColor.RED, .6), MyColor.get(MyColor.GREEN, .3),
      MyColor.get(MyColor.BROWN, .6), MyColor.get(MyColor.MAGENTA, .6),
      MyColor.get(MyColor.ORANGE, .5), };

  public static Color colorForLabel(String label) {
    return dc[MyMath.mod(indexForLabel(label), dc.length)];
  }

  public static int indexForLabel(String label) {
    int i = 0;
    if (label != null && label.length() > 0) {
      i = Character.toUpperCase(label.charAt(0)) - 'A';
    }
    return Math.max(0, i);
  }

}

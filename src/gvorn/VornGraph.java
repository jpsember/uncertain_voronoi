package gvorn;

import java.awt.*;
import java.util.*;
import testbed.*;
import base.*;

public class VornGraph extends Graph implements Renderable {

  public void setAppearance(Color c, int stroke, int markType, boolean labels) {
    if (c == null)
      c = this.color;
    if (c == null)
      c = MyColor.cDARKGREEN;
    this.color = c;
    this.stroke = stroke;
    this.markType = markType;
    this.labels = labels;
  }

  private Color color = MyColor.cDARKGREEN;
  private int stroke = -1;
  private int markType = -1;
  private boolean labels = false;

  public boolean isInf(int node) {
    return node == infNode;
  }

  public static Hyperbola[] addCrossings(Hyperbola[] edges) {
    final boolean db = true;
    final double EPS = 1e-4;
    DArray e = new DArray(edges);
    DArray ipt = new DArray();
    for (int i = 0; i < e.size(); i++) {
      for (int j = i + 1; j < e.size(); j++) {
        Hyperbola hi = (Hyperbola) e.get(i);
        Hyperbola hj = (Hyperbola) e.get(j);

        if (hi.getData(0) == hj.getData(0) && hi.getData(1) == hj.getData(1))
          continue;
        if (db && T.update())
          T.msg("testing intersections between " + hi + " and " + hj);
        Hyperbola.findIntersections(hi, hj, ipt);
        for (int k = 0; k < ipt.size(); k++) {
          FPoint2 pt = ipt.getFPoint2(k);
          double ti = hi.calcParameter(pt);
          double tj = hj.calcParameter(pt);

          if (db && T.update())
            T.msg(" intersection at " + pt + ", ti=" + ti + ", tj=" + tj);
          //          if (!hi.containsPoint(ti))
          //            continue;
          //          if (!hj.containsPoint(tj))
          //            continue;

          boolean inBoth = (hi.containsPoint(ti - EPS) || hi.containsPoint(ti
              + EPS))
              && (hj.containsPoint(tj - EPS) || hj.containsPoint(tj + EPS));
          if (inBoth) {

            if (hi.containsPoint(ti - EPS) && hi.containsPoint(ti + EPS)) {
              if (db && T.update())
                T.msg("  splitting in two");
              // split into two
              Hyperbola h = new Hyperbola(hi);
              h.clipAllBut(hi.minParameter(), ti);
              e.set(i, h);
              h = new Hyperbola(hi);
              h.clipAllBut(ti, hi.maxParameter());
              e.add(h);
            }

            if (hj.containsPoint(tj - EPS) && hj.containsPoint(tj + EPS)) {
              // split into two
              //            Hyperbola
              Hyperbola h = new Hyperbola(hj);
              h.clipAllBut(hj.minParameter(), tj);
              e.set(j, h);
              h = new Hyperbola(hj);
              h.clipAllBut(tj, hj.maxParameter());
              e.add(h);
            }
          }
        }
      }
    }
    //    if (db)
    //      Streams.out.println("added crossings, original=" + edges.length
    //          + "\n new=" + e.size());
    edges = (Hyperbola[]) e.toArray(Hyperbola.class);
    return edges;
  }
  /**
   * Construct a graph from a Voronoi diagram.
   * Each node's data is the FPoint2 containing its location.
   * The first node is the point at infinity.
   * To access edge data, use getHyperbola() and getSite() methods.
   * @param edges : Hyperbolas
   * @return Graph
   */
  public VornGraph(Hyperbola[] edges) {

    //    setAppearance(null, -1, -1, false);

    final boolean db = false;

    Map map = new HashMap();

    if (db)
      Streams.out.println("buildGraph");

    Map locMap = new HashMap();
    infNode = newNode();
    for (int i = 0; i < edges.length; i++) {

      Hyperbola h = edges[i];
      map.put(h, h);

      if (db)
        Streams.out.println(" edge #" + i + ": " + h + " isEmpty="
            + h.isEmpty());
      if (h.isEmpty())
        continue;

      if (db)
        Streams.out
            .println("t=" + h.minParameter() + " .. " + h.maxParameter());

      EdgeData[] edgeDataSets = new EdgeData[2];
      int[] nodeIds = new int[2];
      for (int endpoint = 0; endpoint < 2; endpoint++) {

        double t = (endpoint == 0) ? h.minParameter() : h.maxParameter();
        int nodeId = infNode;
        if (t != ((endpoint == 0) ? Hyperbola.CLIP_MIN : Hyperbola.CLIP_MAX)) {
          FPoint2 pt = h.calcPoint(t);
          nodeId = findPoint(locMap, pt);

          if (nodeId < 0) {
            nodeId = newNode(pt);
            storePoint(locMap, pt, nodeId);
            if (db)
              Streams.out.println(" ep " + endpoint + " t=" + Tools.f(t)
                  + " pt=" + pt + " node=" + nodeId);
          }
        }
        nodeIds[endpoint] = nodeId;

        EdgeData ed = new EdgeData(h, endpoint != 0);
        edgeDataSets[endpoint] = ed;
      }

      // if there is already an edge between these nodes, with the
      // same data, skip.
      int ti = 0;
      if (nodeIds[0] > nodeIds[1])
        ti = 1;
      String key = nodeIds[ti] + "/" + nodeIds[ti ^ 1] + "/"
          + edgeDataSets[ti].getData(0) + "/" + edgeDataSets[ti ^ 1].getData(0);
      if (db)
        Streams.out.println(" key=" + key);

      if (map.containsKey(key))
        continue;
      map.put(key, key);

      addEdgesBetween(nodeIds[0], nodeIds[1], edgeDataSets[0], edgeDataSets[1]);
      if (db)
        Streams.out.println("...adding edge from " + nodeIds[0] + " to "
            + nodeIds[1]);

    }

    if (db) {
      for (Iterator it = getNodeList().iterator(); it.hasNext();) {
        int id = ((Integer) it.next()).intValue();
        Streams.out.println("node " + id);
        for (int j = 0; j < nCount(id); j++) {
          Hyperbola h = getHyperbola(id, j);

          Object obj = getSite(id, j, 0);

          Object obj2 = getSite(id, j, 1);
          Streams.out.println(" edge #" + j + " h=" + h + ", right=" + obj
              + " left=" + obj2);
        }
      }
    }

    // sort edges of non-infinite vertices
    {
      for (Iterator it = getNodeList().iterator(); it.hasNext();) {
        int id = ((Integer) it.next()).intValue();
        if (isInf(id))
          continue;

        sortEdges(id, new Comparator() {

          public int compare(Object o1, Object o2) {
            Object[] a1 = (Object[]) o1;
            Object[] a2 = (Object[]) o2;
            //            VornGraph g = (Graph) a1[0];
            int node = ((Integer) a1[1]).intValue();
            int e1 = ((Integer) a1[2]).intValue();
            int e2 = ((Integer) a2[2]).intValue();

            Hyperbola h1 = getHyperbola(node, e1);
            Hyperbola h2 = getHyperbola(node, e2);
            FPoint2 origin = getVertLoc(node);

            FPoint2 p1 = nearVert(h1, origin);
            FPoint2 p2 = nearVert(h2, origin);

            double an1 = MyMath.polarAngle(origin, p1);
            double an2 = MyMath.polarAngle(origin, p2);
            return (int) Math.signum(an1 - an2);
          }

        });
      }
    }

    if (db)
      Streams.out.println(this);
  }

  private static FPoint2 nearVert(Hyperbola h, FPoint2 origin) {
    double t = h.calcParameter(origin);
    double m = (h.minParameter() + h.maxParameter()) * .5;
    double scl = Math.abs(m - t);
    final double OFFSET = 1e-7;
    if (scl < OFFSET)
      Tools.warn("scale too low"); //: " + scl);
    if (t < m)
      t += OFFSET;
    else
      t -= OFFSET;
    return h.calcPoint(t);
  }

  public FPoint2 getVertLoc(int node) {
    return (FPoint2) nodeData(node);
  }

  /**
   * Get hyperbola associated with an edge
   * @param g : graph 
   * @param node : node
   * @param edgeNumber : edge number within node
   * @return Hyperbola associated with edge
   */
  public Hyperbola getHyperbola(int node, int edgeNumber) {
    EdgeData ed = (EdgeData) edgeData(node, edgeNumber);
    return ed.h;
  }

  /**
   * Get site to one side of an edge
   * @param g : graph
   * @param node : node
   * @param edgeNumber : edge number within node
   * @param side : 0 for site to right, 1 for site to left
   * @return Object 
   */
  public Object getSite(int node, int edgeNumber, int side) {
    EdgeData ed = (EdgeData) edgeData(node, edgeNumber);
    return ed.getData(side);
  }
  /**
   * Get EdCircle to one side of an edge
   * @param g : graph
   * @param node : node
   * @param edgeNumber : edge number within node
   * @param side : 0 for site to right, 1 for site to left
   * @return Object 
   */
  public EdDisc getCircle(int node, int edgeNumber, int side) {
    return (EdDisc) getSite(node, edgeNumber, side);
  }

  /**
   * See if a point is within a map, or one very close to it
   * @param map : map to test
   * @param pt : point to look for
   * @return id associated with point, if found; or -1
   */
  private static int findPoint(Map map, FPoint2 pt) {
    Integer id = (Integer) map.get(getPointKey(pt));
    return (id == null) ? -1 : id.intValue();
  }

  private static FPoint2 getPointKey(FPoint2 pt) {
    final double SCALE = 1e-2;
    FPoint2 pt2 = new FPoint2(pt.x - MyMath.mod(pt.x, SCALE), pt.y
        - MyMath.mod(pt.y, SCALE));
    return pt2;
  }

  private static void storePoint(Map map, FPoint2 pt, int id) {
    map.put(getPointKey(pt), new Integer(id));
  }
  private static class EdgeData {
    public String toString() {
      StringBuilder sb = new StringBuilder();
      sb.append("EdgeData");
      sb.append(" h=" + h.toString());
      sb.append("  flip=" + flip);
      return sb.toString();
    }
    public EdgeData(Hyperbola h, boolean flip) {
      this.h = h;
      this.flip = flip;
    }
    public Object getData(int side) {
      return h.getData(side ^ (flip ? 1 : 0));
    }
    public Hyperbola h;
    boolean flip;
  }

  public void plot(Color c, int stroke, boolean plotVertices, boolean labels) {
    setAppearance(c, stroke, plotVertices ? Globals.MARK_DISC : -1, labels);
    render(null, -1, -1);
  }

  public boolean specialPlot;

  public void render(Color c, int stroke, int markType) {
    if (c != null)
      this.color = c;
    if (stroke >= 0)
      this.stroke = stroke;
    if (markType >= 0)
      this.markType = markType;

    V.pushColor(this.color);
    if (this.stroke >= 0)
      V.pushStroke(this.stroke);
    //    if (c == null)
    //      c = MyColor.cRED;
    //    V.pushColor(c);
    if (false)
      render(true, markType > 0, markType > 0);

    for (Iterator it = getNodeList().iterator(); it.hasNext();) {
      int id = ((Integer) it.next()).intValue();
      for (int n = 0; n < nCount(id); n++) {
        int id2 = neighbor(id, n);
        if (id2 < id)
          continue;
        Hyperbola hi = getHyperbola(id, n);

          hi.render(null, -1, -1);

        if (this.markType >= 0) {
          //          V.pushScale(.4);
          for (int i = 0; i < 2; i++) {
            double t = (i == 0) ? hi.minParameter() : hi.maxParameter();
            if (t >= -2000 && t <= 2000) {
              FPoint2 pt = hi.calcPoint(t);
              V.mark(pt, this.markType);
            }
          }
          //          V.popScale();
        }

        if (labels) {
          V.pushScale(.5);
          {
            double t0 = hi.minParameter();
            double t1 = hi.maxParameter();
            double t = 0;
            if (t1 - t0 < 5)
              t = (t0 + t1) * .5;
            else {
              if (t < t0)
                t = t0 + 2;
              if (t > t1)
                t = t1 - 2;
            }
            FPoint2 pt = hi.calcPoint(t);
            V.draw(hi.toString(), pt, Globals.TX_BGND);
          }
          for (int cp = 0; cp < 2; cp++) {
            Object ca = hi.getData(Hyperbola.MINCLIPPED + cp);
            if (ca != null) {
              String s2 = ca.toString();
              final double INSET = 2;

              double m0 = hi.minParameter(), m1 = hi.maxParameter();

              double t = cp == 0 ? m0 + INSET : m1 - INSET;
              t = MyMath.clamp(t, m0, m1);
              V.draw(s2, hi.calcPoint(t), Globals.TX_BGND);
            }
          }
          V.popScale();
        }
      }
    }
    if (this.stroke >= 0)
      V.popStroke();
    V.popColor();
  }

  /**
   * @deprecated
   * @param plotEdges
   * @param plotVertices
   * @param labels
   */
  public void render(boolean plotEdges, boolean plotVertices, boolean labels) {

    for (Iterator it = getNodeList().iterator(); it.hasNext();) {
      int id = ((Integer) it.next()).intValue();
      for (int n = 0; n < nCount(id); n++) {
        int id2 = neighbor(id, n);
        if (id2 < id)
          continue;
        Hyperbola hi = getHyperbola(id, n);
        if (plotEdges) {
          hi.render(null, -1, -1);
        }

        if (plotVertices) {
          V.pushScale(.4);
          for (int i = 0; i < 2; i++) {
            double t = (i == 0) ? hi.minParameter() : hi.maxParameter();
            if (t >= -2000 && t <= 2000) {
              FPoint2 pt = hi.calcPoint(t);
              V.mark(pt, TestBed.MARK_CIRCLE);
            }
          }
          V.popScale();
        }

        if (labels) {
          V.pushScale(.5);
          {
            double t0 = hi.minParameter();
            double t1 = hi.maxParameter();
            double t = 0;
            if (t1 - t0 < 5)
              t = (t0 + t1) * .5;
            else {
              if (t < t0)
                t = t0 + 2;
              if (t > t1)
                t = t1 - 2;
            }
            FPoint2 pt = hi.calcPoint(t);
            V.draw(hi.toString(), pt, Globals.TX_BGND);
          }
          for (int cp = 0; cp < 2; cp++) {
            Object ca = hi.getData(Hyperbola.MINCLIPPED + cp);
            if (ca != null) {
              String s2 = ca.toString();
              final double INSET = 2;

              double m0 = hi.minParameter(), m1 = hi.maxParameter();

              double t = cp == 0 ? m0 + INSET : m1 - INSET;
              t = MyMath.clamp(t, m0, m1);
              V.draw(s2, hi.calcPoint(t), Globals.TX_BGND);
            }
          }
          V.popScale();
        }
      }
    }
  }

  public void printEdges() {
    // vp vp = TestBed.view();
    //    vp.pushColor(MyColor.get(MyColor.BROWN));
    V.pushScale(.7);

    StringBuilder sb = new StringBuilder();
    sb.append("Edges:\n");

    for (Iterator it = getNodeList().iterator(); it.hasNext();) {
      int id = ((Integer) it.next()).intValue();
      for (int n = 0; n < nCount(id); n++) {
        int id2 = neighbor(id, n);
        if (id2 < id)
          continue;
        Hyperbola hi = getHyperbola(id, n);
        sb.append(toString(hi));
        sb.append('\n');

      }
    }
    V.draw(sb.toString(), 95, 95, Globals.TX_CLAMP | 40);
    V.popScale();
    //    vp.popColor();
  }

  private static String toString(Hyperbola hi) {
    StringBuilder sb = new StringBuilder();
    EdDisc c0 = (EdDisc) hi.getData(0), c1 = (EdDisc) hi.getData(1);
    sb.append('<');
    sb.append(c0.getLabel());
    sb.append(',');
    sb.append(c1.getLabel());
    sb.append('>');
    sb.append("  ");
    if (hi.isEmpty())
      sb.append('X');
    else {
      double min = hi.minParameter();
      double max = hi.maxParameter();
      if (min == Hyperbola.CLIP_MIN)
        sb.append("  <<< ");
      else
        sb.append(Tools.f(min));
      sb.append(" ");
      if (max == Hyperbola.CLIP_MAX)
        sb.append("  >>> ");
      else
        sb.append(Tools.f(max));
    }
    return sb.toString();
  }

  /**
   * Find edge from node that has a particular site to one side
   * @param site : site to search for
   * @param currNode : node whose edges we're examining
   * @param side : side of edge to look at 0,left 1,right
   * @return index of edge, or -1
   */
  public int findNextEdgeBorderingSite(EdDisc site, int currNode, int side) {
    int ne = -1;
    for (int j = 0; j < nCount(currNode); j++) {
      if (getSite(currNode, j, side) != site)
        continue;
      ne = j;
      break;
    }
    return ne;

  }
  private int infNode;

  private static class MyComparator implements Comparator {
    public MyComparator(EdDisc site, VornGraph g) {
      this.site = site;
      this.g = g;
    }
    private EdDisc site;
    private VornGraph g;
    public int compare(Object arg0, Object arg1) {
      int i0 = ((Integer) arg0).intValue();
      int i1 = ((Integer) arg1).intValue();
      int n0 = i0 >> 16;
      int n1 = i1 >> 16;
      int e0 = i0 & 0xffff;
      int e1 = i1 & 0xffff;

      Hyperbola h0 = g.getHyperbola(n0, e0);
      Hyperbola h1 = g.getHyperbola(n1, e1);

      FPoint2 p0 = h0.calcPoint((h0.maxParameter() + h0.maxParameter()) * .5);
      FPoint2 p1 = h1.calcPoint((h1.maxParameter() + h1.maxParameter()) * .5);

      double a0 = MyMath.polarAngle(site.getOrigin(), p0);
      double a1 = MyMath.polarAngle(site.getOrigin(), p1);
      return (int) Math.signum(a0 - a1);
    }
  }

  /**
   * Find list of edges bordering a site (to the left).
   * 
   * @param site 
   * @return DArray of ints, where each int has node id in upper 16 bits,
   *   edge number in lower 16 bits
   */
  public DArray findEdgesBorderingSite(EdDisc site) {
    boolean db = true;

    if (db && T.update())
      T.msg("findEdgesBorderingSite " + site.getLabel() + "\n" + this);

    DArray a = new DArray();

    DArray nl = getNodeList();
    for (Iterator it = nl.iterator(); it.hasNext();) {
      int id = ((Integer) it.next()).intValue();

      for (int j = 0; j < nCount(id); j++) {

        if (db && T.update())
          T.msg(" id=" + id + " j=" + j + " site="
              + getCircle(id, j, 1).getLabel());
        if (getSite(id, j, 1) != site)
          continue;
        a.addInt((id << 16) + j);
        if (db && T.update())
          T.msg(" added it");
      }
    }

    // sort array by polar angle around site
    a.sort(new MyComparator(site, this));
    if (db && T.update())
      T.msg(" after sorting: " + a);
    return a;
  }
}

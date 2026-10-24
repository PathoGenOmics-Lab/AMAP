/* AMAP website: the drawings. Every gene drawn here is a real one, from AMAP
   run on Mycobacterium tuberculosis CDC1551 (NC_002755.2) against the H37Rv
   reference it ships with, unless a figure says otherwise; the coordinates are
   that run's. tools/check_site.py puts each scenario through the pipeline's
   own code and fails if a drawing keeps something the code would not.
   Plain JavaScript, no build step. */
(function () {
  "use strict";

  function mq(query) {
    return !!(window.matchMedia && window.matchMedia(query).matches);
  }

  var reduceMotion = mq("(prefers-reduced-motion: reduce)");

  function esc(text) {
    return String(text)
      .replace(/&/g, "&amp;")
      .replace(/</g, "&lt;")
      .replace(/>/g, "&gt;")
      .replace(/"/g, "&quot;");
  }

  function fmt(n) {
    return String(n).replace(/\B(?=(\d{3})+(?!\d))/g, ",");
  }

  // A tooltip's HTML as one line a screen reader can say.
  function plain(html) {
    // A line break of the tip reads as a new sentence before a capital and
    // as a comma before anything else.
    return String(html)
      .replace(/\s*<br\s*\/?>\s*(<[^>]+>)*\s*(.)/gi, function (_m, _tags, next) {
        return (/[A-Z]/.test(next) ? ". " : ", ") + next;
      })
      .replace(/<[^>]+>/g, "")
      .replace(/&lt;/g, "<")
      .replace(/&gt;/g, ">")
      .replace(/&quot;/g, '"')
      .replace(/&amp;/g, "&")
      .replace(/(\d)-(\d)/g, "$1 to $2")
      .replace(/\(\+\)/g, "plus strand")
      .replace(/\(−\)/g, "minus strand")
      .replace(/([.:;,])\s*\./g, "$1")
      .replace(/\s+/g, " ")
      .trim();
  }

  /* ── The track renderer ──

     A drawing is rows of genes on one coordinate axis:
       { from, to, w, rows: [{ label, cls, phase, genes: [gene] }],
         ribbons: [{ a: [row, gene], b: [row, gene], kind, phase }],
         marks: [{ row, s, e, label }], edges: [{ row, at, side, label }] }
     A gene is { s, e, st, name, cls, tip, open, gap, dim, id }, 1-based and
     inclusive like everything AMAP writes. Drawing a spec leaves on it the
     height it took (_h) and the kinds of mark it used (_used), for a legend
     that names exactly those. */

  var figureCount = 0;

  // open: the ends a contig cut, drawn ragged. flat: the end an alignment
  // stops short of, drawn square.
  function arrow(x1, x2, y, h, strand, open, flat) {
    var tip = Math.min(8, Math.max(0, (x2 - x1) * 0.45));
    var mid = y + h / 2;
    var has = function (s, end) {
      return !!s && s.indexOf(end) !== -1;
    };
    // A partial gene has no arrowhead at its cut end, and a ragged edge there.
    var zig = function (x, dir) {
      var d = 3 * dir;
      return "L" + (x + d) + "," + (y + h * 0.25) + "L" + x + "," + mid + "L" + (x + d) + "," + (y + h * 0.75);
    };
    if (strand === "-") {
      var rightOpen = has(open, "5'"); // 5' end is on the right
      var leftOpen = has(open, "3'");
      var p = "M" + x2 + "," + y;
      p += leftOpen ? "H" + x1 + zig(x1, 1) + "L" + x1 + "," + (y + h) : has(flat, "3'") ? "H" + x1 + "V" + (y + h) : "H" + (x1 + tip) + "L" + x1 + "," + mid + "L" + (x1 + tip) + "," + (y + h);
      p += "H" + x2;
      p += rightOpen ? zig(x2, -1) + "Z" : "Z";
      return p;
    }
    var endOpen = has(open, "3'");
    var startOpen = has(open, "5'");
    var q = "M" + x1 + "," + y;
    q += endOpen ? "H" + x2 + zig(x2, -1) + "L" + x2 + "," + (y + h) : has(flat, "3'") ? "H" + x2 + "V" + (y + h) : "H" + (x2 - tip) + "L" + x2 + "," + mid + "L" + (x2 - tip) + "," + (y + h);
    q += "H" + x1;
    q += startOpen ? zig(x1, 1) + "Z" : "Z";
    return q;
  }

  function textWidth(text, size) {
    return String(text).length * size * 0.58;
  }

  // Which steps of a walk-through an element shows in, as attributes.
  function steps(o) {
    return (o.from ? ' data-from="' + o.from + '"' : "") + (o.until ? ' data-until="' + o.until + '"' : "");
  }

  function r1(v) {
    return Math.round(v * 10) / 10;
  }

  function draw(spec, width) {
    var id = "f" + ++figureCount;
    var used = {};
    // Drawn at the width it is shown at, so its text is always the page's size.
    var W = Math.max(280, Math.round(width || spec.w || 640));
    var labelW = spec.labelW == null ? 104 : spec.labelW;
    var x0 = labelW;
    var x1 = W - (spec.padRight == null ? 8 : spec.padRight);
    var geneH = spec.geneH || 16;
    var rowGap = spec.rowGap || 30;
    var top = spec.top == null ? 10 : spec.top;
    var span = spec.to - spec.from;
    var live = !spec.inert;
    function X(p) {
      return x0 + ((p - spec.from) / span) * (x1 - x0);
    }
    // Genes of a row that overlap by more than a few codons go on lanes of
    // their own, one under the other (Rv0590 and Rv0590A share 52 bp).
    var LANE = geneH + 4;
    var lanes = spec.rows.map(function (row) {
      var ends = [];
      var order = row.genes.map(function (g, k) {
        return k;
      }).sort(function (a, b) {
        return row.genes[a].s - row.genes[b].s;
      });
      var lane = [];
      order.forEach(function (k) {
        var g = row.genes[k];
        var l = 0;
        while (l < ends.length && ends[l] - g.s > 30) l += 1;
        ends[l] = g.e;
        lane[k] = l;
      });
      return { of: lane, n: Math.max(1, ends.length) };
    });
    // Where each row's genes sit; a row can ask for more room above it.
    var ys = [];
    var y = top;
    spec.rows.forEach(function (row, i) {
      if (row.overlay && ys.length) {
        // Drawn over the row before it: the same genes at another moment.
        ys.push(ys[ys.length - 1]);
        y += Math.max(0, lanes[i].n - lanes[i - 1].n) * LANE;
        return;
      }
      y += row.space || 0;
      ys.push(y);
      y += geneH + (lanes[i].n - 1) * LANE + rowGap;
    });
    function gy0(r, k) {
      return ys[r] + (lanes[r].of[k] || 0) * LANE;
    }
    var height = y - rowGap + (spec.axis ? 30 : 8);
    spec._h = Math.round(height);
    var out = [];
    // A drawing is one stop for the keyboard (a group of named genes); a small
    // one inside a link is a picture, and only the link is a stop.
    out.push(
      '<svg class="trk" viewBox="0 0 ' + W + " " + Math.round(height) + '" ' +
        (live ? 'role="group" aria-roledescription="gene track" aria-label="' + esc(spec.alt || "") + '"' : 'aria-hidden="true" focusable="false"') + ">"
    );
    out.push(
      "<defs>" +
        '<pattern id="hatch-' + id + '" width="6" height="6" patternUnits="userSpaceOnUse" patternTransform="rotate(45)">' +
        '<line class="hatch-line" x1="0" y1="0" x2="0" y2="6"/></pattern>' +
        '<pattern id="dh-' + id + '" width="5" height="5" patternUnits="userSpaceOnUse" patternTransform="rotate(45)">' +
        '<line class="diff-hatch" x1="0" y1="0" x2="0" y2="5"/></pattern>' +
        '<pattern id="dhr-' + id + '" width="5" height="5" patternUnits="userSpaceOnUse" patternTransform="rotate(45)">' +
        '<line class="diff-hatch diff-hatch-ref" x1="0" y1="0" x2="0" y2="5"/></pattern>' +
        "</defs>"
    );
    var masks = [];

    // Ribbons first, under the genes: straight where the two agree, slanted
    // where their coordinates differ. Their sides are drawn in full colour,
    // solid where they agree and dashed where they do not, so the difference
    // shows by shape as well as by hue.
    (spec.ribbons || []).forEach(function (r) {
      var ga = spec.rows[r.a[0]].genes[r.a[1]];
      var gb = spec.rows[r.b[0]].genes[r.b[1]];
      var ya = gy0(r.a[0], r.a[1]) + geneH + 1;
      var yb = gy0(r.b[0], r.b[1]) - 1;
      var clip = function (p) {
        return Math.max(spec.from, Math.min(spec.to, p));
      };
      var ax1 = r1(X(clip(ga.s) - 0.5));
      var ax2 = r1(X(clip(ga.e) + 0.5));
      var bx2 = r1(X(clip(gb.e) + 0.5));
      var bx1 = r1(X(clip(gb.s) - 0.5));
      used[r.kind] = 1;
      out.push(
        '<g class="rib-g' + (r.phase ? " ph-" + r.phase : "") + '"' + steps(r) + (r.k ? ' data-k="' + r.k + '"' : "") + ">" +
          '<polygon class="rib rib-' + r.kind + '" points="' + [ax1, ya, ax2, ya, bx2, yb, bx1, yb].join(",") + '"/>' +
          '<path class="rib-edge rib-edge-' + r.kind + '" d="M' + ax1 + "," + ya + "L" + bx1 + "," + yb + "M" + ax2 + "," + ya + "L" + bx2 + "," + yb + '"/>' +
          "</g>"
      );
    });

    (spec.marks || []).forEach(function (m) {
      var mx1 = X(m.s - 0.5);
      var mx2 = X(m.e + 0.5);
      var my1 = ys[m.rows[0]] - 4;
      var my2 = ys[m.rows[1]] + geneH + (lanes[m.rows[1]].n - 1) * LANE + 4;
      used[m.kind === "ref" ? "diff-ref" : "diff"] = 1;
      out.push('<g class="mark' + (m.phase ? " ph-" + m.phase : "") + '"' + steps(m) + ">");
      out.push(
        '<rect class="diff' + (m.kind ? " diff-" + m.kind : "") + '" style="fill:url(#' + (m.kind === "ref" ? "dhr-" : "dh-") + id + ')" x="' + mx1.toFixed(1) + '" y="' + my1 +
          '" width="' + Math.max(2, mx2 - mx1).toFixed(1) + '" height="' + (my2 - my1) + '" rx="2"/>'
      );
      if (m.label) {
        var lx = m.labelAt === "left" ? mx1 - 4 : m.labelAt === "right" ? mx2 + 4 : (mx1 + mx2) / 2;
        var anchor = m.labelAt === "left" ? "end" : m.labelAt === "right" ? "start" : "middle";
        // A label that would run off either side is pinned to that side.
        var lw = textWidth(m.label, 10.5);
        var left = anchor === "end" ? lx - lw : anchor === "middle" ? lx - lw / 2 : lx;
        if (left + lw > W - 2) {
          lx = W - 2;
          anchor = "end";
        } else if (left < x0) {
          lx = x0;
          anchor = "start";
        }
        out.push('<text class="diff-label' + (m.kind ? " diff-label-" + m.kind : "") + '" x="' + lx.toFixed(1) + '" y="' + (my1 - 4) + '" text-anchor="' + anchor + '">' + esc(m.label) + "</text>");
      }
      out.push("</g>");
    });

    var tabStop = !live;
    spec.rows.forEach(function (row, i) {
      var gy = ys[i];
      // The label and the line stay while a row is swapped for a later
      // version of itself (an overlay), so they only take its first step.
      // Until then they show faintly, so the reader sees where it will come.
      out.push('<g class="row-frame' + (row.phase ? " fr-" + row.phase : "") + '"' + (row.from ? ' data-ghost-from="' + row.from + '"' : "") + ">");
      if (labelW > 0 && !row.overlay) {
        var lines = String(row.label).split("\n");
        lines.forEach(function (line, k) {
          out.push(
            '<text class="row-label" x="0" y="' + (gy + geneH / 2 + 4 + (k - (lines.length - 1) / 2) * 12).toFixed(1) + '">' + esc(line) + "</text>"
          );
        });
      }
      if (!row.overlay) {
        out.push('<line class="row-line" x1="' + x0 + '" x2="' + x1 + '" y1="' + (gy + geneH / 2) + '" y2="' + (gy + geneH / 2) + '"/>');
      }
      out.push("</g>");
      out.push('<g class="row' + (row.phase ? " ph-" + row.phase : "") + '"' + steps(row) + ">");
      (row.empty ? [row.empty] : []).forEach(function (note) {
        out.push('<text class="gsub" x="' + ((x0 + x1) / 2).toFixed(1) + '" y="' + (gy + geneH / 2 - 5) + '" text-anchor="middle">' + esc(note) + "</text>");
      });
      // The keyboard starts at the leftmost gene of the first row.
      var first = -1;
      if (!tabStop && row.genes.length) {
        first = 0;
        row.genes.forEach(function (g, k) {
          if (g.s < row.genes[first].s) first = k;
        });
        tabStop = true;
      }
      // Wide genes first, so a narrow one is drawn over its neighbours and
      // its target is never under theirs.
      var order = row.genes.map(function (g, k) {
        return k;
      }).sort(function (a, b) {
        return row.genes[b].e - row.genes[b].s - (row.genes[a].e - row.genes[a].s);
      });
      order.forEach(function (k) {
        var g = row.genes[k];
        var gy = gy0(i, k);
        // A gene that runs past the window is cut at its edge, left open there.
        var cutL = g.s < spec.from;
        var cutR = g.e > spec.to;
        var gx1 = X(Math.max(g.s, spec.from) - 0.5);
        var gx2 = X(Math.min(g.e, spec.to) + 0.5);
        var open = g.open || "";
        if (cutL) open += g.st === "-" ? "3'" : "5'";
        if (cutR) open += g.st === "-" ? "5'" : "3'";
        var cls = g.cls || row.cls;
        var klass = "g g-" + cls;
        var fill = cls === "pseudo" ? ' style="fill:url(#hatch-' + id + ')"' : "";
        // A call with nothing aligned is drawn as a call until the calls are
        // compared, then dashed, and grey if dropped. A figure that only
        // shows the end (final) never shows the first.
        var c1 = spec.phased && g.phaseCls && g.phaseCls.indexOf("c1-call") !== -1;
        if (!(c1 && spec.final)) used[cls === "drop" && g.lk ? g.lk : cls] = 1;
        if (g.short) used.short = 1;
        if (g.partial) used.partial = 1;
        if (c1) used[g.phaseCls.indexOf("ph-dim") !== -1 ? "drop" : "c1"] = 1;
        var attrs = live
          ? ' role="img" tabindex="' + (k === first ? 0 : -1) + '" aria-label="' + esc(plain(g.tip || g.name || "")) + '"' + (g.tip ? ' data-tip="' + esc(g.tip) + '"' : "")
          : "";
        out.push('<g class="hot' + (g.phaseCls ? " " + g.phaseCls : "") + '"' + attrs + (g.k ? ' data-k="' + g.k + '"' : "") + ">");
        if (live) {
          // A target at least 14 px wide and 24 tall, whatever the gene's size.
          var hw = Math.max(14, gx2 - gx1 + 4);
          var hh = Math.max(24, geneH + 8);
          out.push('<rect class="hit" x="' + ((gx1 + gx2) / 2 - hw / 2).toFixed(1) + '" y="' + (gy + geneH / 2 - hh / 2).toFixed(1) + '" width="' + hw.toFixed(1) + '" height="' + hh + '" rx="4"/>');
        }
        var maskAttr = "";
        if (g.short && gx2 - gx1 > 6) {
          // An alignment short of its whole protein fades out on the side it
          // lacks: its end (3'), or for a piece of a cut gene its start (5').
          // The 3' end is on the right of a plus-strand gene.
          var mid = id + "-" + i + "-" + k;
          var fadeW = Math.min(16, (gx2 - gx1) * 0.4);
          var right = (g.short === "3'") === (g.st !== "-");
          var fx1 = right ? gx2 - fadeW : gx1 + fadeW;
          var fx2 = right ? gx2 : gx1;
          masks.push(
            '<linearGradient id="sg-' + mid + '" gradientUnits="userSpaceOnUse" x1="' + fx1.toFixed(1) + '" x2="' + fx2.toFixed(1) + '" y1="0" y2="0">' +
              '<stop offset="0" stop-color="#fff"/><stop offset="1" stop-color="#262626"/></linearGradient>' +
              '<mask id="sm-' + mid + '" maskUnits="userSpaceOnUse" x="' + (gx1 - 4).toFixed(1) + '" y="' + (gy - 4) + '" width="' + (gx2 - gx1 + 8).toFixed(1) + '" height="' + (geneH + 8) + '">' +
              '<rect x="' + (gx1 - 4).toFixed(1) + '" y="' + (gy - 4) + '" width="' + (gx2 - gx1 + 8).toFixed(1) + '" height="' + (geneH + 8) + '" fill="url(#sg-' + mid + ')"/></mask>'
          );
          maskAttr = ' mask="url(#sm-' + mid + ')"';
        }
        out.push('<path class="' + klass + (cutL || cutR ? " g-cut" : "") + (g.short ? " g-short" : "") + '"' + fill + maskAttr + ' d="' + arrow(gx1, gx2, gy, geneH, g.st, open, g.short || "") + '"/>');
        if (g.gap) {
          var gpx = X(g.gap);
          out.push('<line class="gap-mark" x1="' + gpx.toFixed(1) + '" x2="' + gpx.toFixed(1) + '" y1="' + (gy - 5) + '" y2="' + (gy + geneH + 5) + '"/>');
          if (g.gapLabel) {
            out.push('<text class="diff-label diff-label-ref" x="' + gpx.toFixed(1) + '" y="' + (gy + geneH + 16) + '" text-anchor="middle">' + esc(g.gapLabel) + "</text>");
          }
        }
        // A notch over each end AMAP moved away from Bakta's call.
        (g.moved || []).forEach(function (at) {
          var nx = X(at);
          used.moved = 1;
          out.push('<path class="moved" d="M' + (nx - 4).toFixed(1) + "," + (gy - 7) + "h8l-4,5z" + '"/>');
        });
        if (g.name) {
          var size = 11.5;
          var tw = textWidth(g.name, size);
          if (tw + 12 < gx2 - gx1) {
            var nx = (gx1 + gx2) / 2;
            if (g.gap) {
              // Put the name in the wider of the two halves the gap leaves.
              var gpx2 = X(g.gap);
              nx = gpx2 - gx1 > gx2 - gpx2 ? (gx1 + gpx2) / 2 : (gpx2 + gx2) / 2;
            }
            out.push(
              '<text class="gname' + (cls === "amap" ? " gname-on" : "") + '" x="' + nx.toFixed(1) + '" y="' + (gy + geneH / 2 + 4) + '" text-anchor="middle">' + esc(g.name) + "</text>"
            );
          } else if (g.nameAbove) {
            // Too narrow for its name, and the gene the figure is about: the
            // name goes over it, kept inside the drawing.
            var ax = Math.max(x0 + tw / 2, Math.min(W - 2 - tw / 2, (gx1 + gx2) / 2));
            out.push('<text class="gname gname-above" x="' + ax.toFixed(1) + '" y="' + (gy - 5) + '" text-anchor="middle">' + esc(g.name) + "</text>");
          }
        }
        out.push("</g>");
      });
      out.push("</g>");
    });

    (spec.edges || []).forEach(function (e) {
      var ex = X(e.at);
      var ey1 = ys[0] - 6;
      var ey2 = y - rowGap + 4;
      out.push("<g" + steps(e) + ">");
      out.push('<line class="edge" x1="' + ex.toFixed(1) + '" x2="' + ex.toFixed(1) + '" y1="' + ey1 + '" y2="' + ey2 + '"/>');
      out.push('<text class="edge-label" x="' + (ex + (e.side === "left" ? -5 : 5)).toFixed(1) + '" y="' + (ey1 - 4) + '" text-anchor="' + (e.side === "left" ? "end" : "start") + '">' + esc(e.label) + "</text>");
      out.push("</g>");
    });

    if (spec.axis) {
      var ay = y - rowGap + 10;
      out.push('<g class="axis" aria-hidden="true"><line x1="' + x0 + '" x2="' + x1 + '" y1="' + ay + '" y2="' + ay + '"/>');
      var step = spec.tick || niceStep(span);
      var label = spec.tickLabel || fmt;
      // Fewer ticks on a narrow drawing, until their labels fit between them.
      var widest = function (st) {
        var w = 0;
        for (var u = Math.ceil(spec.from / st) * st; u <= spec.to; u += st) w = Math.max(w, textWidth(label(u), 10) + 14);
        return w;
      };
      var grow = [2, 2.5, 2];
      for (var g = 0; widest(step) > (step / span) * (x1 - x0) && g < 12; g++) step *= grow[g % 3];
      var first0 = Math.ceil(spec.from / step) * step;
      for (var t = first0; t <= spec.to; t += step) {
        var tx = X(t);
        out.push('<line x1="' + tx.toFixed(1) + '" x2="' + tx.toFixed(1) + '" y1="' + ay + '" y2="' + (ay + 4) + '"/>');
        out.push('<text x="' + tx.toFixed(1) + '" y="' + (ay + 15) + '" text-anchor="middle">' + esc(label(t)) + "</text>");
      }
      out.push("</g>");
    }
    if (masks.length) out.splice(2, 0, "<defs>" + masks.join("") + "</defs>");
    out.push("</svg>");
    spec._used = used;
    return out.join("");
  }

  // Whole positions on a short stretch; megabases only where they still differ.
  function tickLabel(span) {
    return function (t) {
      return t >= 1e6 && span > 20000 ? (t / 1e6).toFixed(3) + " Mb" : fmt(t);
    };
  }

  function niceStep(span) {
    var raw = span / 5;
    var pow = Math.pow(10, Math.floor(Math.log10(raw)));
    var n = raw / pow;
    return (n < 1.5 ? 1 : n < 3.5 ? 2 : n < 7.5 ? 5 : 10) * pow;
  }

  /* ── What can be seen ──

     A layer of the hero waits for its phase, and an element of a story for
     its step; until then it is hidden, and takes neither the pointer nor the
     keyboard. */

  function isShown(el) {
    if (el.closest(".step-off")) return false;
    var svg = el.closest("svg");
    var phase = svg ? +svg.getAttribute("data-phase") : 0;
    if (phase) {
      if (phase < 2 && el.closest(".ph-ref")) return false;
      if (phase < 3 && el.closest(".ph-case")) return false;
      if (phase < 4 && el.closest(".ph-amap")) return false;
    }
    return true;
  }

  function stops(svg) {
    return Array.prototype.slice.call(svg.querySelectorAll(".hot[tabindex], .stn[tabindex]"));
  }

  // Keeps the drawing's one tab stop on something that can be seen.
  function refocus(svg) {
    if (!svg) return;
    var all = stops(svg);
    var current = all.filter(function (el) {
      return el.getAttribute("tabindex") === "0";
    })[0];
    if (current && isShown(current)) return;
    var next = all.filter(isShown).sort(function (a, b) {
      var ra = a.closest(".row");
      var rb = b.closest(".row");
      if (ra !== rb && ra && rb) return ra.compareDocumentPosition(rb) & 4 ? -1 : 1;
      return box(a).x - box(b).x;
    })[0];
    all.forEach(function (el) {
      el.setAttribute("tabindex", el === next ? "0" : "-1");
    });
  }

  function box(el) {
    var shape = el.querySelector(".g, circle, .net-line, .dot, .hub, .end") || el;
    try {
      return shape.getBBox();
    } catch (e) {
      return { x: 0, y: 0, width: 0, height: 0 };
    }
  }

  /* ── The keyboard ──

     Tab reaches a drawing once. Inside it, Left and Right move along a row,
     Home and End go to its ends, and Up and Down go to the gene the row above
     or below aligns with this one: a Bakta call, the protein placed on it,
     the gene AMAP keeps. */

  var KEYS = { ArrowLeft: 1, ArrowRight: 1, ArrowUp: 1, ArrowDown: 1, Home: 1, End: 1 };

  function byX(list) {
    return list.sort(function (a, b) {
      var ba = box(a);
      var bb = box(b);
      return ba.x - bb.x || ba.y - bb.y;
    });
  }

  function partner(from, row) {
    var a = box(from);
    var best = null;
    var score = -Infinity;
    row.forEach(function (el) {
      var b = box(el);
      var shared = Math.min(a.x + a.width, b.x + b.width) - Math.max(a.x, b.x);
      var s = shared > 0 ? shared : -Math.abs(a.x + a.width / 2 - (b.x + b.width / 2));
      if (s > score) {
        score = s;
        best = el;
      }
    });
    return best;
  }

  function moveFrom(el, key, svg) {
    if (svg.classList.contains("trk")) {
      var rows = Array.prototype.slice.call(svg.querySelectorAll("g.row")).filter(isShown).map(function (row) {
        return byX(Array.prototype.slice.call(row.querySelectorAll(".hot[tabindex]")).filter(isShown));
      }).filter(function (row) {
        return row.length;
      });
      var ri = -1;
      var gi = -1;
      rows.forEach(function (row, i) {
        var k = row.indexOf(el);
        if (k !== -1) {
          ri = i;
          gi = k;
        }
      });
      if (ri === -1) return rows.length ? rows[0][0] : null;
      var row = rows[ri];
      if (key === "ArrowRight") return row[Math.min(row.length - 1, gi + 1)];
      if (key === "ArrowLeft") return row[Math.max(0, gi - 1)];
      if (key === "Home") return row[0];
      if (key === "End") return row[row.length - 1];
      if (key === "ArrowDown") return ri + 1 < rows.length ? partner(el, rows[ri + 1]) : el;
      if (key === "ArrowUp") return ri > 0 ? partner(el, rows[ri - 1]) : el;
      return el;
    }
    var items = stops(svg).filter(isShown);
    if (svg.classList.contains("net")) {
      // The proteins, then the lines between them.
      items = items.filter(function (it) {
        return it.classList.contains("net-node");
      }).concat(items.filter(function (it) {
        return !it.classList.contains("net-node");
      }));
    }
    if (svg.classList.contains("metro")) {
      items = byX(items);
      var here = box(el);
      if (key === "ArrowUp" || key === "ArrowDown") {
        var down = key === "ArrowDown";
        var pick = null;
        var dist = Infinity;
        items.forEach(function (other) {
          var b = box(other);
          var dy = b.y + b.height / 2 - (here.y + here.height / 2);
          if (down ? dy < 6 : dy > -6) return;
          var d = Math.abs(dy) + 2 * Math.abs(b.x + b.width / 2 - (here.x + here.width / 2));
          if (d < dist) {
            dist = d;
            pick = other;
          }
        });
        return pick || el;
      }
    }
    var at = items.indexOf(el);
    if (key === "Home") return items[0];
    if (key === "End") return items[items.length - 1];
    if (key === "ArrowRight" || key === "ArrowDown") return items[Math.min(items.length - 1, at + 1)];
    return items[Math.max(0, at - 1)];
  }

  // One listener per figure, on the element that holds its drawing.
  function roving(container) {
    if (container._roving) return;
    container._roving = true;
    container.addEventListener("keydown", function (event) {
      if (!KEYS[event.key] || event.altKey || event.ctrlKey || event.metaKey) return;
      var el = event.target.closest && event.target.closest(".hot[tabindex], .stn[tabindex]");
      var svg = el && el.closest("svg");
      if (!svg || !container.contains(svg)) return;
      event.preventDefault();
      event.stopPropagation();
      var next = moveFrom(el, event.key, svg);
      if (next && next !== el) {
        stops(svg).forEach(function (other) {
          other.setAttribute("tabindex", other === next ? "0" : "-1");
        });
        next.focus();
      }
    });
    container.addEventListener("focusin", function (event) {
      var el = event.target.closest && event.target.closest(".hot[tabindex], .stn[tabindex]");
      if (!el) return;
      // A gene picked with the pointer becomes the drawing's tab stop.
      stops(el.closest("svg")).forEach(function (other) {
        other.setAttribute("tabindex", other === el ? "0" : "-1");
      });
    });
  }

  /* ── Tooltips ──

     With a mouse the tip floats over the drawing, above the gene with its
     caret on it, so the rows under the gene stay in sight. On a touch
     screen, or a drawing too narrow for a tip beside its rows, it is a line
     of text under the drawing instead, which never covers them. Escape hides
     it. */

  // Under this width a floating tip (up to 290 px) would cover the drawing.
  var NARROW = 420;

  var openTips = [];

  document.addEventListener("keydown", function (event) {
    if (event.key === "Escape") {
      openTips.forEach(function (hide) {
        hide(true);
      });
    }
  });

  function attachTips(container, opts) {
    opts = opts || {};
    var what = opts.what || "a gene for its coordinates.";
    var tip = container.querySelector(":scope > .demo-tip");
    if (!tip) {
      tip = document.createElement("div");
      tip.className = "demo-tip";
      tip.hidden = true;
      container.appendChild(tip);
    }
    var readout = container.querySelector(":scope > .fig-readout");
    if (!readout) {
      readout = document.createElement("p");
      readout.className = "fig-readout";
      readout.setAttribute("aria-live", "polite");
      container.appendChild(readout);
    }
    var hint = '<span class="say-hint"><span class="say-hover">Point at</span><span class="say-touch">Tap</span> ' + esc(what) + "</span>";
    readout.innerHTML = hint;
    function inline() {
      return mq("(hover: none)") || container.clientWidth < NARROW;
    }
    // A mouse moving from gene to gene must not make the page under the
    // drawing jump, so there the line keeps the height of the figure's
    // longest tip. A tap is deliberate: the line keeps three lines' room and
    // grows for a longer tip.
    var fitted = -1;
    function fit() {
      var w = readout.offsetWidth;
      if (w === fitted) return;
      fitted = w;
      readout.style.minHeight = "";
      if (!w || mq("(hover: none)")) return;
      var probe = readout.cloneNode(false);
      probe.removeAttribute("aria-live");
      probe.setAttribute("aria-hidden", "true");
      probe.style.cssText = "position:absolute;visibility:hidden;pointer-events:none;left:0;top:0;min-height:0;width:" + w + "px";
      container.appendChild(probe);
      var tallest = 0;
      container.querySelectorAll(".hot[data-tip]").forEach(function (h) {
        probe.innerHTML = h.getAttribute("data-tip");
        tallest = Math.max(tallest, probe.offsetHeight);
      });
      probe.remove();
      if (tallest) readout.style.minHeight = tallest + "px";
    }
    function sync() {
      var on = inline();
      container.classList.toggle("tips-inline", on);
      if (on) fit();
      else fitted = -1;
    }
    sync();
    if ("ResizeObserver" in window) new ResizeObserver(sync).observe(container);
    var shownTarget = null;
    // A tip stays while the pointer is on it, and for a moment after the
    // pointer leaves a gene, so it can be read and pointed at (WCAG 1.4.13).
    var tipTimer = 0;
    function show(target) {
      var text = target.getAttribute("data-tip");
      if (!text) return;
      clearTimeout(tipTimer);
      sync();
      container.querySelectorAll(".hot.on").forEach(function (el) {
        el.classList.remove("on");
      });
      target.classList.add("on");
      shownTarget = target;
      if (inline()) {
        tip.hidden = true;
        readout.innerHTML = text;
        return;
      }
      tip.innerHTML = text;
      tip.hidden = false;
      var bounds = container.getBoundingClientRect();
      var r = (target.querySelector(".g, circle, .net-line, .net-pill") || target).getBoundingClientRect();
      var centre = r.left - bounds.left + r.width / 2;
      var left = centre - tip.offsetWidth / 2;
      left = Math.max(0, Math.min(left, bounds.width - tip.offsetWidth));
      // Above the gene, over the rows above it or past the drawing's top, so
      // the rows under it stay in sight; below only where the window ends.
      var topPos = r.top - bounds.top - tip.offsetHeight - 10;
      var below = bounds.top + topPos < 72;
      if (below) topPos = r.bottom - bounds.top + 10;
      tip.classList.toggle("below", below);
      tip.style.left = left + "px";
      tip.style.top = topPos + "px";
      tip.style.setProperty("--caret", Math.max(10, Math.min(tip.offsetWidth - 10, centre - left)) + "px");
    }
    function hide(clear) {
      tip.hidden = true;
      shownTarget = null;
      container.querySelectorAll(".hot.on").forEach(function (el) {
        el.classList.remove("on");
      });
      if (clear) readout.innerHTML = hint;
    }
    openTips.push(hide);
    // The gene nearest a tap that missed every one, in the same row and
    // within 16 px of it.
    function nearest(event) {
      var best = null;
      var gap = 17;
      container.querySelectorAll(".hot[data-tip]").forEach(function (h) {
        if (!isShown(h)) return;
        var r = (h.querySelector(".g, circle") || h).getBoundingClientRect();
        if (!r.width || event.clientY < r.top - 8 || event.clientY > r.bottom + 8) return;
        var d = event.clientX < r.left ? r.left - event.clientX : event.clientX > r.right ? event.clientX - r.right : 0;
        if (d < gap) {
          gap = d;
          best = h;
        }
      });
      return best;
    }
    container.addEventListener("mouseover", function (event) {
      if (tip.contains(event.target)) {
        clearTimeout(tipTimer);
        return;
      }
      var hot = event.target.closest(".hot");
      if (hot && container.contains(hot) && isShown(hot)) show(hot);
      // Off every gene, the tip goes, unless it is the focused gene's.
      else if (!inline() && shownTarget && document.activeElement !== shownTarget) {
        clearTimeout(tipTimer);
        tipTimer = setTimeout(hide, 200);
      }
    });
    container.addEventListener("mouseleave", function () {
      if (!inline()) hide();
    });
    container.addEventListener("focusin", function (event) {
      var hot = event.target.closest(".hot");
      // A focused gene's name is read already; the line under the drawing
      // would only say it again.
      readout.setAttribute("aria-live", "off");
      if (hot) show(hot);
    });
    container.addEventListener("focusout", function () {
      if (!inline()) hide();
    });
    container.addEventListener("click", function (event) {
      if (tip.contains(event.target)) return;
      // A tap says what it found, there and only there.
      if (event.detail) readout.setAttribute("aria-live", "polite");
      var hot = event.target.closest(".hot");
      if (!hot || !isShown(hot)) hot = nearest(event);
      if (hot) show(hot);
      else hide(true);
    });
  }

  /* ── Legends ──

     Under each drawing, a chip for every kind of mark it uses and no other,
     drawn with the same fills as the drawing. */

  var LEGEND = [
    ["bakta", "bakta", "Bakta's call"],
    ["ref", "ref", null],
    ["short", "ref short", "alignment short of the whole protein"],
    ["amap", "amap", "what AMAP keeps"],
    ["partial", "amap partial", "partial: cut by a contig end"],
    ["pseudo", "pseudo", "pseudogene"],
    ["c1", "c1", "nothing aligns: case 1"],
    ["drop", "drop", "dropped"],
    ["bad", "drop", "a call Bakta is known to make wrongly"],
    ["agree", "band-agree", "same coordinates"],
    ["conflict", "band-conflict", "coordinates differ"],
    ["moved", "moved", "an end AMAP moved"],
    ["diff", "diff", "a marked stretch"],
    ["diff-ref", "diff-ref", "a reference gene with no call"]
  ];

  function legend(used, reference, leave) {
    return LEGEND.filter(function (l) {
      return used && used[l[0]] && !(leave && leave.indexOf(l[0]) !== -1);
    }).map(function (l) {
      var text = l[2] || (reference || "reference") + " protein, aligned";
      return '<span data-key="' + l[0] + '"><i class="chip ' + l[1] + '"></i>' + esc(text) + "</span>";
    }).join("");
  }

  // Draws spec into stage at the stage's width, and redraws when it changes.
  function render(stage, specFn, after) {
    var last = 0;
    function go() {
      var width = stage.clientWidth;
      if (!width || Math.abs(width - last) < 8) return;
      last = width;
      // The drawing goes in a holder of its own, so the tooltip next to it
      // survives a redraw.
      var holder = stage.querySelector(":scope > .fig-svg");
      if (!holder) {
        holder = document.createElement("div");
        holder.className = "fig-svg";
        stage.insertBefore(holder, stage.firstChild);
      }
      var spec = specFn();
      holder.innerHTML = draw(spec, width);
      if (after) after(holder.querySelector("svg"), spec);
    }
    go();
    roving(stage);
    if ("ResizeObserver" in window) {
      new ResizeObserver(go).observe(stage);
    } else {
      window.addEventListener("resize", go);
    }
    return function () {
      last = 0;
      go();
    };
  }

  /* ── The real genes ──

     The data are in figure-data.js: windows of a real run, with every call
     and placement in them. What follows only decides how each one is drawn
     and worded. */

  var DATA = window.AMAP_FIGURE_DATA || { windows: {} };

  // Where each drawing sits and what it points out.
  var LOOK = {
    hero: { tick: 2000, kbTicks: true, dropMark: true, alt: "The first 12.4 kb of M. tuberculosis CDC1551: twelve Bakta calls, eight H37Rv proteins aligned, and the eight genes AMAP keeps." },
    case2: { from: 1960, to: 3300, tick: 500, alt: "dnaN: Bakta's call and the aligned H37Rv protein start and end at the same bases; AMAP keeps Bakta's call." },
    case3: { from: 23700, to: 25600, tick: 500, alt: "fhaA: the same ends as the aligned H37Rv protein, which has six residues more; AMAP keeps Bakta's call.", gaps: { Rv0020c: { at: 24706, label: "6 residues" } } },
    "case4-start": { from: -60, to: 1640, tick: 500, alt: "dnaA: Bakta's call starts at base 34, the aligned H37Rv protein at base 1; AMAP keeps the reference's start.", marks: [{ s: 1, e: 33, rows: [0, 1], label: "33 bp", labelAt: "right" }] },
    "case4-split": { from: 9850, to: 10950, tick: 250, alt: "Rv0007: two Bakta calls where one H37Rv protein aligns end to end; AMAP keeps one gene." },
    "case4-longer": { from: 2267500, to: 2272800, tick: 1000, alt: "Rv2024c: Bakta's call is 3,303 bp longer than an alignment that stops before the end of its protein; AMAP keeps Bakta's call, named Rv2024c.", marks: [{ s: 2267753, e: 2271055, rows: [0, 1], label: "3,303 bp" }] },
    "case1-kept": { from: 1979300, to: 1985800, tick: 2000, alt: "Four Bakta calls in a stretch CDC1551 has and H37Rv does not; AMAP keeps all four." },
    "case1-dropped": { from: 11330, to: 11620, tick: 100, alt: "A Bakta call covered end to end by one of Bakta's known false positives; AMAP drops it." },
    "small-protein": { from: 3056150, to: 3058450, tick: 500, alt: "E. coli K-12: the exact search places IbsC, 19 residues, between fau and serA, where Bakta called nothing.", marks: [{ s: 3056890, e: 3056949, rows: ["bakta", "ref"], label: "IbsC, by the exact search", kind: "ref" }] },
    "paralogs-1": { from: 400150, to: 402300, tick: 500, alt: "CDC1551 near 400 kb: Rv0336 keeps its own locus, and Rv0515's placement there is dropped." },
    "paralogs-2": { from: 607300, to: 609950, tick: 500, alt: "CDC1551 near 608 kb: Rv0515 keeps its own locus, and Rv0336's placement there is dropped." },
    "contig-end": { from: 100250, to: 102146, tick: 500, alt: "CDC1551 in 40 contigs: Rv1239c runs to the end of contig_10 with no stop codon and this piece is marked partial at its 3' end.", edges: [{ at: 102146.5, label: "end of contig_10", side: "left" }] },
    "contig-end-2": { from: 3700, to: 4591, tick: 250, alt: "CDC1551 in 40 contigs: the other piece of Rv1239c, its last 144 residues, at the very end of contig_9, kept over a short Bakta call and marked partial at its 5' end.", edges: [{ at: 4591.5, label: "end of contig_9", side: "left" }] },
    frameshift: { from: 688900, to: 690500, tick: 500, alt: "M. bovis: Rv0590's frame runs on through Rv0590A; AMAP keeps it at its extent as a pseudogene, and Rv0590A as a gene of its own.", marks: [{ s: 690104, e: 690306, rows: ["amap", "amap"], label: "67 residues past its protein" }] },
    "copies-1": { from: 616700, to: 621600, tick: 1000, labelIds: true, alt: "E. coli O157:H7 Sakai: the first of two whole RhsD copies, kept as b0497." },
    "copies-2": { from: 2045700, to: 2050700, tick: 1000, labelIds: true, alt: "The second RhsD copy, 1.4 Mb away, kept as b0497_2." },
    missing: { from: 33100, to: 34420, tick: 250, alt: "Rv0031: an H37Rv protein aligns between two Bakta calls where Bakta called nothing; AMAP puts the gene back.", marks: [{ s: 33526, e: 33738, rows: [0, 1], label: "no call", kind: "ref" }] }
  };

  var CASE_SAYS = {
    1: "Case 1: no reference protein aligns here.",
    2: "Case 2: same start, end, strand and length as the aligned protein.",
    3: "Case 3: same start, end and strand, but the alignment's length differs.",
    4: "Case 4: overlaps an aligned protein, with other coordinates."
  };

  function span(g) {
    return fmt(g.s) + "-" + fmt(g.e) + " (" + (g.st === "-" ? "−" : "+") + ")";
  }

  function short(tag) {
    return String(tag).replace(/^.*_(\d{5})$/, "$1");
  }

  function overlap(a, b) {
    return Math.max(0, Math.min(a.e, b.e) - Math.max(a.s, b.s) + 1);
  }

  // Why a case 1 call went, in the words of the rule that took it out.
  function dropReason(call, win) {
    var len = call.e - call.s + 1;
    var covered = 0;
    (win.bad || []).forEach(function (f) {
      covered += overlap(call, f);
    });
    if (covered / len >= 0.8) return "one of the calls Bakta makes wrongly on " + refName(win) + " itself covers it";
    var best = 0;
    (win.refMatch || []).forEach(function (m) {
      if (m.tag === call.tag && m.num > 0) best = Math.max(best, m.match / m.num);
    });
    return "its protein matches the " + refName(win) + " genome at " + Math.round(best * 100) + "% identity, so it is not new";
  }

  // The reference a window was corrected against.
  function refName(win) {
    return /coli/.test(win.genome) ? "K-12" : "H37Rv";
  }

  // Where a kept gene's ends lie against the one Bakta call it overlaps on
  // its strand: the ends AMAP moved, and by how much.
  function moves(a, win) {
    var calls = win.bakta.filter(function (c) {
      return c.st === a.st && overlap(c, a) > 30;
    });
    if (calls.length !== 1) return { at: [], words: [] };
    var c = calls[0];
    var minus = a.st === "-";
    var out = { at: [], words: [] };
    // Positive: the gene starts earlier, or ends later, than the call.
    var start = minus ? a.e - c.e : c.s - a.s;
    var end = minus ? c.s - a.s : a.e - c.e;
    if (start) {
      out.at.push(minus ? a.e : a.s);
      out.words.push("starts " + fmt(Math.abs(start)) + " bp " + (start > 0 ? "earlier" : "later"));
    }
    if (end) {
      out.at.push(minus ? a.s : a.e);
      out.words.push("ends " + fmt(Math.abs(end)) + " bp " + (end > 0 ? "later" : "earlier"));
    }
    if (out.words.length) out.words = [out.words.join(" and ") + " than Bakta's call " + short(c.tag)];
    return out;
  }

  function rowsFor(key) {
    var win = DATA.windows[key];
    var look = LOOK[key] || {};
    var reference = refName(win);
    var names = {};
    win.ref.forEach(function (r) {
      names[r.prot] = r.name;
    });
    var kept = {};
    win.amap.forEach(function (a) {
      kept[a.id] = a;
    });
    var same = function (a, b) {
      return a.s === b.s && a.e === b.e && a.st === b.st;
    };
    var bakta = win.bakta.map(function (c) {
      var gone = c.case === 1 && !Object.keys(kept).some(function (id) {
        return id === c.tag;
      });
      var tip = "<b>Bakta call " + short(c.tag) + (c.name ? " (" + c.name + ")" : "") + "</b>" + (c.product ? ", " + esc(c.product) : "") +
        "<br>" + span(c) + ", " + fmt(c.e - c.s + 1) + " bp<br>" + CASE_SAYS[c.case];
      if (c.case === 1) tip += gone ? " Dropped: " + dropReason(c, win) + "." : " Kept: nothing rules it out.";
      return {
        s: c.s, e: c.e, st: c.st, name: c.name, tip: tip, phaseCls: c.case === 1 ? (gone ? "ph-dim c1-call" : "c1-call") : "",
        k: c.case === 2 || c.case === 3 ? "agree" : gone ? "drop" : ""
      };
    });
    var ref = win.ref.map(function (r) {
      // The piece of a gene a contig end cut (the second pass's, in the one
      // window that has one) holds the protein's last residues and lacks its
      // start; any other alignment that is not the whole protein stops short
      // of its end.
      var whole = r.complete ? "the whole protein, to its stop" : r.piece ? "a piece: the protein's last residues, without its start" : "stops short of its end: not the whole protein";
      var gap = look.gaps && look.gaps[r.prot];
      return {
        s: r.s, e: r.e, st: r.st, name: r.prot, gap: gap ? gap.at : null, gapLabel: gap ? gap.label : null, short: r.complete ? "" : r.piece ? "5'" : "3'",
        k: win.bakta.some(function (c) {
          return same(c, r);
        }) ? "agree" : win.amap.some(function (a) {
          return same(a, r);
        }) ? "fixed" : "",
        tip: "<b>" + r.prot + (r.name !== r.prot ? " (" + r.name + ")" : "") + "</b>, " + r.len + " residues<br>aligned at " + span(r) +
          ": " + whole + (r.num !== r.e - r.s + 1 ? "; " + fmt(r.num) + " nt aligned on " + fmt(r.e - r.s + 1) + " bp" : "")
      };
    });
    function keptGene(a) {
      var fromCall = win.bakta.some(function (c) {
        return same(c, a);
      });
      var fromRef = win.ref.some(function (r) {
        return same(r, a);
      });
      var how = fromCall && fromRef ? "Bakta's call, as it was" : fromCall ? "Bakta's call" : fromRef ? "the aligned protein's coordinates" : "extended to its stop codon";
      var moved = fromCall ? { at: [], words: [] } : moves(a, win);
      var marks = (a.pseudo ? "<br>Pseudogene: " + esc(a.note) : "") + (a.partial ? "<br>Partial at its " + a.partial + " end" : "") +
        (moved.words.length ? "<br>It " + moved.words[0] + "." : "");
      var name = a.rvs ? a.rvs.split(",")[0] : "";
      return {
        s: a.s, e: a.e, st: a.st, name: name, cls: a.pseudo ? "pseudo" : null, open: a.partial || "", partial: !!a.partial, moved: moved.at,
        k: fromCall && fromRef ? "agree" : fromRef ? "fixed" : "", id: a.id,
        tip: "<b>" + a.id + "</b>" + (a.rvs ? ", named " + a.rvs : ", with no reference name") + "<br>" + span(a) + ": " + how + marks
      };
    }
    return {
      win: win, look: look, reference: reference, bakta: bakta, ref: ref,
      amap: win.amap.map(keptGene),
      corrected: (win.corrected || win.amap).map(function (c) {
        return keptGene({ id: c.id, rvs: c.rvs, s: c.s, e: c.e, st: c.st });
      }),
      bad: (win.bad || []).map(function (f) {
        return { s: f.s, e: f.e, st: f.st, lk: "bad", tip: "A call Bakta makes on " + reference + " itself that the curated annotation does not have, aligned here at " + span(f) + "." };
      })
    };
  }

  // The genes a figure is about (those under its marks) keep their name even
  // where they are too narrow for it: it goes over them.
  function nameSubjects(spec) {
    (spec.marks || []).forEach(function (m) {
      var lo = Math.min(m.rows[0], m.rows[1]);
      var hi = Math.max(m.rows[0], m.rows[1]);
      spec.rows.forEach(function (row, i) {
        row.genes.forEach(function (g) {
          if (!g.name) return;
          // In the rows the mark spans, a gene it covers at least half of
          // (or that covers half of it); in AMAP's row, the gene it frames.
          var inside = i >= lo && i <= hi && overlap(g, m) >= Math.min(g.e - g.s + 1, m.e - m.s + 1) / 2;
          var framed = row.cls === "amap" && g.s >= m.s && g.e <= m.e;
          if (inside || framed) g.nameAbove = true;
        });
      });
    });
    return spec;
  }

  function windowSpec(key) {
    var parts = rowsFor(key);
    var win = parts.win;
    var look = parts.look;
    var bakta = parts.bakta;
    var ref = parts.ref;
    var amap = parts.amap;
    var rows = [
      { label: "Bakta", cls: "bakta", genes: bakta, empty: bakta.length ? null : "no call" },
      { label: parts.reference + ",\naligned", cls: "ref", genes: ref, empty: ref.length ? null : "no " + parts.reference + " protein aligns here", phase: key === "hero" ? "ref" : null }
    ];
    if (win.placed) {
      // Placements the paralog step dropped: on a locus another protein owns.
      var stays = {};
      win.ref.forEach(function (r) {
        stays[r.prot + ":" + r.s] = true;
      });
      var dropped = win.placed.filter(function (q) {
        return !stays[q.prot + ":" + q.s];
      });
      if (dropped.length) {
        rows.push({
          label: "Dropped:\nanother's locus", cls: "drop", genes: dropped.map(function (q) {
            return { s: q.s, e: q.e, st: q.st, name: q.prot, tip: "<b>" + q.prot + "</b> placed at " + span(q) + ": not where " + q.prot + " belongs, and on a locus another protein owns. Dropped." };
          })
        });
      }
    }
    if (parts.bad.length && key !== "hero") {
      rows.push({ label: "Known false\npositives", cls: "drop", genes: parts.bad });
    }
    if (look.labelIds) {
      amap = amap.map(function (g, i) {
        g.name = win.amap[i].id.replace(win.base + "_", "");
        return g;
      });
    }
    rows.push({ label: "AMAP", cls: "amap", genes: amap, empty: amap.length ? null : "nothing kept", phase: key === "hero" ? "amap" : null, space: key === "hero" ? 6 : 0 });
    // A band from each call to each aligned protein it overlaps on its strand:
    // straight where the coordinates agree, slanted where they do not.
    var ribbons = [];
    win.bakta.forEach(function (c, i) {
      win.ref.forEach(function (r, j) {
        if (r.st !== c.st || overlap(c, r) <= 30) return;
        var agree = r.s === c.s && r.e === c.e;
        ribbons.push({ a: [0, i], b: [1, j], kind: agree ? "agree" : "conflict", k: agree ? "agree" : "fixed", phase: key === "hero" ? "case" : null });
      });
    });
    var marks = (look.marks || []).map(function (m) {
      var named = { bakta: 0, ref: 1, amap: rows.length - 1 };
      var copy = {};
      Object.keys(m).forEach(function (k) {
        copy[k] = m[k];
      });
      copy.rows = m.rows.map(function (r) {
        return typeof r === "string" ? named[r] : r;
      });
      return copy;
    });
    if (look.dropMark) {
      // The calls the hero drops, marked once AMAP has kept what it keeps.
      var gone = win.bakta.filter(function (c, i) {
        return bakta[i].k === "drop";
      });
      if (gone.length) {
        marks.push({
          s: Math.min.apply(null, gone.map(function (c) {
            return c.s;
          })),
          e: Math.max.apply(null, gone.map(function (c) {
            return c.e;
          })),
          rows: [0, 0], label: gone.length + " dropped", phase: "amap"
        });
      }
    }
    var hero = key === "hero";
    return nameSubjects({
      from: look.from != null ? look.from : win.from,
      to: look.to != null ? look.to : win.to,
      labelW: hero ? 66 : 64,
      geneH: hero ? 20 : 18,
      rowGap: hero ? 30 : 22,
      top: hero ? 18 : 18,
      alt: look.alt || "",
      refName: parts.reference,
      phased: true,
      rows: rows,
      ribbons: ribbons,
      marks: marks,
      edges: look.edges || [],
      axis: true,
      tick: look.tick,
      tickLabel: look.kbTicks ? function (t) {
        return t / 1000 + " kb";
      } : tickLabel(win.to - win.from)
    });
  }

  /* ── The hero: the first 12.4 kb of CDC1551, in four phases ──

     It plays once, when it first comes into view, and holds still while the
     reader points at it. Four segments show where it is and jump to a phase;
     the button pauses it, or plays it again. */

  var PHASES = [
    { label: "What Bakta called" },
    { label: "Where the H37Rv proteins align" },
    { label: "Compared, call by call" },
    { label: "What AMAP keeps" }
  ];

  // How long each phase stays before the next; the comparison is the densest.
  var DWELL = [2100, 2100, 3200];

  function countUp(el, n, animate) {
    if (!animate || reduceMotion || n < 2) {
      el.textContent = String(n);
      return;
    }
    var t0 = 0;
    function frame(t) {
      if (!t0) t0 = t;
      var f = Math.min(1, (t - t0) / 600);
      el.textContent = String(Math.round(f * n));
      if (f < 1) requestAnimationFrame(frame);
    }
    el.textContent = "0";
    requestAnimationFrame(frame);
  }

  function mountHero(fig) {
    var stage = fig.querySelector("[data-hero-stage]");
    var stepEl = fig.querySelector("[data-hero-step]");
    var button = fig.querySelector("[data-hero-replay]");
    var win = DATA.windows.hero;
    if (!win || !stage) return;
    var phase = 4;
    var svg = null;
    var lit = null;
    render(stage, function () {
      return windowSpec("hero");
    }, function (drawn, spec) {
      svg = drawn;
      svg.setAttribute("data-phase", String(phase));
      if (lit) svg.setAttribute("data-light", lit);
      refocus(svg);
      var key = fig.querySelector(".legend");
      if (key && !key._built) {
        key._built = true;
        // Its one marked stretch names itself ("3 dropped").
        key.innerHTML = legend(spec._used, spec.refName, ["diff"]);
      }
    });
    attachTips(stage);

    // The tally counts the window's own calls: those that agree (case 2 and 3),
    // the genes that took the reference's coordinates, and the calls dropped.
    var keptIds = {};
    win.amap.forEach(function (a) {
      keptIds[a.id] = true;
    });
    var counts = {
      agree: win.bakta.filter(function (c) {
        return c.case === 2 || c.case === 3;
      }).length,
      fixed: win.amap.filter(function (a) {
        return !win.bakta.some(function (c) {
          return c.s === a.s && c.e === a.e && c.st === a.st;
        });
      }).length,
      drop: win.bakta.filter(function (c) {
        return c.case === 1 && !keptIds[c.tag];
      }).length
    };
    var tallies = {};
    var items = [];
    // The count pressed to stay lit, or null.
    var pinned = null;
    function show(key) {
      lit = key;
      items.forEach(function (it) {
        it.el.classList.toggle("lit", it.key === key);
        it.el.setAttribute("aria-pressed", it.key === pinned ? "true" : "false");
      });
      if (svg) {
        if (lit) svg.setAttribute("data-light", lit);
        else svg.removeAttribute("data-light");
      }
    }
    fig.querySelectorAll("[data-tally]").forEach(function (el) {
      tallies[el.getAttribute("data-tally")] = el;
      el._on = null;
      // Pointing at a count lights the genes and bands it counts; pressing
      // it keeps them lit until it is pressed again.
      var item = el.parentNode;
      var key = el.getAttribute("data-tally");
      items.push({ el: item, key: key });
      item.setAttribute("tabindex", "0");
      item.setAttribute("role", "button");
      item.setAttribute("aria-pressed", "false");
      var toggle = function () {
        pinned = pinned === key ? null : key;
        show(pinned || key);
      };
      item.addEventListener("pointerenter", function () {
        show(key);
      });
      item.addEventListener("pointerleave", function () {
        show(pinned);
      });
      item.addEventListener("focus", function () {
        show(key);
      });
      item.addEventListener("blur", function () {
        show(pinned);
      });
      item.addEventListener("click", toggle);
      item.addEventListener("keydown", function (event) {
        if (event.key !== "Enter" && event.key !== " ") return;
        event.preventDefault();
        toggle();
      });
    });

    // The four segments, built into the step line, and its words.
    var segs = [];
    var words = null;
    if (stepEl) {
      stepEl.setAttribute("aria-live", "off");
      stepEl.textContent = "";
      var group = document.createElement("span");
      group.className = "hero-segs";
      group.setAttribute("role", "group");
      group.setAttribute("aria-label", "The four phases");
      PHASES.forEach(function (p, i) {
        var b = document.createElement("button");
        b.type = "button";
        b.className = "hero-seg";
        b.setAttribute("data-go", String(i + 1));
        b.setAttribute("aria-label", "Phase " + (i + 1) + " of 4: " + p.label);
        b.innerHTML = "<i></i>";
        group.appendChild(b);
        segs.push(b);
      });
      words = document.createElement("span");
      words.className = "hero-phase";
      stepEl.appendChild(group);
      stepEl.appendChild(words);
    }

    // The line keeps the height of its longest phase, so a phase whose words
    // wrap at this width does not push the drawing down while it plays.
    var stepWidth = 0;
    function fitStep() {
      if (!stepEl || !words) return;
      var w = stepEl.offsetWidth;
      if (!w || Math.abs(w - stepWidth) < 2) return;
      stepWidth = w;
      var probe = stepEl.cloneNode(true);
      probe.removeAttribute("aria-live");
      probe.removeAttribute("data-hero-step");
      probe.setAttribute("aria-hidden", "true");
      probe.style.cssText = "position:absolute;visibility:hidden;pointer-events:none;left:0;top:0;min-height:0;width:" + w + "px";
      stepEl.parentNode.appendChild(probe);
      var said = probe.querySelector(".hero-phase");
      var tallest = 0;
      PHASES.forEach(function (p, i) {
        said.textContent = i + 1 + " of 4 · " + p.label;
        tallest = Math.max(tallest, probe.offsetHeight);
      });
      probe.remove();
      stepEl.style.minHeight = tallest + "px";
    }
    fitStep();
    if (stepEl && "ResizeObserver" in window) new ResizeObserver(fitStep).observe(stepEl);

    var auto = false;
    var userPaused = false;
    var hover = false;
    var focus = false;
    var timer = null;
    var started = 0;
    var left = 0;
    var played = false;
    var inView = false;

    // After the reader's first move, and not before, the phases are spoken.
    function spoke() {
      if (stepEl) stepEl.setAttribute("aria-live", "polite");
    }

    function setPhase(n) {
      phase = n;
      if (svg) {
        svg.setAttribute("data-phase", String(n));
        refocus(svg);
      }
      if (words) words.textContent = n + " of 4 · " + PHASES[n - 1].label;
      segs.forEach(function (b, i) {
        b.classList.toggle("done", i + 1 < n);
        b.classList.remove("now");
        if (i + 1 === n) {
          // Reflow, so the fill starts again from empty.
          void b.offsetWidth;
          b.classList.add("now");
          b.style.setProperty("--dwell", (DWELL[i] || 0) + "ms");
          b.setAttribute("aria-current", "step");
        } else {
          b.removeAttribute("aria-current");
        }
      });
      // 'agree' once the calls are compared; the rest once AMAP has kept.
      Object.keys(tallies).forEach(function (k) {
        var on = k === "agree" ? n >= 3 : n >= 4;
        var el = tallies[k];
        el.parentNode.classList.toggle("t-off", !on);
        // A count not yet made shows no number, so nothing reads as a count
        // that is not true yet.
        if (on && el._on !== true) countUp(el, counts[k], el._on !== null);
        else if (!on) el.textContent = "\u2013";
        if (on) el.removeAttribute("aria-hidden");
        else el.setAttribute("aria-hidden", "true");
        el._on = on;
      });
    }

    function held() {
      return hover || focus || userPaused;
    }
    function arm(ms) {
      left = ms;
      if (held()) {
        fig.classList.add("paused");
        return;
      }
      started = Date.now();
      timer = setTimeout(next, ms);
    }
    function next() {
      timer = null;
      setPhase(phase + 1);
      if (phase < 4) arm(DWELL[phase - 1]);
      else finish();
    }
    function finish() {
      auto = false;
      userPaused = false;
      fig.classList.remove("running", "paused");
      label();
    }
    function freeze() {
      if (timer) {
        clearTimeout(timer);
        timer = null;
        left = Math.max(0, left - (Date.now() - started));
      }
      if (auto) fig.classList.add("paused");
    }
    function thaw() {
      if (!auto || held()) return;
      fig.classList.remove("paused");
      if (!timer) {
        started = Date.now();
        timer = setTimeout(next, left);
      }
    }
    function start(from) {
      clearTimeout(timer);
      timer = null;
      auto = true;
      played = true;
      userPaused = false;
      fig.classList.remove("paused");
      fig.classList.add("running");
      setPhase(from);
      arm(DWELL[from - 1]);
      label();
    }
    function stopAuto() {
      clearTimeout(timer);
      timer = null;
      auto = false;
      userPaused = false;
      fig.classList.remove("running", "paused");
      label();
    }
    function label() {
      if (!button) return;
      var playing = auto && !userPaused;
      // Both words stay in the button and one shows, so it keeps the width
      // of the longer and the line beside it never rewraps.
      if (!button.querySelector(".say-2")) button.innerHTML = '<span class="say-1">Play</span><span class="say-2">Pause</span>';
      button.classList.toggle("said-2", playing);
      button.setAttribute("aria-label", playing ? "Pause the drawing" : "Play the drawing");
    }

    if (button) {
      button.addEventListener("click", function () {
        spoke();
        if (auto && !userPaused) {
          userPaused = true;
          freeze();
        } else if (auto) {
          userPaused = false;
          thaw();
        } else {
          start(phase < 4 ? phase : 1);
        }
        label();
      });
    }
    segs.forEach(function (b) {
      b.addEventListener("click", function () {
        spoke();
        played = true;
        stopAuto();
        setPhase(+b.getAttribute("data-go"));
      });
    });
    // Still while the reader points at the drawing or moves through it.
    stage.addEventListener("pointerenter", function (event) {
      if (event.pointerType !== "mouse") return;
      hover = true;
      freeze();
    });
    stage.addEventListener("pointerleave", function (event) {
      if (event.pointerType !== "mouse") return;
      hover = false;
      thaw();
    });
    stage.addEventListener("focusin", function () {
      focus = true;
      freeze();
    });
    stage.addEventListener("focusout", function (event) {
      if (event.relatedTarget && stage.contains(event.relatedTarget)) return;
      focus = false;
      thaw();
    });

    // It plays when it first comes into view, however late. A page that is
    // not being shown (a background tab) shows the last phase meanwhile.
    function maybeStart() {
      if (!played && inView && !document.hidden) start(1);
    }
    if ("IntersectionObserver" in window && !reduceMotion) {
      setPhase(document.hidden ? 4 : 1);
      label();
      var io = new IntersectionObserver(function (entries) {
        inView = entries[entries.length - 1].isIntersecting;
        maybeStart();
        if (played) io.disconnect();
      }, { threshold: 0.35 });
      io.observe(fig);
      document.addEventListener("visibilitychange", function () {
        if (played) return;
        if (document.hidden) {
          setPhase(4);
        } else {
          setPhase(1);
          maybeStart();
        }
      });
    } else {
      setPhase(4);
      label();
    }
  }

  /* ── One drawing per case ── */

  // The legend under a figure of the docs: a chip for each kind of mark its
  // drawing uses.
  function legendAfter(stage) {
    var el = document.createElement("div");
    el.className = "legend fig-legend";
    stage.parentNode.insertBefore(el, stage.nextSibling);
    return el;
  }

  function mountCase(el) {
    var names = el.getAttribute("data-fig").split(" ");
    var stage = el.querySelector("[data-fig-stage]") || el;
    var switcher = el.querySelector("[data-fig-switch]");
    var caption = el.querySelector("[data-fig-caption]");
    var current = names[0];
    // Every drawing the switcher can show, so the figure keeps the height of
    // the tallest and nothing under it moves on a switch.
    var variants = switcher ? Array.prototype.map.call(switcher.querySelectorAll("button[data-show]"), function (b) {
      return b.getAttribute("data-show");
    }).filter(function (n) {
      return DATA.windows[n];
    }) : [];
    var key = el.tagName === "FIGURE" && stage !== el ? legendAfter(stage) : null;
    var redraw = render(stage, function () {
      var spec = windowSpec(current);
      spec.final = true;
      return spec;
    }, function (svg, spec) {
      // A figure shows the end of the story: what AMAP kept and dropped.
      svg.setAttribute("data-phase", "4");
      if (variants.length > 1) {
        var w = svg.viewBox.baseVal.width;
        var tallest = 0;
        variants.forEach(function (n) {
          var s = windowSpec(n);
          draw(s, w);
          tallest = Math.max(tallest, s._h);
        });
        svg.parentNode.style.minHeight = (tallest * svg.getBoundingClientRect().width) / w + "px";
      }
      if (key) key.innerHTML = legend(spec._used, spec.refName);
    });
    attachTips(stage);
    function show(name) {
      if (!DATA.windows[name]) return;
      current = name;
      redraw();
      if (switcher) {
        switcher.querySelectorAll("button").forEach(function (b) {
          b.setAttribute("aria-pressed", b.getAttribute("data-show") === name ? "true" : "false");
        });
      }
      if (caption) {
        el.querySelectorAll("[data-for]").forEach(function (p) {
          p.hidden = p.getAttribute("data-for") !== name;
        });
      }
    }
    if (switcher) {
      switcher.addEventListener("click", function (event) {
        var b = event.target.closest("button[data-show]");
        if (b) show(b.getAttribute("data-show"));
      });
    }
    show(names[0]);
  }

  /* ── The metro map ──

     The stations are main.nf's processes in the order its workflow calls
     them; check_site.py compares the two. Lines are what travels: the sample
     (main), the reference, the eggNOG database, the duplication analysis and
     the pangenome. */

  var STATIONS = {
    SHEET: {
      name: "Sample sheet", proc: "--tsv", x: 70, y: 170, kind: "end", label: "below",
      text: "Two tab-separated columns, sampleId and fasta: one assembled genome per row.",
      reads: "your assemblies (FASTA)", writes: "", flags: "--tsv", doc: "getting-started/inputs.html"
    },
    VALIDATE_INPUTS: {
      name: "Check inputs", proc: "VALIDATE_INPUTS", x: 160, y: 170, label: "above",
      text: "Checks the sample sheet and every FASTA before anything long starts. With --validate_only the run stops here.",
      reads: "the sample sheet", writes: "pipeline_info/validation.log", flags: "--validate_only", doc: "getting-started/inputs.html", exit: "pipeline_info/", exitSide: "below"
    },
    PATHOTYPR: {
      name: "Lineage", proc: "PATHOTYPR", x: 250, y: 170, label: "above",
      text: "pathotypr assigns each genome its MTBC lineage and renames it (L4_G26689). MTBC only.",
      reads: "each genome", writes: "PATHOTYPR/", flags: "--skip_typing turns it off", doc: "how-it-works/workflow.html#2-lineage-typing", exit: "PATHOTYPR/", exitSide: "above", typing: true
    },
    BAKTA_ANNOTATOR: {
      name: "Bakta", proc: "BAKTA_ANNOTATOR", x: 345, y: 170, label: "above",
      text: "Bakta calls the genes: the annotation AMAP then corrects, published untouched for comparison.",
      reads: "each genome, --genus and --species", writes: "BAKTA/", flags: "--genus, --species, --bakta_db", doc: "how-it-works/workflow.html#3-annotation-with-bakta", exit: "BAKTA/", exitSide: "below"
    },
    MINIPROT: {
      name: "miniprot", proc: "MINIPROT", x: 440, y: 170, label: "below", hub: true,
      text: "miniprot aligns the reference proteome onto each genome, plus the false-positive set, plus Bakta's proteins back onto the reference genome; then, for each placement a contig end cut, what it lacks, to find the gene's other piece.",
      reads: "each genome, Bakta's proteins, the reference", writes: "nothing published: the alignments go on to Classify", flags: "", doc: "how-it-works/workflow.html#4-protein-alignment-with-miniprot"
    },
    TYPE_CASES: {
      name: "Cases", proc: "TYPE_CASES", x: 535, y: 170, label: "above", hub: true,
      text: "Every Bakta call sorted into case 1 to 4, plus the reference genes Bakta missed; conflicts resolved and missing genes put back.",
      reads: "the alignments, Bakta's calls, the genome", writes: "CDS_CASES/", flags: "", doc: "how-it-works/correction-algorithm.html", exit: "CDS_CASES/", exitSide: "below"
    },
    FIX_CDS: {
      name: "Correction", proc: "FIX_CDS", x: 625, y: 170, label: "below", hub: true,
      text: "Translates every gene: extends to a stop codon, flags pseudogenes, marks genes a contig end cut, keeps the real copies of a gene.",
      reads: "the corrected coordinates, the genome", writes: "AMAP/: *.AMAP.tsv, .faa, .fna", flags: "", doc: "how-it-works/correction-algorithm.html#7-stop-codon-check", exit: "AMAP/", exitSide: "above"
    },
    EGGNOG_DB: {
      name: "eggNOG database", proc: "EGGNOG_DB", x: 735, y: 60, kind: "end", label: "left",
      text: "Builds the eggNOG database on the first run: eggnog.db, 41.4 GB unpacked from a 6.8 GB download, and a bacteria-only Diamond index of about 3 GB. It reuses them after.",
      reads: "", writes: "the folder --eggnog_db points at", flags: "--eggnog_db; skipped with --skip_functional_annotation", doc: "how-it-works/workflow.html#1-eggnog-database-setup", functional: true
    },
    ANNOT_EGGNOGMAP_BATCH: {
      name: "eggNOG-mapper", proc: "ANNOT_EGGNOGMAP_BATCH", x: 735, y: 170, label: "below", hub: true, hubUp: 9, hubDown: 27,
      text: "eggNOG-mapper annotates the corrected proteins of all samples in one go: COG, KEGG, GO, EC, PFAM, CAZy.",
      reads: "every sample's corrected proteins", writes: "nothing published: the annotations go on", flags: "--skip_functional_annotation stops before this", doc: "how-it-works/workflow.html#7-functional-annotation", functional: true
    },
    GET_OUTPUTS: {
      name: "GFF3, GenBank", proc: "GET_OUTPUTS", x: 860, y: 170, label: "above",
      text: "Writes the corrected annotation: GFF3, GenBank, the master table, the clean proteome and the QC table.",
      reads: "the corrected genes and their functions", writes: "AMAP/: *-AMAP.gff3, .gbk, mastertable, .clean.faa, .qc.tsv", flags: "--skip_mycobrowser", doc: "user-guide/outputs.html#the-key-files", exit: "AMAP/", exitSide: "below", functional: true
    },
    GET_COUNT: {
      name: "Counts", proc: "GET_COUNT", x: 955, y: 170, label: "above",
      text: "Counts the COG and Mycobrowser categories of each sample.",
      reads: "the master table", writes: "AMAP/: *-AMAP.COG.tsv, *-AMAP.MycoBrowser.tsv", flags: "", doc: "user-guide/outputs.html#the-full-tree", functional: true
    },
    GET_DUPLICATIONS: {
      name: "Duplications", proc: "GET_DUPLICATIONS", x: 890, y: 250, label: "above", lineOnly: "dup",
      text: "Each proteome against itself with Diamond: duplicated genes grouped into families.",
      reads: "the corrected proteins and their functions", writes: "DUPLICATIONS/{sample}/", flags: "", doc: "how-it-works/duplications.html", exit: "DUPLICATIONS/", exitSide: "below", functional: true
    },
    GET_ALLDUP: {
      name: "All genomes", proc: "GET_ALLDUP", x: 980, y: 250, label: "above", labelDx: 10, lineOnly: "dup",
      text: "All samples' duplicated genes against each other, as communities with network metrics.",
      reads: "every sample's duplicated genes", writes: "DUPLICATIONS/homologous_communities.tsv, .html", flags: "", doc: "how-it-works/duplications.html", functional: true
    },
    PANGENE: {
      name: "Pangenome", proc: "PANGENE", x: 930, y: 320, label: "below", lineOnly: "pan",
      text: "The reference proteome and the genes no reference gene matches, aligned with miniprot to every genome, its own included, then one pangene graph for the run and a presence/absence matrix with a column per sample.",
      reads: "the reference proteome, the corrected genes, and every genome: as PATHOTYPR renamed it, or with --skip_typing as the sample sheet gives it, under its sample ID", writes: "PANGENE/pangenome.gfa, .bubbles.txt, .Rtab, .stat", flags: "", doc: "how-it-works/workflow.html#9-pangenome", exit: "PANGENE/", exitSide: "right", functional: true
    },
    FINAL_REPORT: {
      name: "Report", proc: "FINAL_REPORT", x: 1050, y: 170, kind: "hub-end", label: "above",
      text: "One self-contained HTML file from everything above: quality, cases, missing genes, functions, duplications, pangenome.",
      reads: "every summary (REPORT_INPUTS gathers them)", writes: "final_report.html", flags: "", doc: "user-guide/report.html", exit: "final_report.html", exitSide: "below", functional: true
    },
    AMAP_REPORT: {
      name: "Interactive report", proc: "AMAP_REPORT", x: 1050, y: 100, kind: "hub-end", label: "above-end",
      text: "One self-contained HTML file organised by question: how each genome annotated, a map of the reference genes in every genome, one gene in every genome, what the correction changed, the pseudogenes by lineage. It needs only the correction, so it is written with or without eggNOG, whose steps add their own sections.",
      reads: "every genome's corrected genes, case tables and Bakta folder; the master tables, the pangenome and the communities when the run has them", writes: "amap_report.html", flags: "", doc: "user-guide/amap-report.html", exit: "amap_report.html", exitSide: "above"
    },
    REFERENCE: {
      name: "Reference", proc: "BUNDLED_REFERENCE", procCustom: "PREPARE_REFERENCE", x: 395, y: 60, kind: "end", label: "left",
      text: "The H37Rv proteome and its false-positive set, bundled in the container.",
      textCustom: "Your --ref_annotation and --ref_genome: the proteome extracted, and Bakta's false positives on your reference worked out by annotating it.",
      reads: "", writes: "nothing published", flags: "--ref_annotation, --ref_genome", doc: "user-guide/custom-reference.html"
    }
  };

  // The order main.nf's workflow calls them in, for check_site.py.
  var METRO_ORDER = [
    "VALIDATE_INPUTS", "EGGNOG_DB", "PATHOTYPR", "BAKTA_ANNOTATOR", "REFERENCE", "MINIPROT", "TYPE_CASES", "FIX_CDS",
    "ANNOT_EGGNOGMAP_BATCH", "GET_OUTPUTS", "GET_DUPLICATIONS", "PANGENE", "GET_COUNT", "GET_ALLDUP", "FINAL_REPORT",
    "AMAP_REPORT"
  ];

  // The sample's line runs on past Correction only to the steps that need
  // eggNOG, so it has a stretch of its own that a run without them leaves
  // out.
  var LINES = [
    { key: "main", d: "M70,170 H625" },
    { key: "fn", d: "M625,170 H1050" },
    { key: "ref", d: "M395,60 H400 Q410,60 410,70 V151 Q410,161 420,161 H650", endTick: [650, 161] },
    { key: "db", d: "M735,60 V160" },
    { key: "dup", d: "M735,179 H800 Q820,179 826,200 L836,230 Q842,250 862,250 H1000 Q1020,250 1028,230 L1040,195 Q1046,179 1050,179" },
    { key: "pan", d: "M735,188 H790 Q805,188 809,210 L815,295 Q818,320 840,320 H1010 Q1034,320 1040,295 L1048,195 Q1050,188 1050,181" },
    { key: "pan-feed", d: "M250,178 V330 Q250,350 270,350 H790 Q812,350 812,330", dUntyped: "M70,178 Q70,187 58,187 H24 Q14,187 14,197 V330 Q14,350 34,350 H790 Q812,350 812,330", feed: true },
    // The interactive report reads the correction's outputs directly, and the
    // functional steps' when they run (the dashed stretch).
    { key: "rep", d: "M625,152 V118 Q625,100 643,100 H1050" },
    { key: "rep-feed", d: "M1050,158 V114", dash: true }
  ];

  // The map is drawn once; the switches and the selection then only change
  // classes, so lines and stations fade rather than blink.
  function metroSvg(state) {
    var out = [];
    out.push('<svg class="metro" viewBox="0 0 1100 380" role="group" aria-label="The AMAP pipeline as a metro map: stations are the steps, lines are what travels between them.">');
    LINES.forEach(function (l) {
      if (l.feed) {
        // Every genome, for the pangenome: renamed by Lineage, or from the
        // sample sheet when typing is off (PANGENE_GENOME names them).
        out.push('<path class="ml ml-pan" data-line="' + l.key + '" d="' + l.d + '" style="stroke-width:3;stroke-dasharray:2 6"/>');
        out.push('<text class="exit" data-line="' + l.key + '" data-feed x="270" y="343">the renamed genomes, from Lineage</text>');
        return;
      }
      out.push('<path class="ml ml-' + l.key + '" data-line="' + l.key + '" d="' + l.d + '"' + (l.dash ? ' style="stroke-width:3;stroke-dasharray:2 6"' : "") + "/>");
      if (l.endTick) {
        out.push('<line class="ml ml-' + l.key + '" data-line="' + l.key + '" x1="' + l.endTick[0] + '" x2="' + l.endTick[0] + '" y1="' + (l.endTick[1] - 6) + '" y2="' + (l.endTick[1] + 6) + '" style="stroke-width:4"/>');
      }
    });
    Object.keys(STATIONS).forEach(function (key, n) {
      var s = STATIONS[key];
      out.push('<g class="stn" data-stn="' + key + '" tabindex="' + (n ? -1 : 0) + '" role="button" aria-pressed="false">');
      if (s.hub) {
        // A hub spans every track that runs through it.
        var top = s.y - (s.hubUp == null ? 18 : s.hubUp);
        var bottom = s.y + (s.hubDown == null ? 9 : s.hubDown);
        out.push('<rect class="hub" x="' + (s.x - 8) + '" y="' + top + '" width="16" height="' + (bottom - top) + '" rx="8"/>');
      } else if (s.kind === "end") {
        out.push('<circle class="end" cx="' + s.x + '" cy="' + s.y + '" r="8"/>');
        out.push('<circle class="dot" cx="' + s.x + '" cy="' + s.y + '" r="3.5" style="stroke-width:0;fill:var(--bg-elev)"/>');
      } else if (s.kind === "hub-end") {
        out.push('<circle class="hub" cx="' + s.x + '" cy="' + s.y + '" r="12"/>');
        out.push('<circle class="end" cx="' + s.x + '" cy="' + s.y + '" r="5"/>');
      } else {
        out.push('<circle class="dot" cx="' + s.x + '" cy="' + s.y + '" r="8"/>');
      }
      // labelDx nudges a name clear of its neighbour's.
      var lx = s.x + (s.labelDx || 0);
      var anchor = "middle";
      var ny;
      var py;
      if (s.label === "above") {
        ny = s.y - 34;
        py = s.y - 21;
      } else if (s.label === "above-end") {
        // Above, ending at the station: at the map's right edge
        lx = s.x + 12;
        anchor = "end";
        ny = s.y - 34;
        py = s.y - 21;
      } else if (s.label === "below") {
        ny = s.y + 32;
        py = s.y + 45;
      } else if (s.label === "left") {
        lx = s.x - 16;
        anchor = "end";
        ny = s.y - 2;
        py = s.y + 11;
      } else {
        lx = s.x + 16;
        anchor = "start";
        ny = s.y - 2;
        py = s.y + 11;
      }
      if (s.hub && s.label === "above") {
        ny -= 6;
        py -= 6;
      }
      if (s.hub && s.label === "below" && s.hubDown) {
        ny += s.hubDown - 9;
        py += s.hubDown - 9;
      }
      out.push('<text class="name" x="' + lx + '" y="' + ny + '" text-anchor="' + anchor + '">' + esc(s.name) + "</text>");
      out.push('<text class="proc" x="' + lx + '" y="' + py + '" text-anchor="' + anchor + '">' + esc(s.proc) + "</text>");
      out.push("</g>");
    });
    out.push("</svg>");
    return out.join("");
  }

  function paintMetro(svg, state) {
    svg.querySelectorAll("[data-line]").forEach(function (el) {
      el.classList.toggle("off", lineOff(el.getAttribute("data-line"), state));
    });
    var feed = LINES.filter(function (l) {
      return l.feed;
    })[0];
    var untyped = state.skipTyping;
    svg.querySelector('path[data-line="pan-feed"]').setAttribute("d", untyped ? feed.dUntyped : feed.d);
    var feedText = svg.querySelector("[data-feed]");
    feedText.setAttribute("x", untyped ? "34" : "270");
    feedText.textContent = untyped ? "the genomes of the sample sheet, under their sample IDs" : "the renamed genomes, from Lineage";
    svg.querySelectorAll(".stn").forEach(function (g) {
      var key = g.getAttribute("data-stn");
      var s = STATIONS[key];
      var proc = state.custom && s.procCustom ? s.procCustom : s.proc;
      g.classList.toggle("off", stationOff(key, state));
      g.classList.toggle("skipped", key === "PATHOTYPR" && state.skipTyping);
      g.classList.toggle("on", state.selected === key);
      g.setAttribute("aria-pressed", state.selected === key ? "true" : "false");
      g.setAttribute("aria-label", s.name + ", " + proc + (stationOff(key, state) ? ", not run" : ""));
      g.querySelector("text.proc").textContent = proc;
    });
  }

  function lineOff(key, state) {
    if (key === "main") return false;
    if (key === "ref") return false;
    if (state.skipFunctional && (key === "fn" || key === "db" || key === "dup" || key === "pan" || key === "pan-feed" || key === "rep-feed")) return true;
    return false;
  }

  function stationOff(key, state) {
    var s = STATIONS[key];
    if (state.skipFunctional && s.functional) return true;
    if (state.skipTyping && key === "PATHOTYPR") return true;
    return false;
  }

  function writes(state) {
    var dirs = [];
    if (!state.skipTyping) dirs.push("PATHOTYPR/");
    dirs.push("BAKTA/", "CDS_CASES/", "AMAP/");
    if (!state.skipFunctional) dirs.push("DUPLICATIONS/", "PANGENE/", "final_report.html");
    dirs.push("amap_report.html", "pipeline_info/");
    return dirs;
  }

  function infoHtml(key, state) {
    if (!key) {
      return '<p class="hint">Pick a station to see what it does, what it reads and where its files go. The switches above show what a run leaves out.</p>';
    }
    var s = STATIONS[key];
    var proc = state.custom && s.procCustom ? s.procCustom : s.proc;
    var text = state.custom && s.textCustom ? s.textCustom : s.text;
    var status = "";
    if (stationOff(key, state)) status = '<p class="note">Not run with the switches above.</p>';
    if (key === "PANGENE" && state.skipTyping && !state.skipFunctional) status = '<p class="note">With <code>--skip_typing</code> it searches the genomes of the sample sheet, each under its sample ID.</p>';
    return (
      "<div><h3>" + esc(s.name) + "<code>" + esc(proc) + "</code></h3><p>" + esc(text) + "</p>" + status +
      (s.doc ? '<p class="note"><a href="' + s.doc + '">More on this step</a></p>' : "") + "</div>" +
      "<dl>" +
      (s.reads ? "<dt>Reads</dt><dd>" + esc(s.reads) + "</dd>" : "") +
      (s.writes ? "<dt>Writes</dt><dd>" + esc(s.writes) + "</dd>" : "") +
      (s.flags ? "<dt>Options</dt><dd>" + esc(s.flags) + "</dd>" : "") +
      "</dl>"
    );
  }

  function mountMetro(el) {
    var stage = el.querySelector("[data-metro-stage]");
    var info = el.querySelector("[data-metro-info]");
    var writesEl = el.querySelector("[data-metro-writes]");
    // The line chips of the legend, where the page has them, by their line.
    var chips = Array.prototype.map.call(el.querySelectorAll(".metro-lines span"), function (span) {
      var i = span.querySelector("i");
      var m = /--m-(\w+)/.exec((i && i.getAttribute("style")) || "");
      return { span: span, line: m ? m[1] : "" };
    });
    var state = { skipTyping: false, skipFunctional: false, custom: false, selected: null };
    stage.innerHTML = metroSvg(state);
    var svg = stage.querySelector("svg");
    roving(stage);

    // Once a station is picked, the panel under the map keeps the height of
    // the tallest station's text, so picking another never moves what
    // follows. Before that it holds only its short hint.
    function fit() {
      if (!info || !info.offsetWidth) return;
      info.style.minHeight = "";
      if (!state.selected) return;
      var probe = info.cloneNode(false);
      probe.removeAttribute("data-metro-info");
      probe.removeAttribute("aria-live");
      probe.setAttribute("aria-hidden", "true");
      probe.style.cssText = "position:absolute;visibility:hidden;pointer-events:none;left:0;top:0;min-height:0;width:" + info.offsetWidth + "px";
      info.parentNode.appendChild(probe);
      var tallest = 0;
      Object.keys(STATIONS).forEach(function (key) {
        probe.innerHTML = infoHtml(key, state);
        tallest = Math.max(tallest, probe.offsetHeight);
      });
      probe.remove();
      var kept = parseFloat(getComputedStyle(info).minHeight) || 0;
      info.style.minHeight = Math.max(kept, tallest) + "px";
    }
    var fitWidth = 0;
    if (info && "ResizeObserver" in window) {
      new ResizeObserver(function () {
        if (Math.abs(info.offsetWidth - fitWidth) < 2) return;
        fitWidth = info.offsetWidth;
        fit();
      }).observe(info);
    }

    // On a narrow screen the map scrolls sideways; a fade at the edge says
    // there is more of it.
    function edges() {
      var more = stage.scrollWidth - stage.clientWidth;
      stage.classList.toggle("more-right", more > 2 && stage.scrollLeft < more - 2);
      stage.classList.toggle("more-left", more > 2 && stage.scrollLeft > 2);
    }
    stage.addEventListener("scroll", edges, { passive: true });
    if ("ResizeObserver" in window) new ResizeObserver(edges).observe(stage);

    function paint() {
      paintMetro(svg, state);
      chips.forEach(function (c) {
        c.span.classList.toggle("off", state.skipFunctional && (c.line === "db" || c.line === "dup" || c.line === "pan"));
      });
      if (writesEl) {
        writesEl.innerHTML = "This run writes " + writes(state).map(function (d) {
          return "<code>" + d + "</code>";
        }).join(" ") + ".";
      }
      describe();
      edges();
    }
    // The panel is a live region: it is only rewritten when what it says
    // changes, so a switch that changes nothing there is not read out again.
    function describe() {
      if (!info) return;
      var html = infoHtml(state.selected, state);
      if (html !== info._said) {
        info._said = html;
        info.innerHTML = html;
      }
    }
    function pick(g) {
      var key = g.getAttribute("data-stn");
      if (key === state.selected) return;
      var first = !state.selected;
      state.selected = key;
      if (first) fit();
      paint();
      // The station picked comes to the middle of a map that scrolls.
      if (stage.scrollWidth > stage.clientWidth + 2) {
        var r = g.getBoundingClientRect();
        var sr = stage.getBoundingClientRect();
        var by = r.left + r.width / 2 - (sr.left + sr.width / 2);
        if (Math.abs(by) > 8) stage.scrollBy({ left: by, behavior: reduceMotion ? "auto" : "smooth" });
      }
    }
    el.addEventListener("change", function (event) {
      var t = event.target;
      if (t.matches("[data-metro-opt]")) {
        state[t.getAttribute("data-metro-opt")] = t.checked;
        paint();
        fit();
      }
    });
    stage.addEventListener("click", function (event) {
      var g = event.target.closest(".stn");
      if (g) pick(g);
    });
    stage.addEventListener("keydown", function (event) {
      if (event.key !== "Enter" && event.key !== " ") return;
      var g = event.target.closest(".stn");
      if (!g) return;
      event.preventDefault();
      pick(g);
    });
    paint();
    fit();
  }

  /* ── The overview: what goes into a decision, and what comes out ── */

  // The hero's 12.4 kb of CDC1551 again, with the three things each decision
  // reads: Bakta's calls, the reference proteins aligned, and the calls Bakta
  // makes wrongly on the reference, aligned here.
  function inputsSpec() {
    var parts = rowsFor("hero");
    return {
      from: parts.win.from, to: parts.win.to, labelW: 84, geneH: 18, rowGap: 22, top: 18,
      alt: "The first 12.4 kb of M. tuberculosis CDC1551: Bakta's twelve calls, the eight H37Rv proteins aligned, and two of Bakta's known false positives.",
      refName: parts.reference,
      rows: [
        { label: "Bakta's\ncalls", cls: "bakta", genes: parts.bakta },
        { label: parts.reference + "\nproteins", cls: "ref", genes: parts.ref },
        { label: "Known false\npositives", cls: "drop", genes: parts.bad }
      ],
      // Where the false positives lie, and the call Bakta makes there.
      marks: parts.bad.length ? [{
        s: Math.min.apply(null, parts.bad.map(function (f) {
          return f.s;
        })),
        e: Math.max.apply(null, parts.bad.map(function (f) {
          return f.e;
        })),
        rows: [0, 2], label: "a call Bakta is known to make wrongly", labelAt: "left"
      }] : [],
      axis: true, tick: 2000,
      tickLabel: function (t) {
        return t / 1000 + " kb";
      }
    };
  }

  // A small drawing of one window: its rows without their names, no axis.
  // It sits inside a link, so it is a picture: no tab stops, no tips.
  function miniSpec(key) {
    var spec = windowSpec(key);
    spec.labelW = 0;
    spec.axis = false;
    spec.geneH = 14;
    spec.rowGap = 12;
    spec.top = spec.marks.length ? 16 : 4;
    spec.padRight = 2;
    spec.inert = true;
    spec.marks.forEach(function (m) {
      m.label = "";
    });
    spec.rows.forEach(function (row) {
      row.genes.forEach(function (g) {
        g.gapLabel = null;
        g.nameAbove = false;
        g.moved = [];
      });
    });
    spec.edges = [];
    return spec;
  }

  function mountInputs(el) {
    var stage = el.querySelector("[data-fig-stage]");
    if (!DATA.windows.hero) return;
    var key = legendAfter(stage);
    render(stage, inputsSpec, function (svg, spec) {
      key.innerHTML = legend(spec._used, spec.refName);
    });
    attachTips(stage);
  }

  function mountOutcomes(el) {
    el.querySelectorAll("[data-mini]").forEach(function (card) {
      var key = card.getAttribute("data-mini");
      var stage = card.querySelector("[data-fig-stage]");
      if (!DATA.windows[key]) return;
      render(stage, function () {
        return miniSpec(key);
      }, function (svg) {
        svg.setAttribute("data-phase", "4");
      });
    });
  }

  /* ── Duplications ──

     A family of proteins and the pair of homologous communities that holds
     it, from a real run of GET_DUPLICATIONS and GET_ALLDUP (dup-data.js);
     check_site.py reruns both scripts on the hits and fails if the links,
     the family, the communities or their numbers come out different. */

  var DUP = window.AMAP_DUP_DATA || null;

  function pct(x) {
    return (Math.round(x * 10) / 10).toFixed(1) + "%";
  }

  // On a line, the share rounded down, so a link never reads as a miss.
  function floorPct(x) {
    return Math.floor(x + 1e-9) + "%";
  }

  // A network: nodes at x, y between 0 and 1, links between them, and a soft
  // outline around the nodes of a family. One tab stop, like a gene track:
  // the arrows move through the proteins, then the lines.
  function netSvg(o, width) {
    var W = Math.max(260, Math.round(width));
    var H = Math.round(W < 520 ? W * (o.tall || 0.92) : Math.min(o.maxH || 380, W * (o.ratio || 0.56)));
    var pad = 38;
    var X = function (x) {
      return pad + x * (W - 2 * pad);
    };
    var Y = function (y) {
      return pad - 6 + y * (H - 2 * pad);
    };
    var out = ['<svg class="net" viewBox="0 0 ' + W + " " + H + '" role="group" aria-roledescription="network" aria-label="' + esc(o.alt) + '">'];
    if (o.hull && o.hull.length) {
      var pts = o.hull.map(function (i) {
        return X(o.nodes[i].x).toFixed(1) + "," + Y(o.nodes[i].y).toFixed(1);
      });
      out.push('<polygon class="net-hull" points="' + pts.join(" ") + '"/>');
    }
    // Every node, its name and every label so far, so that the next label
    // finds a stretch of its line where it covers none of them.
    var boxes = o.nodes.map(function (n) {
      return { x: X(n.x) - n.r - 2, y: Y(n.y) - n.r - 2, w: 2 * n.r + 4, h: 2 * n.r + 20 };
    });
    var clear = function (b) {
      return boxes.every(function (c) {
        return b.x + b.w < c.x || c.x + c.w < b.x || b.y + b.h < c.y || c.y + c.h < b.y;
      });
    };
    var item = function (tip) {
      return ' role="img" tabindex="-1" aria-label="' + esc(plain(tip)) + '" data-tip="' + esc(tip) + '"';
    };
    o.links.forEach(function (l) {
      var a = o.nodes[l.a];
      var b = o.nodes[l.b];
      out.push(
        '<g class="hot net-link ' + l.cls + '"' + item(l.tip) + ">" +
          '<line class="net-hit" x1="' + X(a.x).toFixed(1) + '" y1="' + Y(a.y).toFixed(1) + '" x2="' + X(b.x).toFixed(1) + '" y2="' + Y(b.y).toFixed(1) + '"/>' +
          '<line class="net-line" x1="' + X(a.x).toFixed(1) + '" y1="' + Y(a.y).toFixed(1) + '" x2="' + X(b.x).toFixed(1) + '" y2="' + Y(b.y).toFixed(1) + '"/></g>'
      );
    });
    // The labels over every line, the links' first, each where it fits.
    var rank = { kept: 0, near: 1, far: 2, ghost: 3 };
    o.links.slice().sort(function (x, y) {
      return rank[x.cls] - rank[y.cls];
    }).forEach(function (l) {
      if (!l.label) return;
      var a = o.nodes[l.a];
      var b = o.nodes[l.b];
      var lw = textWidth(l.label, 10.5) + 10;
      var spot = null;
      [0.5, 0.36, 0.64, 0.26, 0.74, 0.44, 0.56].some(function (f) {
        var mx = X(a.x) + f * (X(b.x) - X(a.x));
        var my = Y(a.y) + f * (Y(b.y) - Y(a.y));
        var box = { x: mx - lw / 2, y: my - 8, w: lw, h: 16 };
        if (clear(box)) spot = box;
        return spot;
      });
      if (!spot) return;
      boxes.push(spot);
      // The label repeats its line, which says it all to a screen reader.
      out.push(
        '<g class="hot net-link ' + l.cls + '" aria-hidden="true" data-tip="' + esc(l.tip) + '">' +
          // A finger's target: 24 px tall, though the pill is drawn 16.
          '<rect class="net-pill-hit" x="' + spot.x.toFixed(1) + '" y="' + (spot.y - 4).toFixed(1) + '" width="' + Math.max(lw, 24).toFixed(1) + '" height="24"/>' +
          '<rect class="net-pill" x="' + spot.x.toFixed(1) + '" y="' + spot.y.toFixed(1) + '" width="' + lw.toFixed(1) + '" height="16" rx="8"/>' +
          '<text class="net-pct" x="' + (spot.x + lw / 2).toFixed(1) + '" y="' + (spot.y + 11.5).toFixed(1) + '" text-anchor="middle">' + esc(l.label) + "</text></g>"
      );
    });
    o.nodes.forEach(function (n, i) {
      var cx = X(n.x);
      var cy = Y(n.y);
      // The name goes under the node, or over it at the bottom of the frame.
      var below = cy + n.r + 14 < H - 2;
      out.push('<g class="hot net-node ' + n.cls + '"' + item(n.tip).replace('tabindex="-1"', 'tabindex="' + (i ? -1 : 0) + '"') + ">");
      out.push('<circle class="net-ring" cx="' + cx.toFixed(1) + '" cy="' + cy.toFixed(1) + '" r="' + (n.r + 4).toFixed(1) + '"/>');
      out.push('<circle cx="' + cx.toFixed(1) + '" cy="' + cy.toFixed(1) + '" r="' + n.r.toFixed(1) + '"/>');
      out.push('<text class="net-name" x="' + cx.toFixed(1) + '" y="' + (below ? cy + n.r + 13 : cy - n.r - 6).toFixed(1) + '" text-anchor="middle">' + esc(n.label) + "</text>");
      out.push("</g>");
    });
    out.push("</svg>");
    return out.join("");
  }

  // Draws into stage at its width, and again when that changes.
  function onWidth(stage, paint) {
    var holder = document.createElement("div");
    holder.className = "fig-svg";
    stage.insertBefore(holder, stage.firstChild);
    var last = 0;
    function go() {
      var w = stage.clientWidth;
      if (!w || Math.abs(w - last) < 8) return;
      last = w;
      holder.innerHTML = paint(w);
    }
    go();
    if ("ResizeObserver" in window) new ResizeObserver(go).observe(stage);
    else window.addEventListener("resize", go);
    attachTips(stage, { what: "a protein or a line for its numbers." });
    roving(stage);
  }

  function radius(aa, longest) {
    return 6 + 10 * Math.sqrt(aa / longest);
  }

  function hullOf(points) {
    // Andrew's monotone chain, over indices into points.
    var idx = points.map(function (p, i) {
      return i;
    }).sort(function (a, b) {
      return points[a].x - points[b].x || points[a].y - points[b].y;
    });
    var cross = function (o, a, b) {
      return (points[a].x - points[o].x) * (points[b].y - points[o].y) - (points[a].y - points[o].y) * (points[b].x - points[o].x);
    };
    var lower = [];
    var upper = [];
    idx.forEach(function (i) {
      while (lower.length >= 2 && cross(lower[lower.length - 2], lower[lower.length - 1], i) <= 0) lower.pop();
      lower.push(i);
    });
    idx.slice().reverse().forEach(function (i) {
      while (upper.length >= 2 && cross(upper[upper.length - 2], upper[upper.length - 1], i) <= 0) upper.pop();
      upper.push(i);
    });
    return lower.slice(0, -1).concat(upper.slice(0, -1));
  }

  function proteinTip(n, genome) {
    return "<b>" + esc(n.id) + "</b>" + (n.gene ? " (" + esc(n.gene) + ")" : "") + (n.product ? ", " + esc(n.product.charAt(0).toLowerCase() + n.product.slice(1)) : "") +
      "<br>" + fmt(n.aa) + " residues" + (genome ? ", " + esc(genome) : "");
  }

  function mountDupFamily(el) {
    var fam = DUP && DUP.family;
    if (!fam) return;
    var stage = el.querySelector("[data-fig-stage]");
    var longest = Math.max.apply(null, fam.nodes.map(function (n) {
      return n.aa;
    }));
    var members = [];
    var nodes = fam.nodes.map(function (n, i) {
      if (n.member) members.push(i);
      return { x: n.x, y: n.y, r: radius(n.aa, longest), label: n.id, cls: n.member ? "in" : "out", tip: proteinTip(n) + (n.member ? "" : "<br>Not in the family: none of its hits is a link") };
    });
    // Each protein outside the family is drawn with its best hit only.
    var best = {};
    fam.pairs.forEach(function (p) {
      [p.a, p.b].forEach(function (i) {
        if (!fam.nodes[i].member && (!best[i] || p.cover > best[i].cover)) best[i] = p;
      });
    });
    var links = fam.pairs.filter(function (p) {
      var a = fam.nodes[p.a];
      var b = fam.nodes[p.b];
      return (a.member && b.member) || best[a.member ? p.b : p.a] === p;
    }).map(function (p) {
      var a = fam.nodes[p.a];
      var b = fam.nodes[p.b];
      var inside = a.member && b.member;
      var cls = p.kept ? "kept" : inside ? "near" : "far";
      return {
        a: p.a, b: p.b, cls: cls, label: floorPct(p.cover),
        tip: "<b>" + esc(a.id) + " and " + esc(b.id) + "</b><br>" + fmt(p.nident) + " identical residues: " + pct((100 * p.nident) / p.len[0]) + " of " + esc(a.id) +
          ", " + pct((100 * p.nident) / p.len[1]) + " of " + esc(b.id) + "<br>" + (p.kept ? "A link: the family holds both." : inside ? "No link, though both are in the family through others." : "No link: the best hit " + esc(a.member ? b.id : a.id) + " has.")
      };
    });
    // Links that do not count under the ones that do.
    var order = { far: 0, near: 1, kept: 2 };
    links.sort(function (x, y) {
      return order[x.cls] - order[y.cls];
    });
    var hull = hullOf(members.map(function (i) {
      return nodes[i];
    })).map(function (k) {
      return members[k];
    });
    onWidth(stage, function (w) {
      return netSvg({
        nodes: nodes, links: links, hull: hull, ratio: 0.62, maxH: 420, tall: 1.05,
        alt: fam.label + ": " + members.length + " proteins make one duplication family; the lines are the pairs Diamond matched, with the share of the shorter protein their identical residues cover."
      }, w);
    });
  }

  function mountDupCommunities(el) {
    var com = DUP && DUP.communities;
    if (!com) return;
    var panels = el.querySelectorAll("[data-net]");
    var sides = com.sides;
    var longest = 0;
    sides.forEach(function (s) {
      s.nodes.forEach(function (n) {
        longest = Math.max(longest, n.aa);
      });
    });
    // Each side's links, and the other side's as faint ghosts: a link one
    // genome has and the other does not.
    var keyOf = function (side, e) {
      var a = side.nodes[e[0]].twin;
      var b = side.nodes[e[1]].twin;
      return a < b ? a + "|" + b : b + "|" + a;
    };
    var has = sides.map(function (s) {
      var set = {};
      s.edges.forEach(function (e) {
        set[keyOf(s, e)] = true;
      });
      return set;
    });
    panels.forEach(function (panel, k) {
      var side = sides[k];
      var other = sides[1 - k];
      var stage = panel.querySelector("[data-fig-stage]");
      var head = panel.querySelector("[data-net-head]");
      if (head) head.innerHTML = "<b>" + esc(side.label) + "</b><span>" + esc(side.community) + ", " + side.nodes.length + " proteins, " + side.edges.length + " links</span>";
      var byTwin = {};
      side.nodes.forEach(function (n, i) {
        byTwin[n.twin] = i;
      });
      var nodes = side.nodes.map(function (n) {
        return { x: n.x, y: n.y, r: radius(n.aa, longest), label: n.id, cls: "in", tip: proteinTip(n, side.label) };
      });
      var links = side.edges.map(function (e) {
        return { a: e[0], b: e[1], cls: "kept", label: "", tip: "<b>" + esc(side.nodes[e[0]].id) + " and " + esc(side.nodes[e[1]].id) + "</b><br>A link in " + esc(side.label) + "." };
      });
      other.edges.forEach(function (e) {
        var key = keyOf(other, e);
        if (has[k][key]) return;
        var a = byTwin[other.nodes[e[0]].twin];
        var b = byTwin[other.nodes[e[1]].twin];
        if (a == null || b == null) return;
        links.unshift({ a: a, b: b, cls: "ghost", label: "", tip: "<b>" + esc(side.nodes[a].id) + " and " + esc(side.nodes[b].id) + "</b><br>A link in " + esc(other.label) + ", not in " + esc(side.label) + "." });
      });
      var hull = hullOf(nodes);
      onWidth(stage, function (w) {
        return netSvg({
          nodes: nodes, links: links, hull: hull,
          alt: side.label + ": " + side.community + ", " + side.nodes.length + " proteins and " + side.edges.length + " links."
        }, w);
      });
    });
    var bridge = el.querySelector("[data-net-bridge]");
    if (bridge) bridge.textContent = com.cross + " hits between the two: a homologous pair";
    var table = el.querySelector("[data-net-metrics]");
    if (table) {
      var rows = [
        ["Size", "size", 0], ["AvgDegree", "avg_degree", 2], ["Density", "density", 2],
        ["AvgClustering", "avg_clustering", 2], ["AvgPathLength", "avg_path_length", 2]
      ];
      table.innerHTML =
        "<thead><tr><th>Column</th>" + sides.map(function (s, k) {
          return '<th class="num">' + esc(s.genome) + '<span class="th-sub">columns ending in ' + (k + 1) + "</span></th>";
        }).join("") + "</tr></thead><tbody>" +
        rows.map(function (r) {
          return "<tr><td><code>" + r[0] + "</code></td>" + sides.map(function (s) {
            return '<td class="num">' + s.metrics[r[1]].toFixed(r[2]) + "</td>";
          }).join("") + "</tr>";
        }).join("") + "</tbody>";
    }
  }

  /* ── Mount ── */

  document.querySelectorAll("[data-hero]").forEach(mountHero);
  document.querySelectorAll("[data-fig]").forEach(mountCase);
  document.querySelectorAll("[data-metro]").forEach(mountMetro);
  document.querySelectorAll("[data-inputs]").forEach(mountInputs);
  document.querySelectorAll("[data-outcomes]").forEach(mountOutcomes);
  document.querySelectorAll('[data-dup="family"]').forEach(mountDupFamily);
  document.querySelectorAll('[data-dup="communities"]').forEach(mountDupCommunities);

  window.AMAPFigures = {
    draw: draw, render: render, tickLabel: tickLabel, attachTips: attachTips, rowsFor: rowsFor, windowSpec: windowSpec, fmt: fmt, esc: esc,
    legend: legend, refocus: refocus, isShown: isShown, nameSubjects: nameSubjects, overlap: overlap,
    stations: STATIONS, metroOrder: METRO_ORDER
  };
})();

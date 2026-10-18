/* AMAP website: the drawings. Every gene drawn here is a real one, from AMAP
   run on Mycobacterium tuberculosis CDC1551 (NC_002755.2) against the H37Rv
   reference it ships with, unless a figure says otherwise; the coordinates are
   that run's. tools/check_site.py puts each scenario through the pipeline's
   own code and fails if a drawing keeps something the code would not.
   Plain JavaScript, no build step. */
(function () {
  "use strict";

  var reduceMotion = window.matchMedia && window.matchMedia("(prefers-reduced-motion: reduce)").matches;

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

  /* ── The track renderer ──

     A drawing is rows of genes on one coordinate axis:
       { from, to, w, rows: [{ label, cls, phase, genes: [gene] }],
         ribbons: [{ a: [row, gene], b: [row, gene], kind, phase }],
         marks: [{ row, s, e, label }], edges: [{ row, at, side, label }] }
     A gene is { s, e, st, name, cls, tip, open, gap, dim, id }, 1-based and
     inclusive like everything AMAP writes. */

  var figureCount = 0;

  function arrow(x1, x2, y, h, strand, open) {
    var tip = Math.min(8, Math.max(0, (x2 - x1) * 0.45));
    var mid = y + h / 2;
    // A partial gene has no arrowhead at its cut end, and a ragged edge there.
    var zig = function (x, dir) {
      var d = 3 * dir;
      return "L" + (x + d) + "," + (y + h * 0.25) + "L" + x + "," + mid + "L" + (x + d) + "," + (y + h * 0.75);
    };
    if (strand === "-") {
      var rightOpen = open && open.indexOf("5'") !== -1; // 5' end is on the right
      var leftOpen = open && open.indexOf("3'") !== -1;
      var p = "M" + x2 + "," + y;
      p += leftOpen ? "H" + x1 + zig(x1, 1) + "L" + x1 + "," + (y + h) : "H" + (x1 + tip) + "L" + x1 + "," + mid + "L" + (x1 + tip) + "," + (y + h);
      p += "H" + x2;
      p += rightOpen ? zig(x2, -1) + "Z" : "Z";
      return p;
    }
    var endOpen = open && open.indexOf("3'") !== -1;
    var startOpen = open && open.indexOf("5'") !== -1;
    var q = "M" + x1 + "," + y;
    q += endOpen ? "H" + x2 + zig(x2, -1) + "L" + x2 + "," + (y + h) : "H" + (x2 - tip) + "L" + x2 + "," + mid + "L" + (x2 - tip) + "," + (y + h);
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

  function draw(spec, width) {
    var id = "f" + ++figureCount;
    // Drawn at the width it is shown at, so its text is always the page's size.
    var W = Math.max(280, Math.round(width || spec.w || 640));
    var labelW = spec.labelW == null ? 104 : spec.labelW;
    var x0 = labelW;
    var x1 = W - (spec.padRight == null ? 8 : spec.padRight);
    var geneH = spec.geneH || 16;
    var rowGap = spec.rowGap || 30;
    var top = spec.top == null ? 10 : spec.top;
    var span = spec.to - spec.from;
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
    var out = [];
    out.push(
      '<svg class="trk" viewBox="0 0 ' + W + " " + Math.round(height) + '" role="img" aria-label="' + esc(spec.alt || "") + '">'
    );
    out.push(
      '<defs><pattern id="hatch-' + id + '" width="6" height="6" patternUnits="userSpaceOnUse" patternTransform="rotate(45)">' +
        '<line class="hatch-line" x1="0" y1="0" x2="0" y2="6"/></pattern></defs>'
    );

    // Ribbons first, under the genes: straight where the two agree, slanted
    // where their coordinates differ.
    (spec.ribbons || []).forEach(function (r) {
      var ga = spec.rows[r.a[0]].genes[r.a[1]];
      var gb = spec.rows[r.b[0]].genes[r.b[1]];
      var ya = gy0(r.a[0], r.a[1]) + geneH + 1;
      var yb = gy0(r.b[0], r.b[1]) - 1;
      var clip = function (p) {
        return Math.max(spec.from, Math.min(spec.to, p));
      };
      var pts = [X(clip(ga.s) - 0.5), ya, X(clip(ga.e) + 0.5), ya, X(clip(gb.e) + 0.5), yb, X(clip(gb.s) - 0.5), yb];
      out.push(
        '<polygon class="rib rib-' + r.kind + (r.phase ? " ph-" + r.phase : "") + '"' + steps(r) + ' points="' + pts.map(function (v) {
          return Math.round(v * 10) / 10;
        }).join(",") + '"/>'
      );
    });

    (spec.marks || []).forEach(function (m) {
      var mx1 = X(m.s - 0.5);
      var mx2 = X(m.e + 0.5);
      var my1 = ys[m.rows[0]] - 4;
      var my2 = ys[m.rows[1]] + geneH + (lanes[m.rows[1]].n - 1) * LANE + 4;
      out.push('<g class="' + (m.phase ? "ph-" + m.phase : "") + '"' + steps(m) + ">");
      out.push('<rect class="diff' + (m.kind ? " diff-" + m.kind : "") + '" x="' + mx1.toFixed(1) + '" y="' + my1 + '" width="' + Math.max(2, mx2 - mx1).toFixed(1) + '" height="' + (my2 - my1) + '" rx="2"/>');
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

    spec.rows.forEach(function (row, i) {
      var gy = ys[i];
      // The label and the line stay while a row is swapped for a later
      // version of itself (an overlay), so they only take its first step.
      out.push('<g class="row-frame' + (row.phase ? " ph-" + row.phase : "") + '"' + steps({ from: row.from }) + ">");
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
      row.genes.forEach(function (g, k) {
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
        var group = '<g class="hot' + (g.phaseCls ? " " + g.phaseCls : "") + '" tabindex="0"' + (g.tip ? ' data-tip="' + esc(g.tip) + '"' : "") + ">";
        out.push(group);
        out.push('<path class="' + klass + (cutL || cutR ? " g-cut" : "") + '"' + fill + ' d="' + arrow(gx1, gx2, gy, geneH, g.st, open) + '"/>');
        if (g.gap) {
          var gpx = X(g.gap);
          out.push('<line class="gap-mark" x1="' + gpx.toFixed(1) + '" x2="' + gpx.toFixed(1) + '" y1="' + (gy - 5) + '" y2="' + (gy + geneH + 5) + '"/>');
          if (g.gapLabel) {
            out.push('<text class="diff-label diff-label-ref" x="' + gpx.toFixed(1) + '" y="' + (gy + geneH + 16) + '" text-anchor="middle">' + esc(g.gapLabel) + "</text>");
          }
        }
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
            out.push('<text class="gname" x="' + ((gx1 + gx2) / 2).toFixed(1) + '" y="' + (gy - 5) + '" text-anchor="middle">' + esc(g.name) + "</text>");
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
      out.push('<g class="axis"><line x1="' + x0 + '" x2="' + x1 + '" y1="' + ay + '" y2="' + ay + '"/>');
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
      var first = Math.ceil(spec.from / step) * step;
      for (var t = first; t <= spec.to; t += step) {
        var tx = X(t);
        out.push('<line x1="' + tx.toFixed(1) + '" x2="' + tx.toFixed(1) + '" y1="' + ay + '" y2="' + (ay + 4) + '"/>');
        out.push('<text x="' + tx.toFixed(1) + '" y="' + (ay + 15) + '" text-anchor="middle">' + esc(label(t)) + "</text>");
      }
      out.push("</g>");
    }
    out.push("</svg>");
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

  // A tooltip for whatever in a drawing carries data-tip.
  function attachTips(container) {
    var tip = container.querySelector(".demo-tip");
    if (!tip) {
      tip = document.createElement("div");
      tip.className = "demo-tip";
      tip.hidden = true;
      container.appendChild(tip);
    }
    function show(target) {
      var text = target.getAttribute("data-tip");
      if (!text) return;
      container.querySelectorAll(".hot.on").forEach(function (el) {
        el.classList.remove("on");
      });
      target.classList.add("on");
      tip.innerHTML = text;
      tip.hidden = false;
      var box = container.getBoundingClientRect();
      var r = target.getBoundingClientRect();
      var left = r.left - box.left + r.width / 2 - tip.offsetWidth / 2;
      left = Math.max(0, Math.min(left, box.width - tip.offsetWidth));
      var topPos = r.top - box.top - tip.offsetHeight - 8;
      if (topPos < 0) topPos = r.bottom - box.top + 8;
      tip.style.left = left + "px";
      tip.style.top = topPos + "px";
    }
    function hide() {
      tip.hidden = true;
      container.querySelectorAll(".hot.on").forEach(function (el) {
        el.classList.remove("on");
      });
    }
    container.addEventListener("mouseover", function (event) {
      var hot = event.target.closest(".hot");
      if (hot && container.contains(hot)) show(hot);
    });
    container.addEventListener("mouseleave", hide);
    container.addEventListener("focusin", function (event) {
      var hot = event.target.closest(".hot");
      if (hot) show(hot);
    });
    container.addEventListener("focusout", hide);
    container.addEventListener("click", function (event) {
      var hot = event.target.closest(".hot");
      if (hot) show(hot);
      else hide();
    });
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
      holder.innerHTML = draw(specFn(), width);
      if (after) after(holder.querySelector("svg"));
    }
    go();
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
    hero: { tick: 2000, kbTicks: true, alt: "The first 12.4 kb of M. tuberculosis CDC1551: twelve Bakta calls, eight H37Rv proteins aligned, and the eight genes AMAP keeps." },
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
    var bakta = win.bakta.map(function (c) {
      var gone = c.case === 1 && !Object.keys(kept).some(function (id) {
        return id === c.tag;
      });
      var tip = "<b>Bakta call " + short(c.tag) + (c.name ? " (" + c.name + ")" : "") + "</b>" + (c.product ? ", " + esc(c.product) : "") +
        "<br>" + span(c) + ", " + fmt(c.e - c.s + 1) + " bp<br>" + CASE_SAYS[c.case];
      if (c.case === 1) tip += gone ? " Dropped: " + dropReason(c, win) + "." : " Kept: nothing rules it out.";
      return { s: c.s, e: c.e, st: c.st, name: c.name, tip: tip, phaseCls: c.case === 1 ? (gone ? "ph-dim c1-call" : "c1-call") : "" };
    });
    var ref = win.ref.map(function (r) {
      var whole = r.complete ? "the whole protein, to its stop" : "not the whole protein";
      var gap = look.gaps && look.gaps[r.prot];
      return {
        s: r.s, e: r.e, st: r.st, name: r.prot, gap: gap ? gap.at : null, gapLabel: gap ? gap.label : null,
        tip: "<b>" + r.prot + (r.name !== r.prot ? " (" + r.name + ")" : "") + "</b>, " + r.len + " residues<br>aligned at " + span(r) +
          ": " + whole + (r.num !== r.e - r.s + 1 ? "; " + fmt(r.num) + " nt aligned on " + fmt(r.e - r.s + 1) + " bp" : "")
      };
    });
    function keptGene(a) {
      var fromCall = win.bakta.some(function (c) {
        return c.s === a.s && c.e === a.e && c.st === a.st;
      });
      var fromRef = win.ref.some(function (r) {
        return r.s === a.s && r.e === a.e && r.st === a.st;
      });
      var how = fromCall && fromRef ? "Bakta's call, as it was" : fromCall ? "Bakta's call" : fromRef ? "the aligned protein's coordinates" : "extended to its stop codon";
      var marks = (a.pseudo ? "<br>Pseudogene: " + esc(a.note) : "") + (a.partial ? "<br>Partial at its " + a.partial + " end" : "");
      var name = a.rvs ? a.rvs.split(",")[0] : "";
      return {
        s: a.s, e: a.e, st: a.st, name: name, cls: a.pseudo ? "pseudo" : null, open: a.partial || "",
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
        return { s: f.s, e: f.e, st: f.st, tip: "A call Bakta makes on " + reference + " itself that the curated annotation does not have, aligned here at " + span(f) + "." };
      })
    };
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
        ribbons.push({ a: [0, i], b: [1, j], kind: r.s === c.s && r.e === c.e ? "agree" : "conflict", phase: key === "hero" ? "case" : null });
      });
    });
    var hero = key === "hero";
    return {
      from: look.from != null ? look.from : win.from,
      to: look.to != null ? look.to : win.to,
      labelW: hero ? 66 : 64,
      geneH: hero ? 20 : 18,
      rowGap: hero ? 30 : 22,
      top: hero ? 6 : 18,
      alt: look.alt || "",
      rows: rows,
      ribbons: ribbons,
      marks: (look.marks || []).map(function (m) {
        var named = { bakta: 0, ref: 1, amap: rows.length - 1 };
        var copy = {};
        Object.keys(m).forEach(function (k) {
          copy[k] = m[k];
        });
        copy.rows = m.rows.map(function (r) {
          return typeof r === "string" ? named[r] : r;
        });
        return copy;
      }),
      edges: look.edges || [],
      axis: true,
      tick: look.tick,
      tickLabel: look.kbTicks ? function (t) {
        return t / 1000 + " kb";
      } : tickLabel(win.to - win.from)
    };
  }

  var PHASES = [
    { label: "What Bakta called", tally: false },
    { label: "Where the H37Rv proteins align", tally: false },
    { label: "Compared, call by call", tally: true },
    { label: "What AMAP keeps", tally: true }
  ];

  function mountHero(fig) {
    var stage = fig.querySelector("[data-hero-stage]");
    var stepLabel = fig.querySelector("[data-hero-step]");
    var replay = fig.querySelector("[data-hero-replay]");
    var win = DATA.windows.hero;
    if (!win) return;
    var phase = 4;
    var svg = null;
    render(stage, function () {
      return windowSpec("hero");
    }, function (drawn) {
      svg = drawn;
      svg.setAttribute("data-phase", String(phase));
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
    fig.querySelectorAll("[data-tally]").forEach(function (el) {
      tallies[el.getAttribute("data-tally")] = el;
    });
    var timer = null;
    function setPhase(n) {
      phase = n;
      if (svg) svg.setAttribute("data-phase", String(n));
      stepLabel.textContent = n + " of 4 · " + PHASES[n - 1].label;
      var on = PHASES[n - 1].tally;
      Object.keys(tallies).forEach(function (k) {
        tallies[k].textContent = String(counts[k]);
        tallies[k].parentNode.classList.toggle("t-off", !on);
      });
    }
    function play() {
      clearTimeout(timer);
      if (reduceMotion) {
        setPhase(4);
        return;
      }
      var n = 1;
      setPhase(n);
      (function next() {
        timer = setTimeout(function () {
          n += 1;
          setPhase(n);
          if (n < 4) next();
        }, 2100);
      })();
    }
    replay.addEventListener("click", play);
    // Start when the drawing comes into view, so the reader sees it happen.
    if ("IntersectionObserver" in window && !reduceMotion) {
      setPhase(1);
      var seen = false;
      var io = new IntersectionObserver(function (entries) {
        if (!seen && entries[0].isIntersecting) {
          seen = true;
          io.disconnect();
          play();
        }
      }, { threshold: 0.35 });
      io.observe(fig);
      // A page that is never shown (a hidden tab) still ends up drawn.
      setTimeout(function () {
        if (!seen) {
          seen = true;
          io.disconnect();
          setPhase(4);
        }
      }, 9000);
    } else {
      setPhase(4);
    }
  }

  /* ── One drawing per case ── */

  function mountCase(el) {
    var names = el.getAttribute("data-fig").split(" ");
    var stage = el.querySelector("[data-fig-stage]") || el;
    var switcher = el.querySelector("[data-fig-switch]");
    var caption = el.querySelector("[data-fig-caption]");
    var current = names[0];
    var redraw = render(stage, function () {
      return windowSpec(current);
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
      text: "Builds the eggNOG database on the first run, a bacteria-only index of about 3 GB plus the metadata it downloads, and reuses it after.",
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
      name: "All genomes", proc: "GET_ALLDUP", x: 980, y: 250, label: "above", lineOnly: "dup",
      text: "All samples' duplicated genes against each other, as communities with network metrics.",
      reads: "every sample's duplicated genes", writes: "DUPLICATIONS/homologous_communities.tsv, .html", flags: "", doc: "how-it-works/duplications.html", functional: true
    },
    PANGENE: {
      name: "Pangenome", proc: "PANGENE", x: 930, y: 320, label: "below", lineOnly: "pan",
      text: "Each sample's corrected proteins aligned with miniprot onto every other sample's genome, then a pangene graph and a presence/absence matrix.",
      reads: "the corrected proteins, and every genome: as PATHOTYPR renamed it, or with --skip_typing as the sample sheet gives it, under its sample ID", writes: "PANGENE/{sample}/", flags: "", doc: "how-it-works/workflow.html#9-pangenome", exit: "PANGENE/", exitSide: "right", functional: true
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

  var LINES = [
    { key: "main", d: "M70,170 H1050" },
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

  function metroSvg(state) {
    var out = [];
    out.push('<svg class="metro" viewBox="0 0 1100 380" role="img" aria-label="The AMAP pipeline as a metro map: stations are the steps, lines are what travels between them.">');
    LINES.forEach(function (l) {
      var off = lineOff(l.key, state);
      if (l.feed) {
        // Every genome, for the pangenome: renamed by Lineage, or from the
        // sample sheet when typing is off (PANGENE_GENOME names them).
        var untyped = state.skipTyping;
        out.push('<path class="ml ml-pan' + (off ? " off" : "") + '" d="' + (untyped ? l.dUntyped : l.d) + '" style="stroke-width:3;stroke-dasharray:2 6"/>');
        out.push('<text class="exit' + (off ? " off" : "") + '" x="' + (untyped ? 34 : 270) + '" y="343">' + (untyped ? "the genomes of the sample sheet, under their sample IDs" : "the renamed genomes, from Lineage") + "</text>");
        return;
      }
      out.push('<path class="ml ml-' + l.key + (off ? " off" : "") + '" d="' + l.d + '"' + (l.dash ? ' style="stroke-width:3;stroke-dasharray:2 6"' : "") + "/>");
      if (l.endTick) {
        out.push('<line class="ml ml-' + l.key + (off ? " off" : "") + '" x1="' + l.endTick[0] + '" x2="' + l.endTick[0] + '" y1="' + (l.endTick[1] - 6) + '" y2="' + (l.endTick[1] + 6) + '" style="stroke-width:4"/>');
      }
    });
    Object.keys(STATIONS).forEach(function (key) {
      var s = STATIONS[key];
      var off = stationOff(key, state);
      var skipped = key === "PATHOTYPR" && state.skipTyping;
      var cls = "stn" + (off ? " off" : "") + (skipped ? " skipped" : "") + (state.selected === key ? " on" : "");
      out.push('<g class="' + cls + '" data-stn="' + key + '" tabindex="0" role="button" aria-label="' + esc(s.name + ", " + (state.custom && s.procCustom ? s.procCustom : s.proc)) + '">');
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
      var proc = state.custom && s.procCustom ? s.procCustom : s.proc;
      var lx = s.x;
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
      out.push('<text class="proc" x="' + lx + '" y="' + py + '" text-anchor="' + anchor + '">' + esc(proc) + "</text>");
      out.push("</g>");
    });
    out.push("</svg>");
    return out.join("");
  }

  function lineOff(key, state) {
    if (key === "main") return false;
    if (key === "ref") return false;
    if (state.skipFunctional && (key === "db" || key === "dup" || key === "pan" || key === "pan-feed" || key === "rep-feed")) return true;
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

  function mountMetro(el) {
    var stage = el.querySelector("[data-metro-stage]");
    var info = el.querySelector("[data-metro-info]");
    var writesEl = el.querySelector("[data-metro-writes]");
    var state = { skipTyping: false, skipFunctional: false, custom: false, selected: null };
    function render() {
      stage.innerHTML = metroSvg(state);
      if (writesEl) {
        writesEl.innerHTML = "This run writes " + writes(state).map(function (d) {
          return "<code>" + d + "</code>";
        }).join(" ") + ".";
      }
      describe();
    }
    function describe() {
      if (!info) return;
      var key = state.selected;
      if (!key) {
        info.innerHTML = '<p class="hint">Pick a station to see what it does, what it reads and where its files go. The switches above show what a run leaves out.</p>';
        return;
      }
      var s = STATIONS[key];
      var proc = state.custom && s.procCustom ? s.procCustom : s.proc;
      var text = state.custom && s.textCustom ? s.textCustom : s.text;
      var status = "";
      if (stationOff(key, state)) status = '<p class="note">Not run with the switches above.</p>';
      if (key === "PANGENE" && state.skipTyping && !state.skipFunctional) status = '<p class="note">With <code>--skip_typing</code> it aligns onto the genomes of the sample sheet, each under its sample ID.</p>';
      info.innerHTML =
        "<div><h3>" + esc(s.name) + "<code>" + esc(proc) + "</code></h3><p>" + esc(text) + "</p>" + status +
        (s.doc ? '<p class="note"><a href="' + s.doc + '">More on this step</a></p>' : "") + "</div>" +
        "<dl>" +
        (s.reads ? "<dt>Reads</dt><dd>" + esc(s.reads) + "</dd>" : "") +
        (s.writes ? "<dt>Writes</dt><dd>" + esc(s.writes) + "</dd>" : "") +
        (s.flags ? "<dt>Options</dt><dd>" + esc(s.flags) + "</dd>" : "") +
        "</dl>";
    }
    el.addEventListener("change", function (event) {
      var t = event.target;
      if (t.matches("[data-metro-opt]")) {
        state[t.getAttribute("data-metro-opt")] = t.checked;
        render();
      }
    });
    stage.addEventListener("click", function (event) {
      var g = event.target.closest(".stn");
      if (!g) return;
      state.selected = g.getAttribute("data-stn");
      render();
    });
    stage.addEventListener("keydown", function (event) {
      if (event.key !== "Enter" && event.key !== " ") return;
      var g = event.target.closest(".stn");
      if (!g) return;
      event.preventDefault();
      state.selected = g.getAttribute("data-stn");
      render();
      var again = stage.querySelector('[data-stn="' + state.selected + '"]');
      if (again) again.focus();
    });
    render();
  }

  /* ── The overview: what goes into a decision, and what comes out ── */

  // The hero's 12.4 kb of CDC1551 again, with the three things each decision
  // reads: Bakta's calls, the reference proteins aligned, and the calls Bakta
  // makes wrongly on the reference, aligned here.
  function inputsSpec() {
    var parts = rowsFor("hero");
    return {
      from: parts.win.from, to: parts.win.to, labelW: 84, geneH: 18, rowGap: 22, top: 8,
      alt: "The first 12.4 kb of M. tuberculosis CDC1551: Bakta's twelve calls, the eight H37Rv proteins aligned, and two of Bakta's known false positives.",
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
  function miniSpec(key) {
    var spec = windowSpec(key);
    spec.labelW = 0;
    spec.axis = false;
    spec.geneH = 14;
    spec.rowGap = 12;
    spec.top = spec.marks.length ? 16 : 4;
    spec.padRight = 2;
    spec.marks.forEach(function (m) {
      m.label = "";
    });
    spec.rows.forEach(function (row) {
      row.genes.forEach(function (g) {
        g.gapLabel = null;
      });
    });
    spec.edges = [];
    return spec;
  }

  function mountInputs(el) {
    var stage = el.querySelector("[data-fig-stage]");
    if (!DATA.windows.hero) return;
    render(stage, inputsSpec);
    attachTips(stage);
  }

  function mountOutcomes(el) {
    el.querySelectorAll("[data-mini]").forEach(function (card) {
      var key = card.getAttribute("data-mini");
      var stage = card.querySelector("[data-fig-stage]");
      if (!DATA.windows[key]) return;
      render(stage, function () {
        return miniSpec(key);
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
  // outline around the nodes of a family.
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
    var out = ['<svg class="net" viewBox="0 0 ' + W + " " + H + '" role="img" aria-label="' + esc(o.alt) + '">'];
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
    o.links.forEach(function (l) {
      var a = o.nodes[l.a];
      var b = o.nodes[l.b];
      out.push(
        '<g class="hot net-link ' + l.cls + '" tabindex="0" data-tip="' + esc(l.tip) + '">' +
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
      out.push(
        '<g class="hot net-link ' + l.cls + '" data-tip="' + esc(l.tip) + '">' +
          '<rect class="net-pill" x="' + spot.x.toFixed(1) + '" y="' + spot.y.toFixed(1) + '" width="' + lw.toFixed(1) + '" height="16" rx="8"/>' +
          '<text class="net-pct" x="' + (spot.x + lw / 2).toFixed(1) + '" y="' + (spot.y + 11.5).toFixed(1) + '" text-anchor="middle">' + esc(l.label) + "</text></g>"
      );
    });
    o.nodes.forEach(function (n) {
      var cx = X(n.x);
      var cy = Y(n.y);
      // The name goes under the node, or over it at the bottom of the frame.
      var below = cy + n.r + 14 < H - 2;
      out.push('<g class="hot net-node ' + n.cls + '" tabindex="0" data-tip="' + esc(n.tip) + '">');
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
    attachTips(stage);
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

  window.AMAPFigures = { draw: draw, render: render, tickLabel: tickLabel, attachTips: attachTips, rowsFor: rowsFor, windowSpec: windowSpec, fmt: fmt, esc: esc, stations: STATIONS, metroOrder: METRO_ORDER };
})();

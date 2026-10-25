/* How it works: one story at a time, through the seven steps AMAP takes.
   Every gene drawn is a real one, from the windows in figure-data.js, which
   tools/check_site.py runs through bin/; the words of each step are the rules
   of docs/how-it-works/correction-algorithm.md applied to those genes.
   Plain JavaScript, no build step. */
(function () {
  "use strict";

  var F = window.AMAPFigures;
  var DATA = window.AMAP_FIGURE_DATA;
  if (!F || !DATA) return;
  var reduceMotion = window.matchMedia && window.matchMedia("(prefers-reduced-motion: reduce)").matches;

  // The numbers the stories give for each rule. tools/check_site.py holds
  // them against the code; this block is JSON.
  var RULES = {"single_bp": 550, "multiple_bp": 500, "false_positive_pct": 80, "identity_pct": 80, "account_cover": 0.5, "synteny_bp": 10000, "anchors": 5, "claim_cover": 0.5, "exact_residues": 8, "copy_length": 0.9, "stop_search_bp": 300, "past_reference_residues": 10};

  var STEPS = [
    "Bakta calls the genes",
    "The reference proteins are placed",
    "Each call is compared",
    "Disagreements are settled",
    "Every gene is translated",
    "Copies are sorted out",
    "Written to the annotation"
  ];

  var NOTHING = {
    6: "One gene for each reference protein here: nothing to sort out.",
    copies: "Two genes carry the same protein only where a genome has it twice."
  };

  /* The stories. Each names the windows it draws ("panels"), says what each
     step does to them ("text", one per step) and ends on the lines AMAP
     writes ("out"). */

  var STORIES = [
    {
      key: "agree", group: "The five outcomes, in eight stories", label: "Agrees", tag: "Case 2",
      panels: [{ key: "case2", from: 1960, to: 3300, tick: 500 }],
      text: [
        "Bakta calls dnaN from 2,052 to 3,260 on the plus strand.",
        "miniprot aligns H37Rv's Rv0002, dnaN, onto CDC1551: the whole protein, from 2,052 to 3,260, stop codon included.",
        "Same start, same end, same strand and the same length: case 2, a perfect match.",
        "Nothing to settle. The call stays as Bakta made it, and is tied to Rv0002.",
        "It ends in a stop codon: nothing to change.",
        NOTHING[6],
        "It keeps Bakta's locus tag and name, with Rv0002 in its Dbxref."
      ],
      out: [
        "L4_CDC1551  AMAP  CDS  2052  3260  .  +  0  ID=L4_CDC1551_00002;Parent=gene-L4_CDC1551_00002;Name=dnaN;product=DNA polymerase III (beta chain) DnaN (DNA nucleotidyltransferase);Dbxref=AMAP_ids:Rv0002"
      ]
    },
    {
      key: "same-ends", group: "The five outcomes, in eight stories", label: "Same ends, other length", tag: "Case 3",
      panels: [{ key: "case3", from: 23700, to: 25600, tick: 500 }],
      text: [
        "Bakta calls fhaA from 23,861 to 25,426 on the minus strand: 1,566 bp.",
        "H37Rv's Rv0020c, fhaA, aligns over exactly the same stretch, the whole protein to its stop codon. Its protein is six residues longer: miniprot aligns 1,584 nt of it on these 1,566 bp, with a gap where CDC1551 lacks them.",
        "Same start, same end, same strand, but not the same length: 1,566 bp called, 1,584 nt aligned. That is case 3, an imperfect match.",
        "Nothing to settle. A case 3 call is reliable enough: it stays as Bakta made it, tied to Rv0020c.",
        "It ends in a stop codon: nothing to change.",
        NOTHING[6],
        "It keeps Bakta's locus tag, name and product, with Rv0020c in its Dbxref."
      ],
      out: [
        "L4_CDC1551  AMAP  CDS  23861  25426  .  -  0  ID=L4_CDC1551_00025;Parent=gene-L4_CDC1551_00025;Name=fhaA;product=Conserved protein with FHA domain%2C FhaA;Dbxref=AMAP_ids:Rv0020c"
      ]
    },
    {
      key: "late-start", group: "The five outcomes, in eight stories", label: "Starts late", tag: "Case 4",
      panels: [{ key: "case4-start", from: -60, to: 1640, tick: 500, marks: [{ s: 1, e: 33, rows: [0, 1], label: "33 bp", labelAt: "right", from: 3 }] }],
      text: [
        "Bakta calls dnaA from base 34 to 1,524.",
        "H37Rv's Rv0001, dnaA, aligns end to end from base 1 to 1,524.",
        "Same strand, same end, but a start 33 bp apart: case 4, a conflicting overlap. One call against one alignment.",
        "With a single overlap the reference wins, unless the call is longer by more than 550 bp and the alignment is not complete in its frame. Here the alignment is the longer of the two: AMAP keeps 1 to 1,524.",
        "It ends in a stop codon: nothing to change.",
        NOTHING[6],
        "No Bakta call has these coordinates, so it is written under the reference protein's name."
      ],
      out: ["L4_CDC1551  AMAP  CDS  1  1524  .  +  0  ID=L4_CDC1551_Rv0001;Parent=gene-L4_CDC1551_Rv0001;Name=;product=;Dbxref=AMAP_ids:Rv0001"]
    },
    {
      key: "split", group: "The five outcomes, in eight stories", label: "One gene, two calls", tag: "Case 4",
      panels: [{ key: "case4-split", from: 9850, to: 10950, tick: 250 }],
      text: [
        "Bakta calls two short genes, from 9,963 to 10,160 and from 10,523 to 10,828, both hypothetical proteins.",
        "H37Rv's Rv0007 aligns end to end over both, from 9,914 to 10,828.",
        "Both calls overlap the one alignment with other coordinates: case 4, with several calls on one alignment.",
        "The call that overlaps the alignment most is weighed against it. The reference wins if it is at least as long, if the call is longer by less than 500 bp, if the alignment is complete in the call's frame, or if the call spans other reference genes too. Rv0007 is longer than either call: one gene goes in, 9,914 to 10,828.",
        "It ends in a stop codon: nothing to change.",
        NOTHING[6],
        "Written under Rv0007, as no Bakta call has its coordinates."
      ],
      out: ["L4_CDC1551  AMAP  CDS  9914  10828  .  +  0  ID=L4_CDC1551_Rv0007;Parent=gene-L4_CDC1551_Rv0007;Name=;product=;Dbxref=AMAP_ids:Rv0007"]
    },
    {
      key: "longer", group: "The five outcomes, in eight stories", label: "Much longer", tag: "Case 4",
      panels: [{ key: "case4-longer", from: 2267500, to: 2272800, tick: 1000, marks: [{ s: 2267753, e: 2271055, rows: [0, 1], label: "3,303 bp", from: 3 }] }],
      text: [
        "Bakta calls a 4,821 bp site-specific DNA methyltransferase on the minus strand.",
        "H37Rv's Rv2024c, 516 residues, aligns over 1,518 bp at one end of it: residues 1 to 506, short of the protein's end, with no stop codon.",
        "One call against one alignment, with other coordinates: case 4.",
        "A call longer by more than 550 bp wins over an alignment that is not complete in its frame: such an alignment stopped short, and dropping a whole gene for it would be the worse mistake. Bakta's call, 3,303 bp longer, is kept.",
        "It ends in a stop codon: nothing to change.",
        NOTHING[6],
        "It keeps Bakta's locus tag and product, and takes Rv2024c's name in its Dbxref: no other gene carries it, and the alignment spans more than half of the protein."
      ],
      out: ["L4_CDC1551  AMAP  CDS  2267753  2272573  .  -  0  ID=L4_CDC1551_02139;Parent=gene-L4_CDC1551_02139;Name=;product=site-specific DNA-methyltransferase (adenine-specific);Dbxref=AMAP_ids:Rv2024c"]
    },
    {
      key: "new", group: "The five outcomes, in eight stories", label: "Not in the reference", tag: "Case 1",
      panels: [{ key: "case1-kept", from: 1979300, to: 1985800, tick: 2000 }],
      text: [
        "Bakta calls four genes in a row: a mannosyltransferase, a sulfite oxidase-like oxidoreductase, an MMPL transporter and a hypothetical protein.",
        "No H37Rv protein aligns anywhere along them.",
        "A call with no alignment over it on its strand is case 1, potentially novel.",
        "A case 1 call goes if Bakta's known false positives cover 80% of it, or if its protein matches the H37Rv genome at 80% identity or more. None is covered, three of the proteins do not align to H37Rv at all and the fourth matches at 60%: all four are kept.",
        "They end in stop codons: nothing to change.",
        "No reference protein to share: nothing to sort out.",
        "They keep Bakta's locus tags and products, with no reference name."
      ],
      out: [
        "L4_CDC1551  AMAP  CDS  1980912  1982015  .  +  0  ID=L4_CDC1551_01842;Parent=gene-L4_CDC1551_01842;Name=;product=Sulfite oxidase-like oxidoreductase;Dbxref=AMAP_ids:",
        "L4_CDC1551  AMAP  CDS  1982216  1985053  .  +  0  ID=L4_CDC1551_01843;Parent=gene-L4_CDC1551_01843;Name=;product=Membrane transport protein MMPL domain-containing protein;Dbxref=AMAP_ids:"
      ]
    },
    {
      key: "false", group: "The five outcomes, in eight stories", label: "A known false call", tag: "Case 1",
      panels: [{ key: "case1-dropped", from: 11330, to: 11620, tick: 100 }],
      text: [
        "Bakta calls a 108 bp hypothetical protein on the minus strand.",
        "No H37Rv protein aligns here. Bakta's known false positives are aligned too: calls Bakta makes on H37Rv itself that its curated annotation does not have. One covers this call end to end.",
        "Case 1, potentially novel.",
        "Covered 100% by a known false positive, more than the 80% that rules a call out: it is dropped.",
        "Nothing left here to translate.",
        "Nothing to sort out.",
        "Nothing is written here."
      ],
      out: []
    },
    {
      key: "missed", group: "The five outcomes, in eight stories", label: "Missed by Bakta", tag: "Missing",
      panels: [{ key: "missing", from: 33100, to: 34420, tick: 250, marks: [{ s: 33526, e: 33738, rows: [0, 1], label: "no call", kind: "ref", from: 3 }] }],
      text: [
        "Bakta calls a gene from 33,168 to 33,497 and bioF2 from 34,239, and nothing between them.",
        "H37Rv's Rv0031, 71 residues, aligns between the two, end to end, from 33,526 to 33,738.",
        "An alignment with no call over it on its strand: a missing reference gene. The two calls match Rv0030 and Rv0032 exactly, case 2.",
        "A protein is missing when no gene accounts for any of its placements: none sits on it, in the same frame, covering half of it. So Rv0031 goes back in at its best placement.",
        "It ends in a stop codon: nothing to change.",
        NOTHING[6],
        "Written under Rv0031."
      ],
      out: ["L4_CDC1551  AMAP  CDS  33526  33738  .  +  0  ID=L4_CDC1551_Rv0031;Parent=gene-L4_CDC1551_Rv0031;Name=;product=;Dbxref=AMAP_ids:Rv0031"]
    },
    {
      key: "paralog", group: "The harder ones", label: "On its paralog's locus", tag: "Paralogs",
      panels: [
        { key: "paralogs-1", from: 400150, to: 402300, tick: 500, title: "Near 400 kb" },
        { key: "paralogs-2", from: 607300, to: 609950, tick: 500, title: "Near 608 kb" }
      ],
      text: [
        "Two loci of CDC1551, about 208 kb apart. Bakta calls a REP13E12 repeat protein at each, from 400,321 and from 607,987.",
        "H37Rv's Rv0336 and Rv0515 are paralogs, and miniprot places each at both loci, as it reports every placement within 5% of a protein's best. The proteins placed only once around each locus, five on each side, say where it sits in H37Rv: the first at about 400 kb, where Rv0336 belongs, the second at about 607 kb, where Rv0515 does. Each protein's placement on the other's locus is dropped.",
        "Each REP13E12 call now overlaps one alignment, starting 69 bp after it: case 4 at both loci. The calls around them match their proteins exactly, case 2.",
        "The reference wins both: Rv0336 at 400,252 and Rv0515 at 607,918.",
        "Both end in stop codons: nothing to change.",
        NOTHING[6],
        "Each locus is written under its own protein."
      ],
      out: [
        "L4_CDC1551  AMAP  CDS  400252  401763  .  +  0  ID=L4_CDC1551_Rv0336;Parent=gene-L4_CDC1551_Rv0336;Name=;product=;Dbxref=AMAP_ids:Rv0336",
        "L4_CDC1551  AMAP  CDS  607918  609429  .  +  0  ID=L4_CDC1551_Rv0515;Parent=gene-L4_CDC1551_Rv0515;Name=;product=;Dbxref=AMAP_ids:Rv0515"
      ]
    },
    {
      key: "small", group: "The harder ones", label: "Too small for miniprot", tag: "Exact search",
      panels: [{ key: "small-protein", from: 3056150, to: 3058450, tick: 500, marks: [{ s: 3056890, e: 3056949, rows: [0, 1], label: "IbsC, by the exact search", kind: "ref", from: 2 }] }],
      text: [
        "E. coli K-12 against its own annotation. Bakta calls fau and serA, and nothing in the 388 bp between them.",
        "miniprot places fau and serA. IbsC, a toxic peptide of 19 residues, is too small to seed an alignment of its own: miniprot placed it only at 3,194,766, on its neighbour IbsD's locus, where 3 of its residues differ. Every protein of 8 residues or more without an exact, complete placement is looked for in the six-frame translation of the genome, a start codon, its residues, then a stop: IbsC is here, from 3,056,890 to 3,056,949 on the minus strand.",
        "fau: Bakta's call starts 54 bp before the alignment, case 4. serA matches exactly, case 2. IbsC has no call over it: a missing reference gene.",
        "fau takes the reference's start, and IbsC goes in where the exact search found it.",
        "All three end in stop codons: nothing to change.",
        NOTHING[6],
        "IbsC is written under its own name, b4665."
      ],
      out: ["K12self  AMAP  CDS  3056890  3056949  .  -  0  ID=K12self_b4665;Parent=gene-K12self_b4665;Name=;product=;Dbxref=AMAP_ids:b4665"]
    },
    {
      key: "contig-end", group: "The harder ones", label: "Cut by a contig end", tag: "Partial",
      panels: [
        { key: "contig-end", from: 100250, to: 102146, tick: 500, title: "The end of contig_10", edges: [{ at: 102146.5, label: "end of contig_10", side: "left", from: 1 }],
          marks: [{ s: 102120, e: 102146, rows: [3, 3], label: "partial: 3' end", labelAt: "left", from: 5 }] },
        { key: "contig-end-2", from: 3700, to: 4591, tick: 250, title: "The end of contig_9", edges: [{ at: 4591.5, label: "end of contig_9", side: "left", from: 1 }],
          marks: [{ s: 4565, e: 4591, rows: [3, 3], label: "partial: 5' end", labelAt: "left", from: 5 }] }
      ],
      text: [
        "CDC1551 cut into 40 contigs. Near the end of contig_10 Bakta calls Rv1240's gene and two short genes on the minus strand, and nothing on the plus strand. At the end of contig_9 it calls sugC and a short hypothetical protein.",
        "H37Rv's Rv1239c, 367 residues, aligns on the plus strand from 101,478 to the very end of contig_10, 102,146: its first 223 residues. Its last 144, which miniprot does not report beside the larger piece, are aligned again on their own, and land at the very end of contig_9: the other piece.",
        "Rv1240's call and sugC's match exactly, case 2. The two short calls on contig_10 have nothing aligned, case 1. The short call on contig_9 overlaps the piece with other coordinates, case 4. On contig_10 Rv1239c has no call over it on its strand: a missing reference gene.",
        "Rv1239c goes back in at its placement on contig_10, and on contig_9 the piece wins over the shorter call, as the reference does. The two short calls on contig_10 go: their proteins match the H37Rv genome at 98 and 99% identity, so they are not new.",
        "Each piece's frame reaches its contig's end with no stop codon, from within 30 bp of it: the assembly cut the gene. The piece on contig_10 is marked partial at its 3' end, the one on contig_9 at its 5' end. Neither is a pseudogene.",
        "Two genes carry Rv1239c, but they are the pieces of one cut gene, not copies, and both stay. Their IDs must differ: one is written Rv1239c and the other Rv1239c_2, both with Rv1239c in their Dbxref.",
        "The cut ends are written the way NCBI writes them: partial=true with end_range in the GFF3, and a > before the cut end in the GenBank location."
      ],
      out: [
        "contig_9   AMAP  CDS  4160  4591  .  -  0  ID=L4_CDC1551d_Rv1239c;Parent=gene-L4_CDC1551d_Rv1239c;Name=;product=;Dbxref=AMAP_ids:Rv1239c;partial=true;end_range=4591,.",
        "contig_10  AMAP  CDS  101478  102146  .  +  0  ID=L4_CDC1551d_Rv1239c_2;Parent=gene-L4_CDC1551d_Rv1239c_2;Name=;product=;Dbxref=AMAP_ids:Rv1239c;partial=true;end_range=102146,.",
        "     CDS             complement(4160..>4591)",
        "     CDS             101478..>102146"
      ]
    },
    {
      key: "frameshift", group: "The harder ones", label: "Runs into the next gene", tag: "Pseudogene",
      panels: [{ key: "frameshift", from: 688900, to: 690500, tick: 500, marks: [{ s: 690104, e: 690306, rows: [3, 3], label: "67 residues past its protein", from: 5 }] }],
      text: [
        "M. bovis AF2122/97. Bakta calls one gene from 689,275 to 690,306, between mce2A's and mce2C's.",
        "Two H37Rv proteins align inside it: Rv0590 end to end, with one base too many for its frame, and Rv0590A after it.",
        "The call overlaps both alignments with other coordinates: case 4.",
        "A call over more than one alignment is several reference genes read as one, and the reference wins: Rv0590 and Rv0590A go in as two genes.",
        "Rv0590 has no stop codon where its protein ends. The next stop in its frame is 67 residues past the end of the protein, at 690,306, the end of Rv0590A and of Bakta's call. An extension may reach 10 residues past the protein. So Rv0590 keeps its extent and is flagged a pseudogene, with the reason.",
        NOTHING[6],
        "The pseudogene is written the way NCBI writes one: pseudo=true and the reason in Note, /pseudo and /note in the GenBank, with no translation."
      ],
      out: [
        "La1_Mbovis  AMAP  CDS  689275  690103  .  +  0  ID=La1_Mbovis_Rv0590;Parent=gene-La1_Mbovis_Rv0590;Name=;product=;Dbxref=AMAP_ids:Rv0590;pseudo=true;Note=frameshift or lost stop codon: its frame runs on 67 residues past the end of the reference protein",
        "La1_Mbovis  AMAP  CDS  690052  690306  .  +  0  ID=La1_Mbovis_Rv0590A;Parent=gene-La1_Mbovis_Rv0590A;Name=;product=;Dbxref=AMAP_ids:Rv0590A"
      ]
    },
    {
      key: "copies", group: "The harder ones", label: "Two real copies", tag: "Copies",
      labelIds: true,
      panels: [
        { key: "copies-1", from: 616700, to: 621600, tick: 1000, title: "At 617 kb" },
        { key: "copies-2", from: 2045700, to: 2050700, tick: 1000, title: "At 2,046 kb" }
      ],
      text: [
        "E. coli O157:H7 Sakai, against K-12. Bakta calls a 4.2 kb RhsD gene at 617,058 and another at 2,046,050, 1.4 Mb apart.",
        "K-12's RhsD (b0497) aligns at both, over most of each gene but short of its end.",
        "Each RhsD call overlaps its alignment with other coordinates, case 4. The other calls have nothing aligned, case 1: the three at the second locus are kept, and the short one at the first goes, as one of the calls Bakta makes on K-12 itself covers it.",
        "Each call is longer than its alignment by less than 550 bp: the reference wins at both.",
        "Neither alignment reaches a stop codon. Each gene is extended to the next stop in its frame, which completes the protein and is where Bakta's calls end.",
        "Two genes now carry b0497. A genome can have a gene twice: the second stays, as it is at least 90% as long as the first, clear of it and on no other reference gene's locus. The copy later in the genome is named b0497_2, and the duplication analysis counts it.",
        "Both copies are written, each with an ID of its own."
      ],
      out: [
        "Sakai  AMAP  CDS  617058  621254  .  +  0  ID=Sakai_b0497;Parent=gene-Sakai_b0497;Name=;product=;Dbxref=AMAP_ids:b0497",
        "Sakai  AMAP  CDS  2046050  2050252  .  +  0  ID=Sakai_b0497_2;Parent=gene-Sakai_b0497_2;Name=;product=;Dbxref=AMAP_ids:b0497"
      ]
    }
  ];

  /* ── Drawing a story ── */

  function sameGenes(a, b) {
    var key = function (g) {
      return [g.s, g.e, g.st, g.cls || "", g.open || ""].join(":");
    };
    return a.map(key).join("|") === b.map(key).join("|");
  }

  function overlap(a, b) {
    return Math.max(0, Math.min(a.e, b.e) - Math.max(a.s, b.s) + 1);
  }

  function copyOf(g, name) {
    var c = {};
    Object.keys(g).forEach(function (k) {
      c[k] = g[k];
    });
    c.name = name;
    return c;
  }

  function panelSpec(story, panel) {
    var parts = F.rowsFor(panel.key);
    var win = parts.win;
    var reference = parts.reference;
    var rows = [{ label: "Bakta", cls: "bakta", genes: parts.bakta, empty: parts.bakta.length ? null : "no call", from: 1 }];
    var refRow = rows.length;
    rows.push({ label: reference + ",\naligned", cls: "ref", genes: parts.ref, empty: parts.ref.length ? null : "no " + reference + " protein aligns here", from: 2 });
    // What the paralog step and the exact search changed: placements the
    // first dropped, and the ones the second added.
    if (win.placed) {
      var kept = {};
      win.ref.forEach(function (r) {
        kept[r.prot + ":" + r.s] = true;
      });
      var dropped = win.placed.filter(function (p) {
        return !kept[p.prot + ":" + p.s];
      });
      if (dropped.length) {
        rows.push({
          label: "Dropped:\nanother's locus", cls: "drop", from: 2, genes: dropped.map(function (p) {
            return { s: p.s, e: p.e, st: p.st, name: p.prot, tip: "<b>" + p.prot + "</b> placed at " + F.fmt(p.s) + "-" + F.fmt(p.e) + ": not where " + p.prot + " belongs, and on a locus another protein owns. Dropped." };
          })
        });
      }
    }
    if (parts.bad.length) rows.push({ label: "Known false\npositives", cls: "drop", genes: parts.bad, from: 2 });
    // With labelIds the genes go by their IDs. The copy's _2 comes from the
    // step that sorts out copies, so until then it goes by the protein's
    // name alone.
    var idOf = function (id) {
      return id.replace(win.base + "_", "");
    };
    var corrected = story.labelIds ? parts.corrected.map(function (g) {
      return copyOf(g, idOf(g.id));
    }) : parts.corrected;
    var final = parts.amap;
    var renamed = null;
    if (story.labelIds) {
      renamed = parts.amap.map(function (g, i) {
        return copyOf(g, idOf(win.amap[i].id));
      });
      final = renamed.map(function (g, i) {
        return copyOf(g, win.amap[i].rvs ? g.name.replace(/_\d+$/, "") : g.name);
      });
      if (renamed.map(function (g) {
        return g.name;
      }).join("|") === final.map(function (g) {
        return g.name;
      }).join("|")) renamed = null;
    }
    var changed = !sameGenes(corrected, final);
    var amapRow = rows.length;
    var layers = [corrected];
    if (changed) layers.push(final);
    if (renamed) layers.push(renamed);
    // Each layer of what AMAP keeps from its own step: the classification's
    // genes at 4, the translated ones at 5, the copies named at 6.
    var steps = [4, changed ? 5 : null, renamed ? 6 : null].filter(function (n) {
      return n;
    });
    layers.forEach(function (genes, i) {
      rows.push({
        label: "AMAP", cls: "amap", genes: genes, overlay: i > 0, empty: i === 0 && !genes.length ? "nothing kept" : null,
        from: steps[i], until: i + 1 < layers.length ? steps[i + 1] - 1 : null
      });
    });
    var ribbons = [];
    win.bakta.forEach(function (c, i) {
      win.ref.forEach(function (r, j) {
        // Neighbours share a codon or two at their ends; a band for that
        // would cross the drawing for nothing.
        if (r.st !== c.st || overlap(c, r) <= 30) return;
        ribbons.push({ a: [0, i], b: [refRow, j], kind: r.s === c.s && r.e === c.e ? "agree" : "conflict", from: 3 });
      });
    });
    var marks = (panel.marks || []).map(function (m) {
      var copy = {};
      Object.keys(m).forEach(function (k) {
        copy[k] = m[k];
      });
      // Marks name rows by their place among Bakta, the reference and AMAP.
      copy.rows = m.rows.map(function (r) {
        return r === 3 ? amapRow : r === 1 ? refRow : r;
      });
      return copy;
    });
    return F.nameSubjects({
      from: panel.from, to: panel.to, labelW: 96, geneH: 18, rowGap: 24, top: 20,
      alt: story.label + ": " + story.text[0],
      refName: reference, phased: true,
      rows: rows, ribbons: ribbons, marks: marks, edges: panel.edges || [],
      axis: true, tick: panel.tick,
      tickLabel: F.tickLabel(panel.to - panel.from)
    });
  }

  // Which steps change the drawing of a story, or the lines it ends on: the
  // others go by quicker when it plays, and say so in the list of steps.
  function changesOf(story) {
    var ch = [false, true, false, false, false, false, false, false];
    var note = function (o) {
      if (o.from) ch[o.from] = true;
      if (o.until && o.until < 7) ch[o.until + 1] = true;
    };
    story.panels.forEach(function (panel) {
      var spec = panelSpec(story, panel);
      spec.rows.forEach(note);
      spec.ribbons.forEach(note);
      spec.marks.forEach(note);
      spec.edges.forEach(note);
      spec.rows.forEach(function (row) {
        row.genes.forEach(function (g) {
          // Calls with nothing aligned are dashed at 3, and those dropped
          // fade at 4.
          if (/c1-call/.test(g.phaseCls || "")) ch[3] = true;
          if (/ph-dim/.test(g.phaseCls || "")) ch[4] = true;
        });
      });
    });
    if (story.out.length) ch[7] = true;
    return ch;
  }

  // A step's text, with the species named in italics as everywhere else,
  // and the words that stay together as the docs keep them (F.tie).
  function setText(el, text) {
    el.textContent = "";
    String(text).split(/\b([A-Z]\. (?:tuberculosis|bovis|coli))\b/).forEach(function (part, i) {
      if (!part) return;
      if (i % 2) {
        var em = document.createElement("i");
        em.textContent = part.replace(" ", "\u00a0");
        el.appendChild(em);
      } else {
        el.appendChild(document.createTextNode(F.tie(part)));
      }
    });
  }

  function words(text) {
    return String(text).split(/\s+/).filter(Boolean).length;
  }

  // How long Play stays on a step: time to read it, between 3 and 22 s, and
  // 2 s for a step that changes nothing.
  function dwell(story, n, changes) {
    var ch = changes || changesOf(story);
    if (!ch[n]) return 2000;
    return Math.max(3000, Math.min(22000, 1500 + 300 * words(story.text[n - 1])));
  }

  /* ── The page ── */

  var root = document.querySelector("[data-how]");
  if (!root) return;
  var picker = root.querySelector("[data-stories]");
  var stepList = root.querySelector("[data-steps]");
  var card = root.querySelector(".how-stage-card");
  var stage = root.querySelector("[data-stage]");
  var textBox = root.querySelector(".how-text");
  var kicker = root.querySelector("[data-kicker]");
  var title = root.querySelector("[data-step-title]");
  var text = root.querySelector("[data-step-text]");
  var outBox = root.querySelector("[data-out]");
  var progress = root.querySelector("[data-progress]");
  var prev = root.querySelector("[data-prev]");
  var next = root.querySelector("[data-next]");
  var playBtn = root.querySelector("[data-play]");

  // Each button holds both of its words and shows one, so it keeps the size
  // of the longer and nothing beside or under it moves when the word changes.
  function say(btn, words, second) {
    if (!btn.querySelector(".say-2")) {
      btn.innerHTML = '<span class="say-1">' + F.esc(words[0]) + '</span><span class="say-2">' + F.esc(words[1]) + "</span>";
    }
    btn.classList.toggle("said-2", !!second);
  }
  var genomeLabel = root.querySelector("[data-genome]");
  var legendBox = root.querySelector(".how-legend");

  function find(key) {
    return STORIES.filter(function (s) {
      return s.key === key;
    })[0];
  }

  // The story a page opens on, without a story in its address.
  var story = find("agree") || STORIES[0];
  var step = 1;
  var timer = null;
  var changes = changesOf(story);

  // The page carries the chips and the steps, so it does not move when this
  // runs; a story the page lacks still gets a chip.
  STORIES.forEach(function (s) {
    if (picker.querySelector('[data-story="' + s.key + '"]')) return;
    var group = null;
    picker.querySelectorAll(".how-group").forEach(function (g) {
      if (g.getAttribute("data-group") === s.group) group = g.querySelector(".seg");
    });
    if (!group) {
      var g = document.createElement("div");
      g.className = "how-group";
      g.setAttribute("data-group", s.group);
      g.innerHTML = '<p class="q-title">' + F.esc(s.group) + '</p><div class="seg"></div>';
      picker.appendChild(g);
      group = g.querySelector(".seg");
    }
    var a = document.createElement("a");
    a.href = "how-it-works.html#" + s.key;
    a.setAttribute("data-story", s.key);
    a.innerHTML = F.esc(s.label) + ' <span class="how-tag">' + F.esc(s.tag) + "</span>";
    group.appendChild(a);
  });
  if (!stepList.querySelector("[data-step]")) {
    STEPS.forEach(function (name, i) {
      var li = document.createElement("li");
      li.innerHTML = '<button type="button" data-step="' + (i + 1) + '"><span class="how-n">' + (i + 1) + "</span>" + F.esc(name) + "</button>";
      stepList.appendChild(li);
    });
  }

  function address() {
    try {
      history.replaceState(null, "", "#" + story.key + (step > 1 ? "/" + step : ""));
    } catch (e) {
      /* file:// in a strict browser */
    }
  }

  function flash(el) {
    if (reduceMotion) return;
    el.classList.remove("just-on");
    void el.getBoundingClientRect();
    el.classList.add("just-on");
    setTimeout(function () {
      el.classList.remove("just-on");
    }, 750);
  }

  function showStep(moved) {
    stage.querySelectorAll("[data-from]").forEach(function (el) {
      var from = +el.getAttribute("data-from");
      var until = +el.getAttribute("data-until") || 99;
      var off = step < from || step > until;
      var was = el.classList.contains("step-off");
      el.classList.toggle("step-off", off);
      if (moved && was && !off) flash(el);
    });
    stage.querySelectorAll("[data-ghost-from]").forEach(function (el) {
      var ghost = step < +el.getAttribute("data-ghost-from");
      el.classList.toggle("ghost", ghost);
      // A row still to come is named when it arrives, not before.
      if (ghost) el.setAttribute("aria-hidden", "true");
      else el.removeAttribute("aria-hidden");
    });
    // The calls with nothing aligned are dashed once compared (3), and the
    // ones dropped fade once the disagreements are settled (4).
    stage.querySelectorAll("svg").forEach(function (svg) {
      svg.setAttribute("data-phase", String(Math.min(step, 4)));
      F.refocus(svg);
    });
    stepList.querySelectorAll("[data-step]").forEach(function (b) {
      var n = +b.getAttribute("data-step");
      b.setAttribute("aria-current", n === step ? "step" : "false");
      b.classList.toggle("done", n < step);
    });
    kicker.textContent = "Step " + step + " of 7 · " + story.label;
    title.textContent = STEPS[step - 1];
    setText(text, story.text[step - 1]);
    outBox.hidden = step < 7 || !story.out.length;
    if (!timer) {
      progress.style.transition = "";
      progress.style.width = (step / 7) * 100 + "%";
    }
    prev.setAttribute("aria-disabled", step === 1 ? "true" : "false");
    say(next, ["Next", "Start again"], step === 7);
  }

  // The text under the drawing keeps the height of the story's longest step,
  // so nothing below it moves while stepping. On a phone the step 7 file
  // excerpt is left out of that height: it would open every story on a
  // screenful of empty space, and the controls sit above the text, so what
  // it pushes down is only the page under it.
  function fitText() {
    var w = textBox.offsetWidth;
    if (!w) return;
    var probe = textBox.cloneNode(true);
    probe.removeAttribute("aria-live");
    probe.setAttribute("aria-hidden", "true");
    probe.style.cssText = "position:absolute;visibility:hidden;pointer-events:none;left:0;top:0;min-height:0;width:" + w + "px";
    card.appendChild(probe);
    var parts = ["[data-kicker]", "[data-step-title]", "[data-step-text]"].map(function (sel) {
      return probe.querySelector(sel);
    });
    var out = probe.querySelector("[data-out]");
    out.querySelector("pre").textContent = story.out.join("\n");
    var tallest = 0;
    for (var n = 1; n <= 7; n++) {
      parts[0].textContent = "Step " + n + " of 7 · " + story.label;
      parts[1].textContent = STEPS[n - 1];
      setText(parts[2], story.text[n - 1]);
      out.hidden = n < 7 || !story.out.length || w < 600;
      tallest = Math.max(tallest, probe.offsetHeight);
    }
    probe.remove();
    textBox.style.minHeight = tallest + "px";
  }

  // center: bring the story's chip into the middle of its row, when the
  // reader picked the story or a link named it; on a first visit the rows
  // start at their first chip.
  function drawStory(center) {
    stop();
    stage.innerHTML = "";
    var used = {};
    story.panels.forEach(function (panel) {
      var box = document.createElement("div");
      box.className = "how-panel";
      if (panel.title) {
        var h = document.createElement("p");
        h.className = "how-panel-title";
        h.textContent = panel.title;
        box.appendChild(h);
      }
      var holder = document.createElement("div");
      holder.className = "how-draw";
      box.appendChild(holder);
      stage.appendChild(box);
      F.render(holder, function () {
        return panelSpec(story, panel);
      }, function (svg, spec) {
        Object.keys(spec._used || {}).forEach(function (k) {
          used[k] = 1;
        });
        showStep();
      });
      F.attachTips(holder);
    });
    var win = DATA.windows[story.panels[0].key];
    // Species in italics, as everywhere else on the site.
    genomeLabel.innerHTML = F.esc(win.genome).replace(/^([A-Z]\.) ([a-z]+)/, "<i>$1&nbsp;$2</i>") + (/coli/.test(win.genome) ? ", against K-12" : ", against H37Rv");
    if (legendBox) legendBox.innerHTML = F.legend(used, /coli/.test(win.genome) ? "K-12" : "H37Rv");
    outBox.querySelector("pre").textContent = story.out.join("\n");
    picker.querySelectorAll("[data-story]").forEach(function (b) {
      if (b.getAttribute("data-story") === story.key) b.setAttribute("aria-current", "true");
      else b.removeAttribute("aria-current");
    });
    // Steps that change nothing in this story show hollow.
    changes = changesOf(story);
    stepList.querySelectorAll("[data-step]").forEach(function (b) {
      var n = +b.getAttribute("data-step");
      var same = !changes[n];
      b.classList.toggle("same", same);
      var note = b.querySelector(".how-same");
      if (note) note.remove();
      // A hollow step is named with its note, as one phrase.
      if (same) b.setAttribute("aria-label", n + " " + STEPS[n - 1] + ", no change in this story");
      else b.removeAttribute("aria-label");
      if (same) b.title = "No change in this story";
      else b.removeAttribute("title");
    });
    showStep();
    fitText();
    // The picked story's chip in view, in a row that scrolls sideways.
    var chip = center && picker.querySelector('[data-story="' + story.key + '"]');
    var row = chip && chip.parentNode;
    if (row && row.scrollWidth > row.clientWidth + 2) {
      row.scrollLeft = Math.max(0, chip.offsetLeft - row.offsetLeft - (row.clientWidth - chip.offsetWidth) / 2);
    }
  }

  // The drawing in view when a story is picked or opened by a link, when it
  // is not already.
  function reveal() {
    var r = card.getBoundingClientRect();
    if (r.top < 64 || r.top + Math.min(r.height, 320) > window.innerHeight) {
      card.scrollIntoView({ block: "start", behavior: reduceMotion ? "auto" : "smooth" });
    }
  }

  function stop() {
    clearTimeout(timer);
    timer = null;
    say(playBtn, ["Play", "Pause"], false);
    card.classList.remove("playing");
    progress.style.transition = "";
    progress.style.width = (step / 7) * 100 + "%";
  }

  function go(n) {
    var before = step;
    step = Math.max(1, Math.min(7, n));
    showStep(step !== before);
    address();
  }

  // While it plays, the bar runs on to the next step in the time the step
  // is given.
  function tick() {
    if (step >= 7) {
      stop();
      return;
    }
    var ms = dwell(story, step, changes);
    progress.style.transition = "none";
    progress.style.width = (step / 7) * 100 + "%";
    void progress.offsetWidth;
    progress.style.transition = "width " + ms + "ms linear";
    progress.style.width = ((step + 1) / 7) * 100 + "%";
    timer = setTimeout(function () {
      go(step + 1);
      tick();
    }, ms);
  }

  picker.addEventListener("click", function (event) {
    var b = event.target.closest("[data-story]");
    if (!b || event.button || event.metaKey || event.ctrlKey || event.shiftKey || event.altKey) return;
    event.preventDefault();
    stop();
    story = find(b.getAttribute("data-story")) || story;
    step = 1;
    drawStory(true);
    address();
    reveal();
  });
  stepList.addEventListener("click", function (event) {
    var b = event.target.closest("[data-step]");
    if (!b) return;
    stop();
    go(+b.getAttribute("data-step"));
  });
  prev.addEventListener("click", function () {
    if (step === 1) return;
    stop();
    go(step - 1);
  });
  next.addEventListener("click", function () {
    stop();
    go(step === 7 ? 1 : step + 1);
  });
  playBtn.addEventListener("click", function () {
    if (timer) {
      stop();
      return;
    }
    if (step === 7) go(1);
    say(playBtn, ["Play", "Pause"], true);
    card.classList.add("playing");
    tick();
  });
  root.addEventListener("keydown", function (event) {
    if (/input|textarea|select/i.test(event.target.tagName) || event.altKey || event.ctrlKey || event.metaKey) return;
    // Inside a drawing the arrows move from gene to gene.
    if (event.target.closest && event.target.closest("svg")) return;
    if (event.key === "ArrowRight") {
      stop();
      go(step + 1);
    } else if (event.key === "ArrowLeft") {
      stop();
      go(step - 1);
    }
  });
  // Which keys do what, while the keyboard is in the walk-through.
  root.addEventListener("focusin", function (event) {
    root.setAttribute("data-keys", event.target.closest && event.target.closest("svg") ? "genes" : "steps");
  });
  root.addEventListener("focusout", function (event) {
    if (!event.relatedTarget || !root.contains(event.relatedTarget)) root.removeAttribute("data-keys");
  });
  if ("ResizeObserver" in window) {
    var fitted = 0;
    new ResizeObserver(function () {
      if (Math.abs(textBox.offsetWidth - fitted) < 2) return;
      fitted = textBox.offsetWidth;
      fitText();
    }).observe(textBox);
  }

  // A story named in the address opens first, so a link can point at one,
  // and #key/4 at its fourth step.
  function wanted() {
    var m = /^#([\w-]+)(?:\/([1-7]))?$/.exec(location.hash || "");
    var found = m && find(m[1]);
    return found ? { story: found, step: +(m[2] || 1) } : null;
  }
  var want = wanted();
  if (want) {
    story = want.story;
    step = want.step;
  }
  drawStory(!!want);
  if (want) reveal();
  window.addEventListener("hashchange", function () {
    var w = wanted();
    if (!w || (w.story === story && w.step === step)) return;
    stop();
    var redraw = w.story !== story;
    story = w.story;
    step = w.step;
    if (redraw) drawStory(true);
    else showStep(true);
    reveal();
  });

  window.AMAPHow = {
    stories: STORIES, rules: RULES,
    dwell: function (key, n) {
      var s = find(key);
      return s ? dwell(s, n) : null;
    }
  };
})();

/* AMAP website: theme, tabs, copy buttons, the command builder, the embedded
   report and the documentation pages' search and contents. Plain JavaScript,
   no build step. The drawings are in figures.js. */
(function () {
  "use strict";

  /* ── Theme ── */

  // Dark unless the reader chose light.
  var THEME_KEY = "amap.theme";

  function currentTheme() {
    return document.documentElement.dataset.theme === "light" ? "light" : "dark";
  }

  function remember(theme) {
    try {
      localStorage.setItem(THEME_KEY, theme);
    } catch (e) {
      /* private mode: the choice lasts for this page only */
    }
  }

  document.querySelectorAll("[data-theme-toggle]").forEach(function (button) {
    button.addEventListener("click", function () {
      var next = currentTheme() === "dark" ? "light" : "dark";
      document.documentElement.dataset.theme = next;
      remember(next);
    });
  });

  /* ── Tabs ── */

  document.querySelectorAll('[role="tablist"]').forEach(function (list) {
    var tabs = Array.prototype.slice.call(list.querySelectorAll('[role="tab"]'));
    function select(tab) {
      tabs.forEach(function (t) {
        var on = t === tab;
        t.setAttribute("aria-selected", on ? "true" : "false");
        t.tabIndex = on ? 0 : -1;
        var panel = document.getElementById(t.getAttribute("aria-controls"));
        if (panel) panel.hidden = !on;
      });
    }
    tabs.forEach(function (tab, i) {
      tab.tabIndex = tab.getAttribute("aria-selected") === "true" ? 0 : -1;
      tab.addEventListener("click", function () {
        select(tab);
      });
      tab.addEventListener("keydown", function (event) {
        var step = event.key === "ArrowRight" ? 1 : event.key === "ArrowLeft" ? -1 : 0;
        if (!step) return;
        event.preventDefault();
        var next = tabs[(i + step + tabs.length) % tabs.length];
        select(next);
        next.focus();
      });
    });
  });

  /* ── Copy buttons ── */

  function plainText(node) {
    var clone = node.cloneNode(true);
    clone.querySelectorAll(".p").forEach(function (prompt) {
      prompt.remove();
    });
    return clone.textContent.replace(/\s+$/, "") + "\n";
  }

  function copyText(text) {
    if (navigator.clipboard && window.isSecureContext) {
      return navigator.clipboard.writeText(text);
    }
    return new Promise(function (resolve, reject) {
      var area = document.createElement("textarea");
      area.value = text;
      area.setAttribute("readonly", "");
      area.style.position = "fixed";
      area.style.opacity = "0";
      document.body.appendChild(area);
      area.select();
      try {
        document.execCommand("copy") ? resolve() : reject();
      } catch (e) {
        reject(e);
      }
      area.remove();
    });
  }

  function addCopy(button, source) {
    button.addEventListener("click", function () {
      copyText(plainText(source)).then(
        function () {
          button.textContent = "Copied";
          button.classList.add("done");
          setTimeout(function () {
            button.textContent = "Copy";
            button.classList.remove("done");
          }, 1600);
        },
        function () {
          button.textContent = "Select and copy";
        }
      );
    });
  }

  document.querySelectorAll("[data-copy-from]").forEach(function (button) {
    var source = document.querySelector(button.getAttribute("data-copy-from"));
    if (source) addCopy(button, source);
  });

  // Every other code block gets a button too. A one-line block wraps rather
  // than run under the button; a block of lines keeps its layout.
  function decorateCode(root) {
    root.querySelectorAll(".code").forEach(function (block) {
      var pre = block.querySelector("pre");
      if (pre && pre.textContent.trim().indexOf("\n") === -1) block.classList.add("one-line");
      if (block.querySelector(".copy") || !pre) return;
      var button = document.createElement("button");
      button.type = "button";
      button.className = "copy";
      button.textContent = "Copy";
      block.appendChild(button);
      addCopy(button, pre);
    });
  }
  decorateCode(document);

  /* ── Command builder ── */

  // Step 3 of the front page, and the quickstart's copy of it: two commands
  // from one set of choices. The check is the run with --validate_only, so it
  // looks at the databases and the reference the run will use. Every option
  // the out.push lines below write is checked against main.nf and
  // nextflow.config by tools/check_site.py.
  if (document.querySelector("[data-builder]")) {
    // The boot script in the page's head sets this before the first paint;
    // a page without it still gets its choices shown.
    document.documentElement.classList.add("js");
  }

  // The site's root, from this script's own address: the notes link to the
  // docs from the front page and from a page a folder down.
  var builderRoot = (function () {
    var me = document.querySelector('script[src$="assets/site.js"]');
    return me ? me.getAttribute("src").replace("assets/site.js", "") : "";
  })();

  document.querySelectorAll("[data-builder]").forEach(function (builder) {
    var checkEl = builder.querySelector("[data-check]");
    var commandEl = builder.querySelector("[data-command]");
    var explainEl = builder.querySelector("[data-explain]");
    var noteEl = builder.querySelector("[data-builder-note]");
    var choice = { engine: "docker", organism: "mtbc" };
    // Whether "Only the corrected genes" was ticked by the laptop option, not
    // by the reader: unticking the laptop then unticks it too.
    var onlyByLaptop = false;
    var said = null;

    var PROFILES = {
      docker: { profile: "docker", says: "Every tool runs inside the AMAP image, through Docker.", db: "~/databases/eggnog" },
      singularity: { profile: "singularity", says: "The same image, through Singularity or Apptainer.", db: "~/databases/eggnog" },
      slurm: { profile: "singularity,slurm", says: "Each task is a SLURM job, run in the image.", db: "/shared/databases/eggnog" }
    };

    function box(name) {
      return builder.querySelector('[data-opt="' + name + '"]');
    }

    function opt(name) {
      var b = box(name);
      return b && b.checked && !b.disabled;
    }

    function esc(text) {
      return String(text).replace(/&/g, "&amp;").replace(/</g, "&lt;").replace(/>/g, "&gt;");
    }

    // Each line is [flags as they are typed, what they do].
    function lines() {
      var p = PROFILES[choice.engine];
      var out = [];
      out.push([["-profile", p.profile], p.says]);
      out.push([["--tsv", "samples.tsv"], "Your sample sheet: sampleId and fasta, tab separated."]);
      if (choice.organism === "other") {
        out.push([["--skip_typing", "--skip_mycobrowser"], "Lineage typing and Mycobrowser categories are for the tuberculosis complex."]);
        out.push([["--genus", "Escherichia", "--species", '"Escherichia coli"'], "Passed to Bakta, which uses them to name genes."]);
        out.push([["--ref_annotation", "ref.gbk", "--ref_genome", "ref.fna"], "Your trusted reference: its annotation and the genome it describes, the same assembly."]);
      }
      if (opt("small")) out.push([["--max_cpus", "4", "--max_memory", "8.GB"], "Ceilings: every task is clamped to them."]);
      if (opt("only")) out.push([["--skip_functional_annotation"], "Stops after the correction, with the corrected genes as a table and sequences: no eggNOG, GFF3 or GenBank, duplications, pangenome or final_report.html; amap_report.html is still written."]);
      if (opt("db")) out.push([["--eggnog_db", p.db], "Builds the eggNOG database there once, and reuses it."]);
      return out;
    }

    // A command as the page writes it: a span per line, with its newline in
    // it, and a span per word, a flag in .f and a value in .v, so that a
    // narrow screen wraps between words, under the line it continues.
    function commandHtml(all) {
      var rows = ['<span class="p">$ </span><span class="t">nextflow</span> <span class="t">run</span> <span class="t">PathoGenOmics-Lab/AMAP</span>'];
      all.forEach(function (line) {
        rows.push(
          "    " +
            line[0]
              .map(function (part) {
                return '<span class="' + (/^-/.test(part) ? "f" : "v") + '">' + esc(part) + "</span>";
              })
              .join(" ")
        );
      });
      return rows
        .map(function (row, k) {
          return '<span class="ln">' + row + (k < rows.length - 1 ? " \\\n" : "") + "</span>";
        })
        .join("");
    }

    function explainHtml(all) {
      return all
        .map(function (line) {
          var flags = line[0].filter(function (part) {
            return /^-/.test(part);
          });
          return (
            '<li><span class="flags">' +
            flags
              .map(function (flag) {
                return "<code>" + esc(flag) + "</code>";
              })
              .join(" ") +
            "</span><span>" + esc(line[1]) + "</span></li>"
          );
        })
        .join("");
    }

    function notes() {
      var out = [];
      var small = opt("small");
      if (choice.engine === "slurm") {
        out.push("On a cluster the slurm profile sets the ceilings, so a laptop's do not apply.");
        out.push(
          "Bakta's jobs ask for <code>--qos=medium</code>; if your cluster has none, see " +
            '<a href="' + builderRoot + 'deployment/hpc.html#fitting-your-site">HPC clusters</a>.'
        );
      }
      if (choice.engine !== "docker") {
        out.push(
          "The singularity profile keeps the converted image in <code>singularity_cache/</code> under <code>--outdir</code>, " +
            "and that wins over <code>NXF_SINGULARITY_CACHEDIR</code>. To keep one shared copy, put " +
            "<code>singularity.cacheDir = '/shared/singularity_cache'</code> in a config passed with <code>-c</code>."
        );
      }
      if (small && choice.engine === "docker") {
        out.push("8 GB is what Docker Desktop usually gives: if <code>docker info</code> reports more, raise <code>--max_memory</code> to it.");
      } else if (small) {
        out.push("Set <code>--max_memory</code> to the memory the machine has free.");
      }
      if (small && opt("only") && onlyByLaptop) {
        out.push("<b>Only the corrected genes</b> is ticked for you: eggNOG loads its whole database into memory and asks for 50 GB. Untick it on a machine that has that much.");
      } else if (small && !opt("only")) {
        out.push("eggNOG loads its whole database into memory and asks for 50 GB, which a laptop does not have: tick <b>Only the corrected genes</b>, or run it on a larger machine.");
      }
      if (choice.organism === "other") {
        out.push(
          "The annotation and the genome must be the same assembly. How many proteins AMAP took from it is in the run's work folder: " +
            "<code>grep -h 'Curated proteins' work/*/*/.command.out</code>."
        );
      }
      return out;
    }

    function render() {
      var small = box("small");
      if (small) small.disabled = choice.engine === "slurm";
      var db = box("db");
      if (db) db.disabled = opt("only");
      var all = lines();
      commandEl.innerHTML = commandHtml(all);
      if (checkEl) checkEl.innerHTML = commandHtml(all.concat([[["--validate_only"], ""]]));
      explainEl.innerHTML = explainHtml(all);
      if (noteEl) {
        var list = notes();
        noteEl.innerHTML = list
          .map(function (note) {
            return "<p>" + note + "</p>";
          })
          .join("");
        noteEl.hidden = !list.length;
      }
      // What changed, for a screen reader (announce() is the page's status
      // line): each line by its first flag, so a new profile reads as one.
      var now = {};
      all.forEach(function (line) {
        now[line[0][0]] = line[0].join(" ");
      });
      if (said) {
        var parts = [];
        Object.keys(now).forEach(function (flag) {
          if (!(flag in said)) parts.push("adds " + now[flag]);
          else if (said[flag] !== now[flag]) parts.push(now[flag]);
        });
        Object.keys(said).forEach(function (flag) {
          if (!(flag in now)) parts.push("drops " + said[flag]);
        });
        if (parts.length && typeof announce === "function") announce("Command updated: " + parts.join("; "));
      }
      said = now;
    }

    // Each question is a radio group: one tab stop, and the arrow keys move
    // the choice along it.
    function pick(button, focus) {
      var group = button.closest("[data-choice]");
      choice[group.getAttribute("data-choice")] = button.getAttribute("data-value");
      group.querySelectorAll('[role="radio"]').forEach(function (b) {
        var on = b === button;
        b.setAttribute("aria-checked", on ? "true" : "false");
        b.tabIndex = on ? 0 : -1;
      });
      if (focus) button.focus();
      var only = box("only");
      if (choice.engine === "slurm" && onlyByLaptop && only) {
        only.checked = false;
        onlyByLaptop = false;
      }
      render();
    }

    builder.addEventListener("click", function (event) {
      var button = event.target.closest('[data-choice] [role="radio"]');
      if (button) pick(button, false);
    });

    builder.addEventListener("keydown", function (event) {
      var button = event.target.closest('[data-choice] [role="radio"]');
      if (!button) return;
      var radios = Array.prototype.slice.call(button.closest("[data-choice]").querySelectorAll('[role="radio"]'));
      var at = radios.indexOf(button);
      var to = { ArrowRight: at + 1, ArrowDown: at + 1, ArrowLeft: at - 1, ArrowUp: at - 1, Home: 0, End: radios.length - 1 }[event.key];
      if (to === undefined) return;
      event.preventDefault();
      pick(radios[(to + radios.length) % radios.length], true);
    });

    builder.addEventListener("change", function (event) {
      var only = box("only");
      if (event.target === box("small") && only) {
        // A laptop cannot hold eggNOG's database: the run stops before it,
        // unless the reader asks otherwise.
        if (event.target.checked && !only.checked) {
          only.checked = true;
          onlyByLaptop = true;
        } else if (!event.target.checked && onlyByLaptop) {
          only.checked = false;
          onlyByLaptop = false;
        }
      } else if (event.target === only) {
        onlyByLaptop = false;
      }
      render();
    });
    render();
  });

  /* ── Embedded report ── */

  // The example report is 1.9 MB of scripts and plots, so the page shows a
  // picture of it until the reader asks for it. Once open it takes the mouse
  // and the keyboard for its own plots, and gives them back when the pointer
  // or focus leaves it, a touch lands outside it, Escape is pressed or Done
  // is tapped. Until then it is inert: none of its links is a tab stop.
  document.querySelectorAll("[data-report-stage]").forEach(function (stage) {
    var cover = stage.querySelector(".report-cover");
    var frame = stage.querySelector("iframe");
    var done = stage.querySelector(".report-done");
    if (!cover || !frame) return;
    frame.inert = true;

    function rest(refocus) {
      if (!stage.classList.contains("live")) return;
      stage.classList.remove("live");
      frame.inert = true;
      if (done) done.hidden = true;
      if (refocus) cover.focus();
    }

    function onEscape(event) {
      if (event.key === "Escape") rest(true);
    }

    function open() {
      if (!frame.getAttribute("src") && frame.getAttribute("data-src")) {
        frame.addEventListener(
          "load",
          function () {
            stage.classList.add("loaded");
            var label = cover.querySelector("span");
            if (label) label.textContent = "Explore the example report";
            // Escape inside the report, where the page can reach it (not
            // from a file:// address, which Chrome keeps apart).
            try {
              frame.contentWindow.document.addEventListener("keydown", onEscape);
            } catch (e) {
              /* another origin: Escape works from the page only */
            }
          },
          { once: true }
        );
        frame.setAttribute("src", frame.getAttribute("data-src"));
      }
      stage.classList.add("live");
      frame.inert = false;
      if (done) done.hidden = false;
      try {
        frame.focus();
      } catch (e) {
        /* nothing to focus yet */
      }
    }

    cover.addEventListener("click", open);
    if (done) {
      done.addEventListener("click", function () {
        rest(true);
      });
    }
    stage.addEventListener("mouseleave", function () {
      rest(false);
    });
    stage.addEventListener("focusout", function (event) {
      if (event.relatedTarget && !stage.contains(event.relatedTarget)) rest(false);
    });
    document.addEventListener(
      "touchstart",
      function (event) {
        if (!stage.contains(event.target)) rest(false);
      },
      { passive: true }
    );
    document.addEventListener("keydown", function (event) {
      if (stage.classList.contains("live")) onEscape(event);
    });
  });

  /* ── Documentation: search, page contents, menu ── */

  var searchBox = document.querySelector("[data-docs-search]");
  var resultsBox = document.querySelector("[data-docs-results]");
  if (searchBox && resultsBox && window.AMAP_DOCS_INDEX) {
    // The index lists URLs from the docs root; this page may sit a folder down.
    var indexScript = document.querySelector('script[src$="assets/search-index.js"]');
    var root = indexScript ? indexScript.getAttribute("src").replace("assets/search-index.js", "") : "";
    var active = -1;

    var escapeHtml = function (text) {
      return text.replace(/&/g, "&amp;").replace(/</g, "&lt;").replace(/>/g, "&gt;").replace(/"/g, "&quot;");
    };

    // Words match where a word starts, so "skip" finds --skip_typing and not
    // the middle of another word.
    var startsWord = function (w, flags) {
      return new RegExp("(^|[^a-z0-9])(" + w.replace(/[.*+?^${}()|[\]\\]/g, "\\$&") + ")", flags);
    };

    var highlight = function (text, words) {
      var safe = escapeHtml(text);
      words.forEach(function (w) {
        safe = safe.replace(startsWord(escapeHtml(w), "ig"), "$1<mark>$2</mark>");
      });
      return safe;
    };

    var snippet = function (text, word) {
      var at = text.toLowerCase().indexOf(word);
      if (at < 0) return text.slice(0, 140);
      var start = Math.max(0, at - 60);
      return (start > 0 ? "…" : "") + text.slice(start, start + 150) + (start + 150 < text.length ? "…" : "");
    };

    var search = function (query) {
      var words = query.toLowerCase().split(/\s+/).filter(function (w) {
        return w.length > 1;
      });
      if (!words.length) return [words, []];
      var hits = [];
      window.AMAP_DOCS_INDEX.forEach(function (entry) {
        var page = entry.page.toLowerCase();
        var section = entry.section.toLowerCase();
        var text = entry.text.toLowerCase();
        var score = 0;
        for (var k = 0; k < words.length; k++) {
          var pattern = startsWord(words[k], "i");
          var inSection = pattern.test(section);
          var inPage = pattern.test(page);
          var inText = pattern.test(text);
          if (!inSection && !inPage && !inText) return;
          score += (inSection ? 5 : 0) + (inPage ? 3 : 0) + (inText ? 1 : 0);
        }
        hits.push({ entry: entry, score: score });
      });
      hits.sort(function (a, b) {
        return b.score - a.score;
      });
      return [words, hits.slice(0, 8)];
    };

    var renderResults = function () {
      var found = search(searchBox.value.trim());
      var words = found[0];
      var hits = found[1];
      active = -1;
      if (!words.length) {
        resultsBox.hidden = true;
        return;
      }
      resultsBox.hidden = false;
      if (!hits.length) {
        resultsBox.innerHTML = '<div class="none">Nothing matches. Try another word.</div>';
        return;
      }
      resultsBox.innerHTML = hits
        .map(function (hit) {
          var e = hit.entry;
          var title = e.section ? e.page + " › " + e.section : e.page;
          return (
            '<a href="' + escapeHtml(root + e.url) + '"><b>' + highlight(title, words) + "</b><small>" +
            highlight(snippet(e.text, words[0]), words) + "</small></a>"
          );
        })
        .join("");
    };

    var move = function (step) {
      var links = resultsBox.querySelectorAll("a");
      if (!links.length) return;
      active = (active + step + links.length) % links.length;
      links.forEach(function (link, k) {
        link.classList.toggle("active", k === active);
      });
      links[active].scrollIntoView({ block: "nearest" });
    };

    searchBox.addEventListener("input", renderResults);
    searchBox.addEventListener("keydown", function (event) {
      if (event.key === "ArrowDown") {
        event.preventDefault();
        move(1);
      } else if (event.key === "ArrowUp") {
        event.preventDefault();
        move(-1);
      } else if (event.key === "Enter") {
        var links = resultsBox.querySelectorAll("a");
        var pick = links[active >= 0 ? active : 0];
        if (pick) window.location.href = pick.getAttribute("href");
      } else if (event.key === "Escape") {
        searchBox.value = "";
        resultsBox.hidden = true;
      }
    });
    document.addEventListener("click", function (event) {
      if (!event.target.closest(".docs-search")) resultsBox.hidden = true;
    });
    // "/" jumps to the search box, as on most documentation sites.
    document.addEventListener("keydown", function (event) {
      if (event.key === "/" && !/input|textarea|select/i.test(document.activeElement.tagName)) {
        event.preventDefault();
        searchBox.focus();
      }
    });
  }

  // The page contents follow the reader down the page.
  var tocLinks = Array.prototype.slice.call(document.querySelectorAll(".docs-toc a"));
  if (tocLinks.length) {
    var targets = tocLinks
      .map(function (link) {
        return document.getElementById(link.getAttribute("href").slice(1));
      })
      .filter(Boolean);
    var ticking = false;
    var markCurrent = function () {
      ticking = false;
      var current = 0;
      targets.forEach(function (heading, k) {
        if (heading.getBoundingClientRect().top < 120) current = k;
      });
      tocLinks.forEach(function (link, k) {
        link.classList.toggle("active", k === current);
      });
    };
    window.addEventListener(
      "scroll",
      function () {
        if (!ticking) {
          ticking = true;
          requestAnimationFrame(markCurrent);
        }
      },
      { passive: true }
    );
    markCurrent();
  }

  // On a phone the page list starts folded, so the page itself comes first.
  var docsMenu = document.querySelector(".docs-menu");
  if (docsMenu && window.matchMedia && window.matchMedia("(max-width: 860px)").matches) {
    docsMenu.open = false;
  }
})();

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

  document.querySelectorAll("[data-builder]").forEach(function (builder) {
    var commandEl = builder.querySelector("[data-command]");
    var explainEl = builder.querySelector("[data-explain]");
    var noteEl = builder.querySelector("[data-builder-note]");
    var choice = { engine: "docker", organism: "mtbc" };

    var PROFILES = {
      docker: { profile: "docker", says: "Every tool runs inside the AMAP image, through Docker.", db: "~/databases/eggnog" },
      singularity: { profile: "singularity", says: "The same image, through Singularity or Apptainer.", db: "~/databases/eggnog" },
      slurm: { profile: "singularity,slurm", says: "Each task is a SLURM job, run in the image.", db: "/shared/databases/eggnog" }
    };

    function opt(name) {
      var box = builder.querySelector('[data-opt="' + name + '"]');
      return box && box.checked && !box.disabled;
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
      if (opt("small")) out.push([["--max_cpus", "4", "--max_memory", "12.GB"], "Ceilings: every task is clamped to them."]);
      if (opt("only")) out.push([["--skip_functional_annotation"], "Stops after the correction, with the corrected genes as a table and sequences: no eggNOG, GFF3 or GenBank, duplications, pangenome or report."]);
      if (opt("db")) out.push([["--eggnog_db", p.db], "Builds the eggNOG database there once, and reuses it."]);
      return out;
    }

    function render() {
      var all = lines();
      var text = '<span class="p">$ </span>nextflow run PathoGenOmics-Lab/AMAP';
      all.forEach(function (line) {
        var flags = line[0];
        text += " \\\n    ";
        text += flags
          .map(function (part) {
            return /^-/.test(part) ? '<span class="f">' + esc(part) + "</span>" : '<span class="v">' + esc(part) + "</span>";
          })
          .join(" ");
      });
      commandEl.innerHTML = text;
      explainEl.innerHTML = all
        .map(function (line) {
          return "<li><code>" + esc(line[0].filter(function (part) {
            return /^-/.test(part);
          }).join(" ")) + "</code><span>" + esc(line[1]) + "</span></li>";
        })
        .join("");
      var db = builder.querySelector('[data-opt="db"]');
      if (db) db.disabled = opt("only");
      if (noteEl) {
        var notes = [];
        if (choice.engine === "docker" && opt("small")) {
          notes.push("On macOS and Windows, Docker Desktop has its own memory limit, and it is the one that counts: keep <code>--max_memory</code> at or below it (Settings, Resources, Memory).");
        }
        if (choice.organism === "other") {
          notes.push("The annotation and the genome must be the same assembly. AMAP prints how many proteins it extracted: check it looks like a genome's worth.");
        }
        noteEl.innerHTML = notes.join(" ");
        noteEl.hidden = !notes.length;
      }
    }

    builder.addEventListener("click", function (event) {
      var button = event.target.closest("[data-choice] button");
      if (!button) return;
      var group = button.closest("[data-choice]");
      choice[group.getAttribute("data-choice")] = button.getAttribute("data-value");
      group.querySelectorAll("button").forEach(function (b) {
        b.setAttribute("aria-pressed", b === button ? "true" : "false");
      });
      render();
    });
    builder.addEventListener("change", render);
    render();
  });

  /* ── Embedded report ── */

  // The report takes the mouse wheel for its own plots, so it only takes the
  // mouse after a click, and gives it back when the pointer leaves.
  document.querySelectorAll("[data-report-stage]").forEach(function (stage) {
    var cover = stage.querySelector(".report-cover");
    var frame = stage.querySelector("iframe");
    if (!cover || !frame) return;
    cover.addEventListener("click", function () {
      stage.classList.add("live");
      try {
        frame.focus();
      } catch (e) {
        /* nothing to focus yet */
      }
    });
    stage.addEventListener("mouseleave", function () {
      stage.classList.remove("live");
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

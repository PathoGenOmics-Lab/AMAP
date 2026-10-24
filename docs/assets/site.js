/* AMAP website: theme, tabs, copy buttons, the command builder, the embedded
   report, the search on every page, and the documentation pages' contents.
   Plain JavaScript, no build step. The drawings are in figures.js. */
(function () {
  "use strict";

  // The site's root, read from this script's own address: the search index and
  // the pages it lists sit there, whatever folder this page is in (the 404 page
  // is served from any path).
  var SITE_SCRIPT = document.currentScript || document.querySelector('script[src*="assets/site.js"]');
  var ROOT = SITE_SCRIPT ? SITE_SCRIPT.src.replace(/assets\/site\.js(?:[?#].*)?$/, "") : "";

  function reducedMotion() {
    return !!(window.matchMedia && window.matchMedia("(prefers-reduced-motion: reduce)").matches);
  }

  // Storage can be missing or refused (private windows, blocked cookies); the
  // page works without it.
  function sessionGet(key) {
    try {
      return JSON.parse(sessionStorage.getItem(key));
    } catch (e) {
      return null;
    }
  }

  function sessionSet(key, value) {
    try {
      if (value === null) sessionStorage.removeItem(key);
      else sessionStorage.setItem(key, JSON.stringify(value));
    } catch (e) {
      /* nothing to remember with */
    }
  }

  /* ── Theme ── */

  // Dark unless the reader chose light. The button is a toggle for the light
  // theme: pressed when light is on, and its icon shows what a press gives you.
  // The browser's own bar (theme-color) follows the page.
  var THEME_KEY = "amap.theme";
  var THEME_COLOR = { dark: "#0a0e0d", light: "#f5f7f6" };
  var MOON =
    '<svg class="i-moon" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round" stroke-linejoin="round" aria-hidden="true">' +
    '<path d="M20.4 14.6A8.5 8.5 0 0 1 9.4 3.6a8.5 8.5 0 1 0 11 11z"/></svg>';

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

  var themeToggles = document.querySelectorAll("[data-theme-toggle]");

  function syncTheme() {
    var light = currentTheme() === "light";
    themeToggles.forEach(function (button) {
      button.setAttribute("aria-pressed", light ? "true" : "false");
      button.title = light ? "Switch to the dark theme" : "Switch to the light theme";
    });
    var meta = document.querySelector('meta[name="theme-color"]');
    if (meta) meta.setAttribute("content", THEME_COLOR[light ? "light" : "dark"]);
  }

  themeToggles.forEach(function (button) {
    var sun = button.querySelector("svg");
    if (sun) sun.classList.add("i-sun");
    button.insertAdjacentHTML("beforeend", MOON);
    button.setAttribute("aria-label", "Light theme");
    button.addEventListener("click", function () {
      var next = currentTheme() === "dark" ? "light" : "dark";
      var apply = function () {
        document.documentElement.dataset.theme = next;
        syncTheme();
      };
      // A short cross-fade where the browser can draw one, none under reduced motion.
      if (document.startViewTransition && !reducedMotion()) document.startViewTransition(apply);
      else apply();
      remember(next);
    });
  });
  syncTheme();
  // Anything else that sets the theme (another tab, a test) keeps the button
  // and the browser bar in step.
  if (window.MutationObserver) {
    new MutationObserver(syncTheme).observe(document.documentElement, { attributes: true, attributeFilter: ["data-theme"] });
  }

  /* ── Status ── */

  // One line that screen readers read out, for changes that happen away from
  // focus: a copy, the command builder's new command, the number of results.
  var statusLine = document.createElement("div");
  statusLine.className = "visually-hidden";
  statusLine.setAttribute("role", "status");
  document.body.insertBefore(statusLine, document.body.firstChild);
  var statusTimer = 0;

  function announce(text, delay) {
    clearTimeout(statusTimer);
    statusLine.textContent = "";
    // Emptied first, so the same words twice are still read twice.
    statusTimer = setTimeout(function () {
      statusLine.textContent = text;
    }, delay || 60);
  }

  /* ── Tabs ── */

  // A file excerpt too wide for a phone scrolls; it opens on what it is
  // there to show, its first highlighted mark, with a little before it.
  function revealMark(panel) {
    var pre = panel && panel.querySelector("pre");
    var mark = pre && pre.querySelector("mark");
    if (!mark || pre.scrollWidth <= pre.clientWidth + 1) return;
    var left = mark.getBoundingClientRect().left - pre.getBoundingClientRect().left + pre.scrollLeft;
    pre.scrollLeft = Math.max(0, left - pre.clientWidth * 0.3);
  }

  document.querySelectorAll('[role="tablist"]').forEach(function (list) {
    var tabs = Array.prototype.slice.call(list.querySelectorAll('[role="tab"]'));
    function select(tab) {
      tabs.forEach(function (t) {
        var on = t === tab;
        t.setAttribute("aria-selected", on ? "true" : "false");
        t.tabIndex = on ? 0 : -1;
        var panel = document.getElementById(t.getAttribute("aria-controls"));
        if (panel) panel.hidden = !on;
        if (on) revealMark(panel);
      });
    }
    tabs.forEach(function (tab) {
      if (tab.getAttribute("aria-selected") === "true") revealMark(document.getElementById(tab.getAttribute("aria-controls")));
    });
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

  /* ── Sideways scrollers ── */

  // On a phone the case cards, the story chips and the metro map scroll
  // sideways. What the keyboard focuses there comes into view: a case's
  // whole card, anything else just enough. Not on a screen wide enough for
  // them, where nothing scrolls.
  ["cases", "how-group", "metro-scroll"].forEach(function (name) {
    document.querySelectorAll("." + name).forEach(function (box) {
      box.addEventListener("focusin", function (event) {
        var scroller = name === "how-group" ? event.target.closest(".seg") : box;
        if (!scroller || scroller.scrollWidth <= scroller.clientWidth + 1) return;
        var card = name === "cases" && event.target.closest(".case");
        (card || event.target).scrollIntoView({ block: "nearest", inline: card ? "start" : "nearest" });
      });
    });
  });

  /* ── Scroll regions ── */

  // A table, a command or a program's output too wide for the screen
  // scrolls sideways. While it does, it is a tab stop with a name, so the
  // keyboard can scroll it; while it fits, it is neither.
  function scrollRegion(el) {
    var scrolls = el.scrollWidth > el.clientWidth + 1;
    var was = el.hasAttribute("data-scroll-region");
    if (scrolls === was) return;
    if (scrolls) {
      el.setAttribute("data-scroll-region", "");
      el.setAttribute("tabindex", "0");
      el.setAttribute("role", "region");
      if (!el.hasAttribute("aria-labelledby") && !el.hasAttribute("aria-label")) {
        var code = el.closest(".code");
        var title = code && code.querySelector(".code-title");
        var kind = el.classList.contains("table-wrap")
          ? "Table"
          : title
            ? title.textContent.trim()
            : el.closest(".cmd, .builder") || el.querySelector(".ln")
              ? "Command"
              : el.closest(".excerpt")
                ? "File excerpt"
                : el.classList.contains("out") || el.closest(".out")
                  ? "Output"
                  : "Text";
        el.setAttribute("aria-label", kind + ", scrolls sideways");
        el.setAttribute("data-scroll-label", "");
      }
    } else {
      el.removeAttribute("data-scroll-region");
      el.removeAttribute("tabindex");
      el.removeAttribute("role");
      if (el.hasAttribute("data-scroll-label")) {
        el.removeAttribute("aria-label");
        el.removeAttribute("data-scroll-label");
      }
    }
  }
  var scrollers = document.querySelectorAll(".table-wrap, .code pre, .out, .excerpt pre");
  if ("ResizeObserver" in window) {
    var watchScroll = new ResizeObserver(function (records) {
      records.forEach(function (r) {
        scrollRegion(r.target);
      });
    });
    scrollers.forEach(function (el) {
      watchScroll.observe(el);
    });
  } else {
    scrollers.forEach(scrollRegion);
  }

  /* ── Copy buttons ── */

  function plainText(node) {
    var clone = node.cloneNode(true);
    clone.querySelectorAll(".p").forEach(function (prompt) {
      prompt.remove();
    });
    return clone.textContent.replace(/\s+$/, "") + "\n";
  }

  // The clipboard API when the page may use it, the old way when it may not
  // (a page opened from a clone, a browser that refuses without a gesture).
  function copyText(text) {
    var fallback = function () {
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
    };
    if (navigator.clipboard && window.isSecureContext) {
      return navigator.clipboard.writeText(text).catch(fallback);
    }
    return fallback();
  }

  function addCopy(button, source) {
    var label = button.textContent;
    button.addEventListener("click", function () {
      copyText(plainText(source)).then(
        function () {
          button.textContent = "Copied";
          button.classList.add("done");
          announce("Copied to the clipboard");
          setTimeout(function () {
            button.textContent = label;
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

  // A short note beside what was copied (after a heading's last word, or
  // above a flag in a table). It sits over the page, in the page's own
  // coordinates, so nothing moves and it scrolls with the text.
  function noteBeside(target, text, above) {
    var range = document.createRange();
    range.selectNodeContents(target);
    var rects = range.getClientRects();
    var box = rects.length && !above ? rects[rects.length - 1] : target.getBoundingClientRect();
    var note = document.createElement("span");
    note.className = "copied-note" + (above ? " above" : "");
    note.setAttribute("aria-hidden", "true");
    note.textContent = text;
    note.style.left = Math.round((above ? box.left : box.right + 10) + window.scrollX) + "px";
    note.style.top = Math.round(box.top + window.scrollY + (above ? -6 : box.height / 2)) + "px";
    document.body.appendChild(note);
    // A note that would run past the window's edge (a long heading on a
    // phone) sits above the line instead, inside the edge, so the page never
    // widens under it.
    var room = document.documentElement.clientWidth - 8;
    if (note.getBoundingClientRect().right > room) {
      note.classList.add("above");
      note.style.left = Math.round(Math.max(8, room - note.offsetWidth) + window.scrollX) + "px";
      note.style.top = Math.round(box.top + window.scrollY - 6) + "px";
    }
    setTimeout(function () {
      note.remove();
    }, 1500);
  }

  // A heading's "#" copies the link to that heading, and still moves the
  // address bar there. An option's flag in the Parameters tables copies the
  // flag.
  document.querySelectorAll('tr[id^="opt-"] > :first-child code').forEach(function (code) {
    code.classList.add("copyable");
    code.title = "Click to copy";
  });

  document.addEventListener("click", function (event) {
    if (event.button || event.metaKey || event.ctrlKey || event.shiftKey || event.altKey) return;
    var target = event.target.closest ? event.target : event.target.parentElement;
    if (!target) return;
    var anchor = target.closest("a.anchor");
    if (anchor) {
      var hash = anchor.getAttribute("href") || "";
      if (hash.charAt(0) !== "#") return;
      var heading = anchor.closest("h1, h2, h3, h4, h5, h6") || document.getElementById(decodeURIComponent(hash.slice(1)));
      copyText(location.href.split("#")[0] + hash).then(function () {
        noteBeside(heading || anchor, "Link copied");
        announce("Link copied");
      }, function () {});
      return;
    }
    var flag = target.closest(".copyable");
    if (flag && !String(window.getSelection ? window.getSelection() : "")) {
      var text = flag.textContent.trim();
      copyText(text).then(function () {
        noteBeside(flag, "Copied", true);
        announce("Copied " + text);
      }, function () {});
    }
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
    // A continued line keeps its indent with its first word and its closing
    // backslash with its last, so neither sits on a line of its own on a
    // phone; what Copy takes is the same text.
    function commandHtml(all) {
      var rows = [['<span class="p">$ </span><span class="t">nextflow</span>', '<span class="t">run</span>', '<span class="t">PathoGenOmics-Lab/AMAP</span>']];
      all.forEach(function (line) {
        rows.push(
          line[0].map(function (part) {
            return '<span class="' + (/^-/.test(part) ? "f" : "v") + '">' + esc(part) + "</span>";
          })
        );
      });
      return rows
        .map(function (row, k) {
          var parts = row.slice();
          var more = k < rows.length - 1;
          if (k > 0) parts[0] = '<span class="nw">    ' + parts[0] + "</span>";
          if (more) parts[parts.length - 1] = '<span class="nw">' + parts[parts.length - 1] + " \\</span>";
          // An indented line that wraps hangs four characters further in,
          // as the docs' commands do (--i is its indent).
          return '<span class="ln"' + (k > 0 ? ' style="--i:4"' : "") + ">" + parts.join(" ") + (more ? "\n" : "") + "</span>";
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
                return '<code class="nw">' + esc(flag) + "</code>";
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
          "Bakta's jobs ask for <code class='nw'>--qos=medium</code>; if your cluster has none, see " +
            '<a href="' + builderRoot + 'deployment/hpc.html#fitting-your-site">HPC clusters</a>.'
        );
      }
      if (choice.engine !== "docker") {
        out.push(
          "The singularity profile keeps the converted image in <code class='nw'>singularity_cache/</code> under <code class='nw'>--outdir</code>, " +
            "and that wins over <code class='nw'>NXF_SINGULARITY_CACHEDIR</code>. To keep one shared copy, put " +
            "<code>singularity.cacheDir = '/shared/singularity_cache'</code> in a config passed with <code class='nw'>-c</code>."
        );
      }
      if (small && choice.engine === "docker") {
        out.push("8 GB is what Docker Desktop usually gives: if <code class='nw'>docker info</code> reports more, raise <code class='nw'>--max_memory</code> to it.");
      } else if (small) {
        out.push("Set <code class='nw'>--max_memory</code> to the memory the machine has free.");
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

  // The example report is 4.5 MB of scripts and plots, so the page shows a
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

  // An old anchor (a heading renamed since, or GitHub's anchor for a
  // changelog heading) leads to the heading's current one, which goes in the
  // address so the contents and the arrival mark follow it.
  function followAlias() {
    var id = decodeURIComponent(location.hash.slice(1));
    var alias = id && document.getElementById(id);
    if (!alias || !alias.classList.contains("anchor-alias")) return;
    var to = alias.getAttribute("data-to");
    var target = document.getElementById(to);
    if (!target) return;
    try {
      history.replaceState(history.state, "", "#" + to);
    } catch (e) {}
    target.scrollIntoView();
  }
  followAlias();
  window.addEventListener("hashchange", followAlias);

  var header = document.querySelector("header.top");
  var HEADER_GAP = 61; // the sticky header's height, its border included

  function escapeHtml(text) {
    return String(text).replace(/&/g, "&amp;").replace(/</g, "&lt;").replace(/>/g, "&gt;").replace(/"/g, "&quot;");
  }

  function escapeRe(text) {
    return text.replace(/[.*+?^${}()|[\]\\]/g, "\\$&");
  }

  // ── The index ──
  // One file, assets/search-index.js, serves every page. It loads the first
  // time a reader reaches for search (focus, pointer, shortcut), not with the
  // page; a page that already has it uses it as it is.

  var indexState = window.AMAP_DOCS_INDEX ? "ready" : "idle";
  var indexWaiting = [];

  function indexEntries() {
    var entries = window.AMAP_DOCS_INDEX || [];
    if (!entries.prepared) {
      entries.forEach(function (e) {
        e._page = (e.page || "").toLowerCase();
        e._section = (e.section || "").toLowerCase();
        e._text = (e.text || "").toLowerCase();
        e._code = (e.code || "").toLowerCase();
        e._base = (e.url || "").split("#")[0];
        e._kind = e.kind || (/^about\/changelog/.test(e.url) ? "changelog" : e.section ? "section" : "page");
      });
      entries.prepared = true;
    }
    return entries;
  }

  function loadIndex(then) {
    if (window.AMAP_DOCS_INDEX) indexState = "ready";
    if (indexState === "ready" || indexState === "failed") {
      if (then) then();
      return;
    }
    if (then && indexWaiting.indexOf(then) < 0) indexWaiting.push(then);
    if (indexState === "loading") return;
    indexState = "loading";
    var script = document.createElement("script");
    script.src = ROOT + "assets/search-index.js";
    script.async = true;
    var settle = function () {
      indexState = window.AMAP_DOCS_INDEX ? "ready" : "failed";
      var waiting = indexWaiting;
      indexWaiting = [];
      waiting.forEach(function (fn) {
        fn();
      });
    };
    script.onload = settle;
    script.onerror = settle;
    document.head.appendChild(script);
  }

  // ── Matching and ranking ──

  // Words people type for what the docs call something else.
  var SYNONYMS = {
    oom: ["memory"],
    db: ["database"],
    cluster: ["hpc", "slurm"],
    mamba: ["conda"],
    gbk: ["genbank"],
    vanish: ["missing", "ignored"],
    vanished: ["missing", "ignored"],
    disappear: ["missing", "ignored"],
    disappeared: ["missing", "ignored"],
    dropped: ["missing", "ignored"],
    proteome: ["faa", "proteins"],
    laptop: ["docker desktop"]
  };

  // Words a question in plain words is made of and every page has: "how to
  // run on slurm" searches for run and slurm. A query of nothing else keeps
  // them.
  var STOP = {};
  ("a an the to of on in is are be do does did how what where why which when who can i my me we you use using " +
    "much many mean means get with for it this that there").split(" ").forEach(function (w) {
    STOP[w] = true;
  });

  // A query word matches where a word starts, so "skip" finds --skip_typing and
  // not the middle of another word. A lone number matches only as a whole
  // number ("case 4", not 40). A long word also matches inside other words
  // ("mastertable" in AMAPmastertable.tsv).
  function term(word, weight) {
    var core = escapeRe(word);
    var tail = /^[0-9]+$/.test(word) ? "(?![0-9]|[.,][0-9])" : "";
    return {
      word: word,
      weight: weight,
      start: new RegExp("(^|[^a-z0-9])(" + core + ")" + tail, "g"),
      inside: word.length >= 5 ? new RegExp("()(" + core + ")", "g") : null
    };
  }

  function queryWords(query) {
    var words = query
      .toLowerCase()
      .split(/\s+/)
      .map(function (w) {
        return w.replace(/^[^a-z0-9-]+|[^a-z0-9]+$/g, "");
      })
      .filter(function (w) {
        return w.length > 1 || /^[0-9]$/.test(w);
      });
    var content = words.filter(function (w) {
      return !STOP[w];
    });
    return content.length ? content : words;
  }

  function termGroups(words) {
    return words.map(function (w) {
      return [term(w, 1)].concat(
        (SYNONYMS[w] || []).map(function (s) {
          return term(s, 0.8);
        })
      );
    });
  }

  function count(re, text) {
    if (!text) return 0;
    re.lastIndex = 0;
    var n = 0;
    while (re.exec(text) && n < 60) n++;
    return n;
  }

  function hits(t, text) {
    var n = count(t.start, text);
    return n || !t.inside ? n : count(t.inside, text);
  }

  // One query word against one entry: the best of the word and its synonyms.
  function scoreTerm(group, e) {
    var best = null;
    group.forEach(function (t) {
      var page = hits(t, e._page);
      var section = hits(t, e._section);
      var text = hits(t, e._text);
      var code = hits(t, e._code);
      if (!page && !section && !text && !code) return;
      var s =
        1 +
        (page ? 5 : 0) +
        (section ? 3 : 0) +
        (e._section === t.word || e._page === t.word ? 3 : 0) +
        Math.log(1 + text) +
        0.5 * Math.log(1 + code);
      s *= t.weight;
      if (!best || s > best.score) best = { score: s, tf: page + section + text + code };
    });
    return best;
  }

  // Words typed next to each other score more where they stand together.
  function phraseBonus(words, e) {
    var bonus = 0;
    for (var k = 0; k + 1 < words.length; k++) {
      var tail = /^[0-9]+$/.test(words[k + 1]) ? "(?![0-9]|[.,][0-9])" : "";
      var re = new RegExp("(^|[^a-z0-9])" + escapeRe(words[k]) + "[^a-z0-9]{1,3}" + escapeRe(words[k + 1]) + tail, "g");
      if (count(re, e._section)) bonus += 4;
      var inText = count(re, e._text);
      if (inText) bonus += 2 + Math.log(1 + inText);
    }
    return bonus;
  }

  var vocabulary = null;

  function buildVocabulary() {
    if (vocabulary) return vocabulary;
    vocabulary = {};
    indexEntries().forEach(function (e) {
      (e._page + " " + e._section + " " + e._text).replace(/[a-z][a-z0-9]{2,}/g, function (w) {
        vocabulary[w] = (vocabulary[w] || 0) + 1;
        return w;
      });
    });
    return vocabulary;
  }

  // Edits between two words, a swap of neighbours counting as one.
  function distance(a, b, limit) {
    if (Math.abs(a.length - b.length) > limit) return limit + 1;
    var prev2 = null;
    var prev = [];
    for (var j = 0; j <= b.length; j++) prev.push(j);
    for (var i = 1; i <= a.length; i++) {
      var row = [i];
      var low = i;
      for (var k = 1; k <= b.length; k++) {
        var cost = a.charAt(i - 1) === b.charAt(k - 1) ? 0 : 1;
        var v = Math.min(prev[k] + 1, row[k - 1] + 1, prev[k - 1] + cost);
        if (prev2 && i > 1 && k > 1 && a.charAt(i - 1) === b.charAt(k - 2) && a.charAt(i - 2) === b.charAt(k - 1)) {
          v = Math.min(v, prev2[k - 2] + 1);
        }
        row.push(v);
        if (v < low) low = v;
      }
      if (low > limit) return limit + 1;
      prev2 = prev;
      prev = row;
    }
    return prev[b.length];
  }

  // With nothing found: the nearest words the docs do use.
  function suggestion(words) {
    var vocab = buildVocabulary();
    var changed = false;
    var fixed = words.map(function (w) {
      var group = termGroups([w])[0];
      var known = Object.keys(vocab).some(function (v) {
        return group.some(function (t) {
          return v.indexOf(t.word) === 0 || (t.inside && v.indexOf(t.word) >= 0);
        });
      });
      if (known || w.length < 4 || !/^[a-z]+$/.test(w)) return w;
      var limit = w.length <= 5 ? 1 : 2;
      var best = null;
      var bestD = limit + 1;
      Object.keys(vocab).forEach(function (v) {
        var d = distance(w, v, limit);
        if (d < bestD || (d === bestD && best && vocab[v] > vocab[best])) {
          best = v;
          bestD = d;
        }
      });
      if (best && bestD <= limit) {
        changed = true;
        return best;
      }
      return w;
    });
    return changed ? fixed.join(" ") : null;
  }

  // Every entry that has every word, grouped by page. A page ranks by its best
  // entry plus how much of the page is about the words; Changelog entries
  // count half, so the page that explains a thing comes before the release
  // that changed it. When no entry has every word of three or more, the
  // entries that miss one word come instead, scored by how many they have.
  function searchDocs(query) {
    var words = queryWords(query);
    var result = { words: words, groups: [], total: 0, suggestion: null, partial: false };
    if (!words.length) return result;
    var groupsOfTerms = termGroups(words);
    var byPage = {};
    var order = [];
    var entries = indexEntries();
    // Each word against each entry, once. A word most entries have ("run")
    // weighs less than one few have ("slurm"), so in a question the word
    // that names the subject decides.
    var per = entries.map(function (e) {
      return groupsOfTerms.map(function (group) {
        return scoreTerm(group, e);
      });
    });
    var weights = groupsOfTerms.map(function (_g, k) {
      var df = 0;
      per.forEach(function (row) {
        if (row[k]) df++;
      });
      return Math.log(1 + entries.length / (1 + df));
    });
    var mean = weights.reduce(function (a, b) {
      return a + b;
    }, 0) / weights.length;
    var allowMissing = 0;
    var matched = function (row) {
      var score = 0;
      var tf = 0;
      var missing = 0;
      for (var k = 0; k < row.length; k++) {
        var s = row[k];
        if (!s) {
          if (++missing > allowMissing) return null;
          continue;
        }
        score += s.score * (mean ? weights[k] / mean : 1);
        tf += s.tf;
      }
      score *= (row.length - missing) / row.length;
      return { score: score, tf: tf };
    };
    var matches = per.map(matched);
    if (words.length >= 3 && !matches.some(Boolean)) {
      allowMissing = 1;
      result.partial = true;
      matches = per.map(matched);
    }
    entries.forEach(function (e, i) {
      var m = matches[i];
      if (!m) return;
      var score = m.score;
      var tf = m.tf;
      score += phraseBonus(words, e);
      if (e._kind === "changelog") score *= 0.5;
      var key = e._base;
      var g = byPage[key];
      if (!g) {
        g = byPage[key] = { page: e.page, url: e._base, top: null, sections: [], best: 0, tf: 0 };
        order.push(g);
      }
      var hit = { entry: e, score: score };
      if (!e.section && e.url === e._base) g.top = hit;
      else g.sections.push(hit);
      g.best = Math.max(g.best, score);
      g.tf += tf;
      result.total++;
    });
    order.forEach(function (g) {
      g.score = g.best + Math.log(1 + g.tf) * (g.top || g.sections.length ? 1 : 0);
      if (/^about\/changelog/.test(g.url)) g.score *= 0.6;
      g.sections.sort(function (a, b) {
        return b.score - a.score;
      });
    });
    order.sort(function (a, b) {
      return b.score - a.score;
    });
    result.groups = order;
    if (!result.total) result.suggestion = suggestion(words);
    return result;
  }

  // Where the words are in a piece of text, as [start, end) ranges.
  function ranges(text, words) {
    var lower = text.toLowerCase();
    var found = [];
    termGroups(words).forEach(function (group) {
      group.forEach(function (t) {
        [t.start, t.inside].forEach(function (re) {
          if (!re) return;
          re.lastIndex = 0;
          var m;
          var n = 0;
          while ((m = re.exec(lower)) && n++ < 40) {
            var at = m.index + m[1].length;
            found.push([at, at + m[2].length]);
          }
        });
      });
    });
    found.sort(function (a, b) {
      return a[0] - b[0];
    });
    var merged = [];
    found.forEach(function (r) {
      var last = merged[merged.length - 1];
      if (last && r[0] <= last[1]) last[1] = Math.max(last[1], r[1]);
      else merged.push(r.slice());
    });
    return merged;
  }

  function highlight(text, words) {
    var out = "";
    var at = 0;
    ranges(text, words).forEach(function (r) {
      out += escapeHtml(text.slice(at, r[0])) + "<mark>" + escapeHtml(text.slice(r[0], r[1])) + "</mark>";
      at = r[1];
    });
    return out + escapeHtml(text.slice(at));
  }

  // A line of prose around the first place the words appear, from the
  // section's text and never from its code.
  function snippet(text, words) {
    if (!text) return "";
    var found = ranges(text, words);
    var at = found.length ? found[0][0] : 0;
    var start = 0;
    if (at > 70) {
      var sentence = Math.max(text.lastIndexOf(". ", at), text.lastIndexOf("? ", at), text.lastIndexOf(": ", at));
      start = sentence >= 0 && at - sentence < 110 ? sentence + 2 : text.indexOf(" ", at - 60) + 1;
    }
    var end = start + 170;
    if (end < text.length) end = text.lastIndexOf(" ", end) > start ? text.lastIndexOf(" ", end) : end;
    else end = text.length;
    return (start > 0 ? "…" : "") + text.slice(start, end) + (end < text.length ? "…" : "");
  }

  // ── The search box: an ARIA combobox ──
  // Focus stays in the box; Up and Down move through the results, Enter opens
  // one, Esc closes. The panel closes when focus leaves the search.

  var searchBoxes = [];

  function bindSearch(input, panel, options) {
    var box = { input: input, panel: panel };
    var listId = panel.id || (panel.id = options.id);
    var container = options.container || input.parentElement;
    var list = [];
    var active = -1;
    var words = [];

    input.setAttribute("role", "combobox");
    input.setAttribute("aria-expanded", "false");
    input.setAttribute("aria-controls", listId);
    input.setAttribute("aria-autocomplete", "list");
    input.setAttribute("autocomplete", "off");
    input.setAttribute("spellcheck", "false");
    input.setAttribute("enterkeyhint", "go");
    panel.setAttribute("role", "listbox");
    panel.setAttribute("aria-label", "Search results");
    // A list that scrolls is a Tab stop of its own in some browsers; this one
    // is reached with the arrows, so Tab goes on past it and closes it.
    panel.setAttribute("tabindex", "-1");

    function show() {
      panel.hidden = false;
      input.setAttribute("aria-expanded", "true");
      if (options.onToggle) options.onToggle(true);
    }

    function hide() {
      panel.hidden = true;
      input.setAttribute("aria-expanded", "false");
      input.removeAttribute("aria-activedescendant");
      active = -1;
      if (options.onToggle) options.onToggle(false);
    }

    function select(k) {
      active = k;
      list.forEach(function (option, j) {
        var on = j === k;
        option.classList.toggle("active", on);
        option.setAttribute("aria-selected", on ? "true" : "false");
      });
      if (k < 0 || !list[k]) {
        input.removeAttribute("aria-activedescendant");
        return;
      }
      input.setAttribute("aria-activedescendant", list[k].id);
      var scroller = panel.scrollHeight > panel.clientHeight ? panel : options.scroller;
      if (scroller) {
        var o = list[k].getBoundingClientRect();
        var s = scroller.getBoundingClientRect();
        var pad = 36;
        if (o.top < s.top + pad) scroller.scrollTop -= s.top + pad - o.top;
        else if (o.bottom > s.bottom - pad) scroller.scrollTop += o.bottom - (s.bottom - pad);
      }
    }

    function option(k, cls, href, inner, extra) {
      return (
        '<a role="option" id="' + listId + "-" + k + '" class="' + cls + '" href="' + escapeHtml(href) + '" tabindex="-1" aria-selected="false"' +
        (extra || "") + ">" + inner + "</a>"
      );
    }

    function render(query) {
      var found = searchDocs(query);
      words = found.words;
      var html = "";
      var k = 0;
      var best = -1;
      var bestScore = -1;
      if (!found.total) {
        html += '<div class="r-count" aria-hidden="true">Nothing matches</div>';
        if (found.suggestion) {
          html += option(k++, "r-suggest", "#", "Did you mean <b>" + escapeHtml(found.suggestion) + "</b>?", ' data-suggest="' + escapeHtml(found.suggestion) + '"');
        } else {
          html += '<div class="none" aria-hidden="true">Try another word, or a flag such as --skip_typing.</div>';
        }
        announce("Nothing matches " + query + (found.suggestion ? ". Did you mean " + found.suggestion + "?" : ""), 350);
      } else {
        var groups = found.groups.slice(0, 8);
        var pages = found.groups.length;
        var line =
          (found.partial ? "No page has every word: " : "") +
          found.total + (found.total === 1 ? " result" : " results") + (found.partial ? " with all but one," : "") +
          " in " + pages + (pages === 1 ? " page" : " pages");
        html += '<div class="r-count" aria-hidden="true">' + line + (pages > groups.length ? ", the best " + groups.length + " shown" : "") + "</div>";
        groups.forEach(function (g, gi) {
          var more = g.sections.length;
          var meta = more ? more + (more === 1 ? " section" : " sections") : "";
          var top = g.top ? snippet(g.top.entry.text, words) : "";
          // Enter opens the best match of the best page.
          if (gi === 0 && g.top && g.top.score > bestScore) {
            best = k;
            bestScore = g.top.score;
          }
          html += option(
            k++,
            "r-page",
            ROOT + g.url,
            '<span class="r-title">' + highlight(g.page, words) + "</span>" +
              (meta ? '<span class="r-meta">' + meta + "</span>" : "") +
              (top ? "<small>" + highlight(top, words) + "</small>" : "")
          );
          g.sections.slice(0, 2).forEach(function (hit) {
            if (gi === 0 && hit.score > bestScore) {
              best = k;
              bestScore = hit.score;
            }
            var e = hit.entry;
            html += option(
              k++,
              "r-sec",
              ROOT + e.url,
              '<span class="r-title">' + highlight(e.section, words) + "</span>" +
                "<small>" + highlight(snippet(e.text, words), words) + "</small>"
            );
          });
        });
        announce(
          (found.partial ? "No page has every word. " : "") +
            found.total + (found.total === 1 ? " result" : " results") + (found.partial ? " with all but one word" : "") +
            " for " + query + ", in " + pages + (pages === 1 ? " page" : " pages"),
          350
        );
      }
      html += '<div class="r-foot" aria-hidden="true"><kbd>↑</kbd> <kbd>↓</kbd> to move, <kbd>Enter</kbd> to open, <kbd>Esc</kbd> to close</div>';
      panel.innerHTML = html;
      list = Array.prototype.slice.call(panel.querySelectorAll('[role="option"]'));
      panel.scrollTop = 0;
      show();
      // The match that Enter opens is marked from the start.
      select(found.total ? best : -1);
    }

    function run() {
      var query = input.value.trim();
      if (options.onQuery) options.onQuery(query);
      if (!queryWords(query).length) {
        list = [];
        hide();
        return;
      }
      if (indexState !== "ready") {
        loadIndex(run);
        if (indexState !== "ready") {
          panel.innerHTML =
            '<div class="none">' + (indexState === "failed" ? "The search index did not load. Reload the page to try again." : "Loading the search index…") + "</div>";
          list = [];
          show();
          return;
        }
      }
      render(query);
    }

    function choose(link, viaKey) {
      if (link.hasAttribute("data-suggest")) {
        input.value = link.getAttribute("data-suggest");
        run();
        input.focus();
        return false;
      }
      sessionSet("amap.searchHit", { url: link.href, words: words });
      var samePage = link.href.split("#")[0] === location.href.split("#")[0];
      if (samePage) {
        hide();
        if (options.onChoose) options.onChoose();
      }
      if (viaKey) location.href = link.href;
      return true;
    }

    input.addEventListener("input", run);
    input.addEventListener("focus", function () {
      loadIndex();
      if (input.value.trim() && panel.hidden) run();
    });
    container.addEventListener("pointerenter", function () {
      loadIndex();
    });
    input.addEventListener("keydown", function (event) {
      if (event.key === "ArrowDown" || event.key === "ArrowUp") {
        event.preventDefault();
        if (panel.hidden) run();
        if (!list.length) return;
        var step = event.key === "ArrowDown" ? 1 : -1;
        select(active < 0 ? (step > 0 ? 0 : list.length - 1) : (active + step + list.length) % list.length);
      } else if (event.key === "Home" && event.altKey) {
        event.preventDefault();
        if (list.length) select(0);
      } else if (event.key === "Enter") {
        var pick = list[active >= 0 ? active : 0];
        if (pick && !panel.hidden) {
          event.preventDefault();
          choose(pick, true);
        }
      } else if (event.key === "Escape" && !options.sheet) {
        if (!panel.hidden) {
          event.preventDefault();
          hide();
        } else if (input.value) {
          event.preventDefault();
          input.value = "";
        }
      }
    });
    // A press inside the panel keeps focus in the box, so the panel is still
    // there when the click lands.
    panel.addEventListener("mousedown", function (event) {
      event.preventDefault();
    });
    panel.addEventListener("click", function (event) {
      var link = event.target.closest('[role="option"]');
      if (!link) return;
      if (!choose(link, false)) event.preventDefault();
    });
    panel.addEventListener("mousemove", function (event) {
      var link = event.target.closest('[role="option"]');
      if (link && list.indexOf(link) !== active) select(list.indexOf(link));
    });
    container.addEventListener("focusout", function (event) {
      if (!event.relatedTarget || !container.contains(event.relatedTarget)) hide();
    });
    document.addEventListener("click", function (event) {
      // The header's Search button hands focus to this box; its click is not
      // a click away.
      if (!container.contains(event.target) && !event.target.closest(".search-btn")) hide();
    });

    box.run = run;
    box.hide = hide;
    searchBoxes.push(box);
    return box;
  }

  // The sidebar box of the documentation pages.
  var sideInput = document.querySelector("[data-docs-search]");
  var sidePanel = document.querySelector("[data-docs-results]");
  if (sideInput && sidePanel) {
    var key = document.createElement("kbd");
    key.className = "search-key";
    key.setAttribute("aria-hidden", "true");
    key.textContent = "/";
    sideInput.insertAdjacentElement("afterend", key);
    sideInput.setAttribute("aria-keyshortcuts", "/ Control+K Meta+K");
    bindSearch(sideInput, sidePanel, { id: "docs-results", container: sideInput.closest(".docs-search") });
  }

  // ── Search from every page: a button in the header, and a sheet under it ──

  function pageName() {
    var title = document.title.split(" | ")[0];
    if (/^AMAP:/.test(title)) return "Front page";
    if (/^How AMAP decides/.test(title)) return "How it works, illustrated";
    return title;
  }

  // The page read before this one, for the sheet's "Last visited".
  var here = { title: pageName(), url: location.href.split("#")[0] };
  var trail = (sessionGet("amap.trail") || []).filter(function (p) {
    return p && p.url && p.url !== here.url;
  });
  var lastVisited = trail[0] || null;
  if (here.title !== "Not here") sessionSet("amap.trail", [here].concat(trail).slice(0, 5));

  var GLASS =
    '<svg viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2.2" stroke-linecap="round" aria-hidden="true"><circle cx="10.5" cy="10.5" r="6.5"/><path d="M15.5 15.5 21 21"/></svg>';
  var QUICK = [
    ["user-guide/parameters.html", "Parameters", "Every option and its default"],
    ["user-guide/outputs.html", "Outputs", "Every file a run writes"],
    ["deployment/troubleshooting.html", "Troubleshooting", "When a run stops or looks wrong"],
    ["getting-started/tutorial.html", "Tutorial", "One genome, start to finish"],
    ["how-it-works.html", "How it works", "The correction, step by step"]
  ];

  var nav = document.querySelector(".top .nav");
  var searchButton = null;
  var sheet = null;
  var sheetBox = null;
  var sheetOpener = null;
  var scrim = null;

  function sidebarVisible() {
    return !!(sideInput && window.matchMedia && window.matchMedia("(min-width: 861px)").matches && sideInput.getClientRects().length);
  }

  function openSearch(opener) {
    loadIndex();
    if (sidebarVisible()) {
      sideInput.focus();
      sideInput.select();
      return;
    }
    if (sheet) openSheet(opener);
  }

  if (nav && header) {
    searchButton = document.createElement("button");
    searchButton.type = "button";
    searchButton.className = "search-btn";
    searchButton.setAttribute("aria-label", "Search the documentation");
    searchButton.setAttribute("aria-keyshortcuts", "/ Control+K Meta+K");
    searchButton.setAttribute("aria-expanded", "false");
    searchButton.setAttribute("aria-controls", "search-sheet");
    searchButton.innerHTML = GLASS + '<span class="search-btn-label">Search</span><kbd class="search-key" aria-hidden="true">/</kbd>';
    var before = nav.querySelector(".gh-link") || nav.querySelector("[data-theme-toggle]");
    nav.insertBefore(searchButton, before);
    searchButton.addEventListener("pointerenter", function () {
      loadIndex();
    });
    searchButton.addEventListener("click", function () {
      if (sheet && !sheet.hidden) closeSheet(true);
      else openSearch(searchButton);
    });

    var quick = QUICK.map(function (q) {
      return '<li><a href="' + escapeHtml(ROOT + q[0]) + '"><b>' + escapeHtml(q[1]) + "</b><small>" + escapeHtml(q[2]) + "</small></a></li>";
    });
    if (lastVisited) {
      quick.push(
        '<li><a class="sheet-last" href="' + escapeHtml(lastVisited.url) + '"><b>' + escapeHtml(lastVisited.title) +
          "</b><small>Last visited</small></a></li>"
      );
    }
    quick.push('<li><a href="' + escapeHtml(ROOT + "documentation.html") + '"><b>All the docs</b><small>Every page, by section</small></a></li>');

    sheet = document.createElement("div");
    sheet.className = "search-sheet";
    sheet.id = "search-sheet";
    sheet.hidden = true;
    sheet.setAttribute("role", "dialog");
    sheet.setAttribute("aria-modal", "true");
    sheet.setAttribute("aria-label", "Search the documentation");
    sheet.innerHTML =
      '<div class="wrap">' +
      '<div class="docs-search sheet-search">' + GLASS +
      '<input type="search" placeholder="Search the docs" aria-label="Search the documentation">' +
      '<button type="button" class="sheet-close" aria-label="Close the search (Esc)"><span class="on-keys">Esc</span><span class="on-touch">Close</span></button>' +
      '<div class="docs-results" hidden></div>' +
      "</div>" +
      '<div class="sheet-empty">' +
      '<nav class="sheet-toc" aria-labelledby="sheet-toc-label" hidden><p class="sheet-label" id="sheet-toc-label">On this page</p><ol></ol></nav>' +
      '<nav aria-labelledby="sheet-go-label"><p class="sheet-label" id="sheet-go-label">Go to</p><ul class="sheet-links">' + quick.join("") + "</ul></nav>" +
      "</div>" +
      "</div>";
    scrim = document.createElement("div");
    scrim.className = "search-scrim";
    scrim.hidden = true;
    header.insertAdjacentElement("afterend", sheet);
    header.insertAdjacentElement("afterend", scrim);

    var sheetInput = sheet.querySelector("input");
    sheetBox = bindSearch(sheetInput, sheet.querySelector(".docs-results"), {
      id: "sheet-results",
      container: sheet,
      scroller: sheet,
      sheet: true,
      onQuery: function (query) {
        sheet.classList.toggle("has-query", queryWords(query).length > 0);
      },
      onChoose: function () {
        closeSheet(false);
      }
    });
    sheet.querySelector(".sheet-close").addEventListener("click", function () {
      closeSheet(true);
    });
    sheet.addEventListener("keydown", function (event) {
      if (event.key === "Escape") {
        event.preventDefault();
        closeSheet(true);
      }
    });
    sheet.addEventListener("focusout", function (event) {
      var to = event.relatedTarget;
      if (to && !sheet.contains(to) && to !== sheetOpener) closeSheet(false);
    });
    sheet.addEventListener("click", function (event) {
      var link = event.target.closest(".sheet-empty a");
      if (link && link.href.split("#")[0] === location.href.split("#")[0]) closeSheet(false);
    });
    scrim.addEventListener("click", function () {
      closeSheet(true);
    });
    document.addEventListener("click", function (event) {
      if (!sheet.hidden && !sheet.contains(event.target) && !event.target.closest(".search-btn, .here-bar")) closeSheet(false);
    });
  }

  function openSheet(opener) {
    if (!sheet) return;
    sheetOpener = opener || document.activeElement;
    sheet.hidden = false;
    scrim.hidden = false;
    [searchButton, hereBar].forEach(function (b) {
      if (b) b.setAttribute("aria-expanded", "true");
    });
    document.documentElement.classList.add("sheet-open");
    // The page under the sheet is out of reach while it is open; the header
    // stays, and focus that goes there closes the sheet.
    Array.prototype.forEach.call(document.body.children, function (el) {
      if (el === sheet || el === scrim || el === header || el === hereBar || el === statusLine || el.inert || /^(SCRIPT|STYLE)$/.test(el.tagName)) return;
      el.inert = true;
      el.setAttribute("data-sheet-inert", "");
    });
    sheet.scrollTop = 0;
    // From the bar that names the section, the contents come first, with no
    // keyboard popping up over them; from anywhere else, the search box.
    var on = sheet.querySelector(".sheet-toc:not([hidden]) a.active");
    if (opener && opener === hereBar && on) {
      on.focus({ preventScroll: true });
      var top = on.getBoundingClientRect().top - sheet.getBoundingClientRect().top;
      if (top > sheet.clientHeight - 60) sheet.scrollTop = top - sheet.clientHeight / 2;
      return;
    }
    var input = sheet.querySelector("input");
    input.focus();
    input.select();
  }

  function closeSheet(restore) {
    if (!sheet || sheet.hidden) return;
    sheet.hidden = true;
    scrim.hidden = true;
    sheetBox.hide();
    [searchButton, hereBar].forEach(function (b) {
      if (b) b.setAttribute("aria-expanded", "false");
    });
    document.documentElement.classList.remove("sheet-open");
    document.querySelectorAll("[data-sheet-inert]").forEach(function (el) {
      el.inert = false;
      el.removeAttribute("data-sheet-inert");
    });
    if (restore && sheetOpener && sheetOpener.focus) sheetOpener.focus();
  }

  // "/" and Ctrl or Cmd+K open the search from anywhere on the site. A key
  // of a single character can be turned off (WCAG 2.1.4), from the footer of
  // every page; this browser remembers it.
  var slashOn = true;
  try {
    slashOn = localStorage.getItem("amap.slash") !== "off";
  } catch (e) {}
  var slashToggle = null;
  function setSlash(on, save) {
    slashOn = on;
    if (save) {
      try {
        if (on) localStorage.removeItem("amap.slash");
        else localStorage.setItem("amap.slash", "off");
      } catch (e) {}
    }
    document.documentElement.classList.toggle("no-slash", !on);
    document.querySelectorAll("[aria-keyshortcuts]").forEach(function (el) {
      el.setAttribute("aria-keyshortcuts", (on ? "/ " : "") + "Control+K Meta+K");
    });
    if (slashToggle) {
      slashToggle.textContent = on ? "The / key opens the search: turn it off" : "The / key does nothing: turn it on";
    }
  }
  var foot = document.querySelector("footer.foot");
  var footBase = foot && foot.querySelector(".foot-base");
  if (foot && !footBase) {
    footBase = document.createElement("div");
    footBase.className = "wrap foot-base";
    foot.appendChild(footBase);
  }
  if (footBase && searchButton) {
    slashToggle = document.createElement("button");
    slashToggle.type = "button";
    slashToggle.className = "slash-toggle";
    slashToggle.addEventListener("click", function () {
      setSlash(!slashOn, true);
      announce(slashOn ? "The / key opens the search" : "The / key no longer opens the search");
    });
    footBase.appendChild(slashToggle);
  }
  setSlash(slashOn, false);

  document.addEventListener("keydown", function (event) {
    var t = event.target;
    var typing = t && (/^(input|textarea|select)$/i.test(t.tagName) || t.isContentEditable);
    var slash = slashOn && event.key === "/" && !typing && !event.ctrlKey && !event.metaKey && !event.altKey;
    var k = (event.key === "k" || event.key === "K") && (event.ctrlKey || event.metaKey) && !event.altKey && !event.shiftKey;
    if (!slash && !k) return;
    event.preventDefault();
    openSearch(document.activeElement && document.activeElement !== document.body ? document.activeElement : searchButton);
  });

  // ── The header keeps to one row ──
  // The stylesheet drops header items at widths measured with the usual
  // fonts. A wider font, larger text or extra letter spacing can still run
  // the row over, so the header then sheds the same items, one step at a
  // time, until it fits (the squeeze-1 to squeeze-6 rules in site.css).
  var SQUEEZE_STEPS = 6;
  function fitHeader() {
    for (var s = 1; s <= SQUEEZE_STEPS; s++) header.classList.remove("squeeze-" + s);
    for (var step = 1; step <= SQUEEZE_STEPS && nav.scrollWidth > nav.clientWidth + 1; step++) {
      header.classList.add("squeeze-" + step);
    }
  }

  if (nav && header) {
    fitHeader();
    var fitQueued = false;
    window.addEventListener("resize", function () {
      if (fitQueued) return;
      fitQueued = true;
      requestAnimationFrame(function () {
        fitQueued = false;
        fitHeader();
      });
    });
  }

  // ── The page contents follow the reader down the page ──

  var tocLinks = Array.prototype.slice.call(document.querySelectorAll(".docs-toc a, .toc-inline a"));
  var tocTargets = [];
  var seenTarget = {};
  tocLinks.forEach(function (link) {
    var id = decodeURIComponent((link.getAttribute("href") || "").slice(1));
    var el = id && !seenTarget[id] ? document.getElementById(id) : null;
    if (!el) return;
    seenTarget[id] = true;
    tocTargets.push({ id: id, el: el, text: link.textContent.trim(), level: link.closest(".toc-3") ? 3 : 2 });
  });

  var hereBar = null;
  var hereLabel = null;
  var hereProgress = null;
  var currentId = null;

  if (tocTargets.length) {
    // The sheet lists this page's sections, so a phone has its contents one
    // tap away (from the header's search, or the bar that names the section).
    if (sheet) {
      var sheetToc = sheet.querySelector(".sheet-toc");
      sheetToc.querySelector("ol").innerHTML = tocTargets
        .filter(function (t) {
          return t.level === 2;
        })
        .map(function (t) {
          return '<li><a href="#' + escapeHtml(t.id) + '">' + escapeHtml(t.text) + "</a></li>";
        })
        .join("");
      sheetToc.hidden = false;
      Array.prototype.push.apply(tocLinks, sheetToc.querySelectorAll("a"));
    }
    var main = document.querySelector("main.prose");
    if (main && header) {
      hereBar = document.createElement("button");
      hereBar.type = "button";
      hereBar.className = "here-bar";
      hereBar.setAttribute("aria-controls", "search-sheet");
      hereBar.setAttribute("aria-expanded", "false");
      hereBar.innerHTML =
        '<span class="here-label"></span>' +
        '<svg viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round" stroke-linejoin="round" aria-hidden="true"><path d="m6 9 6 6 6-6"/></svg>' +
        '<i class="here-progress" aria-hidden="true"></i>';
      hereLabel = hereBar.querySelector(".here-label");
      hereProgress = hereBar.querySelector(".here-progress");
      (sheet || header).insertAdjacentElement("afterend", hereBar);
      document.documentElement.classList.add("has-here-bar");
      hereBar.addEventListener("click", function () {
        openSheet(hereBar);
      });
    }
  }

  // Where an anchored heading lands: the page's scroll padding plus the
  // heading's own scroll margin, read from the stylesheet rather than assumed.
  var landing = [];
  function measureLanding() {
    var pad = parseFloat(getComputedStyle(document.documentElement).scrollPaddingTop) || 0;
    landing = tocTargets.map(function (t) {
      return pad + (parseFloat(getComputedStyle(t.el).scrollMarginTop) || 0);
    });
  }

  function setCurrent(id) {
    currentId = id;
    tocLinks.forEach(function (link) {
      var on = link.getAttribute("href") === "#" + id;
      link.classList.toggle("active", on);
      if (on) link.setAttribute("aria-current", "location");
      else link.removeAttribute("aria-current");
    });
    updateHere();
  }

  function updateHere() {
    if (!hereBar) return;
    var title = document.querySelector("main.prose h1");
    var past = title ? title.getBoundingClientRect().bottom < HEADER_GAP : window.scrollY > 240;
    hereBar.classList.toggle("on", past);
    var current = tocTargets.filter(function (t) {
      return t.id === currentId;
    })[0];
    var text = current ? current.text : pageName();
    if (hereLabel.textContent !== text) {
      hereLabel.textContent = text;
      hereBar.setAttribute("aria-label", text + ": this page's contents, and search");
    }
    var room = document.documentElement.scrollHeight - window.innerHeight;
    var progress = room > 0 ? Math.min(1, Math.max(0, window.scrollY / room)) : 0;
    hereProgress.style.transform = "scaleX(" + progress.toFixed(3) + ")";
  }

  // A click in the contents marks that link at once; the scroll it starts
  // does not move the mark until it ends.
  var locked = false;
  var lockTimer = 0;
  function lockMark() {
    locked = true;
    clearTimeout(lockTimer);
    lockTimer = setTimeout(unlockMark, 700);
  }

  function unlockMark() {
    locked = false;
    clearTimeout(lockTimer);
  }

  if (tocTargets.length) {
    var ticking = false;
    var markCurrent = function () {
      ticking = false;
      updateHere();
      if (locked) return;
      var doc = document.documentElement;
      var scrolls = doc.scrollHeight > window.innerHeight + 2;
      var atBottom = scrolls && window.innerHeight + window.scrollY >= doc.scrollHeight - 2;
      var k = 0;
      if (atBottom) {
        k = tocTargets.length - 1;
      } else {
        tocTargets.forEach(function (t, j) {
          if (t.el.getBoundingClientRect().top < landing[j] + 24) k = j;
        });
      }
      if (tocTargets[k].id !== currentId) setCurrent(tocTargets[k].id);
    };
    measureLanding();
    window.addEventListener(
      "scroll",
      function () {
        if (locked) {
          // While a jump is under way, wait for it to settle.
          clearTimeout(lockTimer);
          lockTimer = setTimeout(unlockMark, 180);
        }
        if (!ticking) {
          ticking = true;
          requestAnimationFrame(markCurrent);
        }
      },
      { passive: true }
    );
    window.addEventListener("scrollend", function () {
      if (locked) unlockMark();
    });
    window.addEventListener("resize", function () {
      measureLanding();
      markCurrent();
    });
    tocLinks.forEach(function (link) {
      link.addEventListener("click", function () {
        var id = decodeURIComponent((link.getAttribute("href") || "").slice(1));
        if (!seenTarget[id]) return;
        lockMark();
        setCurrent(id);
      });
    });
    var arrivedAt = decodeURIComponent(location.hash.slice(1));
    if (arrivedAt && seenTarget[arrivedAt]) {
      setCurrent(arrivedAt);
      lockMark();
    } else {
      markCurrent();
    }
  }

  // ── Arriving at a heading ──
  // The heading a link lands on carries a short accent bar, and the words of
  // the search that brought the reader here are marked in its section.

  function targetOfHash() {
    var id = decodeURIComponent(location.hash.slice(1));
    return id ? document.getElementById(id) : null;
  }

  function arrive(target) {
    if (!target || reducedMotion() || !/^H[1-6]$/.test(target.tagName)) return;
    target.classList.remove("arrive");
    void target.offsetWidth;
    target.classList.add("arrive");
    setTimeout(function () {
      target.classList.remove("arrive");
    }, 1200);
  }

  function headingLevel(node) {
    if (!node || node.nodeType !== 1) return 0;
    var h = /^H[1-6]$/.test(node.tagName) ? node : node.classList.contains("hx") ? node.querySelector("h1, h2, h3, h4, h5, h6") : null;
    return h ? +h.tagName.charAt(1) : 0;
  }

  function markWords(root, words, limit) {
    var walker = document.createTreeWalker(root, NodeFilter.SHOW_TEXT, {
      acceptNode: function (node) {
        return node.parentElement.closest("pre, script, style, svg, button, .anchor, mark, kbd") ? NodeFilter.FILTER_REJECT : NodeFilter.FILTER_ACCEPT;
      }
    });
    var nodes = [];
    while (walker.nextNode()) nodes.push(walker.currentNode);
    var made = 0;
    nodes.forEach(function (node) {
      if (made >= limit) return;
      var found = ranges(node.nodeValue, words);
      if (!found.length) return;
      var frag = document.createDocumentFragment();
      var at = 0;
      found.forEach(function (r) {
        if (made >= limit) return;
        frag.appendChild(document.createTextNode(node.nodeValue.slice(at, r[0])));
        var mark = document.createElement("mark");
        mark.className = "search-hit";
        mark.textContent = node.nodeValue.slice(r[0], r[1]);
        frag.appendChild(mark);
        at = r[1];
        made++;
      });
      frag.appendChild(document.createTextNode(node.nodeValue.slice(at)));
      node.parentNode.replaceChild(frag, node);
    });
  }

  function markSearchArrival() {
    var hit = sessionGet("amap.searchHit");
    if (!hit || !hit.url || !hit.words || !hit.words.length) return;
    if (hit.url.split("#")[0] !== location.href.split("#")[0] || (hit.url.split("#")[1] || "") !== location.hash.slice(1)) return;
    sessionSet("amap.searchHit", null);
    var main = document.querySelector("main");
    var start = targetOfHash() || (main && main.querySelector("h1"));
    if (!start) return;
    var level = headingLevel(start);
    if (!level) {
      markWords(start, hit.words, 40);
      return;
    }
    var node = start.parentElement && start.parentElement.classList.contains("hx") ? start.parentElement : start;
    var section = [node];
    for (var next = node.nextElementSibling; next; next = next.nextElementSibling) {
      var l = headingLevel(next);
      if (l && l <= Math.max(level, 2)) break;
      section.push(next);
    }
    var left = 40;
    section.forEach(function (el) {
      markWords(el, hit.words, left);
      left = 40 - document.querySelectorAll("mark.search-hit").length;
    });
  }

  arrive(targetOfHash());
  markSearchArrival();
  window.addEventListener("hashchange", function () {
    arrive(targetOfHash());
    markSearchArrival();
  });

  // ── The page list keeps the current page in view ──
  // The list scrolls on its own; it opens where it was on the page before,
  // and if the current page is out of sight there, with that page centred.
  var menu = document.querySelector(".docs-menu");
  if (menu) {
    // On a phone the list starts folded, so the page itself comes first. A
    // page that folds it before the first paint (a script right after the
    // list) is left alone; any other page is folded here.
    var foldsItself = menu.nextElementSibling && menu.nextElementSibling.tagName === "SCRIPT";
    if (!foldsItself && window.matchMedia && window.matchMedia("(max-width: 860px)").matches) menu.open = false;
    var currentPage = menu.querySelector('a[aria-current="page"]');
    if (menu.scrollHeight > menu.clientHeight + 1) {
      var saved = sessionGet("amap.menuScroll");
      if (typeof saved === "number") menu.scrollTop = saved;
      if (currentPage) {
        var m = menu.getBoundingClientRect();
        var r = currentPage.getBoundingClientRect();
        if (r.top < m.top + 8 || r.bottom > m.bottom - 8) {
          menu.scrollTop += r.top - m.top - (menu.clientHeight - r.height) / 2;
        }
      }
    }
    var saveMenu = function () {
      sessionSet("amap.menuScroll", menu.scrollTop);
    };
    menu.addEventListener("scroll", saveMenu, { passive: true });
    window.addEventListener("pagehide", saveMenu);
  }

  // ── The header's "you are here" ──
  // On the front page, "Start here" is current while its section is on screen.
  var startLink = document.querySelector('.top .nav a[href="#start"]');
  var startSection = document.getElementById("start");
  if (startLink && startSection && "IntersectionObserver" in window) {
    new IntersectionObserver(
      function (entries) {
        entries.forEach(function (entry) {
          if (entry.isIntersecting) startLink.setAttribute("aria-current", "location");
          else startLink.removeAttribute("aria-current");
        });
      },
      { rootMargin: "-45% 0px -45% 0px" }
    ).observe(startSection);
  }
})();

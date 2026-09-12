(function () {
  "use strict";

  // This report owns both themes; prevent automatic page recolouring.
  if (!document.querySelector('meta[name="darkreader-lock"]')) {
    var themeLock = document.createElement("meta");
    themeLock.name = "darkreader-lock";
    document.head.appendChild(themeLock);
  }

  var themeKey = "dnaepico.theme";
  var theme = "dark";
  try {
    var savedTheme = window.localStorage.getItem(themeKey);
    if (savedTheme === "light" || savedTheme === "dark") theme = savedTheme;
  } catch (error) { /* Local file storage may be unavailable. */ }
  if (window.location.protocol === "file:") {
    var fileTheme = window.name.match(/;dnaepico-theme:(light|dark)/);
    if (fileTheme) theme = fileTheme[1];
  }
  document.documentElement.dataset.dnaTheme = theme;

  function saveFilePreference(key, value) {
    if (window.location.protocol !== "file:" ||
        (window.name && window.name.indexOf("dnaepico-figure-zoom:") !== 0)) return;
    var baseName = window.name.replace(new RegExp(";dnaepico-" + key + ":[^;]*", "g"), "") ||
      "dnaepico-figure-zoom:1";
    window.name = baseName + ";dnaepico-" + key + ":" + value;
  }

  function themeIcon(theme) {
    var namespace = "http://www.w3.org/2000/svg";
    var icon = document.createElementNS(namespace, "svg");
    icon.setAttribute("viewBox", "0 0 24 24");
    icon.setAttribute("aria-hidden", "true");
    icon.setAttribute("focusable", "false");
    var disc = document.createElementNS(namespace, theme === "dark" ? "circle" : "path");
    disc.setAttribute("class", "dnaepico-theme-disc");
    if (theme === "dark") {
      disc.setAttribute("cx", "12");
      disc.setAttribute("cy", "12");
      disc.setAttribute("r", "4");
      var rays = document.createElementNS(namespace, "path");
      rays.setAttribute("d", "M12 2v2M12 20v2M2 12h2M20 12h2M4.93 4.93l1.42 1.42M17.65 17.65l1.42 1.42M4.93 19.07l1.42-1.42M17.65 6.35l1.42-1.42");
      icon.appendChild(rays);
    } else {
      disc.setAttribute("d", "M20.5 13.1A8.5 8.5 0 0 1 10.9 3.5a8.5 8.5 0 1 0 9.6 9.6Z");
    }
    icon.appendChild(disc);
    return icon;
  }

  function initializeTheme() {
    var brand = document.querySelector("#quarto-header .navbar-brand-container");
    if (!brand || document.getElementById("dnaepico-theme-toggle")) return;
    var button = document.createElement("button");
    button.id = "dnaepico-theme-toggle";
    button.type = "button";
    function updateTheme() {
      document.documentElement.dataset.dnaTheme = theme;
      button.replaceChildren(themeIcon(theme));
      var label = "Switch to " + (theme === "dark" ? "light" : "dark") + " mode";
      button.setAttribute("aria-label", label);
      button.title = label;
    }
    button.addEventListener("click", function () {
      theme = theme === "dark" ? "light" : "dark";
      try { window.localStorage.setItem(themeKey, theme); } catch (error) { /* optional */ }
      saveFilePreference("theme", theme);
      updateTheme();
    });
    updateTheme();
    brand.appendChild(button);
  }

  function configureSearchForProtocol() {
    if (window.location.protocol !== "file:") {
      return;
    }

    document.documentElement.classList.add("dnaepico-file-protocol");
    var search = document.getElementById("quarto-search");
    if (search) {
      search.setAttribute("aria-hidden", "true");
      search.remove();
    }
  }

  function cardTitle(control) {
    var card = control.closest(".bslib-card");
    var header = card ? card.querySelector(":scope > .card-header") : null;
    var title = header ? header.textContent.trim() : "report panel";
    return title || "report panel";
  }

  function updateControl(control) {
    var card = control.closest(".bslib-card");
    var action = card && card.dataset.fullScreen === "true" ?
      "Collapse " : "Expand ";
    control.setAttribute("role", "button");
    control.setAttribute("aria-label", action + cardTitle(control));
  }

  function initializeControl(control) {
    if (control.dataset.keyboardExpandInitialized === "true") {
      updateControl(control);
      return;
    }

    control.dataset.keyboardExpandInitialized = "true";
    updateControl(control);
    control.addEventListener("keydown", function (event) {
      if (event.key !== "Enter" && event.key !== " ") {
        return;
      }
      event.preventDefault();
      control.click();
    });
  }

  function initializeAll() {
    initializeTheme();
    configureSearchForProtocol();
    document.querySelectorAll(".bslib-full-screen-enter").forEach(
      initializeControl
    );
  }

  initializeAll();
  if (document.readyState === "loading") {
    document.addEventListener("DOMContentLoaded", initializeAll, {
      once: true
    });
  }

  if (window.MutationObserver) {
    var observer = new MutationObserver(function (changes) {
      configureSearchForProtocol();
      changes.forEach(function (change) {
        if (
          change.type === "attributes" &&
          change.target.matches(".bslib-card")
        ) {
          change.target.querySelectorAll(
            ".bslib-full-screen-enter"
          ).forEach(updateControl);
        }
      });
      initializeAll();
    });
    observer.observe(document.documentElement, {
      attributes: true,
      attributeFilter: ["data-full-screen"],
      childList: true,
      subtree: true
    });
  }
}());

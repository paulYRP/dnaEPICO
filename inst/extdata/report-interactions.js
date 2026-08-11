(function () {
  "use strict";

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

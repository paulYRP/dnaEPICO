(function () {
  "use strict";

  var minimumZoom = 0.25;
  var maximumZoom = 4;
  var zoomStep = 1.15;
  var zoomStorageKey = "dnaepico.figureZoom";
  var windowNamePrefix = "dnaepico-figure-zoom:";

  function clamp(value, minimum, maximum) {
    return Math.min(maximum, Math.max(minimum, value));
  }

  function validZoom(value) {
    if (value === null || value === undefined || value === "") {
      return null;
    }
    var parsed = Number(value);
    if (!Number.isFinite(parsed)) {
      return null;
    }
    return clamp(parsed, minimumZoom, maximumZoom);
  }

  function readStoredZoom() {
    var stored = null;
    try {
      stored = validZoom(window.sessionStorage.getItem(zoomStorageKey));
    } catch (error) {
      stored = null;
    }
    if (stored !== null) {
      return stored;
    }

    if (
      window.location.protocol === "file:" &&
      window.name.indexOf(windowNamePrefix) === 0
    ) {
      stored = validZoom(window.name.slice(windowNamePrefix.length));
    }
    return stored === null ? 1 : stored;
  }

  var sharedZoom = readStoredZoom();

  function storeSharedZoom(value) {
    sharedZoom = clamp(value, minimumZoom, maximumZoom);
    try {
      window.sessionStorage.setItem(zoomStorageKey, String(sharedZoom));
    } catch (error) {
      // Storage can be unavailable for local-file reports.
    }
    if (
      window.location.protocol === "file:" &&
      (!window.name || window.name.indexOf(windowNamePrefix) === 0)
    ) {
      window.name = windowNamePrefix + String(sharedZoom);
    }
    document.dispatchEvent(new CustomEvent(
      "dnaepico:figure-zoom-change",
      { detail: { zoom: sharedZoom } }
    ));
  }

  function parseFigures(controls) {
    var data = controls.querySelector('[data-role="figure-data"]');
    if (!data) {
      return [];
    }

    try {
      return JSON.parse(data.textContent);
    } catch (error) {
      return [];
    }
  }

  function initializeViewer(content) {
    if (content.dataset.figureViewerInitialized === "true") {
      return;
    }

    var browserId = content.dataset.browserId;
    var controls = document.getElementById(browserId + "-controls");
    var card = content.closest(".dnaepico-selected-figure");
    if (!controls || !card) {
      return;
    }

    var figures = parseFigures(controls);
    var select = controls.querySelector('[data-role="figure-select"]');
    var title = content.querySelector('[data-role="figure-title"]');
    var canvas = content.querySelector('[data-role="figure-canvas"]');
    var stage = content.querySelector('[data-role="figure-stage"]');
    var image = content.querySelector('[data-role="figure-image"]');
    var fallback = content.querySelector('[data-role="figure-fallback"]');
    var description = controls.querySelector(
      '[data-role="figure-description"]'
    );
    var download = controls.querySelector('[data-role="figure-download"]');
    var count = controls.querySelector('[data-role="figure-count"]');
    var previous = controls.querySelector('[data-role="figure-prev"]');
    var next = controls.querySelector('[data-role="figure-next"]');
    var zoomStatus = content.querySelector(
      '[data-role="figure-zoom-status"]'
    );

    if (
      !figures.length || !select || !title || !canvas || !stage ||
      !image || !fallback || !download || !count || !previous || !next
    ) {
      return;
    }

    content.dataset.figureViewerInitialized = "true";
    var zoom = sharedZoom;
    var geometryFrame = 0;

    function updateZoomStatus() {
      var percentage = Math.round(zoom * 100);
      var message = percentage === 100 ?
        "Image fitted to the available panel." :
        "Image zoom is " + percentage + "% of the fitted size.";
      if (zoomStatus) {
        zoomStatus.textContent = message;
      }
    }

    function renderGeometry() {
      if (!image.naturalWidth || !image.naturalHeight || image.hidden) {
        return;
      }

      var availableWidth = Math.max(1, canvas.clientWidth);
      var availableHeight = Math.max(1, canvas.clientHeight);
      var fitScale = Math.min(
        availableWidth / image.naturalWidth,
        availableHeight / image.naturalHeight
      );
      var imageWidth = Math.max(1, image.naturalWidth * fitScale * zoom);
      var imageHeight = Math.max(1, image.naturalHeight * fitScale * zoom);

      image.style.width = imageWidth + "px";
      image.style.height = imageHeight + "px";
      stage.style.width = Math.max(availableWidth, imageWidth) + "px";
      stage.style.height = Math.max(availableHeight, imageHeight) + "px";
      updateZoomStatus();
    }

    function scheduleGeometry() {
      if (geometryFrame) {
        window.cancelAnimationFrame(geometryFrame);
      }
      geometryFrame = window.requestAnimationFrame(function () {
        geometryFrame = 0;
        renderGeometry();
      });
    }

    function changeZoom(nextZoom, offsetX, offsetY) {
      var oldWidth = Math.max(1, canvas.scrollWidth);
      var oldHeight = Math.max(1, canvas.scrollHeight);
      var anchorX = (canvas.scrollLeft + offsetX) / oldWidth;
      var anchorY = (canvas.scrollTop + offsetY) / oldHeight;

      zoom = clamp(nextZoom, minimumZoom, maximumZoom);
      storeSharedZoom(zoom);
      renderGeometry();
      window.requestAnimationFrame(function () {
        canvas.scrollLeft = anchorX * canvas.scrollWidth - offsetX;
        canvas.scrollTop = anchorY * canvas.scrollHeight - offsetY;
      });
    }

    function resetZoom() {
      zoom = 1;
      storeSharedZoom(zoom);
      renderGeometry();
      canvas.scrollTo(0, 0);
    }

    function showFigure() {
      var index = clamp(Number(select.value) - 1, 0, figures.length - 1);
      var item = figures[index];
      title.textContent = item.title;
      download.href = item.downloadPath;
      download.download = item.downloadName;
      download.textContent = "Download " + item.downloadName;
      if (description) {
        description.textContent = item.description || "";
      }
      count.textContent = (index + 1) + " of " + figures.length + " figures";
      previous.disabled = index === 0;
      next.disabled = index === figures.length - 1;
      zoom = sharedZoom;

      if (item.browserReady) {
        fallback.hidden = true;
        image.hidden = false;
        image.alt = item.title;
        image.src = item.previewPath;
        if (image.complete && image.naturalWidth) {
          scheduleGeometry();
        }
      } else {
        image.removeAttribute("src");
        image.hidden = true;
        fallback.textContent =
          "Browser preview is unavailable for this TIFF file.";
        fallback.hidden = false;
        updateZoomStatus();
      }
      canvas.scrollTo(0, 0);
    }

    function move(delta) {
      var index = clamp(
        Number(select.value) - 1 + delta,
        0,
        figures.length - 1
      );
      select.value = String(index + 1);
      showFigure();
    }

    select.addEventListener("change", showFigure);
    previous.addEventListener("click", function () {
      move(-1);
    });
    next.addEventListener("click", function () {
      move(1);
    });
    image.addEventListener("load", scheduleGeometry);
    image.addEventListener("error", function () {
      image.hidden = true;
      fallback.textContent = "The browser could not display this figure.";
      fallback.hidden = false;
    });

    canvas.addEventListener("wheel", function (event) {
      if (!event.ctrlKey && !event.metaKey) {
        return;
      }
      event.preventDefault();
      var bounds = canvas.getBoundingClientRect();
      var factor = event.deltaY < 0 ? zoomStep : 1 / zoomStep;
      changeZoom(
        zoom * factor,
        event.clientX - bounds.left,
        event.clientY - bounds.top
      );
    }, { passive: false });

    canvas.addEventListener("keydown", function (event) {
      var offsetX = canvas.clientWidth / 2;
      var offsetY = canvas.clientHeight / 2;
      if (event.key === "+" || event.key === "=") {
        event.preventDefault();
        changeZoom(zoom * zoomStep, offsetX, offsetY);
      } else if (event.key === "-") {
        event.preventDefault();
        changeZoom(zoom / zoomStep, offsetX, offsetY);
      } else if (event.key === "0") {
        event.preventDefault();
        resetZoom();
      }
    });

    canvas.addEventListener("dblclick", resetZoom);
    controls.addEventListener("keydown", function (event) {
      if (event.target.matches("input,select,button,a")) {
        return;
      }
      if (event.key === "ArrowLeft") {
        move(-1);
      } else if (event.key === "ArrowRight") {
        move(1);
      }
    });
    window.addEventListener("resize", scheduleGeometry);
    document.addEventListener("shown.bs.tab", scheduleGeometry);
    document.addEventListener(
      "dnaepico:figure-zoom-change",
      function (event) {
        zoom = validZoom(event.detail && event.detail.zoom) || 1;
        scheduleGeometry();
      }
    );

    if (window.ResizeObserver) {
      var canvasObserver = new ResizeObserver(scheduleGeometry);
      canvasObserver.observe(canvas);
    }

    if (window.MutationObserver) {
      var expansionObserver = new MutationObserver(scheduleGeometry);
      expansionObserver.observe(card, {
        attributes: true,
        attributeFilter: ["data-full-screen"]
      });
    }

    showFigure();
  }

  function initializeAll() {
    document.querySelectorAll(
      ".dnaepico-figure-content[data-browser-id]"
    ).forEach(initializeViewer);
  }

  initializeAll();
  if (document.readyState === "loading") {
    document.addEventListener("DOMContentLoaded", initializeAll, {
      once: true
    });
  }
}());

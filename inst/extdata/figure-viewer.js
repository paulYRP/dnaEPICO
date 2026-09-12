(function () {
  "use strict";

  var minimumZoom = 0.25;
  var maximumZoom = 4;
  var zoomStep = 1.15;

  function clamp(value, minimum, maximum) {
    return Math.min(maximum, Math.max(minimum, value));
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

  function enableImageDownload(link, format) {
    if (!/\.(png|jpe?g)$/i.test(format.href)) return;
    var pending = false;
    link.addEventListener("click", function (event) {
      event.preventDefault();
      event.stopPropagation();
      if (pending) return;
      pending = true;
      link.setAttribute("aria-busy", "true");
      link.setAttribute("aria-disabled", "true");
      link.textContent = "Preparing " + format.label + "...";
      var script = document.createElement("script");
      var timer = window.setTimeout(failed, 30000);

      function finish() {
        window.clearTimeout(timer);
        script.onload = script.onerror = null;
        script.remove();
        pending = false;
        link.removeAttribute("aria-busy");
        link.removeAttribute("aria-disabled");
      }

      function failed() {
        finish();
        link.textContent = "Retry " + format.label + " download";
        link.title = "The image download could not be prepared. Please try again.";
      }

      script.onload = function () {
        try {
          var binary = window.atob(script.dataset.downloadData);
          var bytes = new Uint8Array(binary.length);
          for (var index = 0; index < binary.length; index++) {
            bytes[index] = binary.charCodeAt(index);
          }
          var blob = new Blob([bytes], { type: script.dataset.downloadType });
          var url = URL.createObjectURL(blob);
          var save = document.createElement("a");
          save.href = url;
          save.download = format.name;
          save.hidden = true;
          document.body.appendChild(save);
          save.click();
          save.remove();
          window.setTimeout(function () { URL.revokeObjectURL(url); }, 60000);
          finish();
          link.textContent = "Download " + format.label;
          link.removeAttribute("title");
        } catch (error) {
          failed();
        }
      };
      script.onerror = failed;
      script.src = format.href + ".download.js";
      document.head.appendChild(script);
    });
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
    var zoom = 1;

    function updateZoomStatus() {
      var percentage = Math.round(zoom * 100);
      var message = percentage === 100 ?
        "Image fitted to the available panel." :
        "Image zoom is " + percentage + "% of the fitted size.";
      if (zoomStatus) {
        zoomStatus.textContent = message;
      }
    }

    function applyZoom() {
      // Relative sizes cannot feed rounded client dimensions back into layout.
      canvas.style.setProperty("--dnaepico-stage-size", Math.max(1, zoom) * 100 + "%");
      canvas.style.setProperty("--dnaepico-image-size", Math.min(1, zoom) * 100 + "%");
      updateZoomStatus();
    }

    function updateImageAspect() {
      if (image.naturalWidth && image.naturalHeight && !image.hidden) {
        canvas.style.setProperty("--dnaepico-figure-aspect-ratio",
          image.naturalWidth + " / " + image.naturalHeight);
      }
    }

    function changeZoom(nextZoom, offsetX, offsetY) {
      var oldWidth = Math.max(1, canvas.scrollWidth);
      var oldHeight = Math.max(1, canvas.scrollHeight);
      var anchorX = (canvas.scrollLeft + offsetX) / oldWidth;
      var anchorY = (canvas.scrollTop + offsetY) / oldHeight;

      zoom = clamp(nextZoom, minimumZoom, maximumZoom);
      applyZoom();
      window.requestAnimationFrame(function () {
        canvas.scrollLeft = anchorX * canvas.scrollWidth - offsetX;
        canvas.scrollTop = anchorY * canvas.scrollHeight - offsetY;
      });
    }

    function resetZoom() {
      zoom = 1;
      applyZoom();
      canvas.scrollTo(0, 0);
    }

    function showFigure() {
      var index = clamp(Number(select.value) - 1, 0, figures.length - 1);
      var item = figures[index];
      title.textContent = item.title;
      download.replaceChildren();
      var formats = item.downloads && item.downloads.length ? item.downloads : [{
        href: item.downloadPath, name: item.downloadName, label: "Original"
      }];
      formats.forEach(function (format) {
        var link = document.createElement("a");
        link.className = "btn btn-sm btn-outline-primary dnaepico-download";
        link.href = format.href;
        link.download = format.name;
        link.textContent = "Download " + format.label;
        link.setAttribute("aria-label", "Download " + item.title + " as " + format.label);
        enableImageDownload(link, format);
        download.appendChild(link);
      });
      if (description) {
        description.textContent = item.description || "";
      }
      count.textContent = (index + 1) + " of " + figures.length + " figures";
      previous.disabled = index === 0;
      next.disabled = index === figures.length - 1;
      zoom = 1;
      applyZoom();

      if (item.browserReady) {
        fallback.hidden = true;
        image.hidden = false;
        image.alt = item.title;
        image.src = item.previewPath;
        if (image.complete && image.naturalWidth) {
          updateImageAspect();
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
    image.addEventListener("load", updateImageAspect);
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
        (event.clientX - bounds.left) * canvas.offsetWidth / bounds.width,
        (event.clientY - bounds.top) * canvas.offsetHeight / bounds.height
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

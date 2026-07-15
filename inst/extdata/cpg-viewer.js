(function () {
  "use strict";

  window.dnaEPICOResultChunks = window.dnaEPICOResultChunks || {};
  window.dnaEPICOResultManifests = window.dnaEPICOResultManifests || {};

  function formatCount(value) {
    return Number(value).toLocaleString("en-AU");
  }

  function chunkFile(basePath, number) {
    return basePath + "/chunk-" + String(number).padStart(4, "0") + ".js";
  }

  function chunkKey(manifest, number) {
    return manifest.key + ":" + number;
  }

  function loadChunk(manifest, number) {
    var cacheKey = chunkKey(manifest, number);
    if (window.dnaEPICOResultChunks[cacheKey]) {
      return Promise.resolve(window.dnaEPICOResultChunks[cacheKey].rows);
    }

    return new Promise(function (resolve, reject) {
      var script = document.createElement("script");
      script.src = new URL(chunkFile(manifest.basePath, number), document.baseURI).href;
      script.async = true;
      script.onload = function () {
        script.remove();
        if (window.dnaEPICOResultChunks[cacheKey]) {
          resolve(window.dnaEPICOResultChunks[cacheKey].rows);
        } else {
          reject(new Error("The requested result chunk did not register correctly."));
        }
      };
      script.onerror = function () {
        script.remove();
        reject(new Error("Unable to load result chunk " + number + "."));
      };
      document.head.appendChild(script);
    });
  }

  function releaseChunk(manifest, number) {
    delete window.dnaEPICOResultChunks[chunkKey(manifest, number)];
  }

  async function loadRange(manifest, start, end) {
    var rows = [];
    var position = start;
    while (position < end) {
      var chunkNumber = Math.floor(position / manifest.chunkSize) + 1;
      var chunkStart = (chunkNumber - 1) * manifest.chunkSize;
      var chunkRows = await loadChunk(manifest, chunkNumber);
      var localStart = position - chunkStart;
      var localEnd = Math.min(chunkRows.length, end - chunkStart);
      rows = rows.concat(chunkRows.slice(localStart, localEnd));
      position = chunkStart + localEnd;
      if (localEnd >= chunkRows.length && position < end) {
        position = chunkNumber * manifest.chunkSize;
      }
    }
    return rows;
  }

  async function loadIndexedRows(manifest, indices) {
    var rows = new Array(indices.length);
    var grouped = {};
    indices.forEach(function (globalIndex, outputIndex) {
      var number = Math.floor(globalIndex / manifest.chunkSize) + 1;
      grouped[number] = grouped[number] || [];
      grouped[number].push({
        outputIndex: outputIndex,
        localIndex: globalIndex - ((number - 1) * manifest.chunkSize)
      });
    });

    var chunkNumbers = Object.keys(grouped);
    for (var i = 0; i < chunkNumbers.length; i += 1) {
      var chunkNumber = Number(chunkNumbers[i]);
      var chunkRows = await loadChunk(manifest, chunkNumber);
      grouped[chunkNumber].forEach(function (position) {
        rows[position.outputIndex] = chunkRows[position.localIndex];
      });
    }
    return rows;
  }

  function makeDownloadLink(href, text) {
    var link = document.createElement("a");
    link.href = href;
    link.textContent = text;
    link.className = "btn btn-outline-primary btn-sm";
    link.setAttribute("download", "");
    return link;
  }

  function matchesFilter(value, operator, query) {
    var cell = value == null ? "" : String(value);
    if (operator === "contains") {
      return cell.toLocaleLowerCase().includes(query.toLocaleLowerCase());
    }
    if (operator === "equals") {
      return cell.toLocaleLowerCase() === query.toLocaleLowerCase();
    }

    if (!cell.trim()) {
      return false;
    }
    var numericCell = Number(cell);
    var numericQuery = Number(query);
    if (!Number.isFinite(numericCell) || !Number.isFinite(numericQuery)) {
      return false;
    }
    if (operator === "lt") return numericCell < numericQuery;
    if (operator === "lte") return numericCell <= numericQuery;
    if (operator === "gt") return numericCell > numericQuery;
    if (operator === "gte") return numericCell >= numericQuery;
    return false;
  }

  function initializeViewer(viewer) {
    if (viewer.dataset.initialized === "true") {
      return;
    }
    viewer.dataset.initialized = "true";

    var key = viewer.dataset.resultKey;
    var manifest = window.dnaEPICOResultManifests[key];
    var status = viewer.querySelector('[data-role="status"]');
    if (!manifest) {
      status.textContent = "Result metadata could not be loaded.";
      status.classList.add("dnaepico-viewer-error");
      return;
    }

    var head = viewer.querySelector('[data-role="head"]');
    var body = viewer.querySelector('[data-role="body"]');
    var pageSizeControl = viewer.querySelector('[data-role="page-size"]');
    var pageNumberControls = Array.from(viewer.querySelectorAll('[data-role="page-number"]'));
    var pageCounts = Array.from(viewer.querySelectorAll('[data-role="page-count"]'));
    var identifierSearch = viewer.querySelector('[data-role="cpg-search"]');
    var filterColumn = viewer.querySelector('[data-role="filter-column"]');
    var filterOperator = viewer.querySelector('[data-role="filter-operator"]');
    var filterValue = viewer.querySelector('[data-role="filter-value"]');
    var filterSummary = viewer.querySelector('[data-role="filter-summary"]');
    var buttons = {
      first: Array.from(viewer.querySelectorAll('[data-role="first"]')),
      previous: Array.from(viewer.querySelectorAll('[data-role="previous"]')),
      next: Array.from(viewer.querySelectorAll('[data-role="next"]')),
      last: Array.from(viewer.querySelectorAll('[data-role="last"]')),
      find: viewer.querySelector('[data-role="find-cpg"]'),
      applyFilter: viewer.querySelector('[data-role="apply-filter"]'),
      clearFilter: viewer.querySelector('[data-role="clear-filter"]')
    };
    var itemSingular = manifest.itemSingular || "row";
    var itemPlural = manifest.itemPlural || "rows";
    var currentPage = 1;
    var renderToken = 0;
    var filterToken = 0;
    var filteredIndices = null;

    var headerRow = document.createElement("tr");
    manifest.columns.forEach(function (column, index) {
      var th = document.createElement("th");
      th.scope = "col";
      th.textContent = column;
      headerRow.appendChild(th);

      var option = document.createElement("option");
      option.value = String(index);
      option.textContent = column;
      filterColumn.appendChild(option);
    });
    head.replaceChildren(headerRow);

    var pValueColumn = manifest.columns.findIndex(function (column) {
      return /p[._ -]?value/i.test(column);
    });
    if (pValueColumn >= 0) {
      filterColumn.value = String(pValueColumn);
      filterOperator.value = "lte";
    }

    var downloads = viewer.querySelector('[data-role="downloads"]');
    (manifest.downloads || []).forEach(function (download) {
      downloads.appendChild(makeDownloadLink(download.href, download.label));
    });

    function pageSize() {
      return Number(pageSizeControl.value) || 25;
    }

    function visibleRowCount() {
      return filteredIndices === null ? manifest.totalRows : filteredIndices.length;
    }

    function totalPages() {
      return Math.max(1, Math.ceil(visibleRowCount() / pageSize()));
    }

    function setFilterBusy(busy) {
      buttons.applyFilter.disabled = busy;
      filterColumn.disabled = busy;
      filterOperator.disabled = busy;
      filterValue.disabled = busy;
      buttons.clearFilter.disabled = busy || filteredIndices === null;
    }

    function updateControls() {
      var count = totalPages();
      pageNumberControls.forEach(function (control) {
        control.max = String(count);
        control.value = String(currentPage);
      });
      pageCounts.forEach(function (countElement) {
        countElement.textContent = "of " + formatCount(count);
      });
      buttons.first.concat(buttons.previous).forEach(function (button) {
        button.disabled = currentPage <= 1;
      });
      buttons.next.concat(buttons.last).forEach(function (button) {
        button.disabled = currentPage >= count;
      });
      buttons.clearFilter.disabled = filteredIndices === null;
    }

    async function loadVisibleRange(start, end) {
      if (filteredIndices === null) {
        return loadRange(manifest, start, end);
      }
      return loadIndexedRows(manifest, filteredIndices.slice(start, end));
    }

    async function renderPage(requestedPage) {
      var count = totalPages();
      currentPage = Math.min(Math.max(1, Number(requestedPage) || 1), count);
      updateControls();

      var rowCount = visibleRowCount();
      var start = (currentPage - 1) * pageSize();
      var end = Math.min(start + pageSize(), rowCount);
      var token = ++renderToken;
      status.textContent = "Loading " + itemPlural + "…";
      status.classList.remove("dnaepico-viewer-error");

      try {
        var rows = rowCount ? await loadVisibleRange(start, end) : [];
        if (token !== renderToken) {
          return;
        }
        var fragment = document.createDocumentFragment();
        rows.forEach(function (row) {
          var tr = document.createElement("tr");
          row.forEach(function (value) {
            var td = document.createElement("td");
            td.textContent = value == null ? "" : String(value);
            tr.appendChild(td);
          });
          fragment.appendChild(tr);
        });
        body.replaceChildren(fragment);
        if (rowCount) {
          var label = rowCount === 1 ? itemSingular : itemPlural;
          status.textContent = "Showing " + formatCount(start + 1) + "–" +
            formatCount(end) + " of " + formatCount(rowCount) + " " + label +
            (filteredIndices === null ? "." : " after filtering (" +
              formatCount(manifest.totalRows) + " total)." );
        } else if (filteredIndices === null) {
          status.textContent = "No " + itemPlural + " are available in this table.";
        } else {
          status.textContent = "No " + itemPlural + " match the current filter.";
        }
      } catch (error) {
        if (token !== renderToken) {
          return;
        }
        body.replaceChildren();
        status.textContent = error.message;
        status.classList.add("dnaepico-viewer-error");
      }
    }

    async function applyFilter() {
      var columnIndex = Number(filterColumn.value);
      var query = filterValue.value.trim();
      if (!Number.isInteger(columnIndex) || columnIndex < 0) {
        status.textContent = "Choose a column to filter.";
        return;
      }
      if (!query) {
        status.textContent = "Enter a filter value.";
        return;
      }

      var operator = filterOperator.value;
      var token = ++filterToken;
      var matches = [];
      setFilterBusy(true);
      status.classList.remove("dnaepico-viewer-error");

      try {
        for (var i = 0; i < manifest.chunks.length; i += 1) {
          if (token !== filterToken) return;
          var descriptor = manifest.chunks[i];
          var cacheWasPresent = Boolean(window.dnaEPICOResultChunks[chunkKey(manifest, descriptor.number)]);
          status.textContent = "Filtering " + itemPlural + ": chunk " +
            formatCount(i + 1) + " of " + formatCount(manifest.chunks.length) + "…";
          var rows = await loadChunk(manifest, descriptor.number);
          rows.forEach(function (row, localIndex) {
            if (matchesFilter(row[columnIndex], operator, query)) {
              matches.push(((descriptor.number - 1) * manifest.chunkSize) + localIndex);
            }
          });
          if (!cacheWasPresent) {
            releaseChunk(manifest, descriptor.number);
          }
        }
        if (token !== filterToken) return;
        filteredIndices = matches;
        filterSummary.textContent = formatCount(matches.length) + " matches";
        currentPage = 1;
        await renderPage(1);
      } catch (error) {
        if (token !== filterToken) return;
        status.textContent = error.message;
        status.classList.add("dnaepico-viewer-error");
      } finally {
        if (token === filterToken) setFilterBusy(false);
      }
    }

    function clearFilter() {
      filterToken += 1;
      filteredIndices = null;
      filterValue.value = "";
      filterSummary.textContent = "";
      setFilterBusy(false);
      renderPage(1);
    }

    async function findIdentifier() {
      var query = identifierSearch.value.trim();
      if (!query) {
        status.textContent = "Enter an identifier to find.";
        return;
      }
      var idIndex = manifest.columns.indexOf(manifest.idColumn);
      if (idIndex < 0 || !manifest.chunks.length) {
        status.textContent = "Identifier lookup is unavailable for this table.";
        return;
      }

      var low = 0;
      var high = manifest.chunks.length - 1;
      var candidate = manifest.chunks.length;
      while (low <= high) {
        var middle = Math.floor((low + high) / 2);
        if (manifest.chunks[middle].lastId.localeCompare(query) >= 0) {
          candidate = middle;
          high = middle - 1;
        } else {
          low = middle + 1;
        }
      }

      var candidates = [candidate, candidate + 1].filter(function (index) {
        return index >= 0 && index < manifest.chunks.length;
      });
      for (var i = 0; i < candidates.length; i += 1) {
        var descriptor = manifest.chunks[candidates[i]];
        var rows = await loadChunk(manifest, descriptor.number);
        var localIndex = rows.findIndex(function (row) {
          var value = row[idIndex] == null ? "" : String(row[idIndex]);
          return value === query || value.startsWith(query);
        });
        if (localIndex >= 0) {
          var globalIndex = (descriptor.number - 1) * manifest.chunkSize + localIndex;
          var visibleIndex = filteredIndices === null ? globalIndex : filteredIndices.indexOf(globalIndex);
          identifierSearch.value = String(rows[localIndex][idIndex]);
          if (visibleIndex < 0) {
            status.textContent = "The identifier exists but is excluded by the current filter.";
            return;
          }
          await renderPage(Math.floor(visibleIndex / pageSize()) + 1);
          return;
        }
      }
      status.textContent = "Identifier not found: " + query;
    }

    buttons.first.forEach(function (button) {
      button.addEventListener("click", function () { renderPage(1); });
    });
    buttons.previous.forEach(function (button) {
      button.addEventListener("click", function () { renderPage(currentPage - 1); });
    });
    buttons.next.forEach(function (button) {
      button.addEventListener("click", function () { renderPage(currentPage + 1); });
    });
    buttons.last.forEach(function (button) {
      button.addEventListener("click", function () { renderPage(totalPages()); });
    });
    buttons.find.addEventListener("click", function () {
      findIdentifier().catch(function (error) {
        status.textContent = error.message;
        status.classList.add("dnaepico-viewer-error");
      });
    });
    buttons.applyFilter.addEventListener("click", function () {
      applyFilter();
    });
    buttons.clearFilter.addEventListener("click", clearFilter);
    identifierSearch.addEventListener("keydown", function (event) {
      if (event.key === "Enter") {
        event.preventDefault();
        buttons.find.click();
      }
    });
    filterValue.addEventListener("keydown", function (event) {
      if (event.key === "Enter") {
        event.preventDefault();
        buttons.applyFilter.click();
      }
    });
    pageNumberControls.forEach(function (control) {
      control.addEventListener("change", function () {
        renderPage(control.value);
      });
    });
    pageSizeControl.addEventListener("change", function () {
      renderPage(1);
    });

    renderPage(1);
  }

  function initializeAll() {
    document.querySelectorAll(".dnaepico-cpg-viewer").forEach(initializeViewer);
  }

  if (document.readyState === "loading") {
    document.addEventListener("DOMContentLoaded", initializeAll, { once: true });
  } else {
    initializeAll();
  }
}());

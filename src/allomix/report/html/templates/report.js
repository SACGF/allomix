/*
 * Progressive enhancement only: the timeline timepoints table already arrives
 * sorted server-side, so a reader with JS off loses nothing. Clicking a header
 * sorts by that column, numeric-aware (strips %, commas, and the NA glyph),
 * toggling ascending / descending. Rows are reordered within tbody; any
 * excluded-row styling rides on the row class and is preserved.
 *
 * Inlined verbatim into the report's <script> block (not parsed as a Jinja
 * template), so it may contain any characters.
 */
(function () {
  "use strict";
  function cellValue(row, idx) {
    var cell = row.children[idx];
    return cell ? cell.textContent.trim() : "";
  }
  function asNumber(s) {
    var cleaned = s.replace(/[%,\s]/g, "").replace("—", "");
    if (cleaned === "" || cleaned === "<0.001") {
      return cleaned === "<0.001" ? 0.0005 : NaN;
    }
    var v = parseFloat(cleaned);
    return isNaN(v) ? NaN : v;
  }
  function columnIsNumeric(rows, idx) {
    var seen = 0;
    for (var i = 0; i < rows.length; i++) {
      var t = cellValue(rows[i], idx);
      if (t === "" || t === "—") continue;
      if (isNaN(asNumber(t))) return false;
      seen++;
    }
    return seen > 0;
  }
  function sortTable(table, idx, th) {
    var tbody = table.tBodies[0];
    if (!tbody) return;
    var rows = Array.prototype.slice.call(tbody.rows);
    var numeric = columnIsNumeric(rows, idx);
    var asc = !th.classList.contains("sort-asc");
    var headers = th.parentNode.children;
    for (var h = 0; h < headers.length; h++) {
      headers[h].classList.remove("sort-asc", "sort-desc");
    }
    th.classList.add(asc ? "sort-asc" : "sort-desc");
    rows.sort(function (a, b) {
      var av = cellValue(a, idx), bv = cellValue(b, idx);
      var cmp;
      if (numeric) {
        var an = asNumber(av), bn = asNumber(bv);
        if (isNaN(an)) an = Infinity;
        if (isNaN(bn)) bn = Infinity;
        cmp = an - bn;
      } else {
        cmp = av.localeCompare(bv);
      }
      return asc ? cmp : -cmp;
    });
    for (var r = 0; r < rows.length; r++) tbody.appendChild(rows[r]);
  }
  document.querySelectorAll("table.sortable").forEach(function (table) {
    var ths = table.tHead ? table.tHead.rows[0].children : [];
    for (var i = 0; i < ths.length; i++) {
      (function (idx, th) {
        th.addEventListener("click", function () { sortTable(table, idx, th); });
      })(i, ths[i]);
    }
  });
})();

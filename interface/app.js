(function () {
  "use strict";

  const state = {
    activeTab: "patients",
    cells: [
      { cell_id: "SRR8617653", bam_path: "/path/to/SRR8617653.sorted.bam", cell_type: "ctc" },
    ],
    callers: [
      { cell_id: "SRR8617653", caller: "haplotypecaller", vcf_path: "/path/to/SRR8617653.haplotypecaller.vcf.gz", sccaller_mode: "" },
      { cell_id: "SRR8617653", caller: "deepvariant", vcf_path: "/path/to/SRR8617653.deepvariant.vcf.gz", sccaller_mode: "" },
      { cell_id: "SRR8617653", caller: "sccaller", vcf_path: "/path/to/SRR8617653.sccaller.vcf", sccaller_mode: "external_bulk" },
    ],
  };

  const $ = (id) => document.getElementById(id);
  const fields = [
    "patientId", "refFasta", "outdir", "cellMetadataName", "bamListName", "germlineMode",
    "germlineVcf", "bulkBam", "leukocyteBam", "leukocyteVcf", "monovarScript",
    "sccallerScript", "sccallerHsnp",
  ];

  function esc(value) {
    return String(value == null ? "" : value).replace(/\t/g, " ").replace(/\r?\n/g, " ").trim();
  }

  function q(value) {
    const s = String(value || "");
    if (!s) return "\"\"";
    return `"${s.replace(/(["\\$`])/g, "\\$1")}"`;
  }

  function mode() {
    return document.querySelector("input[name='mode']:checked").value;
  }

  function getForm() {
    const out = {};
    for (const id of fields) out[id] = $(id).value.trim();
    return out;
  }

  function patientTsv() {
    const f = getForm();
    return [
      "patient_id\tref_fasta\tmonovar_bam_list\tcell_metadata\tgermline_mode\tgermline_vcf\tbulk_bam\tleukocyte_bam\tleukocyte_vcf",
      [
        f.patientId,
        f.refFasta,
        f.bamListName || "NA",
        f.cellMetadataName || "NA",
        f.germlineMode,
        f.germlineVcf || "NA",
        f.bulkBam || "NA",
        f.leukocyteBam || "NA",
        f.leukocyteVcf || "NA",
      ].map(esc).join("\t"),
    ].join("\n") + "\n";
  }

  function cellsTsv() {
    const f = getForm();
    return [
      "patient_id\tcell_id\tbam_path\tcell_type",
      ...state.cells.map((row) => [f.patientId, row.cell_id, row.bam_path, row.cell_type].map(esc).join("\t")),
    ].join("\n") + "\n";
  }

  function bamList() {
    return state.cells
      .filter((row) => row.cell_type.toLowerCase() !== "leukocyte")
      .map((row) => esc(row.bam_path))
      .filter(Boolean)
      .join("\n") + "\n";
  }

  function callerTsv() {
    const f = getForm();
    return [
      "patient_id\tcell_id\tcaller\tvcf_path\tsccaller_mode",
      ...state.callers.map((row) => [f.patientId, row.cell_id, row.caller, row.vcf_path, row.sccaller_mode].map(esc).join("\t")),
    ].join("\n") + "\n";
  }

  function commandText() {
    const f = getForm();
    const m = mode();
    const args = [
      "nextflow run main.nf",
      "  -profile conda",
      `  --patients ${q("configs/patients.tsv")}`,
      `  --outdir ${q(f.outdir || "results")}`,
    ];

    if (m === "monovar" || m === "full") {
      args.push("  --run_monovar true");
      args.push(`  --monovar_script ${q(f.monovarScript)}`);
    }
    if (m === "sccaller" || m === "full") {
      args.push("  --run_sccaller true");
      args.push(`  --sccaller_script ${q(f.sccallerScript)}`);
      args.push(`  --sccaller_hsnp_vcf ${q(f.sccallerHsnp)}`);
    }
    if (m === "compare" || m === "full") {
      args.push("  --run_external_callers true");
      args.push(`  --caller_vcfs ${q("configs/caller_vcfs.tsv")}`);
    }
    if (m === "compare" || m === "sccaller" || m === "full") {
      args.push("  --run_comparison true");
    }
    return args.join(" \\\n") + "\n";
  }

  function checks() {
    const f = getForm();
    const warnings = [];
    const m = mode();
    if (!f.patientId) warnings.push("Missing patient ID.");
    if (!f.refFasta) warnings.push("Missing reference FASTA.");
    if (!f.germlineVcf && ["precomputed_vcf", "combined"].includes(f.germlineMode)) warnings.push("precomputed_vcf/combined modes need a germline VCF.");
    if ((m === "monovar" || m === "full") && !f.monovarScript) warnings.push("MonoVar mode needs --monovar_script.");
    if ((m === "sccaller" || m === "full") && !f.sccallerScript) warnings.push("SCcaller mode needs --sccaller_script.");
    if ((m === "sccaller" || m === "full") && !f.sccallerHsnp) warnings.push("SCcaller mode needs --sccaller_hsnp_vcf.");
    if ((m === "sccaller" || m === "full") && !f.bulkBam) warnings.push("SCcaller mode needs bulk_bam in patients.tsv.");
    if ((m === "compare" || m === "full") && state.callers.length === 0) warnings.push("External comparison mode needs at least one caller VCF row.");
    if (state.cells.length === 0 && m !== "compare") warnings.push("MonoVar or SCcaller modes need at least one cell row.");
    const cellIds = new Set(state.cells.map((r) => r.cell_id).filter(Boolean));
    for (const row of state.callers) {
      if (!row.cell_id || !row.caller || !row.vcf_path) warnings.push("Caller rows need cell ID, caller, and VCF path.");
      if (row.cell_id && cellIds.size && !cellIds.has(row.cell_id)) warnings.push(`Caller row cell '${row.cell_id}' is not present in the cell table.`);
    }
    return warnings.length ? warnings.map((w) => `WARN\t${w}`).join("\n") + "\n" : "OK\tNo obvious input problems found.\n";
  }

  function outputFor(tab) {
    return {
      patients: patientTsv,
      cells: cellsTsv,
      bamlist: bamList,
      callers: callerTsv,
      command: commandText,
      warnings: checks,
    }[tab]();
  }

  function renderRows() {
    $("cellsBody").innerHTML = state.cells.map((row, i) => `
      <tr>
        <td><input data-kind="cells" data-index="${i}" data-field="cell_id" value="${row.cell_id}"></td>
        <td><input data-kind="cells" data-index="${i}" data-field="bam_path" value="${row.bam_path}"></td>
        <td><input data-kind="cells" data-index="${i}" data-field="cell_type" value="${row.cell_type}"></td>
        <td><button type="button" data-remove="cells" data-index="${i}">X</button></td>
      </tr>`).join("");

    $("callersBody").innerHTML = state.callers.map((row, i) => `
      <tr>
        <td><input data-kind="callers" data-index="${i}" data-field="cell_id" value="${row.cell_id}"></td>
        <td>
          <select data-kind="callers" data-index="${i}" data-field="caller">
            ${["monovar", "sccaller", "haplotypecaller", "deepvariant"].map((c) => `<option value="${c}" ${row.caller === c ? "selected" : ""}>${c}</option>`).join("")}
          </select>
        </td>
        <td><input data-kind="callers" data-index="${i}" data-field="vcf_path" value="${row.vcf_path}"></td>
        <td>
          <select data-kind="callers" data-index="${i}" data-field="sccaller_mode">
            ${["", "so_true", "external_bulk", "so_true_and_external_bulk", "no_bulk"].map((c) => `<option value="${c}" ${row.sccaller_mode === c ? "selected" : ""}>${c || "default"}</option>`).join("")}
          </select>
        </td>
        <td><button type="button" data-remove="callers" data-index="${i}">X</button></td>
      </tr>`).join("");
  }

  function renderMode() {
    const m = mode();
    $("callerPanel").classList.toggle("hidden", !(m === "compare" || m === "full"));
  }

  function renderOutput() {
    $("outputBox").textContent = outputFor(state.activeTab);
    const warningCount = checks().split("\n").filter((line) => line.startsWith("WARN")).length;
    $("status").textContent = warningCount ? `${warningCount} warning${warningCount === 1 ? "" : "s"}` : "Ready";
    $("status").style.color = warningCount ? "var(--danger)" : "var(--ok)";
  }

  function render() {
    renderRows();
    renderMode();
    renderOutput();
  }

  function download(name, text) {
    const blob = new Blob([text], { type: "text/plain;charset=utf-8" });
    const url = URL.createObjectURL(blob);
    const a = document.createElement("a");
    a.href = url;
    a.download = name;
    a.click();
    URL.revokeObjectURL(url);
  }

  document.addEventListener("input", (event) => {
    const target = event.target;
    if (target.dataset.kind) {
      state[target.dataset.kind][Number(target.dataset.index)][target.dataset.field] = target.value;
    }
    renderOutput();
  });

  document.addEventListener("change", (event) => {
    const target = event.target;
    if (target.name === "mode") render();
    if (target.dataset.kind) {
      state[target.dataset.kind][Number(target.dataset.index)][target.dataset.field] = target.value;
      renderOutput();
    }
  });

  document.addEventListener("click", async (event) => {
    const target = event.target;
    if (target.dataset.tab) {
      state.activeTab = target.dataset.tab;
      document.querySelectorAll(".tabs button").forEach((btn) => btn.classList.toggle("active", btn.dataset.tab === state.activeTab));
      renderOutput();
    }
    if (target.dataset.remove) {
      state[target.dataset.remove].splice(Number(target.dataset.index), 1);
      render();
    }
    if (target.id === "addCell") {
      state.cells.push({ cell_id: "", bam_path: "", cell_type: "ctc" });
      render();
    }
    if (target.id === "addCaller") {
      state.callers.push({ cell_id: "", caller: "haplotypecaller", vcf_path: "", sccaller_mode: "" });
      render();
    }
    if (target.id === "copyCommand") {
      await navigator.clipboard.writeText(commandText());
      $("status").textContent = "Copied";
      $("status").style.color = "var(--ok)";
    }
    if (target.id === "downloadAll") {
      download("patients.tsv", patientTsv());
      download("cells.tsv", cellsTsv());
      download("monovar_bams.txt", bamList());
      download("caller_vcfs.tsv", callerTsv());
      download("cnpup_command.sh", commandText());
    }
  });

  render();
})();

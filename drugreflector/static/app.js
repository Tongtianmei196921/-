const healthSummary = document.getElementById("health-summary");
const healthModels = document.getElementById("health-models");
const healthMode = document.getElementById("health-mode");
const refreshHealthButton = document.getElementById("refresh-health");
const resultsContainer = document.getElementById("results-container");
const messageBox = document.getElementById("message");
const resultCaption = document.getElementById("result-caption");
const downloadCsvButton = document.getElementById("download-csv");
const downloadJsonButton = document.getElementById("download-json");
const geneCountHint = document.getElementById("gene-count-hint");
const vscoresTextArea = document.getElementById("vscores-text");

const tabs = Array.from(document.querySelectorAll(".tab"));
const tabPanels = Array.from(document.querySelectorAll(".tab-panel"));
let latestPayload = null;
let latestDownloadStem = "drugreflector_results";

function setMessage(text, kind = "") {
  messageBox.textContent = text;
  messageBox.className = "message";
  if (kind) {
    messageBox.classList.add(`is-${kind}`);
  }
}

function parseSignatureText(text) {
  const scores = {};
  const lines = text
    .split(/\r?\n/)
    .map((line) => line.trim())
    .filter(Boolean);

  for (const line of lines) {
    const commaParts = line.split(",");
    let gene;
    let scoreText;

    if (commaParts.length >= 2) {
      gene = commaParts[0].trim();
      scoreText = commaParts.slice(1).join(",").trim();
    } else {
      const parts = line.split(/\s+/);
      if (parts.length < 2) {
        throw new Error(`Unable to parse line: "${line}"`);
      }
      gene = parts[0].trim();
      scoreText = parts.slice(1).join(" ").trim();
    }

    const score = Number(scoreText);
    if (!gene || Number.isNaN(score)) {
      throw new Error(`Invalid gene score line: "${line}"`);
    }

    scores[gene] = score;
  }

  if (!Object.keys(scores).length) {
    throw new Error("Please provide at least one gene score.");
  }

  return scores;
}

function countParsedGenes(text) {
  return text
    .split(/\r?\n/)
    .map((line) => line.trim())
    .filter(Boolean).length;
}

function updateGeneCountHint() {
  const count = countParsedGenes(vscoresTextArea.value);
  geneCountHint.textContent = `${count} genes detected in the current text box.`;
}

function setDownloadState(enabled) {
  downloadCsvButton.disabled = !enabled;
  downloadJsonButton.disabled = !enabled;
}

function triggerDownload(filename, content, mimeType) {
  const blob = new Blob([content], { type: mimeType });
  const url = URL.createObjectURL(blob);
  const link = document.createElement("a");
  link.href = url;
  link.download = filename;
  document.body.appendChild(link);
  link.click();
  link.remove();
  URL.revokeObjectURL(url);
}

function flattenResults(results) {
  return Object.entries(results).flatMap(([signature, rows]) =>
    rows.map((row) => ({
      signature,
      compound: row.compound,
      rank: row.rank,
      logit: row.logit,
      prob: row.prob,
    })),
  );
}

function downloadCsv() {
  if (!latestPayload) {
    return;
  }

  const rows = flattenResults(latestPayload.results);
  const header = ["signature", "compound", "rank", "logit", "prob"];
  const csvLines = [
    header.join(","),
    ...rows.map((row) =>
      [row.signature, row.compound, row.rank, row.logit, row.prob]
        .map((value) => `"${String(value).replaceAll('"', '""')}"`)
        .join(","),
    ),
  ];

  triggerDownload(`${latestDownloadStem}.csv`, csvLines.join("\n"), "text/csv;charset=utf-8");
}

function downloadJson() {
  if (!latestPayload) {
    return;
  }

  triggerDownload(
    `${latestDownloadStem}.json`,
    JSON.stringify(latestPayload, null, 2),
    "application/json;charset=utf-8",
  );
}

function renderResults(results, caption) {
  const names = Object.keys(results);
  if (!names.length) {
    resultsContainer.innerHTML = `
      <div class="empty-state">
        <p>No results returned by the model.</p>
      </div>
    `;
    resultCaption.textContent = caption;
    setDownloadState(false);
    return;
  }

  const blocks = names
    .map((name) => {
      const topHit = results[name][0];
      const rows = results[name]
        .map(
          (row) => `
            <tr>
              <td><span class="rank-pill">${row.rank}</span></td>
              <td>${row.compound}</td>
              <td>${row.logit.toFixed(4)}</td>
              <td>${row.prob.toFixed(6)}</td>
            </tr>
          `,
        )
        .join("");

      return `
        <section class="signature-block">
          <div class="signature-header">
            <div>
              <h3>${name}</h3>
              <div class="signature-meta">${results[name].length} ranked compounds shown</div>
            </div>
            ${topHit ? `<div class="top-hit">Top hit: ${topHit.compound}</div>` : ""}
          </div>
          <table>
            <thead>
              <tr>
                <th>Rank</th>
                <th>Compound</th>
                <th>Logit</th>
                <th>Probability</th>
              </tr>
            </thead>
            <tbody>${rows}</tbody>
          </table>
        </section>
      `;
    })
    .join("");

  resultCaption.textContent = caption;
  resultsContainer.innerHTML = blocks;
  setDownloadState(true);
}

async function refreshHealth() {
  healthSummary.textContent = "Checking model status...";
  healthModels.textContent = "Loading...";
  healthMode.textContent = "Checking...";

  try {
    const response = await fetch("/health");
    const payload = await response.json();
    if (payload.status === "ok") {
      healthSummary.textContent = "Online and ready for inference";
      healthModels.textContent = `${payload.n_compounds} compounds in library`;
      healthMode.textContent = "Research ranking from transcriptional signatures";
    } else {
      healthSummary.textContent = "Service started, but model loading failed";
      healthModels.textContent = payload.detail || "Check checkpoint configuration";
      healthMode.textContent = "Model unavailable";
    }
  } catch (error) {
    healthSummary.textContent = "Unable to reach the local API";
    healthModels.textContent = "Make sure the server is running on this machine";
    healthMode.textContent = "Offline";
  }
}

async function submitVscores(event) {
  event.preventDefault();
  setMessage("Running v-score prediction...");

  try {
    const name = document.getElementById("signature-name").value.trim() || "signature";
    const nTop = Number(document.getElementById("vscores-topn").value || 20);
    const text = document.getElementById("vscores-text").value;
    const scores = parseSignatureText(text);

    const response = await fetch("/predict/vscores", {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify({
        n_top: nTop,
        signatures: [{ name, scores }],
      }),
    });

    if (!response.ok) {
      const payload = await response.json();
      throw new Error(payload.detail || "Prediction request failed.");
    }

    const payload = await response.json();
    latestPayload = payload;
    latestDownloadStem = `${name}_top_${payload.n_top}`;
    renderResults(payload.results, `Showing top ${payload.n_top} compounds from a v-score signature.`);
    setMessage("Prediction complete.", "success");
  } catch (error) {
    setMessage(error.message, "error");
  }
}

async function submitH5ad(event) {
  event.preventDefault();
  setMessage("Uploading H5AD file and running prediction...");

  try {
    const fileInput = document.getElementById("h5ad-file");
    if (!fileInput.files.length) {
      throw new Error("Please choose an H5AD file first.");
    }

    const nTop = Number(document.getElementById("h5ad-topn").value || 20);
    const formData = new FormData();
    formData.append("file", fileInput.files[0]);
    formData.append("n_top", String(nTop));

    const response = await fetch("/predict/h5ad", {
      method: "POST",
      body: formData,
    });

    if (!response.ok) {
      const payload = await response.json();
      throw new Error(payload.detail || "H5AD prediction request failed.");
    }

    const payload = await response.json();
    latestPayload = payload;
    latestDownloadStem = `${payload.input_file.replace(/\.[^.]+$/, "")}_top_${payload.n_top}`;
    renderResults(
      payload.results,
      `Showing top ${payload.n_top} compounds from ${payload.input_file} with ${payload.n_obs} observations and ${payload.n_vars} genes.`,
    );
    setMessage("Prediction complete.", "success");
  } catch (error) {
    setMessage(error.message, "error");
  }
}

function activateTab(targetId) {
  tabs.forEach((tab) => {
    tab.classList.toggle("is-active", tab.dataset.target === targetId);
  });
  tabPanels.forEach((panel) => {
    panel.classList.toggle("is-active", panel.id === targetId);
  });
}

function loadExample() {
  document.getElementById("signature-name").value = "mono_to_fcgr3a_demo";
  vscoresTextArea.value = `TP53,1.2
EGFR,-0.4
CDKN1A,0.8
JUN,0.5
FOS,0.3
STAT1,1.1
IRF1,0.9
CXCL8,-0.7`;
  updateGeneCountHint();
}

tabs.forEach((tab) => {
  tab.addEventListener("click", () => activateTab(tab.dataset.target));
});

document.getElementById("vscores-form").addEventListener("submit", submitVscores);
document.getElementById("h5ad-form").addEventListener("submit", submitH5ad);
document.getElementById("load-example").addEventListener("click", loadExample);
refreshHealthButton.addEventListener("click", refreshHealth);
downloadCsvButton.addEventListener("click", downloadCsv);
downloadJsonButton.addEventListener("click", downloadJson);
vscoresTextArea.addEventListener("input", updateGeneCountHint);

setDownloadState(false);
updateGeneCountHint();
refreshHealth();

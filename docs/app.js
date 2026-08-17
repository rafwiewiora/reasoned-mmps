(() => {
  const data = window.REASONED_MMP_DATA;
  const moves = data.moves;
  const compounds = new Map(data.compounds.map((compound) => [compound.chembl_id, compound]));
  const state = { selected: moves[0]?.reasoned_move_id ?? null };

  const $ = (selector) => document.querySelector(selector);
  const escapeHtml = (value = "") => String(value)
    .replaceAll("&", "&amp;")
    .replaceAll("<", "&lt;")
    .replaceAll(">", "&gt;")
    .replaceAll('"', "&quot;")
    .replaceAll("'", "&#039;");
  const humanize = (value = "") => String(value).replaceAll("_", " ").replace(/\b\w/g, (char) => char.toUpperCase());
  const shortClass = (value) => value.includes("retrospective") ? "retrospective" : value.includes("explicit") ? "explicit" : "implied";
  const reasonOf = (move) => move.layer_2_extracted_design_intent;
  const outcomesOf = (move) => move.layer_4_observed_outcomes.comparisons;
  const factsOf = (move) => move.layer_4_observed_outcomes.unpaired_facts ?? [];
  const membersOf = (reason) => reason.member_chembl_ids ?? [reason.child_chembl_id];

  function primaryOutcome(move) {
    const outcomes = outcomesOf(move);
    const intentStatus = move.layer_4_observed_outcomes.stated_intent_outcome;
    if (intentStatus === "not_applicable_retrospective_explanation") return "retrospective";
    if (intentStatus === "indeterminate_direct_endpoint_unavailable") return "indeterminate";
    const direct = outcomes.filter((row) => row.relation_to_stated_intent === "direct_supporting_endpoint");
    if (direct.some((row) => row.classification === "improved")) return "improved";
    if (direct.some((row) => row.classification === "worsened")) return "worsened";
    if (direct.some((row) => row.classification === "comparable")) return "comparable";
    return "indeterminate";
  }

  function renderMetrics() {
    const counts = data.manifest.counts;
    $("#metrics").innerHTML = [
      [counts.papers, "papers"],
      [counts.rationale_episodes, "evidence episodes"],
      [counts.reason_bearing_compounds, "reason-bearing"],
      [counts.resolved_structures, "resolved structures"],
      [counts.unique_outcome_pairs, "unique assay pairs"],
      [counts.outcome_comparisons, "episode links"],
    ].map(([value, label]) => `<div class="metric"><strong>${value}</strong><span>${label}</span></div>`).join("");
  }

  function populateFilters() {
    const classes = [...new Set(moves.map((move) => reasonOf(move).assertion_class))].sort();
    $("#class-filter").insertAdjacentHTML("beforeend", classes.map((value) => `<option value="${escapeHtml(value)}">${escapeHtml(humanize(value))}</option>`).join(""));
  }

  function filteredMoves() {
    const query = $("#search").value.trim().toLowerCase();
    const classFilter = $("#class-filter").value;
    const outcomeFilter = $("#outcome-filter").value;
    return moves.filter((move) => {
      const reason = reasonOf(move);
      const haystack = [reason.source_title, reason.reason_statement, reason.child_paper_label, ...(reason.member_paper_labels ?? []), reason.intended_property.name, reason.evidence.document_id].join(" ").toLowerCase();
      return (!query || haystack.includes(query))
        && (classFilter === "all" || reason.assertion_class === classFilter)
        && (outcomeFilter === "all" || outcomesOf(move).some((row) => row.classification === outcomeFilter));
    });
  }

  function renderList() {
    const visible = filteredMoves();
    $("#result-count").textContent = `${visible.length} / ${moves.length}`;
    if (!visible.some((move) => move.reasoned_move_id === state.selected)) state.selected = visible[0]?.reasoned_move_id ?? null;
    $("#episode-list").innerHTML = visible.length ? visible.map((move) => {
      const reason = reasonOf(move);
      const outcome = primaryOutcome(move);
      const members = membersOf(reason);
      const title = reason.assertion_scope === "named_series" ? `${members.length}-compound series` : reason.child_paper_label;
      return `<button class="episode-card ${move.reasoned_move_id === state.selected ? "active" : ""}" data-id="${escapeHtml(move.reasoned_move_id)}">
        <span class="card-meta"><span>${escapeHtml(reason.evidence.document_id)}</span><span>${escapeHtml(humanize(reason.decision_kind ?? reason.assertion_scope ?? ""))}</span></span>
        <h3>${escapeHtml(title)} · ${escapeHtml(humanize(reason.intended_property.name))}</h3>
        <p>${escapeHtml(reason.reason_statement)}</p>
        <span class="badge-row">
          <span class="badge ${shortClass(reason.assertion_class)}">${escapeHtml(humanize(reason.assertion_class))}</span>
          ${members.length > 1 ? `<span class="badge">${members.length} named members</span>` : ""}
          <span class="badge ${outcome}">${escapeHtml(humanize(outcome))}</span>
        </span>
      </button>`;
    }).join("") : `<div class="empty">No rationale episodes match these filters.</div>`;
    document.querySelectorAll(".episode-card").forEach((card) => card.addEventListener("click", () => {
      state.selected = card.dataset.id;
      renderList();
      renderDetail();
    }));
  }

  function relationshipFor(move) {
    const structural = move.layer_3_inferred_structural_comparison;
    const explicit = structural.author_explicit_relationship;
    if (explicit) {
      return {
        parentId: explicit.parent_chembl_id,
        semantics: explicit.edge_semantics,
        historical: explicit.historical_synthesis_lineage_claim,
        witness: explicit.structural_witness,
        explicit: true,
      };
    }
    const top = structural.top_candidates[0];
    return top ? { parentId: top.parent_chembl_id, semantics: top.edge_semantics, historical: false, witness: top, explicit: false } : null;
  }

  function moleculeCard(compound, role) {
    if (!compound) return `<div class="molecule"><span class="label">${role}</span><div class="empty">Unresolved</div></div>`;
    return `<div class="molecule">
      <span class="label">${escapeHtml(role)}</span>
      <img src="molecules/${escapeHtml(compound.chembl_id)}.svg" alt="2D structure of ${escapeHtml(compound.paper_label)}">
      <div class="molecule-name"><strong>${escapeHtml(compound.paper_label)}</strong><code>${escapeHtml(compound.chembl_id)}</code></div>
    </div>`;
  }

  function formatMeasurement(relation, value, units) {
    return `${relation === "=" ? "" : escapeHtml(relation)}${escapeHtml(value)} ${escapeHtml(units)}`;
  }

  function deltaText(row) {
    if (row.delta_lower == null && row.delta_upper == null) return "bounded / unknown";
    const prefix = row.comparison_scale === "log10" ? "log10 " : "";
    if (row.delta_lower != null && row.delta_upper == null) return `${prefix}≥ ${row.delta_lower}`;
    if (row.delta_lower == null && row.delta_upper != null) return `${prefix}≤ ${row.delta_upper}`;
    return `${prefix}${row.delta_lower}`;
  }

  function renderDetail() {
    const move = moves.find((candidate) => candidate.reasoned_move_id === state.selected);
    if (!move) {
      $("#detail").innerHTML = `<div class="empty">Select a rationale episode.</div>`;
      return;
    }
    const reason = reasonOf(move);
    const child = compounds.get(reason.child_chembl_id);
    const relationship = relationshipFor(move);
    const parent = relationship ? compounds.get(relationship.parentId) : null;
    const witness = relationship?.witness;
    const transformation = witness?.transformation ?? "No strict MMP witness";
    const topCandidates = move.layer_3_inferred_structural_comparison.top_candidates;
    const outcomes = outcomesOf(move);
    const facts = factsOf(move);
    const memberLabels = reason.member_paper_labels ?? [reason.child_paper_label];
    const isSeries = reason.assertion_scope === "named_series";
    const relationshipLabel = relationship?.explicit ? humanize(relationship.semantics) : "Inferred analytic comparator";
    const additionalRelationships = move.layer_3_inferred_structural_comparison.additional_author_relationships ?? [];

    $("#detail").innerHTML = `
      <div class="detail-header">
        <div>
          <p class="detail-kicker">${escapeHtml(reason.evidence.document_id)} · ${reason.evidence.publication_year ?? ""} · ${isSeries ? "series anchor" : "child"} ${escapeHtml(reason.child_paper_label.replace(/\s*\(series anchor\)$/i, ""))}</p>
          <h2>${escapeHtml(reason.source_title ?? reason.evidence.document_id)}</h2>
        </div>
        <a class="evidence-link" target="_blank" rel="noreferrer" href="${escapeHtml(reason.evidence.citation_url)}">Source evidence ↗</a>
      </div>

      <div class="reason-callout">
        <span class="label">Extracted design assertion · ${escapeHtml(humanize(reason.assertion_class))}</span>
        <p>${escapeHtml(reason.reason_statement)}</p>
      </div>

      ${isSeries ? `<div class="series-scope"><span class="label">Named series scope · ${memberLabels.length} compounds</span><p>${escapeHtml(memberLabels.join(", "))}</p><small>One evidence episode; these are named members, not ${memberLabels.length} independent rationales.</small></div>` : ""}

      <div class="flow">
        ${moleculeCard(parent, "Comparator")}
        <div class="arrow" aria-hidden="true">→</div>
        ${moleculeCard(child, isSeries ? "Representative child" : "Reason-bearing child")}
      </div>
      <div class="transform"><span class="label">Structural witness</span><br><code>${escapeHtml(transformation)}</code></div>

      <div class="grid-two">
        <div class="info-card">
          <span class="label">Relationship semantics</span>
          <p><strong>${escapeHtml(relationshipLabel)}</strong></p>
          <p>${relationship?.historical ? "The source supports historical starting-compound lineage." : "No historical synthesis lineage is inferred."}</p>
          ${additionalRelationships.length ? `<p class="scope-note">Also author-linked: ${additionalRelationships.map((link) => `${escapeHtml(link.parent_paper_label ?? link.parent_chembl_id)} → ${escapeHtml(link.child_paper_label ?? link.child_chembl_id)}`).join("; ")}</p>` : ""}
        </div>
        <div class="info-card">
          <span class="label">Declared objective</span>
          <p class="intent-property">${escapeHtml(humanize(reason.intended_property.name))} · ${escapeHtml(reason.intended_property.direction)}</p>
          <p>${escapeHtml(reason.stated_modification)}</p>
        </div>
      </div>

      <section class="outcomes">
        <h3>Assay-matched outcomes</h3>
        ${outcomes.length ? `<div class="table-wrap"><table>
          <thead><tr><th>Endpoint / state</th><th>Measured comparator</th><th>Measured child</th><th>Result</th><th>Intent relationship</th></tr></thead>
          <tbody>${outcomes.map((row) => `<tr>
            <td><strong>${escapeHtml(humanize(row.endpoint))}</strong><br><span class="scope-note">${escapeHtml(row.state)}</span></td>
            <td class="measurement"><strong>${escapeHtml(compounds.get(row.parent_chembl_id)?.paper_label ?? row.parent_chembl_id)}</strong><br>${formatMeasurement(row.parent_relation, row.parent_value, row.units)}</td>
            <td class="measurement"><strong>${escapeHtml(compounds.get(row.child_chembl_id)?.paper_label ?? row.child_chembl_id)}</strong><br>${formatMeasurement(row.child_relation, row.child_value, row.units)}</td>
            <td><span class="badge ${escapeHtml(row.classification)}">${escapeHtml(row.classification)}</span><br><span class="scope-note">Δ ${escapeHtml(deltaText(row))}</span></td>
            <td class="scope-note">${escapeHtml(humanize(row.relation_to_stated_intent))}</td>
          </tr>`).join("")}</tbody>
        </table></div>` : `<div class="empty compact">No matched parent–child assay pair is available for this episode.</div>`}
        ${facts.length ? `<div class="fact-strip"><span class="label">Child-only outcome facts</span>${facts.map((fact) => `<p><strong>${escapeHtml(humanize(fact.endpoint))}</strong> · ${escapeHtml(fact.state)} · ${formatMeasurement(fact.relation, fact.value, fact.units)}</p>`).join("")}</div>` : ""}
      </section>

      <section class="candidates">
        <h3>Inferred candidate ranking</h3>
        ${topCandidates.length ? topCandidates.map((candidate) => `<div class="candidate-row">
          <strong>#${candidate.parent_rank}</strong>
          <span>${escapeHtml(candidate.parent_label)}<br><span class="scope-note">${escapeHtml(humanize(candidate.edge_semantics))}</span></span>
          <code>${escapeHtml(candidate.transformation)}</code>
          <span class="score">${Number(candidate.scores.ranking_score_uncalibrated).toFixed(3)}</span>
        </div>`).join("") : `<div class="empty">No candidate passes the primary MMP rule.</div>`}
      </section>`;
  }

  ["#search", "#class-filter", "#outcome-filter"].forEach((selector) => $(selector).addEventListener("input", () => {
    renderList();
    renderDetail();
  }));
  renderMetrics();
  populateFilters();
  renderList();
  renderDetail();
})();

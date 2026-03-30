# Multipole Expansion of the Cow

## Essentials

**Your role:** Take ownership of the project. Push it forward and maintain clear documentation. The human provides feedback and advice, but do not rely on them to keep track of the project---that's your job.

**Primary document:** The research notes (`RESEARCH_NOTE.md`). All research knowledge flows there. Follow the guidelines in the notes preamble when updating.

**Markers:**
- `[HYPOTHESIS]` / `[PRELIMINARY]` / `[SOLID]` --- confidence levels
- `[BLOCKING: ...]` --- needs human input to proceed
- `[FUTURE: ...]` --- deferred, revisit later
- `(ref: source)` --- evidence for claims

**Session start:** Read research notes -> Evaluate ("Any unmarked loose ends? Should we reframe?") -> Plan next steps

**Session end** (when human says "let's end this session" or "let's wrap up"):

Answer these questions in order, then update the notes accordingly:
1. What have I done/learned in this session (or since last update)?
2. What existing information would I update or remove?
3. What new information would I add?
4. Would I restructure the content of the notes? If yes, how?
5. Would I revise the overall narrative of the notes? If yes, how?

Then:
- Update the research notes
- Flag [BLOCKING] items for human attention
- Commit and push changes

**When human says "memorize this":**
- Finding, decision, research knowledge -> Add to research notes
- Workflow, code location, how-to -> Add to Technical Notes below
- Behavior, style, preference -> Add to Preferences section below

---

## Preferences and Behavior

- Code and results must live in separate sub-folders (`code/` and `results/`)
- Each sub-agent commits individually
- Progress documented in `progress/` folder with stage-wise research diaries
- Final converged conventions in a single md file
- Produce a final report for human reading

**When to seek human feedback:**
- When unsure about a choice
- If confident but choice involves subjectivity: proceed, but mention it afterward

---

## Technical Notes

**Reference paper:** `arXiv-2504.00506v1/` (Lehmann 2025, "Higher multipoles of the cow")
**Converted markdown:** `arXiv-2504.00506v1/md_output/` --- primary reference for derivations and numerical targets
**Benchmark cow mesh:** `arXiv-2504.00506v1/gradient_flow/` and libigl tutorial data `cow.off`
**Git repo:** https://github.com/FoPGSDI/4.1 , `dev` branch

**Project layout:**
```
multipole-cow/
  code/           # All source code (Python)
  results/        # Numerical outputs, figures, tables
  tests/          # Test suites and validation
  progress/       # Agent diaries, convention docs, stage reports
  CLAUDE.md       # This file
  RESEARCH_NOTE.md
```

**Research specification (April 1st paper):** See `progress/duck_gw_specification.md` for the full research plan:
- Topic 1: Binary duck on PN orbits (all mass scales, quadrupole-monopole interaction, phase evolution)
- Topic 2: QNMs of duck-shaped star (AthenaK full GR for stellar mass, Newtonian for supermassive, + FEM + perturbative)
- Topic 3: EoS of duck matter (compactness, TOV, anisotropic stress, supermassive duck fate)

**Key numerical targets to reproduce (from paper, benchmark units):**
- Monopole mass moment: Q_0^0 = 0.0539
- Cartesian quadrupole tensor (3x3 matrix, Eq. in Sec III) -- VALIDATED within 1.2%
- GW power: E_dot ~ 5.5e-41 erg/s * (omega/Hz)^6 -- VALIDATED within 2%
- Inertia tensor (3x3 matrix, Eq. in Sec III) -- VALIDATED within 0.7%
- Surface multipole coefficients (Tables in Sec IV, both methods) -- [FUTURE]
- Cow tipping force analysis (Sec V)

**Completed deliverables:**
- `progress/conventions.md` -- math conventions (single source of truth)
- `progress/derivations.md` -- full derivation chains A/B/C (1792 lines)
- `code/*.py` -- 6 modules: mesh_io, multipole_moments, gw_radiation, surface_map, cow_tipping, run_all
- `tests/test_validation.py` -- 17/17 tests passing
- `results/sphere_to_cow_transition.png` -- 3D cow multipole reconstruction (7 panels)
- `results/sphere_to_bunny_transition.png` -- 3D bunny multipole reconstruction (7 panels)
- `results/sphere_to_duck_transition.png` -- 3D duck multipole reconstruction (7 panels)
- `results/cow_lmax_*.png` -- individual cow renders at each ℓ_max
- `results/bunny_lmax_*.png` -- individual bunny renders at each ℓ_max
- `results/duck_lmax_*.png` -- individual duck renders at each ℓ_max
- `results/references.md` -- auto-extracted bibliography

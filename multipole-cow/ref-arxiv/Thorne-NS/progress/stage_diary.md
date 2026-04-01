# Stage Diary: Thorne & Campolattaro (1967) PDF Conversion

## Paper
- **Title**: Non-Radial Pulsation of General-Relativistic Stellar Models. I. Analytic Analysis for l >= 2
- **Authors**: Kip S. Thorne & Alfonso Campolattaro
- **Journal**: The Astrophysical Journal, Vol. 149, September 1967, pp. 591-611
- **Pages**: 22 (21 content + 1 blank)

## Task Overview
Convert the scanned PDF (1967ApJ...149..591T.pdf) into:
1. A complete LaTeX file (main.tex)
2. A BibTeX file (references.bib)
3. A convention/notation file (conventions.md)
4. Structured markdown files by section (md_output/)

## Stage 1: PDF Splitting
- Split 22-page PDF into individual PNG files at 300 DPI using pdftoppm
- Output: pages/page-01.png through page-22.png

## Stage 2: Page-by-Page Transcription (22 parallel agents)
- Launched 22 agents simultaneously, each reading one page image and producing a .tex file
- Each agent was given context about what equations/sections to expect on their page
- All 22 agents completed successfully
- Output: tex/page01.tex through tex/page22.tex
- Key challenge: 1967 journal scan quality required agents to do careful OCR-like reading
- Flagged uncertainties: eq (C7) on page 19 had several coefficients that were partially illegible

## Stage 3: Merge and BibTeX (2 parallel agents)
- **Merge agent**: Read all 22 per-page tex files, joined split sentences at page boundaries, removed duplicate section headers, added proper LaTeX preamble, replaced thebibliography with \bibliography{references}
- **BibTeX agent**: Converted 26 \bibitem entries from page21.tex into proper BibTeX format
- Output: main.tex (1522 lines), references.bib (214 lines, 26 entries)

## Stage 4: Convention File
- Created conventions.md documenting all notation, symbols, coordinate choices, derivative notation, perturbation variables, boundary conditions, and stability criteria
- Manually written based on full paper reading

## Stage 5: LaTeX-to-Markdown Conversion (10 parallel agents)
- Split main.tex at section boundaries into 11 markdown files (Sec I-VI + App A-E)
- Plus content.md (table of contents) and convention.md
- Each agent applied consistent conversion rules:
  - Section headers to # / ## / ###
  - LaTeX math preserved in $...$ and $$...$$
  - Citations, italics, quotes converted to markdown
  - Footnotes converted to [^N] format
- Output: md_output/ with 14 files total

## Stage 6: Cross-Check and Verification
- Deployed verification agent to compare critical equations against PDF images
- Checked: eq (1), (3b), (7b), (14a), (22) and completeness of all numbered equations

## File Inventory

### Root directory
| File | Description | Lines |
|------|-------------|-------|
| main.tex | Complete merged LaTeX | 1522 |
| references.bib | BibTeX bibliography | 214 |
| conventions.md | Notation reference | ~130 |
| input-Thorne-NS.md | Task specification | 17 |

### pages/
- page-01.png through page-22.png (22 PNG files at 300 DPI)

### tex/
- page01.tex through page22.tex (22 per-page transcriptions)

### md_output/
| File | Content |
|------|---------|
| content.md | Table of contents with links |
| convention.md | Notation & conventions |
| 01_motivation.md | Abstract + Section I |
| 02_arbitrary_small_perturbations.md | Section II (a-d) |
| 03_even_parity_normal_modes.md | Section III (a-c) |
| 04_real_pulsations.md | Section IV (a-b) |
| 05_methods.md | Section V |
| 06_where_do_we_go.md | Section VI + Acknowledgments |
| appendix_a.md | Spherical harmonics & RW gauge |
| appendix_b.md | Odd-parity equations of motion |
| appendix_c.md | Even-parity equations of motion |
| appendix_d.md | Boundary conditions |
| appendix_e.md | Eigenfunctions at r=infinity |

### progress/
- stage_diary.md (this file)
- verification_report.md (cross-check results)

## Known Issues
1. **Eq (C7)**: Several coefficients on page 19 (p. 609) were partially illegible in the scan. The tex file contains inline comments flagging uncertain terms. A higher-resolution scan would be needed for full verification.
2. **Eq (18)**: The asymptotic wave form has some notation (SM ln r) that may need clarification — likely omega*M*ln(r).
3. **Eq (30)**: The radiated power equation transcription from page 13 may have OCR artifacts in the exponent.
4. **Page 12 (eq 25-26)**: The kinetic energy expressions involve labeled factors above each term; the underbrace rendering approximates the original layout.

## Total Agent Usage
- 22 page-transcription agents (parallel)
- 2 merge/bib agents (parallel)
- 10 markdown conversion agents (parallel)
- 1 verification agent
- **Total: 35 sub-agents**

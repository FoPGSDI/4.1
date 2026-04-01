**Central Task**
* Central target: convert `/data/haiyangw/claude/4.1/multipole-cow/ref-arxiv/Thorne-NS/1967ApJ...149..591T.pdf` paper to tex and bib file, and md files by section
* 1. separate the pdf to single pages pdf or png
* 2. for each page, generate a tex file (translated)
* 3. merger the tex file and add the bibtex file
* 4. check and create a convention file
* 5. run skill /latex-to-md
* In editing mode, deploy agents to repeat the task and cross-check.
* Final agent: check and finalize the results.

**Files and Links**
* `/data/haiyangw/claude/4.1/multipole-cow/ref-arxiv/Thorne-NS` contains the source paper and is the working directory

**Commit Principles**
* The progress of each sub-agent should be documented in md files; the final converged convention should be documented in one md file; stage progress should be documented in one md file following the original guideline. These should live in a `progress/` folder.
* We should have stage-wise research diaries and a final report for human reading.
